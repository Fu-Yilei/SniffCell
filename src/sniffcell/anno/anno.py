import os
import json
import pandas as pd
from sniffcell.anno.kmeans import kmeans_cluster_cells
from sniffcell.anno.methyl_matrix import methyl_matrix_from_bam
from sniffcell.anno.filter_bed_based_on_variants import filter_bed_based_on_variants
from sniffcell.anno.vcf_to_df import read_vcf_to_df
from sniffcell.anno.variant_assignment import assign_sv_celltypes
from tqdm import tqdm
import multiprocessing as mp
import numpy as np
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(processName)s] %(levelname)s: %(message)s")


def _write_anno_run_manifest(
    *,
    output_dir: str,
    bam: str,
    vcf: str,
    reference: str,
    bed: str,
    reads_classification: str,
    blocks_classification: str,
    window: int,
    threads: int,
    kanpig_read_names: str | None,
    read_assignment_mode: str,
) -> str:
    manifest_path = os.path.join(output_dir, "anno_run_manifest.json")
    payload = {
        "command": "anno",
        "version": "v1",
        "inputs": {
            "bam": os.path.abspath(bam),
            "vcf": os.path.abspath(vcf),
            "reference": os.path.abspath(reference),
            "bed": os.path.abspath(bed),
            "kanpig_read_names": os.path.abspath(kanpig_read_names) if kanpig_read_names else None,
        },
        "runtime": {
            "window": int(window),
            "threads": int(threads),
            "read_assignment_mode": str(read_assignment_mode),
        },
        "outputs": {
            "output_dir": os.path.abspath(output_dir),
            "reads_classification": os.path.abspath(reads_classification),
            "blocks_classification": os.path.abspath(blocks_classification),
            "sv_assignment": os.path.abspath(os.path.join(output_dir, "sv_assignment.tsv")),
            "sv_assignment_readable": os.path.abspath(os.path.join(output_dir, "sv_assignment_readable.tsv")),
            "sv_assignment_readable_long": os.path.abspath(os.path.join(output_dir, "sv_assignment_readable_long.tsv")),
        },
    }
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
    return manifest_path

def _parse_pipe_values(value) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        return [str(v).strip() for v in value if str(v).strip()]
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, (str, bytes)):
        text = str(value).strip()
        if not text:
            return []
        return [v.strip() for v in text.split("|") if v.strip()]
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    return [v.strip() for v in text.split("|") if v.strip()]


def _resolve_cell_types_and_targets(row: dict, best_group: str) -> tuple[list[str], list[str]]:
    # Prefer explicit global schema produced by find.
    cell_types = _parse_pipe_values(row.get("code_order"))
    has_explicit_schema = bool(cell_types)
    if not cell_types:
        cell_types = _parse_pipe_values(row.get("leaf_order"))
        has_explicit_schema = bool(cell_types)

    # Backward-compatible fallback for legacy find outputs.
    if not cell_types:
        for k, v in row.items():
            if not isinstance(k, str):
                continue
            if not k.startswith("mean_"):
                continue
            if k in ("mean_best_value", "mean_rest_value", "mean_margin"):
                continue
            if pd.isna(v):
                continue
            cell_types.append(k[len("mean_"):])

    if (not has_explicit_schema) and best_group and best_group not in cell_types:
        cell_types.append(best_group)
    cell_types = list(dict.fromkeys(cell_types))

    # Explicit target set for this ctDMR (single cell type or a combination).
    target_cell_types = _parse_pipe_values(row.get("best_group_leaves"))
    if not target_cell_types:
        target_cell_types = [best_group] if best_group else []

    # Keep target order consistent with code_order.
    target_set = set(target_cell_types)
    target_cell_types = [ct for ct in cell_types if ct in target_set]
    if not target_cell_types and best_group in cell_types:
        target_cell_types = [best_group]

    return cell_types, target_cell_types


def _parse_celltype_metric_map(value, cast_type):
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return {}
    text = str(value).strip()
    if not text:
        return {}

    out: dict[str, object] = {}
    for token in text.split(";"):
        token = token.strip()
        if not token or ":" not in token:
            continue
        key, raw = token.rsplit(":", 1)
        key = key.strip()
        raw = raw.strip()
        if not key or not raw:
            continue
        try:
            out[key] = cast_type(raw)
        except (TypeError, ValueError):
            continue
    return out


def _to_finite_float(value) -> float | None:
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(out):
        return None
    return out


def _find_dmr_query_interval(start: int, end: int) -> tuple[int, int]:
    """
    find ctDMR BED rows are G-anchored on the left edge in the atlas index.
    Shift the query interval left by 1 bp so anno re-extracts the left-edge CpG.
    """
    if start <= 0:
        return start, end
    return start - 1, end


def _assign_target_by_reference_mean(
    read_means: np.ndarray,
    *,
    mean_best_value,
    mean_rest_value,
) -> np.ndarray | None:
    mb = _to_finite_float(mean_best_value)
    mr = _to_finite_float(mean_rest_value)
    if mb is None or mr is None:
        return None
    d_best = np.abs(np.asarray(read_means, dtype=float) - mb)
    d_rest = np.abs(np.asarray(read_means, dtype=float) - mr)
    return d_best <= d_rest


def _build_sv_readable_reports(sv_assignment_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_cols = [
        "id",
        "sv_chr",
        "sv_pos",
        "sv_len",
        "vaf",
        "n_supporting",
        "n_overlapped",
        "overlap_pct",
        "majority_pct",
        "classified_celltypes",
        "classified_celltype_count",
        "classified_celltype_counts",
        "classified_celltype_fractions",
        "classification_summary",
        "is_multi_celltype_link",
    ]
    long_cols = [
        "id",
        "sv_chr",
        "sv_pos",
        "sv_len",
        "celltype",
        "rank",
        "supporting_read_count",
        "supporting_read_fraction",
        "n_supporting",
        "n_overlapped",
        "overlap_pct",
    ]
    if sv_assignment_df.empty:
        return pd.DataFrame(columns=summary_cols), pd.DataFrame(columns=long_cols)

    summary_rows = []
    long_rows = []
    for _, row in sv_assignment_df.iterrows():
        linked = _parse_pipe_values(row.get("linked_celltypes", ""))
        count_map = _parse_celltype_metric_map(row.get("linked_celltype_counts", ""), int)
        frac_map = _parse_celltype_metric_map(row.get("linked_celltype_fractions", ""), float)

        for ct in count_map:
            if ct not in linked:
                linked.append(ct)
        for ct in frac_map:
            if ct not in linked:
                linked.append(ct)

        detail_items = []
        for rank, ct in enumerate(linked, start=1):
            ct_count = count_map.get(ct, pd.NA)
            ct_frac = frac_map.get(ct, pd.NA)
            has_count = pd.notna(ct_count)
            has_frac = pd.notna(ct_frac)

            if has_count and has_frac:
                detail_items.append(f"{ct} ({ct_count} reads, {float(ct_frac):.3f})")
            elif has_count:
                detail_items.append(f"{ct} ({ct_count} reads)")
            elif has_frac:
                detail_items.append(f"{ct} ({float(ct_frac):.3f})")
            else:
                detail_items.append(ct)

            long_rows.append(
                {
                    "id": row.get("id", ""),
                    "sv_chr": row.get("sv_chr", ""),
                    "sv_pos": row.get("sv_pos", pd.NA),
                    "sv_len": row.get("sv_len", pd.NA),
                    "celltype": ct,
                    "rank": rank,
                    "supporting_read_count": ct_count,
                    "supporting_read_fraction": ct_frac,
                    "n_supporting": row.get("n_supporting", pd.NA),
                    "n_overlapped": row.get("n_overlapped", pd.NA),
                    "overlap_pct": row.get("overlap_pct", pd.NA),
                }
            )

        summary_rows.append(
            {
                "id": row.get("id", ""),
                "sv_chr": row.get("sv_chr", ""),
                "sv_pos": row.get("sv_pos", pd.NA),
                "sv_len": row.get("sv_len", pd.NA),
                "vaf": row.get("vaf", pd.NA),
                "n_supporting": row.get("n_supporting", pd.NA),
                "n_overlapped": row.get("n_overlapped", pd.NA),
                "overlap_pct": row.get("overlap_pct", pd.NA),
                "majority_pct": row.get("majority_pct", pd.NA),
                "classified_celltypes": "|".join(linked),
                "classified_celltype_count": len(linked),
                "classified_celltype_counts": row.get("linked_celltype_counts", ""),
                "classified_celltype_fractions": row.get("linked_celltype_fractions", ""),
                "classification_summary": "; ".join(detail_items),
                "is_multi_celltype_link": row.get("is_multi_celltype_link", pd.NA),
            }
        )

    summary_df = pd.DataFrame(summary_rows, columns=summary_cols)
    long_df = pd.DataFrame(long_rows, columns=long_cols)
    if not summary_df.empty:
        summary_df["classified_celltypes"] = summary_df["classified_celltypes"].astype("string")
        summary_df["classified_celltype_counts"] = summary_df["classified_celltype_counts"].astype("string")
        summary_df["classified_celltype_fractions"] = summary_df["classified_celltype_fractions"].astype("string")
        summary_df["classification_summary"] = summary_df["classification_summary"].astype("string")
    if not long_df.empty:
        long_df["celltype"] = long_df["celltype"].astype("string")
    return summary_df, long_df

def _one_dmr(args):
    # args: (row_dict, input_bam, reference_fasta)
    logger = logging.getLogger("anno._one_dmr")

    if len(args) >= 4:
        row, input_file, reference, read_assignment_mode = args[:4]
    else:
        row, input_file, reference = args
        read_assignment_mode = "closest_reference_mean"
    read_assignment_mode = str(read_assignment_mode).strip().lower()
    if read_assignment_mode not in {"closest_reference_mean", "kmeans"}:
        raise ValueError("read_assignment_mode must be one of: closest_reference_mean, kmeans")

    chrom = str(row["chr"])
    start = int(row["start"])
    end   = int(row["end"])

    best_group = str(row["best_group"])
    best_dir   = row.get("best_dir", None)
    cell_types, target_cell_types = _resolve_cell_types_and_targets(row, best_group)
    target_set = set(target_cell_types)

    # logger.info(f"[{chrom}:{start}-{end}] best_group={best_group} best_dir={best_dir} "
    #             f"cell_types={cell_types} (n={len(cell_types)})")

    try:
        # load methylation matrix + CpG positions
        query_start, query_end = _find_dmr_query_interval(start, end)
        mm, cpgs = methyl_matrix_from_bam(
            input_file,
            reference,
            chrom=chrom,
            start=query_start,
            end=query_end,
            return_positions=True,
        )
        n_reads_raw = 0 if mm is None else mm.shape[0]
        n_cpgs = len(cpgs)
        if n_cpgs == 0:
            logger.warning(f"[{chrom}:{start}-{end}] no CpGs found; skipping")
            return None

        # drop rows that are entirely NaN across CpGs
        mm = mm.dropna(how="all")
        if mm.empty or mm.shape[0] < 2:
            logger.warning(f"[{chrom}:{start}-{end}] usable_reads={mm.shape[0] if not mm.empty else 0} "
                           f"(raw={n_reads_raw}) < 2; skipping")
            return None

        # CpG bounds from cpgs
        cpgstart = int(cpgs[0])
        cpgend   = int(cpgs[-1])

        # read names from MultiIndex level 0 if present
        if isinstance(mm.index, pd.MultiIndex) and "read_name" in mm.index.names:
            readnames = mm.index.get_level_values("read_name").astype(str).values
        else:
            readnames = (mm.index.astype(str).values if mm.index.dtype == object
                         else np.array([f"read_{i}" for i in range(len(mm))], dtype=str))

        X_imp = mm.astype(float).copy().fillna(mm.astype(float).mean())
        read_mean = X_imp.mean(axis=1).values

        # Assign each read to best_group vs other_group per ctDMR.
        if read_assignment_mode == "closest_reference_mean":
            mask_target = _assign_target_by_reference_mean(
                read_mean,
                mean_best_value=row.get("mean_best_value", np.nan),
                mean_rest_value=row.get("mean_rest_value", np.nan),
            )
            if mask_target is None:
                logger.warning(
                    "[%s:%d-%d] closest_reference_mean unavailable (missing mean_best_value/mean_rest_value); "
                    "falling back to kmeans",
                    chrom,
                    start,
                    end,
                )
                dmr_row = {
                    "best_group": best_group,
                    "best_dir": best_dir,
                    "mean_best_value": row.get("mean_best_value", np.nan),
                    "mean_rest_value": row.get("mean_rest_value", np.nan),
                }
                out = kmeans_cluster_cells(mm, dmr_row=dmr_row)
                mask_target = (out["celltype_or_other"].astype(str).str.lower()
                               == best_group.strip().lower()).values
        else:
            dmr_row = {
                "best_group": best_group,
                "best_dir": best_dir,
                "mean_best_value": row.get("mean_best_value", np.nan),
                "mean_rest_value": row.get("mean_rest_value", np.nan),
            }
            out = kmeans_cluster_cells(mm, dmr_row=dmr_row)
            mask_target = (out["celltype_or_other"].astype(str).str.lower()
                           == best_group.strip().lower()).values

        # Build hierarchical bit codes in a global schema from code_order.
        target_bits = ["1" if ct in target_set else "0" for ct in cell_types]
        if not any(bit == "1" for bit in target_bits) and best_group in cell_types:
            idx = cell_types.index(best_group)
            target_bits[idx] = "1"
            target_set = {best_group}
        other_bits = ["0" if bit == "1" else "1" for bit in target_bits]
        target_code = "".join(target_bits)
        other_code  = "".join(other_bits)
        code_col = np.where(mask_target, target_code, other_code)

        # --- per-read assignments (each read = one row / index) ---
        assign_df = pd.DataFrame({
            "chr": chrom,
            "start": start,
            "end": end,
            "cpgstart": cpgstart,
            "cpgend": cpgend,
            "best_group": best_group,
            "other_group": str(row.get("other_group", "")),
            "is_best_group": mask_target,
            "code_order": "|".join(cell_types),
            "best_group_leaves": "|".join(target_cell_types),
            "other_group_leaves": str(row.get("other_group_leaves", "")),
            "hyper_group_leaves": str(row.get("hyper_group_leaves", "")),
            "hypo_group_leaves": str(row.get("hypo_group_leaves", "")),
            "code": code_col,
        }, index=pd.Index(readnames, name="readname"))

        # per-block means (per-read mean methylation, then avg by target vs other)
        tgt_mean = float(np.nanmean(read_mean[mask_target])) if mask_target.any() else np.nan
        oth_mean = float(np.nanmean(read_mean[~mask_target])) if (~mask_target).any() else np.nan

        # logger.info(f"[{chrom}:{start}-{end}] target_mean={tgt_mean:.4f} other_mean={oth_mean:.4f} "
        #             f"cpg_bounds={cpgstart}-{cpgend}")

        state_payload = {
            "chr": chrom, "start": start, "end": end,
            "cpgstart": cpgstart, "cpgend": cpgend,
        }
        for ct in cell_types:
            state_payload[f"{ct}_methylation"] = tgt_mean if ct in target_set else oth_mean
        state_df = pd.DataFrame([state_payload])

        return assign_df, state_df

    except Exception as e:
        logger.exception(f"[{chrom}:{start}-{end}] failed with error")
        return None

def sv_anno(args):
    logger = logging.getLogger("anno.sv_anno")
    logger.info("Starting SV annotation from pre-annotated reads")
    if args.command == "svanno":
        input_file = args.input
        output_arg = str(args.output)
        # Backward compatibility: treat extension-less output as directory.
        if os.path.isdir(output_arg) or os.path.splitext(output_arg)[1] == "":
            output_dir = output_arg
            sv_assignment_path = os.path.join(output_dir, "sv_assignment.tsv")
        else:
            sv_assignment_path = output_arg
            output_dir = os.path.dirname(output_arg) or "."
    else:
        output_dir = args.output
        input_file = os.path.join(output_dir, "reads_classification.tsv")
        sv_assignment_path = os.path.join(output_dir, "sv_assignment.tsv")

    os.makedirs(output_dir, exist_ok=True)
    if args.kanpig_read_names is not None:
        logger.info(f"Using kanpig read names from: {args.kanpig_read_names}")
    else:
        logger.info("No kanpig read names provided; using Sniffles read names from VCF")
    read_assign_df = pd.read_csv(
        input_file,
        sep="\t",
        index_col=0,
        dtype={
            "code": "string",
            "code_order": "string",
            "best_group": "string",
            "best_group_leaves": "string",
        },
    )

    evidence_mode = str(getattr(args, "evidence_mode", "all_rows")).strip().lower()
    if evidence_mode not in {"all_rows", "per_read"}:
        raise ValueError("evidence_mode must be one of: all_rows, per_read")
    unique_reads_for_overlap = evidence_mode == "per_read"

    min_overlap_pct = float(getattr(args, "min_overlap_pct", 0.0))
    min_agreement_pct = float(getattr(args, "min_agreement_pct", 1.0))
    if not (0.0 <= min_overlap_pct <= 1.0):
        raise ValueError("min_overlap_pct must be in [0, 1]")
    if not (0.0 <= min_agreement_pct <= 1.0):
        raise ValueError("min_agreement_pct must be in [0, 1]")

    logger.info(
        "SV assignment settings: evidence_mode=%s unique_reads_for_overlap=%s min_overlap_pct=%.3f min_agreement_pct=%.3f",
        evidence_mode,
        unique_reads_for_overlap,
        min_overlap_pct,
        min_agreement_pct,
    )

    sv_assignment_df = assign_sv_celltypes(
        read_vcf_to_df(args.vcf, kanpig_read_names=args.kanpig_read_names),
        read_assign_df,
        window=int(getattr(args, "window", 5000)),
        min_overlap_pct=min_overlap_pct,
        min_agreement_pct=min_agreement_pct,
        unique_reads_for_overlap=unique_reads_for_overlap,
    )
    sv_assignment_df.to_csv(sv_assignment_path, sep="\t", index=False)

    readable_report_df, readable_report_long_df = _build_sv_readable_reports(sv_assignment_df)
    readable_report_path = os.path.join(output_dir, "sv_assignment_readable.tsv")
    readable_report_long_path = os.path.join(output_dir, "sv_assignment_readable_long.tsv")
    readable_report_df.to_csv(readable_report_path, sep="\t", index=False)
    readable_report_long_df.to_csv(readable_report_long_path, sep="\t", index=False)
    logger.info(f"Wrote SV assignment report: {sv_assignment_path}")
    logger.info(f"Wrote readable SV report: {readable_report_path}")
    logger.info(f"Wrote long-format readable SV report: {readable_report_long_path}")




def anno_main(args):
    # print(args)
    # return
    logger = logging.getLogger("anno.main")

    bed_input  = args.bed
    base_out   = args.output    # writes <output>.reads.tsv and <output>.blocks.tsv
    input_file = args.input
    reference  = args.reference
    threads    = int(args.threads)
    window     = int(args.window)
    read_assignment_mode = str(getattr(args, "read_assignment_mode", "closest_reference_mean")).strip().lower()
    if read_assignment_mode not in {"closest_reference_mean", "kmeans"}:
        raise ValueError("read_assignment_mode must be one of: closest_reference_mean, kmeans")
    logger.info(
        f"Starting annotation: bed={bed_input} bam={input_file} ref={reference} "
        f"threads={threads} read_assignment_mode={read_assignment_mode} out_base={base_out}"
    )

    # Output paths
    reads_out  = os.path.join(base_out, "reads_classification.tsv")
    blocks_out = os.path.join(base_out, "blocks_classification.tsv")
    manifest_path = _write_anno_run_manifest(
        output_dir=base_out,
        bam=input_file,
        vcf=args.vcf,
        reference=reference,
        bed=bed_input,
        reads_classification=reads_out,
        blocks_classification=blocks_out,
        window=window,
        threads=threads,
        kanpig_read_names=getattr(args, "kanpig_read_names", None),
        read_assignment_mode=read_assignment_mode,
    )
    logger.info("Wrote anno run manifest: %s", manifest_path)

    # Load and (optionally) filter BED
    bed = pd.read_csv(bed_input, sep="\t")
    if not bed.empty and isinstance(bed.columns[0], str) and bed.columns[0].startswith("#"):
        bed.rename(columns={bed.columns[0]: bed.columns[0].lstrip("#")}, inplace=True)
    bed = bed.drop_duplicates(ignore_index=True)
    logger.info(f"Loaded BED with {len(bed)} unique DMR rows")

    sv_df = read_vcf_to_df(args.vcf)
    filtered_bed = filter_bed_based_on_variants(bed, sv_df=sv_df, window=window)

    for col in ["chr", "start", "end", "best_group", "best_dir"]:
        if col not in filtered_bed.columns:
            logger.error(f"BED missing required column: {col}")
            raise ValueError(f"BED missing required column: {col}")

    n_tasks = len(filtered_bed)
    logger.info(f"Filtered BED to {n_tasks} DMRs after variant overlap filtering, window size = {window}")

    tasks = [(dict(row), input_file, reference, read_assignment_mode) for _, row in filtered_bed.iterrows()]

    # --- Prepare outputs: truncate files and reset header flags ---
    # We'll only write headers on the first real chunk for each file.
    open(reads_out,  "w").close()
    open(blocks_out, "w").close()
    reads_header_written  = False
    blocks_header_written = False
    blocks_cols_locked: list[str] | None = None  # we lock schema to the first block we see

    # Stream results and append immediately
    with mp.Pool(threads) as pool:
        for res in tqdm(pool.imap(_one_dmr, tasks, chunksize=1),
                        total=n_tasks, desc="Processing DMRs"):
            if res is None:
                continue

            a_df, s_df = res

            # --- APPEND READS ---
            if a_df is not None and not a_df.empty:
                if not reads_header_written:
                    # first write: include header and index (readname)
                    a_df.to_csv(reads_out, sep="\t", index=True, mode="a", header=True)
                    reads_header_written = True
                else:
                    a_df.to_csv(reads_out, sep="\t", index=True, mode="a", header=False)

            # --- APPEND BLOCKS (variable columns across DMRs) ---
            if s_df is not None and not s_df.empty:
                if not blocks_header_written:
                    # lock the schema to the first encountered block columns
                    blocks_cols_locked = list(s_df.columns)
                    s_df.to_csv(blocks_out, sep="\t", index=False, mode="a", header=True)
                    blocks_header_written = True
                else:
                    # align columns to locked header; drop extras, add missing as NaN
                    assert blocks_cols_locked is not None
                    s_df_aligned = s_df.reindex(columns=blocks_cols_locked)
                    s_df_aligned.to_csv(blocks_out, sep="\t", index=False, mode="a", header=False)

    # If nothing was written, emit empty files with headers to be friendly downstream
    if not reads_header_written:
        empty_reads = pd.DataFrame(
            columns=[
                "chr",
                "start",
                "end",
                "cpgstart",
                "cpgend",
                "best_group",
                "other_group",
                "is_best_group",
                "code_order",
                "best_group_leaves",
                "other_group_leaves",
                "hyper_group_leaves",
                "hypo_group_leaves",
                "code",
            ]
        )
        empty_reads.index.name = "readname"
        empty_reads.to_csv(reads_out, sep="\t", index=True, header=True)
        logger.warning("No per-read assignments generated; wrote empty reads header only")

    if not blocks_header_written:
        pd.DataFrame(columns=["chr","start","end","cpgstart","cpgend"]).to_csv(
            blocks_out, sep="\t", index=False, header=True
        )
        logger.warning("No block states generated; wrote empty blocks header only")
    sv_anno(args)
    logger.info("Annotation complete")
