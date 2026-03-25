import os
import json
from pathlib import Path
import pandas as pd
import pysam
from sniffcell.anno.breakpoint_exclusion import validate_breakpoint_exclusion_frac
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

# ---------------------------------------------------------------------------
# Process-local BAM/FASTA handle cache
# Each worker process opens its own handle once and reuses it across ctDMR
# tasks, avoiding ~90K redundant pysam.AlignmentFile() calls per worker.
# ---------------------------------------------------------------------------
_process_bam_handle: "pysam.AlignmentFile | None" = None
_process_fasta_handle: "pysam.FastaFile | None" = None
_process_bam_path: str = ""
_process_fasta_path: str = ""


def _get_process_bam(path: str) -> "pysam.AlignmentFile":
    global _process_bam_handle, _process_bam_path
    if _process_bam_handle is None or _process_bam_path != path:
        if _process_bam_handle is not None:
            try:
                _process_bam_handle.close()
            except Exception:
                pass
        _process_bam_handle = pysam.AlignmentFile(path, "rb")
        _process_bam_path = path
    return _process_bam_handle


def _get_process_fasta(path: str) -> "pysam.FastaFile":
    global _process_fasta_handle, _process_fasta_path
    if _process_fasta_handle is None or _process_fasta_path != path:
        if _process_fasta_handle is not None:
            try:
                _process_fasta_handle.close()
            except Exception:
                pass
        _process_fasta_handle = pysam.FastaFile(path)
        _process_fasta_path = path
    return _process_fasta_handle


def _write_anno_run_manifest(
    *,
    output_dir: str,
    bam: str,
    variant_input: str,
    variant_source: str,
    reference: str,
    bed: str,
    reads_classification: str,
    blocks_classification: str,
    window: int,
    breakpoint_exclusion_frac: float,
    threads: int,
    kanpig_read_names: str | None,
    read_assignment_mode: str,
) -> str:
    manifest_path = os.path.join(output_dir, "anno_run_manifest.json")
    output_dir_abs = os.path.abspath(output_dir)
    variant_input_abs = os.path.abspath(variant_input)
    payload = {
        "command": "anno",
        "version": "v1",
        "inputs": {
            "bam": os.path.abspath(bam),
            "vcf": (variant_input_abs if variant_source == "vcf" else None),
            "variants": variant_input_abs,
            "variant_source": str(variant_source),
            "reference": os.path.abspath(reference),
            "bed": os.path.abspath(bed),
            "kanpig_read_names": os.path.abspath(kanpig_read_names) if kanpig_read_names else None,
        },
        "runtime": {
            "window": int(window),
            "breakpoint_exclusion_frac": float(breakpoint_exclusion_frac),
            "threads": int(threads),
            "read_assignment_mode": str(read_assignment_mode),
        },
        "outputs": {
            "output_dir": output_dir_abs,
            "reads_classification": os.path.abspath(reads_classification),
            "blocks_classification": os.path.abspath(blocks_classification),
            "variant_assignment": os.path.abspath(os.path.join(output_dir, "variant_assignment.tsv")),
            "variant_assignment_readable": os.path.abspath(os.path.join(output_dir, "variant_assignment_readable.tsv")),
            "variant_assignment_readable_long": os.path.abspath(os.path.join(output_dir, "variant_assignment_readable_long.tsv")),
            "sv_assignment": os.path.abspath(os.path.join(output_dir, "sv_assignment.tsv")),
            "sv_assignment_readable": os.path.abspath(os.path.join(output_dir, "sv_assignment_readable.tsv")),
            "sv_assignment_readable_long": os.path.abspath(os.path.join(output_dir, "sv_assignment_readable_long.tsv")),
        },
    }
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
    return manifest_path


def _detect_variant_source(path_text: str) -> str:
    path = Path(path_text)
    lowered = path.name.lower()
    if lowered.endswith((".vcf", ".vcf.gz", ".bcf")):
        return "vcf"
    if lowered.endswith((".tsv", ".tsv.gz", ".txt", ".txt.gz", ".bed", ".bed.gz", ".csv", ".csv.gz")):
        try:
            cols = set(pd.read_csv(path, sep="\t", nrows=0).columns)
        except Exception:
            cols = set()
        if {"chrom", "start", "end", "variant_id"}.issubset(cols):
            return "harmonized"
    return "vcf"


def _parse_json_read_names(value) -> list[str]:
    try:
        if pd.isna(value):
            return []
    except TypeError:
        pass
    if isinstance(value, (list, tuple, set)):
        return [str(v).strip() for v in value if str(v).strip()]
    text = str(value).strip()
    if not text or text in {".", "NA", "None", "null"}:
        return []
    try:
        parsed = json.loads(text)
    except Exception:
        parsed = None
    if isinstance(parsed, list):
        return [str(v).strip() for v in parsed if str(v).strip()]
    if "|" in text:
        return [tok.strip() for tok in text.split("|") if tok.strip()]
    if "," in text:
        return [tok.strip() for tok in text.split(",") if tok.strip()]
    return [text]


def _normalize_supporting_read_names(read_names: list[str], *, variant_class: str) -> list[str]:
    normalized: list[str] = []
    is_tr = str(variant_class).strip().upper() == "TR"
    for raw_name in read_names:
        name = str(raw_name).strip()
        if not name:
            continue
        if is_tr and "_chr" in name:
            name = name.split("_chr", 1)[0].strip()
        if name:
            normalized.append(name)
    return normalized


def read_harmonized_variants_to_df(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"chrom", "start", "end", "variant_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Harmonized variant TSV missing required columns: {sorted(missing)}")

    out = df.copy()
    for col, default in (
        ("variant_class", "SV"),
        ("variant_subtype", ""),
        ("category", ""),
        ("change_size_bp", pd.NA),
        ("group_a_alt_reads", pd.NA),
        ("group_b_alt_reads", pd.NA),
        ("group_a_read_names", "[]"),
        ("group_b_read_names", "[]"),
    ):
        if col not in out.columns:
            out[col] = default

    out["start"] = pd.to_numeric(out["start"], errors="coerce")
    out["end"] = pd.to_numeric(out["end"], errors="coerce")
    out = out.dropna(subset=["start", "end"]).copy()
    out["start"] = out["start"].astype("int64")
    out["end"] = out["end"].astype("int64")
    bad = out["end"] < (out["start"] + 1)
    if bad.any():
        out.loc[bad, "end"] = out.loc[bad, "start"] + 1

    group_a_lists: list[list[str]] = []
    group_b_lists: list[list[str]] = []
    supporting_lists: list[list[str]] = []
    for _, row in out.iterrows():
        variant_class = str(row.get("variant_class", "SV"))
        a_reads = _normalize_supporting_read_names(
            _parse_json_read_names(row.get("group_a_read_names", "[]")),
            variant_class=variant_class,
        )
        b_reads = _normalize_supporting_read_names(
            _parse_json_read_names(row.get("group_b_read_names", "[]")),
            variant_class=variant_class,
        )
        group_a_lists.append(a_reads)
        group_b_lists.append(b_reads)
        supporting_lists.append(list(dict.fromkeys(a_reads + b_reads)))
    out["group_a_read_names"] = group_a_lists
    out["group_b_read_names"] = group_b_lists
    out["supporting_reads"] = supporting_lists
    out["chr"] = out["chrom"].astype(str)
    out["location"] = out["start"] + 1
    out["ref_start"] = out["start"] + 1
    out["ref_end"] = out["end"]
    out["id"] = out["variant_id"].astype(str)
    out["sv_type"] = out["variant_subtype"].astype("string")
    out["sv_len"] = pd.to_numeric(out["change_size_bp"], errors="coerce").astype("Int64")
    out["vaf"] = pd.Series(pd.NA, index=out.index, dtype="Float64")
    return out


def _load_variant_table(path: str, variant_source: str | None = None) -> pd.DataFrame:
    source = variant_source or _detect_variant_source(path)
    if source == "harmonized":
        return read_harmonized_variants_to_df(path)
    return read_vcf_to_df(path)


def _merge_variant_metadata(assignment_df: pd.DataFrame, variant_df: pd.DataFrame) -> pd.DataFrame:
    extras = [col for col in [
        "id",
        "variant_class",
        "variant_subtype",
        "category",
        "change_size_bp",
        "group_a_alt_reads",
        "group_b_alt_reads",
        "group_a_read_names",
        "group_b_read_names",
    ] if col in variant_df.columns]
    if len(extras) <= 1:
        return assignment_df

    extra_df = variant_df[extras].drop_duplicates(subset=["id"]).copy()
    if "variant_subtype" in extra_df.columns and "sv_type" not in assignment_df.columns:
        extra_df["sv_type"] = extra_df["variant_subtype"]

    merged = assignment_df.merge(extra_df, on="id", how="left")
    if "variant_class" not in merged.columns:
        merged["variant_class"] = "SV"
    else:
        merged["variant_class"] = merged["variant_class"].fillna("SV").astype("string")
    if "variant_subtype" not in merged.columns:
        merged["variant_subtype"] = merged.get("sv_type", pd.Series("", index=merged.index, dtype="string"))
    else:
        merged["variant_subtype"] = merged["variant_subtype"].fillna("").astype("string")
    if "category" not in merged.columns:
        merged["category"] = pd.Series("", index=merged.index, dtype="string")
    else:
        merged["category"] = merged["category"].fillna("").astype("string")
    return merged


def _write_assignment_outputs(
    assignment_df: pd.DataFrame,
    output_dir: str,
) -> tuple[str, str, str]:
    variant_assignment_path = os.path.join(output_dir, "variant_assignment.tsv")
    variant_readable_path = os.path.join(output_dir, "variant_assignment_readable.tsv")
    variant_readable_long_path = os.path.join(output_dir, "variant_assignment_readable_long.tsv")
    sv_assignment_path = os.path.join(output_dir, "sv_assignment.tsv")
    sv_readable_path = os.path.join(output_dir, "sv_assignment_readable.tsv")
    sv_readable_long_path = os.path.join(output_dir, "sv_assignment_readable_long.tsv")

    assignment_df.to_csv(variant_assignment_path, sep="\t", index=False)
    assignment_df.to_csv(sv_assignment_path, sep="\t", index=False)

    readable_report_df, readable_report_long_df = _build_sv_readable_reports(assignment_df)
    readable_report_df.to_csv(variant_readable_path, sep="\t", index=False)
    readable_report_long_df.to_csv(variant_readable_long_path, sep="\t", index=False)
    readable_report_df.to_csv(sv_readable_path, sep="\t", index=False)
    readable_report_long_df.to_csv(sv_readable_long_path, sep="\t", index=False)
    return variant_assignment_path, variant_readable_path, variant_readable_long_path

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
        # load methylation matrix + CpG positions (reuse per-process handles to avoid
        # ~90K redundant BAM/FASTA opens per worker)
        bam_h = _get_process_bam(input_file)
        fasta_h = _get_process_fasta(reference)
        query_start, query_end = _find_dmr_query_interval(start, end)
        mm, cpgs = methyl_matrix_from_bam(
            input_file,
            reference,
            chrom=chrom,
            start=query_start,
            end=query_end,
            return_positions=True,
            bam_handle=bam_h,
            fasta_handle=fasta_h,
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

        # Impute NaN with column means — single array allocation instead of the
        # original three-copy pandas chain (astype→copy→fillna→mean).
        mm_float = mm.to_numpy(dtype=float)
        col_means = np.nanmean(mm_float, axis=0)
        nan_locs = np.isnan(mm_float)
        if nan_locs.any():
            mm_float = mm_float.copy()
            col_idx = np.where(nan_locs)[1]
            mm_float[nan_locs] = col_means[col_idx]
        read_mean = mm_float.mean(axis=1)

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

def _variant_anno_from_reads(args, *, variant_source: str | None = None):
    logger = logging.getLogger("anno.variant_anno_from_reads")
    logger.info("Starting variant annotation from pre-annotated reads")
    if args.command == "svanno":
        input_file = args.input
        output_arg = str(args.output)
        # Backward compatibility: treat extension-less output as directory.
        if os.path.isdir(output_arg) or os.path.splitext(output_arg)[1] == "":
            output_dir = output_arg
        else:
            output_dir = os.path.dirname(output_arg) or "."
    else:
        output_dir = args.output
        input_file = os.path.join(output_dir, "reads_classification.tsv")

    os.makedirs(output_dir, exist_ok=True)
    resolved_variant_source = variant_source or _detect_variant_source(args.vcf)
    if args.kanpig_read_names is not None:
        logger.info("Using kanpig read names from: %s", args.kanpig_read_names)
    else:
        logger.info("No kanpig read names provided; using read names from variant source")
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
    per_read_min_agreement = float(getattr(args, "per_read_min_agreement", 0.66))
    breakpoint_exclusion_frac = validate_breakpoint_exclusion_frac(getattr(args, "breakpoint_exclusion_frac", 0.0))
    if not (0.0 <= min_overlap_pct <= 1.0):
        raise ValueError("min_overlap_pct must be in [0, 1]")
    if not (0.0 <= min_agreement_pct <= 1.0):
        raise ValueError("min_agreement_pct must be in [0, 1]")
    if not (0.0 <= per_read_min_agreement <= 1.0):
        raise ValueError("per_read_min_agreement must be in [0, 1]")

    logger.info(
        "Variant assignment settings: variant_source=%s evidence_mode=%s unique_reads_for_overlap=%s "
        "min_overlap_pct=%.3f min_agreement_pct=%.3f per_read_min_agreement=%.3f "
        "breakpoint_exclusion_frac=%.3f",
        resolved_variant_source,
        evidence_mode,
        unique_reads_for_overlap,
        min_overlap_pct,
        min_agreement_pct,
        per_read_min_agreement,
        breakpoint_exclusion_frac,
    )

    variant_df = _load_variant_table(args.vcf, resolved_variant_source)
    variant_assignment_df = assign_sv_celltypes(
        variant_df,
        read_assign_df,
        window=int(getattr(args, "window", 5000)),
        breakpoint_exclusion_frac=breakpoint_exclusion_frac,
        min_overlap_pct=min_overlap_pct,
        min_agreement_pct=min_agreement_pct,
        unique_reads_for_overlap=unique_reads_for_overlap,
        per_read_min_agreement=per_read_min_agreement,
    )
    variant_assignment_df = _merge_variant_metadata(variant_assignment_df, variant_df)

    variant_assignment_path, readable_report_path, readable_report_long_path = _write_assignment_outputs(
        variant_assignment_df,
        output_dir,
    )
    logger.info("Wrote variant assignment report: %s", variant_assignment_path)
    logger.info("Wrote readable variant report: %s", readable_report_path)
    logger.info("Wrote long-format readable variant report: %s", readable_report_long_path)


def sv_anno(args):
    _variant_anno_from_reads(args)


def _run_annotation_pipeline(args, *, variant_source: str):
    logger = logging.getLogger("anno.main")

    bed_input = args.bed
    base_out = args.output
    input_file = args.input
    reference = args.reference
    threads = int(args.threads)
    window = int(args.window)
    breakpoint_exclusion_frac = validate_breakpoint_exclusion_frac(getattr(args, "breakpoint_exclusion_frac", 0.0))
    read_assignment_mode = str(getattr(args, "read_assignment_mode", "closest_reference_mean")).strip().lower()
    if read_assignment_mode not in {"closest_reference_mean", "kmeans"}:
        raise ValueError("read_assignment_mode must be one of: closest_reference_mean, kmeans")
    logger.info(
        "Starting annotation: variant_source=%s variants=%s bed=%s bam=%s ref=%s threads=%d "
        "read_assignment_mode=%s breakpoint_exclusion_frac=%.3f out_base=%s",
        variant_source,
        args.vcf,
        bed_input,
        input_file,
        reference,
        threads,
        read_assignment_mode,
        breakpoint_exclusion_frac,
        base_out,
    )

    os.makedirs(base_out, exist_ok=True)
    reads_out = os.path.join(base_out, "reads_classification.tsv")
    blocks_out = os.path.join(base_out, "blocks_classification.tsv")
    manifest_path = _write_anno_run_manifest(
        output_dir=base_out,
        bam=input_file,
        variant_input=args.vcf,
        variant_source=variant_source,
        reference=reference,
        bed=bed_input,
        reads_classification=reads_out,
        blocks_classification=blocks_out,
        window=window,
        breakpoint_exclusion_frac=breakpoint_exclusion_frac,
        threads=threads,
        kanpig_read_names=getattr(args, "kanpig_read_names", None),
        read_assignment_mode=read_assignment_mode,
    )
    logger.info("Wrote anno run manifest: %s", manifest_path)

    bed = pd.read_csv(bed_input, sep="\t")
    if not bed.empty and isinstance(bed.columns[0], str) and bed.columns[0].startswith("#"):
        bed.rename(columns={bed.columns[0]: bed.columns[0].lstrip("#")}, inplace=True)
    bed = bed.drop_duplicates(ignore_index=True)
    logger.info("Loaded BED with %d unique DMR rows", len(bed))

    variant_df = _load_variant_table(args.vcf, variant_source)
    filtered_bed = filter_bed_based_on_variants(
        bed,
        sv_df=variant_df,
        window=window,
        breakpoint_exclusion_frac=breakpoint_exclusion_frac,
    )

    for col in ["chr", "start", "end", "best_group", "best_dir"]:
        if col not in filtered_bed.columns:
            logger.error("BED missing required column: %s", col)
            raise ValueError(f"BED missing required column: {col}")

    n_tasks = len(filtered_bed)
    logger.info(
        "Filtered BED to %d DMRs after variant overlap filtering, window size = %d",
        n_tasks,
        window,
    )

    filtered_bed = filtered_bed.sort_values(["chr", "start"], ignore_index=True)
    tasks = [(dict(row), input_file, reference, read_assignment_mode) for _, row in filtered_bed.iterrows()]

    open(reads_out, "w").close()
    open(blocks_out, "w").close()
    reads_header_written = False
    blocks_header_written = False
    blocks_cols_locked: list[str] | None = None

    with mp.Pool(threads) as pool:
        for res in tqdm(pool.imap(_one_dmr, tasks, chunksize=50), total=n_tasks, desc="Processing DMRs"):
            if res is None:
                continue

            a_df, s_df = res

            if a_df is not None and not a_df.empty:
                if not reads_header_written:
                    a_df.to_csv(reads_out, sep="\t", index=True, mode="a", header=True)
                    reads_header_written = True
                else:
                    a_df.to_csv(reads_out, sep="\t", index=True, mode="a", header=False)

            if s_df is not None and not s_df.empty:
                if not blocks_header_written:
                    blocks_cols_locked = list(s_df.columns)
                    s_df.to_csv(blocks_out, sep="\t", index=False, mode="a", header=True)
                    blocks_header_written = True
                else:
                    assert blocks_cols_locked is not None
                    s_df_aligned = s_df.reindex(columns=blocks_cols_locked)
                    s_df_aligned.to_csv(blocks_out, sep="\t", index=False, mode="a", header=False)

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
        pd.DataFrame(columns=["chr", "start", "end", "cpgstart", "cpgend"]).to_csv(
            blocks_out, sep="\t", index=False, header=True
        )
        logger.warning("No block states generated; wrote empty blocks header only")

    _variant_anno_from_reads(args, variant_source=variant_source)
    logger.info("Annotation complete")


def anno_main(args):
    variant_source = _detect_variant_source(args.vcf)
    _run_annotation_pipeline(args, variant_source=variant_source)


def svanno_main(args):
    variant_source = _detect_variant_source(args.vcf)
    full_pipeline_mode = bool(getattr(args, "reference", None)) and bool(getattr(args, "bed", None))
    if full_pipeline_mode:
        _run_annotation_pipeline(args, variant_source=variant_source)
    else:
        _variant_anno_from_reads(args, variant_source=variant_source)
