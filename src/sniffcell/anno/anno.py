from __future__ import annotations

import json
import logging
import os
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pysam

from sniffcell.anno.breakpoint_exclusion import validate_breakpoint_exclusion_frac
from sniffcell.anno.methyl_matrix import methyl_matrix_from_bam
from sniffcell.anno.variant_assignment import assign_sv_celltypes


logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(processName)s] %(levelname)s: %(message)s")


def _parse_json_read_names(value: Any) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        return [str(v).strip() for v in value if str(v).strip()]
    text = "" if value is None else str(value).strip()
    if not text or text in {".", "NA", "None", "null"}:
        return []
    try:
        parsed = json.loads(text)
    except Exception:
        parsed = None
    if isinstance(parsed, list):
        return [str(v).strip() for v in parsed if str(v).strip()]
    return [tok.strip() for tok in re.split(r"[,|;\s]+", text) if tok.strip()]


def _parse_supporting_reads(value: Any) -> list[str]:
    text = "" if value is None else str(value).strip()
    if not text:
        return []
    if text.startswith("@"):
        path = Path(text[1:])
        if not path.exists():
            raise FileNotFoundError(f"Supporting-read file not found: {path}")
        text = path.read_text(encoding="utf-8").strip()
    return list(dict.fromkeys(_parse_json_read_names(text)))


def _parse_variant_location(value: str) -> tuple[str, int, int]:
    text = str(value).strip()
    match = re.fullmatch(r"([^:\s]+):([0-9,]+)(?:-([0-9,]+))?", text)
    if not match:
        raise ValueError("variant location must be formatted as chr:pos or chr:start-end")
    chrom = match.group(1)
    start_1based = int(match.group(2).replace(",", ""))
    end_1based = int(match.group(3).replace(",", "")) if match.group(3) else start_1based
    if start_1based < 1 or end_1based < start_1based:
        raise ValueError("variant location must use positive 1-based coordinates with end >= start")
    return chrom, start_1based - 1, end_1based


def _variant_record(*, variant_name: str, variant_location: str, supporting_reads: list[str]) -> dict[str, Any]:
    chrom, start0, end = _parse_variant_location(variant_location)
    return {
        "chrom": chrom,
        "start": start0,
        "end": end,
        "variant_id": str(variant_name),
        "variant_class": "VAR",
        "variant_subtype": "",
        "category": "",
        "change_size_bp": end - start0,
        "group_a_alt_reads": len(supporting_reads),
        "group_b_alt_reads": 0,
        "group_a_read_names": json.dumps(supporting_reads),
        "group_b_read_names": "[]",
    }


def _read_variant_table(path: str) -> pd.DataFrame:
    return _variant_df_from_records(pd.read_csv(path, sep="\t"))


def _variant_df_from_records(df: pd.DataFrame) -> pd.DataFrame:
    required = {"chrom", "start", "end", "variant_id"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Variant table missing required columns: {sorted(missing)}")
    out = df.copy()
    for col, default in (
        ("variant_class", "VAR"),
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
    out.loc[out["end"] < out["start"] + 1, "end"] = out["start"] + 1
    out["supporting_reads"] = [
        list(dict.fromkeys(_parse_json_read_names(row.get("group_a_read_names", "[]")) + _parse_json_read_names(row.get("group_b_read_names", "[]"))))
        for _, row in out.iterrows()
    ]
    out["chr"] = out["chrom"].astype(str)
    out["location"] = out["start"] + 1
    out["ref_start"] = out["start"] + 1
    out["ref_end"] = out["end"]
    out["id"] = out["variant_id"].astype(str)
    out["sv_type"] = out["variant_subtype"].astype("string")
    out["sv_len"] = pd.to_numeric(out["change_size_bp"], errors="coerce").astype("Int64")
    out["vaf"] = pd.Series(pd.NA, index=out.index, dtype="Float64")
    return out


def _norm_chr(value: Any) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    return text[3:] if text.lower().startswith("chr") else text


def _load_ctdmr_catalog(path_text: str) -> pd.DataFrame:
    path = Path(path_text)
    if not path.exists():
        raise FileNotFoundError(f"ctDMR catalog not found: {path}")
    bed = pd.read_csv(path, sep="\t")
    if not bed.empty and str(bed.columns[0]).startswith("#"):
        bed = bed.rename(columns={bed.columns[0]: str(bed.columns[0]).lstrip("#")})
    required = {"chr", "start", "end", "best_group", "code_order", "mean_best_value", "mean_rest_value"}
    missing = required - set(bed.columns)
    if missing:
        raise ValueError(f"ctDMR catalog missing required columns: {sorted(missing)}")
    bed = bed.drop_duplicates(ignore_index=True).copy()
    bed["start"] = pd.to_numeric(bed["start"], errors="coerce")
    bed["end"] = pd.to_numeric(bed["end"], errors="coerce")
    bed = bed.dropna(subset=["start", "end"]).copy()
    bed[["start", "end"]] = bed[["start", "end"]].astype("int64")
    return bed


def _select_support_read_ctdmrs(
    *,
    catalog: pd.DataFrame,
    variant_row: pd.Series,
    mapping_df: pd.DataFrame,
    breakpoint_exclusion_frac: float,
) -> pd.DataFrame:
    start0 = int(variant_row["ref_start"]) - 1
    end0 = int(variant_row["ref_end"])
    if end0 < start0 + 1:
        end0 = start0 + 1
    var_len = max(1, abs(end0 - start0))
    exclude_bp = int(round(var_len * breakpoint_exclusion_frac))

    mapped = mapping_df[mapping_df.get("is_mapped", False).astype(bool)].copy()
    if mapped.empty:
        return catalog.iloc[0:0].copy()
    mapped["start"] = pd.to_numeric(mapped["start"], errors="coerce")
    mapped["end"] = pd.to_numeric(mapped["end"], errors="coerce")
    mapped = mapped.dropna(subset=["chrom", "start", "end"])
    if mapped.empty:
        return catalog.iloc[0:0].copy()

    selected = []
    for chrom_norm, read_intervals in mapped.groupby(mapped["chrom"].map(_norm_chr), sort=False):
        if not chrom_norm:
            continue
        bed = catalog[catalog["chr"].map(_norm_chr) == chrom_norm].copy()
        if bed.empty:
            continue
        keep = pd.Series(False, index=bed.index)
        for _, read_row in read_intervals.iterrows():
            read_start = int(read_row["start"])
            read_end = int(read_row["end"])
            if read_end < read_start + 1:
                read_end = read_start + 1
            keep |= (bed["start"] < read_end) & (bed["end"] > read_start)
        selected.append(bed.loc[keep])

    if not selected:
        return catalog.iloc[0:0].copy()
    out = pd.concat(selected, ignore_index=True).drop_duplicates()
    variant_chrom = _norm_chr(variant_row["chr"])
    excluded_overlap = (
        (out["chr"].map(_norm_chr) == variant_chrom)
        & (out["start"] <= end0 + exclude_bp)
        & (out["end"] >= start0 - exclude_bp)
    )
    return out.loc[~excluded_overlap].sort_values(["chr", "start", "end"], ignore_index=True)


def _parse_pipe_values(value: Any) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        return [str(v).strip() for v in value if str(v).strip()]
    if value is None or pd.isna(value):
        return []
    return [tok.strip() for tok in str(value).split("|") if tok.strip()]


def _resolve_cell_types_and_targets(row: pd.Series | dict[str, Any]) -> tuple[list[str], list[str]]:
    cell_types = _parse_pipe_values(row.get("code_order", ""))
    best_group = str(row.get("best_group", "")).strip()
    if not cell_types:
        for key, value in row.items():
            if str(key).startswith("mean_") and key not in {"mean_best_value", "mean_rest_value", "mean_margin"} and pd.notna(value):
                cell_types.append(str(key)[len("mean_"):])
    if best_group and best_group not in cell_types:
        cell_types.append(best_group)
    target = _parse_pipe_values(row.get("best_group_leaves", ""))
    if not target and best_group:
        target = [best_group]
    target_set = set(target)
    ordered_target = [ct for ct in cell_types if ct in target_set]
    return list(dict.fromkeys(cell_types)), ordered_target or target


def _to_float(value: Any) -> float | None:
    try:
        if pd.isna(value):
            return None
    except TypeError:
        pass
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def _classify_reads_for_dmr(
    *,
    row: pd.Series,
    bam_path: str,
    reference_path: str,
    support_reads: set[str],
    internal_names: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    chrom = str(row["chr"])
    start = int(row["start"])
    end = int(row["end"])
    query_start = max(0, start - 1)
    mm, cpgs = methyl_matrix_from_bam(
        bam_path,
        reference_path,
        chrom=chrom,
        start=query_start,
        end=end,
        return_positions=True,
        read_name_whitelist=support_reads,
    )
    if mm is None or mm.empty or not cpgs:
        return None
    mm = mm.dropna(how="all")
    mm = mm.dropna(axis=1, how="all")
    if mm.empty:
        return None
    cpgs = [int(cpg) for cpg in mm.columns]

    mm_float = mm.to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        col_means = np.nanmean(mm_float, axis=0)
    nan_locs = np.isnan(mm_float)
    if nan_locs.any():
        mm_float = mm_float.copy()
        col_idx = np.where(nan_locs)[1]
        mm_float[nan_locs] = col_means[col_idx]
    read_mean = np.nanmean(mm_float, axis=1)

    mean_best = _to_float(row.get("mean_best_value"))
    mean_rest = _to_float(row.get("mean_rest_value"))
    if mean_best is None or mean_rest is None:
        return None
    mask_target = np.abs(read_mean - mean_best) <= np.abs(read_mean - mean_rest)

    if isinstance(mm.index, pd.MultiIndex) and "read_name" in mm.index.names:
        readnames = mm.index.get_level_values("read_name").astype(str).to_numpy()
    else:
        readnames = mm.index.astype(str).to_numpy()
    if internal_names:
        readnames = np.array([internal_names.get(read, read) for read in readnames], dtype=str)

    cell_types, target_cell_types = _resolve_cell_types_and_targets(row)
    target_set = set(target_cell_types)
    target_bits = ["1" if ct in target_set else "0" for ct in cell_types]
    if not any(bit == "1" for bit in target_bits) and str(row.get("best_group", "")) in cell_types:
        target_bits[cell_types.index(str(row.get("best_group")))] = "1"
    other_bits = ["0" if bit == "1" else "1" for bit in target_bits]
    target_code = "".join(target_bits)
    other_code = "".join(other_bits)

    assign_df = pd.DataFrame(
        {
            "chr": chrom,
            "start": start,
            "end": end,
            "cpgstart": int(cpgs[0]),
            "cpgend": int(cpgs[-1]),
            "best_group": str(row.get("best_group", "")),
            "other_group": str(row.get("other_group", "")),
            "is_best_group": mask_target,
            "code_order": "|".join(cell_types),
            "best_group_leaves": "|".join(target_cell_types),
            "other_group_leaves": str(row.get("other_group_leaves", "")),
            "hyper_group_leaves": str(row.get("hyper_group_leaves", "")),
            "hypo_group_leaves": str(row.get("hypo_group_leaves", "")),
            "code": np.where(mask_target, target_code, other_code),
        },
        index=pd.Index(readnames, name="readname"),
    )
    state_df = pd.DataFrame(
        [
            {
                "chr": chrom,
                "start": start,
                "end": end,
                "cpgstart": int(cpgs[0]),
                "cpgend": int(cpgs[-1]),
            }
        ]
    )
    return assign_df, state_df


def _support_read_mappings(
    *,
    bam_path: str,
    read_names: list[str],
    variant_name: str,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    window: int = 10000,
    internal_names: dict[str, str] | None = None,
) -> pd.DataFrame:
    wanted = set(read_names)
    found: dict[str, dict[str, Any]] = {}
    query_window = 10000 if window < 0 else int(window)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        if chrom is not None and start is not None and end is not None:
            try:
                query_iter = bam.fetch(str(chrom), max(0, int(start) - query_window), int(end) + query_window)
            except ValueError:
                query_iter = bam.fetch(until_eof=True)
        else:
            query_iter = bam.fetch(until_eof=True)
        for aln in query_iter:
            read_name = str(aln.query_name)
            if read_name not in wanted or read_name in found:
                continue
            if aln.is_unmapped:
                found[read_name] = {"chrom": "", "start": pd.NA, "end": pd.NA, "mapping_quality": pd.NA, "is_reverse": pd.NA, "is_mapped": False}
            else:
                found[read_name] = {
                    "chrom": bam.get_reference_name(aln.reference_id),
                    "start": int(aln.reference_start),
                    "end": int(aln.reference_end),
                    "mapping_quality": int(aln.mapping_quality),
                    "is_reverse": bool(aln.is_reverse),
                    "is_mapped": True,
                }
            if len(found) == len(wanted):
                break
    return pd.DataFrame(
        [
            {
                "variant_id": variant_name,
                "readname": read_name,
                "internal_readname": internal_names.get(read_name, read_name) if internal_names else read_name,
                "bam": os.path.abspath(bam_path),
                **found.get(read_name, {"chrom": "", "start": pd.NA, "end": pd.NA, "mapping_quality": pd.NA, "is_reverse": pd.NA, "is_mapped": False}),
            }
            for read_name in read_names
        ]
    )


def _empty_reads_df() -> pd.DataFrame:
    out = pd.DataFrame(
        columns=[
            "chr", "start", "end", "cpgstart", "cpgend", "best_group", "other_group",
            "is_best_group", "code_order", "best_group_leaves", "other_group_leaves",
            "hyper_group_leaves", "hypo_group_leaves", "code",
        ]
    )
    out.index.name = "readname"
    return out


def _compute_variant_read_classification(
    *,
    variant_row: pd.Series,
    catalog: pd.DataFrame,
    bam_path: str,
    reference_path: str,
    support_reads: list[str],
    window: int,
    breakpoint_exclusion_frac: float,
    internal_names: dict[str, str] | None = None,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mapping_df = _support_read_mappings(
        bam_path=bam_path,
        read_names=support_reads,
        variant_name=str(variant_row["id"]),
        chrom=str(variant_row["chr"]),
        start=int(variant_row["ref_start"]) - 1,
        end=int(variant_row["ref_end"]),
        window=window,
        internal_names=internal_names,
    )
    nearby = _select_support_read_ctdmrs(
        catalog=catalog,
        variant_row=variant_row,
        mapping_df=mapping_df,
        breakpoint_exclusion_frac=breakpoint_exclusion_frac,
    )
    logger.info(
        "%s: selected %d ctDMRs overlapping supporting-read alignments",
        variant_row["id"],
        len(nearby),
    )
    read_tables = []
    block_tables = []
    support_set = set(support_reads)
    for _, dmr in nearby.iterrows():
        result = _classify_reads_for_dmr(
            row=dmr,
            bam_path=bam_path,
            reference_path=reference_path,
            support_reads=support_set,
            internal_names=internal_names,
        )
        if result is None:
            continue
        read_df, block_df = result
        if not read_df.empty:
            read_tables.append(read_df)
        if not block_df.empty:
            block_tables.append(block_df)
    reads_df = pd.concat(read_tables) if read_tables else _empty_reads_df()
    blocks_df = pd.concat(block_tables, ignore_index=True) if block_tables else pd.DataFrame(columns=["chr", "start", "end", "cpgstart", "cpgend"])
    return reads_df, blocks_df, mapping_df


def _readable_reports(assignment_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_cols = ["id", "sv_chr", "sv_pos", "sv_len", "n_supporting", "n_overlapped", "overlap_pct", "majority_pct", "classified_celltypes", "classified_celltype_count", "classified_celltype_counts", "classified_celltype_fractions"]
    long_cols = ["id", "sv_chr", "sv_pos", "celltype", "rank", "supporting_read_count", "supporting_read_fraction", "n_supporting", "n_overlapped", "overlap_pct"]
    if assignment_df.empty:
        return pd.DataFrame(columns=summary_cols), pd.DataFrame(columns=long_cols)
    summary_rows = []
    long_rows = []
    for _, row in assignment_df.iterrows():
        linked = [ct for ct in str(row.get("linked_celltypes", "")).split("|") if ct]
        count_map = _parse_metric_map(row.get("linked_celltype_counts", ""), int)
        frac_map = _parse_metric_map(row.get("linked_celltype_fractions", ""), float)
        for rank, ct in enumerate(linked, start=1):
            long_rows.append({"id": row.get("id", ""), "sv_chr": row.get("sv_chr", ""), "sv_pos": row.get("sv_pos", pd.NA), "celltype": ct, "rank": rank, "supporting_read_count": count_map.get(ct, pd.NA), "supporting_read_fraction": frac_map.get(ct, pd.NA), "n_supporting": row.get("n_supporting", pd.NA), "n_overlapped": row.get("n_overlapped", pd.NA), "overlap_pct": row.get("overlap_pct", pd.NA)})
        summary_rows.append({"id": row.get("id", ""), "sv_chr": row.get("sv_chr", ""), "sv_pos": row.get("sv_pos", pd.NA), "sv_len": row.get("sv_len", pd.NA), "n_supporting": row.get("n_supporting", pd.NA), "n_overlapped": row.get("n_overlapped", pd.NA), "overlap_pct": row.get("overlap_pct", pd.NA), "majority_pct": row.get("majority_pct", pd.NA), "classified_celltypes": "|".join(linked), "classified_celltype_count": len(linked), "classified_celltype_counts": row.get("linked_celltype_counts", ""), "classified_celltype_fractions": row.get("linked_celltype_fractions", "")})
    return pd.DataFrame(summary_rows, columns=summary_cols), pd.DataFrame(long_rows, columns=long_cols)


def _parse_metric_map(value: Any, cast_type) -> dict[str, Any]:
    text = "" if value is None or pd.isna(value) else str(value).strip()
    out = {}
    for token in text.split(";"):
        if ":" not in token:
            continue
        key, raw = token.rsplit(":", 1)
        try:
            out[key.strip()] = cast_type(raw.strip())
        except (TypeError, ValueError):
            continue
    return out


def _write_assignment_outputs(assignment_df: pd.DataFrame, output_dir: str) -> None:
    assignment_df.to_csv(os.path.join(output_dir, "variant_assignment.tsv"), sep="\t", index=False)
    readable, readable_long = _readable_reports(assignment_df)
    readable.to_csv(os.path.join(output_dir, "variant_assignment_readable.tsv"), sep="\t", index=False)
    readable_long.to_csv(os.path.join(output_dir, "variant_assignment_readable_long.tsv"), sep="\t", index=False)


def _assign_from_tables(args, *, variants: pd.DataFrame, reads: pd.DataFrame, output_dir: str) -> None:
    evidence_mode = str(getattr(args, "evidence_mode", "per_read")).strip().lower()
    if evidence_mode not in {"all_rows", "per_read"}:
        raise ValueError("evidence_mode must be one of: all_rows, per_read")
    assignment = assign_sv_celltypes(
        variants,
        reads,
        window=int(getattr(args, "window", 10000)),
        breakpoint_exclusion_frac=validate_breakpoint_exclusion_frac(getattr(args, "breakpoint_exclusion_frac", 0.0)),
        min_overlap_pct=float(getattr(args, "min_overlap_pct", 0.0)),
        min_agreement_pct=float(getattr(args, "min_agreement_pct", 0.0)),
        unique_reads_for_overlap=(evidence_mode == "per_read"),
        per_read_min_agreement=float(getattr(args, "per_read_min_agreement", 0.66)),
        use_read_names=True,
    )
    extras = [col for col in ["id", "variant_class", "variant_subtype", "category", "change_size_bp", "group_a_alt_reads", "group_b_alt_reads", "group_a_read_names", "group_b_read_names"] if col in variants.columns]
    if len(extras) > 1:
        assignment = assignment.merge(variants[extras].drop_duplicates(subset=["id"]), on="id", how="left")
    _write_assignment_outputs(assignment, output_dir)


def _validate_single_args(args) -> None:
    missing = [name for name in ("variant_name", "variant_location", "supporting_reads", "catalog", "input", "reference") if not getattr(args, name, None)]
    if missing:
        rendered = ", ".join("--" + name.replace("_", "-") for name in missing)
        raise ValueError(f"single-variant anno is missing required inputs: {rendered}")


def _run_compact_annotation(args) -> None:
    logger = logging.getLogger("anno.single")
    _validate_single_args(args)
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    support_reads = _parse_supporting_reads(args.supporting_reads)
    if not support_reads:
        raise ValueError("--supporting-reads did not contain any read names")

    variant_path = os.path.join(output_dir, "compact_variant.tsv")
    pd.DataFrame([_variant_record(variant_name=args.variant_name, variant_location=args.variant_location, supporting_reads=support_reads)]).to_csv(variant_path, sep="\t", index=False)
    variants = _read_variant_table(variant_path)
    catalog = _load_ctdmr_catalog(args.catalog)
    reads_df, blocks_df, mapping_df = _compute_variant_read_classification(
        variant_row=variants.iloc[0],
        catalog=catalog,
        bam_path=args.input,
        reference_path=args.reference,
        support_reads=support_reads,
        window=int(args.window),
        breakpoint_exclusion_frac=validate_breakpoint_exclusion_frac(getattr(args, "breakpoint_exclusion_frac", 0.0)),
        logger=logger,
    )
    reads_df.to_csv(os.path.join(output_dir, "reads_classification.tsv"), sep="\t", index=True)
    blocks_df.to_csv(os.path.join(output_dir, "blocks_classification.tsv"), sep="\t", index=False)
    mapping_df.to_csv(os.path.join(output_dir, "support_read_mappings.tsv"), sep="\t", index=False)
    _assign_from_tables(args, variants=variants, reads=reads_df, output_dir=output_dir)
    _write_manifest(output_dir=output_dir, mode="single", inputs={"variant_name": str(args.variant_name), "variant_location": str(args.variant_location), "supporting_read_count": len(support_reads), "bam": os.path.abspath(args.input), "reference": os.path.abspath(args.reference), "catalog": os.path.abspath(args.catalog), "window": int(args.window)})


def _read_batch(path_text: str) -> pd.DataFrame:
    path = Path(path_text)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    df = pd.read_csv(path, sep=sep)
    required = {"variant_name", "variant_location", "supporting_reads", "catalog", "bam"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Batch input missing required columns: {sorted(missing)}")
    return df


def _run_batch_annotation(args) -> None:
    logger = logging.getLogger("anno.batch")
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    batch = _read_batch(args.batch)
    variant_records = []
    read_tables = []
    block_tables = []
    mapping_tables = []
    catalog_cache: dict[str, pd.DataFrame] = {}
    for row_index, row in batch.iterrows():
        variant_name = str(row["variant_name"]).strip()
        support_reads = _parse_supporting_reads(row["supporting_reads"])
        if not variant_name:
            raise ValueError(f"batch row {row_index} has an empty variant_name")
        if not support_reads:
            raise ValueError(f"batch row {row_index} ({variant_name}) has no supporting reads")
        reference = str(row["reference"]).strip() if "reference" in batch.columns and pd.notna(row.get("reference")) else str(getattr(args, "reference", "") or "")
        if not reference:
            raise ValueError("batch anno requires either a reference column or global --reference")
        internal_names = {read: f"{variant_name}::{read}" for read in support_reads}
        internal_reads = [internal_names[read] for read in support_reads]
        variant_records.append(_variant_record(variant_name=variant_name, variant_location=str(row["variant_location"]), supporting_reads=internal_reads))
        variants = _variant_df_from_records(pd.DataFrame([variant_records[-1]]))
        catalog_path = str(row["catalog"])
        if catalog_path not in catalog_cache:
            catalog_cache[catalog_path] = _load_ctdmr_catalog(catalog_path)
        reads_df, blocks_df, mapping_df = _compute_variant_read_classification(
            variant_row=variants.iloc[0],
            catalog=catalog_cache[catalog_path],
            bam_path=str(row["bam"]),
            reference_path=reference,
            support_reads=support_reads,
            window=int(args.window),
            breakpoint_exclusion_frac=validate_breakpoint_exclusion_frac(getattr(args, "breakpoint_exclusion_frac", 0.0)),
            internal_names=internal_names,
            logger=logger,
        )
        read_tables.append(reads_df)
        block_tables.append(blocks_df)
        mapping_tables.append(mapping_df)

    variant_path = os.path.join(output_dir, "batch_variants.tsv")
    pd.DataFrame(variant_records).to_csv(variant_path, sep="\t", index=False)
    variants = _read_variant_table(variant_path)
    reads_all = pd.concat(read_tables) if read_tables else _empty_reads_df()
    blocks_all = pd.concat(block_tables, ignore_index=True) if block_tables else pd.DataFrame(columns=["chr", "start", "end", "cpgstart", "cpgend"])
    mappings_all = pd.concat(mapping_tables, ignore_index=True) if mapping_tables else pd.DataFrame()
    reads_all.to_csv(os.path.join(output_dir, "reads_classification.tsv"), sep="\t", index=True)
    blocks_all.to_csv(os.path.join(output_dir, "blocks_classification.tsv"), sep="\t", index=False)
    mappings_all.to_csv(os.path.join(output_dir, "support_read_mappings.tsv"), sep="\t", index=False)
    _assign_from_tables(args, variants=variants, reads=reads_all, output_dir=output_dir)
    _write_manifest(output_dir=output_dir, mode="batch", inputs={"batch": os.path.abspath(args.batch), "row_count": int(len(batch)), "reference": os.path.abspath(args.reference) if getattr(args, "reference", None) else ""})


def _write_manifest(*, output_dir: str, mode: str, inputs: dict[str, Any]) -> None:
    payload = {
        "command": "anno",
        "mode": mode,
        "inputs": inputs,
        "outputs": {
            "variant_assignment": os.path.abspath(os.path.join(output_dir, "variant_assignment.tsv")),
            "variant_assignment_readable": os.path.abspath(os.path.join(output_dir, "variant_assignment_readable.tsv")),
            "variant_assignment_readable_long": os.path.abspath(os.path.join(output_dir, "variant_assignment_readable_long.tsv")),
            "reads_classification": os.path.abspath(os.path.join(output_dir, "reads_classification.tsv")),
            "blocks_classification": os.path.abspath(os.path.join(output_dir, "blocks_classification.tsv")),
            "support_read_mappings": os.path.abspath(os.path.join(output_dir, "support_read_mappings.tsv")),
        },
    }
    name = "anno_batch_manifest.json" if mode == "batch" else "anno_compact_manifest.json"
    with open(os.path.join(output_dir, name), "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def _empty_reads_df() -> pd.DataFrame:
    out = pd.DataFrame(columns=["chr", "start", "end", "cpgstart", "cpgend", "best_group", "other_group", "is_best_group", "code_order", "best_group_leaves", "other_group_leaves", "hyper_group_leaves", "hypo_group_leaves", "code"])
    out.index.name = "readname"
    return out


def anno_main(args) -> None:
    if getattr(args, "batch", None):
        _run_batch_annotation(args)
        return
    _run_compact_annotation(args)
