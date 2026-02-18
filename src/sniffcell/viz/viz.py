from __future__ import annotations

import json
import logging
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from sniffcell.anno.methyl_matrix import methyl_matrix_from_bam
from sniffcell.anno.variant_assignment import (
    _build_group_leaf_sets,
    _decode_linked_celltypes_from_row,
    _resolve_hierarchy_labels,
)


def _norm_chr(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    return text[3:] if text.lower().startswith("chr") else text


def _parse_support_read_names(value: object) -> list[str]:
    null_tokens = {"", "NA", "NAN", "NONE", "."}

    if isinstance(value, (list, tuple, set, np.ndarray, pd.Series)):
        out: list[str] = []
        for item in value:
            text = str(item).strip()
            if text and text.upper() not in null_tokens:
                out.append(text)
        return out

    if isinstance(value, str):
        text = value.strip()
        if not text or text.upper() in null_tokens:
            return []
        for delim in (",", "|", ";"):
            if delim in text:
                return [t.strip() for t in text.split(delim) if t.strip() and t.strip().upper() not in null_tokens]
        return [text]

    try:
        if pd.isna(value):
            return []
    except TypeError:
        pass

    text = str(value).strip()
    if not text or text.upper() in null_tokens:
        return []
    return [text]


def _load_kanpig_supporting_reads(path: str | None, sv_id: str) -> set[str]:
    if path is None:
        return set()
    mapping = pd.read_csv(path, sep="\t", header=None, names=["sv_id", "read_name"], dtype=str)
    subset = mapping.loc[mapping["sv_id"] == sv_id, "read_name"].dropna().astype(str)
    return {x.strip() for x in subset if x.strip()}


def _first_scalar(value: object) -> object:
    if isinstance(value, (list, tuple)):
        if not value:
            return None
        return value[0]
    return value


def _get_sv_payload(vcf_path: str, sv_id: str) -> dict:
    with pysam.VariantFile(vcf_path) as vf:
        for record in vf.fetch():
            if str(record.id) != str(sv_id):
                continue

            sv_start = int(record.start)
            sv_end = int(record.stop)
            if sv_end <= sv_start:
                svlen_val = _first_scalar(record.info.get("SVLEN", 1))
                try:
                    svlen_abs = abs(int(svlen_val))
                except (TypeError, ValueError):
                    svlen_abs = 1
                sv_end = sv_start + max(1, svlen_abs)

            svtype = str(record.info.get("SVTYPE", "SV"))
            support_names = set(_parse_support_read_names(record.info.get("RNAMES", [])))
            return {
                "id": str(record.id),
                "chrom": str(record.chrom),
                "start": sv_start,
                "end": sv_end,
                "svtype": svtype,
                "svlen": _first_scalar(record.info.get("SVLEN", pd.NA)),
                "supporting_reads": support_names,
            }
    raise ValueError(f"SV ID '{sv_id}' was not found in VCF: {vcf_path}")


def _resolve_output_path(output: str, fmt: str) -> Path:
    out = Path(output)
    if out.suffix.lower() in {".png", ".pdf"}:
        return out
    return out.with_suffix(f".{fmt}")


def _load_anno_manifest(anno_output: str) -> dict:
    manifest_path = Path(anno_output) / "anno_run_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Could not find anno run manifest: {manifest_path}. "
            "Run `sniffcell anno` again (new versions write this file), or pass -i/-v/-r/-b/-a explicitly."
        )
    with manifest_path.open("r", encoding="utf-8") as fh:
        payload = json.load(fh)
    if not isinstance(payload, dict):
        raise ValueError(f"Invalid anno manifest format: {manifest_path}")
    return payload


def _resolve_viz_runtime_inputs(args, logger: logging.Logger) -> dict:
    anno_output = getattr(args, "anno_output", None)
    manifest = {}
    if anno_output:
        manifest = _load_anno_manifest(anno_output)
        logger.info("Loaded anno run manifest from: %s", Path(anno_output) / "anno_run_manifest.json")

    manifest_inputs = manifest.get("inputs", {}) if isinstance(manifest, dict) else {}
    manifest_runtime = manifest.get("runtime", {}) if isinstance(manifest, dict) else {}
    manifest_outputs = manifest.get("outputs", {}) if isinstance(manifest, dict) else {}

    bam_path = args.input or manifest_inputs.get("bam")
    vcf_path = args.vcf or manifest_inputs.get("vcf")
    reference_path = args.reference or manifest_inputs.get("reference")
    bed_path = args.bed or manifest_inputs.get("bed")

    read_assignment = args.read_assignment
    if read_assignment is None:
        read_assignment = manifest_outputs.get("reads_classification")
        if read_assignment is None and anno_output:
            guessed = Path(anno_output) / "reads_classification.tsv"
            if guessed.exists():
                read_assignment = str(guessed)

    if not bam_path or not vcf_path:
        raise ValueError(
            "viz needs BAM and VCF. Provide -i/-v, or provide --anno_output pointing to an anno folder with anno_run_manifest.json."
        )

    if args.output:
        output_path = _resolve_output_path(args.output, args.format)
    elif anno_output:
        output_path = Path(anno_output) / f"{args.sv_id}.viz.{args.format}"
    else:
        output_path = Path.cwd() / f"{args.sv_id}.viz.{args.format}"

    effective_window = int(args.window)
    if anno_output and int(args.window) == 5000 and "window" in manifest_runtime:
        try:
            effective_window = int(manifest_runtime["window"])
        except (TypeError, ValueError):
            pass

    return {
        "bam_path": str(bam_path),
        "vcf_path": str(vcf_path),
        "reference_path": (str(reference_path) if reference_path else None),
        "bed_path": (str(bed_path) if bed_path else None),
        "read_assignment_path": (str(read_assignment) if read_assignment else None),
        "output_path": output_path,
        "window": effective_window,
        "kanpig_read_names": getattr(args, "kanpig_read_names", None),
    }


def _fetch_reads(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    supporting_reads: set[str],
    max_reads: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict] = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.reference_start is None or read.reference_end is None:
                continue
            if read.is_secondary or read.is_supplementary:
                continue

            r_start = max(int(read.reference_start), start)
            r_end = min(int(read.reference_end), end)
            if r_end <= r_start:
                continue

            qname = str(read.query_name)
            rows.append(
                {
                    "read_name": qname,
                    "start": r_start,
                    "end": r_end,
                    "is_supporting": qname in supporting_reads,
                }
            )

    all_reads = pd.DataFrame(rows)
    if all_reads.empty:
        return all_reads, all_reads

    all_reads = all_reads.sort_values(
        ["is_supporting", "start", "end", "read_name"],
        ascending=[False, True, True, True],
        kind="stable",
        ignore_index=True,
    )
    if max_reads <= 0 or len(all_reads) <= max_reads:
        shown = all_reads.copy()
    else:
        supporting_df = all_reads[all_reads["is_supporting"]].copy()
        non_support_df = all_reads[~all_reads["is_supporting"]].copy()
        if len(supporting_df) >= max_reads:
            shown = supporting_df.iloc[:max_reads].copy()
        else:
            n_keep_non = max_reads - len(supporting_df)
            shown = pd.concat([supporting_df, non_support_df.iloc[:n_keep_non]], ignore_index=True)
        shown = shown.sort_values(
            ["is_supporting", "start", "end", "read_name"],
            ascending=[False, True, True, True],
            kind="stable",
            ignore_index=True,
        )
    return shown, all_reads


def _read_ctdmrs(bed_path: str | None, chrom: str, start: int, end: int) -> pd.DataFrame:
    cols = ["chr", "start", "end", "best_group", "best_group_leaves", "best_dir", "name"]
    if bed_path is None:
        return pd.DataFrame(columns=cols + ["label", "chr_norm"])

    dmrs = pd.read_csv(bed_path, sep="\t")
    if dmrs.empty:
        return pd.DataFrame(columns=cols + ["label", "chr_norm"])
    if dmrs.columns[0].startswith("#"):
        dmrs = dmrs.rename(columns={dmrs.columns[0]: dmrs.columns[0].lstrip("#")})
    required = {"chr", "start", "end"}
    if not required.issubset(dmrs.columns):
        missing = sorted(required - set(dmrs.columns))
        raise ValueError(f"ctDMR file missing required columns: {missing}")

    dmrs = dmrs.copy()
    dmrs["start"] = pd.to_numeric(dmrs["start"], errors="coerce")
    dmrs["end"] = pd.to_numeric(dmrs["end"], errors="coerce")
    dmrs = dmrs.dropna(subset=["chr", "start", "end"])
    dmrs["start"] = dmrs["start"].astype(int)
    dmrs["end"] = dmrs["end"].astype(int)
    dmrs = dmrs[dmrs["end"] > dmrs["start"]]

    chrom_norm = _norm_chr(chrom)
    dmrs["chr_norm"] = dmrs["chr"].map(_norm_chr)
    dmrs = dmrs[(dmrs["chr_norm"] == chrom_norm) & (dmrs["start"] < end) & (dmrs["end"] > start)].copy()
    if dmrs.empty:
        return pd.DataFrame(columns=cols + ["label", "chr_norm"])

    for col in cols:
        if col not in dmrs.columns:
            dmrs[col] = ""
        dmrs[col] = dmrs[col].fillna("")

    dmrs["label"] = np.where(
        dmrs["best_group_leaves"].astype(str).str.strip().ne(""),
        dmrs["best_group_leaves"].astype(str),
        np.where(
            dmrs["best_group"].astype(str).str.strip().ne(""),
            dmrs["best_group"].astype(str),
            dmrs["name"].astype(str),
        ),
    )
    dmrs = dmrs.sort_values(["start", "end"], kind="stable", ignore_index=True)
    mean_cols = [c for c in dmrs.columns if isinstance(c, str) and c.startswith("mean_")]
    keep_cols = cols + ["label", "chr_norm"] + [c for c in mean_cols if c not in cols]
    return dmrs[keep_cols]


def _summarize_ctdmr_overlap(
    dmrs: pd.DataFrame,
    reads: pd.DataFrame,
    sv_start: int,
    sv_end: int,
) -> pd.DataFrame:
    cols = [
        "chr",
        "start",
        "end",
        "label",
        "best_group",
        "best_group_leaves",
        "best_dir",
        "overlaps_sv_core",
        "sv_core_overlap_bp",
        "supporting_read_overlap_count",
        "non_supporting_read_overlap_count",
        "read_overlap_count",
    ]
    if dmrs.empty:
        return pd.DataFrame(columns=cols)
    if reads.empty:
        out = dmrs[["chr", "start", "end", "label", "best_group", "best_group_leaves", "best_dir"]].copy()
        out["overlaps_sv_core"] = (out["start"] < sv_end) & (out["end"] > sv_start)
        out["sv_core_overlap_bp"] = np.maximum(0, np.minimum(out["end"], sv_end) - np.maximum(out["start"], sv_start))
        out["supporting_read_overlap_count"] = 0
        out["non_supporting_read_overlap_count"] = 0
        out["read_overlap_count"] = 0
        return out[cols]

    read_start = reads["start"].to_numpy(np.int64)
    read_end = reads["end"].to_numpy(np.int64)
    read_support = reads["is_supporting"].to_numpy(bool)

    rows = []
    for row in dmrs.itertuples(index=False):
        ov = (read_start < int(row.end)) & (read_end > int(row.start))
        n_all = int(np.count_nonzero(ov))
        n_support = int(np.count_nonzero(ov & read_support))
        n_non_support = n_all - n_support
        sv_overlap_bp = int(max(0, min(int(row.end), sv_end) - max(int(row.start), sv_start)))
        rows.append(
            {
                "chr": row.chr,
                "start": int(row.start),
                "end": int(row.end),
                "label": str(row.label),
                "best_group": str(row.best_group),
                "best_group_leaves": str(row.best_group_leaves),
                "best_dir": str(row.best_dir),
                "overlaps_sv_core": sv_overlap_bp > 0,
                "sv_core_overlap_bp": sv_overlap_bp,
                "supporting_read_overlap_count": n_support,
                "non_supporting_read_overlap_count": n_non_support,
                "read_overlap_count": n_all,
            }
        )
    return pd.DataFrame(rows, columns=cols)


def _load_read_assignment_table(path: str | None) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()

    assignment = pd.read_csv(
        path,
        sep="\t",
        index_col=0,
        dtype={
            "code": "string",
            "code_order": "string",
            "best_group": "string",
            "best_group_leaves": "string",
            "other_group": "string",
            "other_group_leaves": "string",
        },
    )
    if assignment.empty:
        assignment["read_name"] = pd.Series(dtype="string")
        return assignment

    assignment = assignment.copy()
    assignment["read_name"] = assignment.index.astype(str)

    for col in ["chr", "start", "end", "code", "code_order", "best_group", "best_group_leaves", "other_group", "other_group_leaves"]:
        if col not in assignment.columns:
            assignment[col] = pd.NA

    if "is_best_group" not in assignment.columns:
        assignment["is_best_group"] = False
    else:
        assignment["is_best_group"] = (
            assignment["is_best_group"]
            .map(lambda x: str(x).strip().lower() in {"1", "true", "t", "yes"})
            .fillna(False)
            .astype(bool)
        )
    return assignment


def _decode_read_assignment_rows(evidence: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "read_name",
        "chr",
        "chr_norm",
        "start",
        "end",
        "assigned_celltypes",
        "assigned_schema",
        "code",
        "code_order",
    ]
    if evidence.empty:
        return pd.DataFrame(columns=cols)

    group_leaf_sets = _build_group_leaf_sets(evidence)
    rows = []
    for row in evidence.itertuples(index=False):
        row_s = pd.Series(row._asdict())
        schema, linked = _decode_linked_celltypes_from_row(row_s)
        if not linked:
            continue
        resolved = _resolve_hierarchy_labels(linked, group_leaf_sets)
        if not resolved:
            continue
        rows.append(
            {
                "read_name": str(row_s.get("read_name", "")),
                "chr": str(row_s.get("chr", "")),
                "chr_norm": _norm_chr(row_s.get("chr", "")),
                "start": int(row_s.get("start")),
                "end": int(row_s.get("end")),
                "assigned_celltypes": "|".join(resolved),
                "assigned_schema": str(schema),
                "code": str(row_s.get("code", "")),
                "code_order": str(row_s.get("code_order", "")),
            }
        )
    return pd.DataFrame(rows, columns=cols)


def _summarize_supporting_read_assignments(
    assignment_df: pd.DataFrame,
    supporting_reads: set[str],
    sv_chrom: str,
    sv_start: int,
    sv_end: int,
    window: int,
    region_start: int,
    region_end: int,
    assignment_available: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_cols = [
        "read_name",
        "is_supporting",
        "assignment_status",
        "is_assigned",
        "assigned_celltypes",
        "assigned_celltype_counts",
        "n_candidate_rows",
        "n_assignment_rows",
        "n_assignment_regions",
    ]
    detail_cols = [
        "read_name",
        "chr",
        "chr_norm",
        "start",
        "end",
        "assigned_celltypes",
        "assigned_schema",
        "code",
        "code_order",
    ]

    support_list = sorted(str(x) for x in supporting_reads if str(x).strip())
    if not support_list:
        return pd.DataFrame(columns=summary_cols), pd.DataFrame(columns=detail_cols)

    if assignment_df.empty:
        status = "unassigned_no_assignment_file" if not assignment_available else "unassigned_no_assignment_rows"
        summary = pd.DataFrame(
            {
                "read_name": support_list,
                "is_supporting": True,
                "assignment_status": status,
                "is_assigned": False,
                "assigned_celltypes": "",
                "assigned_celltype_counts": "",
                "n_candidate_rows": 0,
                "n_assignment_rows": 0,
                "n_assignment_regions": 0,
            }
        )
        return summary, pd.DataFrame(columns=detail_cols)

    evidence = assignment_df.copy()
    evidence["chr_norm"] = evidence["chr"].map(_norm_chr)
    evidence["start"] = pd.to_numeric(evidence["start"], errors="coerce")
    evidence["end"] = pd.to_numeric(evidence["end"], errors="coerce")
    evidence = evidence.dropna(subset=["chr_norm", "start", "end"])
    evidence["start"] = evidence["start"].astype(int)
    evidence["end"] = evidence["end"].astype(int)

    sv_chrom_norm = _norm_chr(sv_chrom)
    same_chr = evidence["chr_norm"].eq(sv_chrom_norm)
    overlap_padded = (int(sv_start) < (evidence["end"] + int(window))) & (int(sv_end) > (evidence["start"] - int(window)))
    overlap_core = (int(sv_start) <= evidence["end"]) & (int(sv_end) >= evidence["start"])
    in_sv_window = same_chr & overlap_padded & (~overlap_core)
    evidence = evidence[in_sv_window].copy()
    # Safety clip to the fetched plotting window.
    evidence = evidence[(evidence["start"] < int(region_end)) & (evidence["end"] > int(region_start))].copy()
    evidence = evidence[evidence["read_name"].isin(support_list)].copy()

    decoded_detail = _decode_read_assignment_rows(evidence)

    raw_counts = evidence.groupby("read_name", sort=False).size().to_dict()
    by_read_counts: dict[str, dict[str, int]] = {}
    by_read_majority: dict[str, str] = {}
    if not decoded_detail.empty:
        per_read_link_counts = (
            decoded_detail
            .groupby(["read_name", "assigned_celltypes"], sort=False)
            .size()
            .rename("count")
            .reset_index()
            .sort_values(["read_name", "count", "assigned_celltypes"], ascending=[True, False, True], kind="stable")
        )
        winners = per_read_link_counts.drop_duplicates(subset=["read_name"], keep="first")
        by_read_majority = {str(r.read_name): str(r.assigned_celltypes) for r in winners.itertuples(index=False)}

        for read_name, read_df in per_read_link_counts.groupby("read_name", sort=False):
            counter: dict[str, int] = defaultdict(int)
            for token, cnt in zip(read_df["assigned_celltypes"], read_df["count"]):
                for ct in [x.strip() for x in str(token).split("|") if x.strip()]:
                    counter[ct] += int(cnt)
            by_read_counts[str(read_name)] = dict(counter)

    rows = []
    for read_name in support_list:
        raw_n = int(raw_counts.get(read_name, 0))
        ct_counter = by_read_counts.get(read_name, {})
        if ct_counter:
            ranked = sorted(ct_counter.items(), key=lambda x: (-x[1], x[0]))
            majority_link = by_read_majority.get(read_name, "")
            majority_celltypes = [ct for ct in str(majority_link).split("|") if ct.strip()]
            celltypes = majority_celltypes if majority_celltypes else [ct for ct, _ in ranked]
            ct_count_str = ";".join(f"{ct}:{cnt}" for ct, cnt in ranked)
            decoded_n = int(len(decoded_detail[decoded_detail["read_name"] == read_name]))
            decoded_regions_n = int(
                decoded_detail[decoded_detail["read_name"] == read_name][["chr_norm", "start", "end"]]
                .drop_duplicates(ignore_index=True)
                .shape[0]
            )
            rows.append(
                {
                    "read_name": read_name,
                    "is_supporting": True,
                    "assignment_status": "assigned",
                    "is_assigned": True,
                    "assigned_celltypes": "|".join(celltypes),
                    "assigned_celltype_counts": ct_count_str,
                    "n_candidate_rows": raw_n,
                    "n_assignment_rows": decoded_n,
                    "n_assignment_regions": decoded_regions_n,
                }
            )
        else:
            status = "unassigned_unresolved_code" if raw_n > 0 else "unassigned_no_overlap_rows"
            rows.append(
                {
                    "read_name": read_name,
                    "is_supporting": True,
                    "assignment_status": status,
                    "is_assigned": False,
                    "assigned_celltypes": "",
                    "assigned_celltype_counts": "",
                    "n_candidate_rows": raw_n,
                    "n_assignment_rows": 0,
                    "n_assignment_regions": 0,
                }
            )

    summary_df = pd.DataFrame(rows, columns=summary_cols)
    summary_df = summary_df.sort_values(["is_assigned", "read_name"], ascending=[False, True], kind="stable", ignore_index=True)
    return summary_df, decoded_detail


def _compute_supporting_read_ctdmr_methylation(
    *,
    sv_id: str,
    bam_path: str,
    reference_path: str | None,
    dmrs: pd.DataFrame,
    support_assignment_df: pd.DataFrame,
    decoded_assignment_df: pd.DataFrame,
    logger: logging.Logger,
) -> pd.DataFrame:
    cols = [
        "sv_id",
        "read_name",
        "assignment_status",
        "is_assigned",
        "assigned_celltypes",
        "dmr_assigned_celltypes",
        "dmr_was_assigned",
        "chr",
        "start",
        "end",
        "label",
        "best_group",
        "best_group_leaves",
        "best_dir",
        "mean_methylation",
        "n_cpg_observed",
        "n_cpg_in_dmr",
    ]
    if support_assignment_df.empty or dmrs.empty:
        return pd.DataFrame(columns=cols)
    if reference_path is None:
        logger.warning(
            "Skipping per-read methylation on ctDMRs because --reference was not provided."
        )
        return pd.DataFrame(columns=cols)

    support_df = support_assignment_df.copy()
    support_df["read_name"] = support_df["read_name"].astype(str)
    support_lookup = support_df.set_index("read_name", drop=False)
    supporting_reads = list(support_lookup.index.unique())

    dmr_assign_map: dict[tuple[str, str, int, int], str] = {}
    if not decoded_assignment_df.empty:
        grouped = decoded_assignment_df.groupby(["read_name", "chr_norm", "start", "end"], sort=False)["assigned_celltypes"]
        for key, values in grouped:
            merged = []
            seen = set()
            for token in values:
                for ct in [x.strip() for x in str(token).split("|") if x.strip()]:
                    if ct not in seen:
                        seen.add(ct)
                        merged.append(ct)
            dmr_assign_map[(str(key[0]), str(key[1]), int(key[2]), int(key[3]))] = "|".join(merged)

    out_rows = []
    for dmr in dmrs.itertuples(index=False):
        dmr_chr = str(dmr.chr)
        dmr_start = int(dmr.start)
        dmr_end = int(dmr.end)
        dmr_chr_norm = _norm_chr(dmr_chr)
        dmr_cpg_count = 0

        read_mean = pd.Series(dtype="float64")
        read_n_obs = pd.Series(dtype="int64")
        try:
            mm, cpgs = methyl_matrix_from_bam(
                bam_path,
                reference_path,
                chrom=dmr_chr,
                start=dmr_start,
                end=dmr_end,
                return_positions=True,
            )
            dmr_cpg_count = int(len(cpgs))
            if not mm.empty:
                if isinstance(mm.index, pd.MultiIndex) and "read_name" in mm.index.names:
                    read_names = mm.index.get_level_values("read_name").astype(str)
                else:
                    read_names = mm.index.astype(str)

                mm2 = mm.copy()
                mm2["_read_name"] = read_names
                value_cols = [c for c in mm2.columns if c != "_read_name"]
                if value_cols:
                    grouped = mm2.groupby("_read_name", sort=False)[value_cols].mean()
                    read_mean = grouped.mean(axis=1, skipna=True)
                    read_n_obs = grouped.notna().sum(axis=1).astype(int)
        except Exception:
            logger.exception(
                "Failed methylation extraction for %s:%d-%d; writing NA values.",
                dmr_chr,
                dmr_start,
                dmr_end,
            )

        for read_name in supporting_reads:
            base = support_lookup.loc[read_name]
            mean_val = read_mean.get(read_name, np.nan)
            n_obs_val = int(read_n_obs.get(read_name, 0)) if read_name in read_n_obs.index else 0

            dmr_assigned = dmr_assign_map.get((read_name, dmr_chr_norm, dmr_start, dmr_end), "")
            out_rows.append(
                {
                    "sv_id": sv_id,
                    "read_name": read_name,
                    "assignment_status": str(base.get("assignment_status", "")),
                    "is_assigned": bool(base.get("is_assigned", False)),
                    "assigned_celltypes": str(base.get("assigned_celltypes", "")),
                    "dmr_assigned_celltypes": dmr_assigned,
                    "dmr_was_assigned": bool(dmr_assigned),
                    "chr": dmr_chr,
                    "start": dmr_start,
                    "end": dmr_end,
                    "label": str(dmr.label),
                    "best_group": str(dmr.best_group),
                    "best_group_leaves": str(dmr.best_group_leaves),
                    "best_dir": str(dmr.best_dir),
                    "mean_methylation": float(mean_val) if pd.notna(mean_val) else np.nan,
                    "n_cpg_observed": n_obs_val,
                    "n_cpg_in_dmr": dmr_cpg_count,
                }
            )

    return pd.DataFrame(out_rows, columns=cols)


def _build_methylation_heatmap_matrix(
    methyl_df: pd.DataFrame,
    max_reads: int = 35,
    max_dmrs: int = 20,
) -> pd.DataFrame:
    if methyl_df.empty:
        return pd.DataFrame()
    m = methyl_df.copy()
    m = m[pd.notna(m["mean_methylation"])].copy()
    if m.empty:
        return pd.DataFrame()

    m["dmr_key"] = m["label"].astype(str)
    duplicated = m["dmr_key"].duplicated(keep=False)
    m.loc[duplicated, "dmr_key"] = (
        m.loc[duplicated, "label"].astype(str)
        + "["
        + m.loc[duplicated, "start"].astype(str)
        + "-"
        + m.loc[duplicated, "end"].astype(str)
        + "]"
    )

    pivot = m.pivot_table(index="read_name", columns="dmr_key", values="mean_methylation", aggfunc="mean")
    if pivot.empty:
        return pd.DataFrame()

    read_order = (
        m[["read_name", "is_assigned"]]
        .drop_duplicates(ignore_index=True)
        .sort_values(["is_assigned", "read_name"], ascending=[False, True], kind="stable")
        ["read_name"]
        .tolist()
    )
    read_order = [r for r in read_order if r in pivot.index]
    if read_order:
        pivot = pivot.reindex(read_order)

    col_order = (
        m.groupby("dmr_key", sort=False)["mean_methylation"]
        .apply(lambda s: int(s.notna().sum()))
        .sort_values(ascending=False, kind="stable")
        .index.tolist()
    )
    col_order = [c for c in col_order if c in pivot.columns]
    if col_order:
        pivot = pivot.reindex(columns=col_order)

    if len(pivot.index) > max_reads:
        pivot = pivot.iloc[:max_reads, :]
    if len(pivot.columns) > max_dmrs:
        pivot = pivot.iloc[:, :max_dmrs]
    return pivot


def _build_dmr_methylation_stats_map(methyl_df: pd.DataFrame) -> dict[tuple[str, int, int], tuple[float, int]]:
    if methyl_df.empty:
        return {}
    stats = (
        methyl_df.groupby(["chr", "start", "end"], sort=False)
        .agg(
            mean_supporting_methylation=("mean_methylation", "mean"),
            supporting_reads_with_signal=("n_cpg_observed", lambda s: int(np.count_nonzero(np.asarray(s) > 0))),
        )
        .reset_index()
    )
    out: dict[tuple[str, int, int], tuple[float, int]] = {}
    for row in stats.itertuples(index=False):
        if pd.notna(row.mean_supporting_methylation):
            out[(str(row.chr), int(row.start), int(row.end))] = (
                float(row.mean_supporting_methylation),
                int(row.supporting_reads_with_signal),
            )
    return out


def _build_read_methylation_map(methyl_df: pd.DataFrame) -> dict[str, float]:
    if methyl_df.empty:
        return {}
    m = methyl_df[pd.notna(methyl_df["mean_methylation"])].copy()
    if m.empty:
        return {}
    stats = (
        m.groupby("read_name", sort=False)["mean_methylation"]
        .mean()
        .astype(float)
    )
    return {str(k): float(v) for k, v in stats.items() if pd.notna(v)}


def _build_reference_celltype_matrix(dmrs: pd.DataFrame) -> pd.DataFrame:
    """
    Build a ctDMR x cell-type matrix from BED mean_* columns.
    Uses all available mean_<celltype> columns (excluding summary means).
    """
    if dmrs.empty:
        return pd.DataFrame()

    skip = {"mean_best_value", "mean_rest_value", "mean_margin"}
    mean_cols = [
        c for c in dmrs.columns
        if isinstance(c, str) and c.startswith("mean_") and c not in skip
    ]
    if not mean_cols:
        return pd.DataFrame()

    mat = dmrs[mean_cols].apply(pd.to_numeric, errors="coerce")
    if mat.empty:
        return pd.DataFrame()

    rename_map = {c: c[len("mean_"):] for c in mean_cols}
    mat = mat.rename(columns=rename_map)

    # Stable row labels with coordinates to avoid collisions.
    row_labels = (
        dmrs["label"].astype(str)
        + " ["
        + dmrs["start"].astype(int).astype(str)
        + "-"
        + dmrs["end"].astype(int).astype(str)
        + "]"
    )
    mat.index = row_labels
    return mat


def _plot_sv_panel(
    *,
    sv: dict,
    shown_reads: pd.DataFrame,
    all_reads: pd.DataFrame,
    dmrs: pd.DataFrame,
    support_assignment_df: pd.DataFrame,
    methyl_df: pd.DataFrame,
    region_start: int,
    region_end: int,
    window: int,
    output_path: Path,
) -> None:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch
    except ImportError as e:
        raise ImportError(
            "sniffcell viz requires matplotlib. Install it with `pip install matplotlib`."
        ) from e

    ref_matrix = _build_reference_celltype_matrix(dmrs)
    has_ref_matrix = not ref_matrix.empty

    fig_height = max(7.5, 6.2 + 0.08 * max(1, len(shown_reads)) + 0.12 * max(1, len(dmrs)))
    if has_ref_matrix:
        fig_height += max(2.4, 0.10 * len(ref_matrix.index))

    if has_ref_matrix:
        fig = plt.figure(figsize=(18, fig_height), constrained_layout=True)
        gs = fig.add_gridspec(3, 1, height_ratios=[4.8, 2.6, max(2.2, 0.16 * len(ref_matrix.index))], hspace=0.24)
        ax_reads = fig.add_subplot(gs[0, 0])
        ax_dmrs = fig.add_subplot(gs[1, 0], sharex=ax_reads)
        ax_ref = fig.add_subplot(gs[2, 0])
    else:
        fig, axes = plt.subplots(
            2,
            1,
            figsize=(15, fig_height),
            sharex=True,
            constrained_layout=True,
            gridspec_kw={"height_ratios": [4.8, 2.6]},
        )
        ax_reads, ax_dmrs = axes
        ax_ref = None

    title_size = 18
    subtitle_size = 13
    axis_label_size = 14
    tick_label_size = 11

    sv_start = int(sv["start"])
    sv_end = int(sv["end"])
    chrom = str(sv["chrom"])

    assign_map: dict[str, bool] = {}
    assign_status_map: dict[str, str] = {}
    if not support_assignment_df.empty:
        assign_map = (
            support_assignment_df.set_index("read_name")["is_assigned"]
            .astype(bool)
            .to_dict()
        )
        assign_status_map = (
            support_assignment_df.set_index("read_name")["assignment_status"]
            .astype(str)
            .to_dict()
        )

    read_meth_map = _build_read_methylation_map(methyl_df)
    dmr_meth_stats_map = _build_dmr_methylation_stats_map(methyl_df)
    meth_threshold = 0.50
    methylated_color = "#d73027"   # IGV-style: methylated -> red
    unmethylated_color = "#4575b4" # IGV-style: unmethylated -> blue

    for ax in [ax_reads, ax_dmrs]:
        ax.axvspan(sv_start, sv_end, color="#3182bd", alpha=0.20, zorder=0)
        ax.axvline(sv_start, color="#225ea8", linewidth=1.2, alpha=0.85)
        ax.axvline(sv_end, color="#225ea8", linewidth=1.2, alpha=0.85)

    plotted_read_methylation = False
    read_track_map: dict[str, tuple[int, int, int, bool]] = {}
    if shown_reads.empty:
        ax_reads.text(0.01, 0.5, "No reads in selected window", transform=ax_reads.transAxes, va="center", fontsize=axis_label_size)
        ax_reads.set_ylim(0, 1)
    else:
        for i, row in enumerate(shown_reads.itertuples(index=False), start=1):
            read_name = str(row.read_name)
            r_start = int(row.start)
            r_end = int(row.end)
            is_supporting = bool(row.is_supporting)
            read_track_map[read_name] = (i, r_start, r_end, is_supporting)

            if bool(row.is_supporting):
                status = assign_status_map.get(read_name, "assigned")
                is_assigned = bool(assign_map.get(read_name, status == "assigned"))
                if is_assigned:
                    color = "#d62728"
                    marker = ">"
                    linestyle = "-"
                elif status.startswith("unassigned"):
                    color = "#ff7f0e"
                    marker = "X"
                    linestyle = "--"
                else:
                    color = "#ff7f0e"
                    marker = "D"
                    linestyle = ":"
                width = 2.6
                alpha = 0.95

                # Assignment pointer marker at read start.
                ax_reads.scatter(
                    [int(row.start)],
                    [i],
                    marker=marker,
                    s=130,
                    c=[color],
                    edgecolors="#111111",
                    linewidths=0.9,
                    zorder=8,
                )
            else:
                color = "#9e9e9e"
                width = 1.3
                alpha = 0.60
                linestyle = "-"
            ax_reads.hlines(i, r_start, r_end, color=color, linewidth=width, alpha=alpha, linestyles=linestyle, zorder=3)

        # Overlay ctDMR-specific methylation directly on supporting-read segments.
        if not methyl_df.empty:
            mdf = methyl_df[
                pd.notna(methyl_df["mean_methylation"])
                & (pd.to_numeric(methyl_df["n_cpg_observed"], errors="coerce").fillna(0) > 0)
            ].copy()
            for mrow in mdf.itertuples(index=False):
                read_name = str(mrow.read_name)
                if read_name not in read_track_map:
                    continue
                y, r_start, r_end, is_supporting = read_track_map[read_name]
                if not is_supporting:
                    continue
                seg_start = max(r_start, int(mrow.start))
                seg_end = min(r_end, int(mrow.end))
                if seg_end <= seg_start:
                    continue
                meth = float(mrow.mean_methylation)
                seg_color = methylated_color if meth >= meth_threshold else unmethylated_color
                # Black outline for visibility, then colored methylation segment on top.
                ax_reads.hlines(
                    y,
                    seg_start,
                    seg_end,
                    color="#111111",
                    linewidth=6.8,
                    alpha=0.95,
                    zorder=6,
                )
                ax_reads.hlines(
                    y,
                    seg_start,
                    seg_end,
                    color=seg_color,
                    linewidth=5.0,
                    alpha=0.92,
                    zorder=7,
                )
                mid_x = 0.5 * (seg_start + seg_end)
                ax_reads.scatter(
                    [mid_x],
                    [y],
                    marker="o",
                    s=26,
                    c=[seg_color],
                    edgecolors="#111111",
                    linewidths=0.35,
                    zorder=9,
                )
                plotted_read_methylation = True

        # Fallback single marker when ctDMR segment overlap is absent.
        if (not plotted_read_methylation) and read_meth_map:
            read_marker_offset = max(1, int(0.004 * max(1, region_end - region_start)))
            for read_name, meth in read_meth_map.items():
                if read_name not in read_track_map:
                    continue
                y, r_start, r_end, is_supporting = read_track_map[read_name]
                if not is_supporting:
                    continue
                marker_x = max(r_start, r_end - read_marker_offset)
                marker_color = methylated_color if float(meth) >= meth_threshold else unmethylated_color
                ax_reads.scatter(
                    [marker_x],
                    [y],
                    marker="o",
                    s=48,
                    c=[marker_color],
                    edgecolors="#000000",
                    linewidths=0.35,
                    zorder=9,
                )
                plotted_read_methylation = True
        ax_reads.set_ylim(0, len(shown_reads) + 1)

    ax_reads.set_ylabel("Reads", fontsize=axis_label_size)
    ax_reads.set_yticks([])
    ax_reads.grid(axis="x", alpha=0.2)
    ax_reads.tick_params(axis="x", labelsize=tick_label_size)

    plotted_dmr_methylation = False
    if dmrs.empty:
        ax_dmrs.text(0.01, 0.5, "No ctDMRs overlap this window", transform=ax_dmrs.transAxes, va="center", fontsize=axis_label_size)
        ax_dmrs.set_ylim(0, 1)
        ax_dmrs.set_yticks([])
    else:
        dmr_y = []
        dmr_labels = []
        for i, row in enumerate(dmrs.itertuples(index=False)):
            y = i + 0.1
            height = 0.8
            color = "#2ca25f" if str(row.best_dir).lower() == "hyper" else "#756bb1"
            ax_dmrs.broken_barh(
                [(int(row.start), int(row.end) - int(row.start))],
                (y, height),
                facecolors=color,
                alpha=0.75,
            )
            dmr_key = (str(row.chr), int(row.start), int(row.end))
            if dmr_key in dmr_meth_stats_map:
                meth_mean, n_signal = dmr_meth_stats_map[dmr_key]
                center_x = 0.5 * (int(row.start) + int(row.end))
                marker_size = 25 + 6 * min(max(n_signal, 1), 12)
                marker_color = methylated_color if float(meth_mean) >= meth_threshold else unmethylated_color
                ax_dmrs.scatter(
                    [center_x],
                    [i + 0.5],
                    marker="s",
                    s=marker_size,
                    c=[marker_color],
                    edgecolors="#111111",
                    linewidths=0.45,
                    zorder=5,
                )
                plotted_dmr_methylation = True
            dmr_y.append(i + 0.5)
            dmr_labels.append(str(row.label))
        ax_dmrs.set_ylim(0, max(1, len(dmrs)))
        if len(dmrs) <= 20:
            ax_dmrs.set_yticks(dmr_y)
            ax_dmrs.set_yticklabels(dmr_labels, fontsize=tick_label_size - 1)
        else:
            ax_dmrs.set_yticks([])

    ax_dmrs.set_ylabel("ctDMRs", fontsize=axis_label_size)
    ax_dmrs.set_xlabel(f"{chrom} coordinate (bp)", fontsize=axis_label_size)
    ax_dmrs.grid(axis="x", alpha=0.2)
    ax_dmrs.tick_params(axis="x", labelsize=tick_label_size)
    ax_dmrs.set_xlim(region_start, region_end)

    n_support_listed = int(len(set(sv["supporting_reads"])))
    n_support_in_window = int(all_reads["is_supporting"].sum()) if not all_reads.empty else 0
    n_support_shown = int(shown_reads["is_supporting"].sum()) if not shown_reads.empty else 0
    n_assigned = int(support_assignment_df["is_assigned"].sum()) if not support_assignment_df.empty else 0
    n_total_support = int(len(support_assignment_df)) if not support_assignment_df.empty else n_support_listed
    n_unassigned = max(0, n_total_support - n_assigned)

    title = f"SV {sv['id']} ({sv['svtype']}) at {chrom}:{sv_start + 1}-{sv_end} | window +/-{window} bp"
    fig.suptitle(title, fontsize=title_size, y=0.995)
    subtitle = (
        f"supporting reads in VCF: {n_support_listed} | in BAM window: {n_support_in_window} | shown: {n_support_shown} | "
        f"assigned: {n_assigned} | unassigned: {n_unassigned}"
    )
    ax_reads.text(0.01, 1.02, subtitle, transform=ax_reads.transAxes, fontsize=subtitle_size, va="bottom")

    legend_handles = [
        Line2D([0], [0], color="#d62728", lw=2.6, linestyle="-", label="Supporting read (assigned)"),
        Line2D([0], [0], color="#ff7f0e", lw=2.6, linestyle="--", label="Supporting read (unassigned)"),
        Line2D([0], [0], color="#9e9e9e", lw=1.3, label="Other read"),
        Line2D([0], [0], linestyle="None", marker=">", markerfacecolor="#d62728", markeredgecolor="#111111", markersize=9, label="Assignment pointer"),
        Line2D([0], [0], color="#6baed6", lw=5.0, label="Read methylation on ctDMR overlap"),
        Line2D([0], [0], linestyle="None", marker="o", markerfacecolor=methylated_color, markeredgecolor="#111111", markersize=7, label="Methylated (red)"),
        Line2D([0], [0], linestyle="None", marker="o", markerfacecolor=unmethylated_color, markeredgecolor="#111111", markersize=7, label="Unmethylated (blue)"),
        Patch(facecolor="#3182bd", alpha=0.20, label="SV interval"),
        Patch(facecolor="#2ca25f", alpha=0.75, label="ctDMR hyper"),
        Patch(facecolor="#756bb1", alpha=0.75, label="ctDMR hypo/other"),
    ]
    ax_reads.legend(handles=legend_handles, loc="upper right", fontsize=tick_label_size, frameon=False)

    if plotted_read_methylation or plotted_dmr_methylation:
        ax_dmrs.text(
            0.01,
            1.02,
            "Methylation color rule: mean >= 0.50 -> red, mean < 0.50 -> blue",
            transform=ax_dmrs.transAxes,
            fontsize=tick_label_size - 1,
            va="bottom",
        )

    if ax_ref is not None:
        ref_vals = ref_matrix.to_numpy(dtype=float)
        ref_ma = np.ma.masked_invalid(ref_vals)
        cmap = plt.cm.bwr.copy()
        cmap.set_bad("#d9d9d9")
        im = ax_ref.imshow(ref_ma, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0.0, vmax=1.0)
        ax_ref.set_title("Cell-type DNA methylation values on shown ctDMRs (from mean_* in BED)", fontsize=axis_label_size)
        ax_ref.set_ylabel("ctDMRs", fontsize=axis_label_size)
        ax_ref.set_xlabel("Cell types", fontsize=axis_label_size)

        # Show all cell types on x-axis.
        ax_ref.set_xticks(np.arange(ref_matrix.shape[1]))
        ax_ref.set_xticklabels(ref_matrix.columns.tolist(), rotation=60, ha="right", fontsize=max(7, tick_label_size - 2))

        # Show ctDMR labels if manageable.
        if ref_matrix.shape[0] <= 30:
            ax_ref.set_yticks(np.arange(ref_matrix.shape[0]))
            ax_ref.set_yticklabels(ref_matrix.index.tolist(), fontsize=max(7, tick_label_size - 2))
        else:
            ax_ref.set_yticks(np.arange(ref_matrix.shape[0]))
            ax_ref.set_yticklabels(
                [str(x).split(" [", 1)[0] for x in ref_matrix.index.tolist()],
                fontsize=max(6, tick_label_size - 3),
            )

        # Overlay numeric values for smaller matrices so all values are explicit.
        n_rows, n_cols = ref_matrix.shape
        if n_rows * n_cols <= 240:
            for i in range(n_rows):
                for j in range(n_cols):
                    v = ref_vals[i, j]
                    if np.isnan(v):
                        continue
                    txt_color = "white" if (v >= 0.65 or v <= 0.35) else "black"
                    ax_ref.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7, color=txt_color)

        cbar = fig.colorbar(im, ax=ax_ref, fraction=0.02, pad=0.01)
        cbar.set_label("DNA methylation (0-1)", fontsize=tick_label_size)
        cbar.ax.tick_params(labelsize=tick_label_size - 1)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def viz_main(args) -> None:
    logger = logging.getLogger("sniffcell.viz")
    if not logger.handlers:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    max_reads = int(args.max_reads)
    if max_reads <= 0:
        raise ValueError("max_reads must be > 0")

    resolved = _resolve_viz_runtime_inputs(args, logger)
    window = int(resolved["window"])
    if window < 0:
        raise ValueError("window must be >= 0")
    output_path = resolved["output_path"]
    logger.info(
        "Resolved viz inputs: bam=%s vcf=%s ref=%s bed=%s read_assignment=%s output=%s window=%d",
        resolved["bam_path"],
        resolved["vcf_path"],
        resolved["reference_path"],
        resolved["bed_path"],
        resolved["read_assignment_path"],
        output_path,
        window,
    )
    sv = _get_sv_payload(resolved["vcf_path"], args.sv_id)

    override_support = _load_kanpig_supporting_reads(resolved["kanpig_read_names"], args.sv_id)
    if override_support:
        sv["supporting_reads"] = override_support

    region_start = max(0, int(sv["start"]) - window)
    region_end = int(sv["end"]) + window

    shown_reads, all_reads = _fetch_reads(
        resolved["bam_path"],
        sv["chrom"],
        region_start,
        region_end,
        supporting_reads=set(sv["supporting_reads"]),
        max_reads=max_reads,
    )
    dmrs = _read_ctdmrs(resolved["bed_path"], sv["chrom"], region_start, region_end)
    overlap_summary = _summarize_ctdmr_overlap(dmrs, all_reads, int(sv["start"]), int(sv["end"]))

    assignment_df = _load_read_assignment_table(resolved["read_assignment_path"])
    support_assignment_df, decoded_assignment_df = _summarize_supporting_read_assignments(
        assignment_df,
        supporting_reads=set(sv["supporting_reads"]),
        sv_chrom=sv["chrom"],
        sv_start=int(sv["start"]),
        sv_end=int(sv["end"]),
        window=window,
        region_start=region_start,
        region_end=region_end,
        assignment_available=(resolved["read_assignment_path"] is not None),
    )
    methyl_df = _compute_supporting_read_ctdmr_methylation(
        sv_id=str(sv["id"]),
        bam_path=resolved["bam_path"],
        reference_path=resolved["reference_path"],
        dmrs=dmrs,
        support_assignment_df=support_assignment_df,
        decoded_assignment_df=decoded_assignment_df,
        logger=logger,
    )

    _plot_sv_panel(
        sv=sv,
        shown_reads=shown_reads,
        all_reads=all_reads,
        dmrs=dmrs,
        support_assignment_df=support_assignment_df,
        methyl_df=methyl_df,
        region_start=region_start,
        region_end=region_end,
        window=window,
        output_path=output_path,
    )

    logger.info("Wrote SV visualization: %s", output_path)
    if bool(getattr(args, "export_tables", False)):
        summary_path = output_path.with_name(f"{output_path.stem}.summary.tsv")
        assign_path = output_path.with_name(f"{output_path.stem}.supporting_reads_assignment.tsv")
        methyl_path = output_path.with_name(f"{output_path.stem}.supporting_reads_ctdmr_methylation.tsv")
        overlap_summary.to_csv(summary_path, sep="\t", index=False)
        support_assignment_df.to_csv(assign_path, sep="\t", index=False)
        methyl_df.to_csv(methyl_path, sep="\t", index=False)
        logger.info("Wrote ctDMR overlap summary: %s", summary_path)
        logger.info("Wrote supporting-read assignment table: %s", assign_path)
        logger.info("Wrote supporting-read ctDMR methylation table: %s", methyl_path)
