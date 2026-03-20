from __future__ import annotations

import argparse
import csv
import json
import logging
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from sniffcell.postprocess.postprocess import (
    _expand_path,
    _infer_sample_id,
    _sanitize_token,
    _write_json,
)


_HAP_RE = re.compile(r"_(hap\d+)_phased-set")
_REGION_RE = re.compile(r"(chr[^_]+_\d+_\d+)")

_CHANGED_ALLELE_JSON_COLS: set[str] = {
    "change_read_names",
    "change_read_lengths",
    "baseline_read_names",
    "baseline_read_lengths",
}

_CHANGED_ALLELE_COLS: list[str] = [
    "chrom", "start", "end", "LocusID",
    "change_allele", "change_celltype", "baseline_celltype",
    "tdb_change_length", "tdb_baseline_length", "change_length_bp",
    "read_mean_change_length", "read_mean_baseline_length",
    "n_change_reads", "n_baseline_reads",
    "change_read_names", "change_read_lengths",
    "baseline_read_names", "baseline_read_lengths",
    "pairing", "pairing_confidence", "max_abs_delta_bp",
]


@dataclass(frozen=True)
class TrPostArgs:
    split_dir: Path
    output_dir: Path
    sample_id: str
    group_a: str
    group_b: str
    sample_a_label: str
    sample_b_label: str
    merged_tdb: Path
    group_a_fasta: Path
    group_b_fasta: Path
    min_expansion_bp: int
    min_reads_per_hap: int
    max_fold: float
    make_plots: bool


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m sniffcell.postprocess.tr_post_processing",
        description=(
            "Summarize pairwise tandem-repeat length differences between two split groups "
            "using merged TDB output and Medaka trimmed reads."
        ),
    )
    parser.add_argument("--split-dir", required=True, help="deconv_requested_group_splits directory")
    parser.add_argument("--groups", required=True, help="Exactly two group names, comma-separated")
    parser.add_argument("--output-dir", default=None, help="Output directory for the TR post-processing report")
    parser.add_argument("--sample-id", default=None, help="Optional sample ID override")
    parser.add_argument("--sample-a-label", default=None, help="Exact TDB sample label for the first group")
    parser.add_argument("--sample-b-label", default=None, help="Exact TDB sample label for the second group")
    parser.add_argument("--merged-tdb", default=None, help="Merged TDB directory. Defaults to <split-dir>/medaka_tandem/<sample>.medaka.tdb")
    parser.add_argument("--group-a-fasta", default=None, help="trimmed_reads.fasta path for the first group")
    parser.add_argument("--group-b-fasta", default=None, help="trimmed_reads.fasta path for the second group")
    parser.add_argument("--min-expansion-bp", type=int, default=10, help="Minimum absolute paired-haplotype delta to keep a locus. Default=10")
    parser.add_argument("--min-reads-per-hap", type=int, default=4, help="Minimum reads per haplotype before trying haplotype cleanup. Default=4")
    parser.add_argument("--max-fold", type=float, default=2.0, help="Drop plotting outliers above this fold over the per-group median. Default=2.0")
    parser.add_argument("--skip-plots", action="store_true", default=False, help="Skip PNG generation")
    return parser


def _resolve_args(raw_args) -> TrPostArgs:
    split_dir = _expand_path(raw_args.split_dir)
    tokens = [x.strip() for x in str(raw_args.groups).split(",") if x.strip()]
    if len(tokens) != 2:
        raise ValueError("--groups must contain exactly two group names")
    sample_id = raw_args.sample_id or _infer_sample_id(split_dir.parent)
    output_dir = (
        _expand_path(raw_args.output_dir)
        if raw_args.output_dir
        else split_dir / "postprocess" / f"tr_post_processing_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
    group_a, group_b = tokens
    sample_a_label = raw_args.sample_a_label or f"{sample_id}.{group_a}"
    sample_b_label = raw_args.sample_b_label or f"{sample_id}.{group_b}"
    merged_tdb = (
        _expand_path(raw_args.merged_tdb)
        if raw_args.merged_tdb
        else split_dir / "medaka_tandem" / f"{_sanitize_token(sample_id)}.medaka.tdb"
    )
    group_a_fasta = (
        _expand_path(raw_args.group_a_fasta)
        if raw_args.group_a_fasta
        else split_dir / "medaka_tandem" / f"{group_a}.medaka" / "trimmed_reads.fasta"
    )
    group_b_fasta = (
        _expand_path(raw_args.group_b_fasta)
        if raw_args.group_b_fasta
        else split_dir / "medaka_tandem" / f"{group_b}.medaka" / "trimmed_reads.fasta"
    )
    return TrPostArgs(
        split_dir=split_dir,
        output_dir=output_dir,
        sample_id=sample_id,
        group_a=group_a,
        group_b=group_b,
        sample_a_label=sample_a_label,
        sample_b_label=sample_b_label,
        merged_tdb=merged_tdb,
        group_a_fasta=group_a_fasta,
        group_b_fasta=group_b_fasta,
        min_expansion_bp=int(raw_args.min_expansion_bp),
        min_reads_per_hap=int(raw_args.min_reads_per_hap),
        max_fold=float(raw_args.max_fold),
        make_plots=not bool(raw_args.skip_plots),
    )


def _write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_empty_outputs(args: TrPostArgs, *, status: str, reason: str) -> dict[str, str]:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    targets_path = args.output_dir / "targets.tsv"
    reads_path = args.output_dir / "read_lengths.tsv"
    summary_path = args.output_dir / "summary.tsv"
    clustered_path = args.output_dir / "summary_clustered.tsv"
    changed_alleles_path = args.output_dir / "changed_alleles.tsv"
    for path, fields in (
        (targets_path, ["LocusID", "chrom", "start", "end"]),
        (reads_path, ["cell_type", "region", "LocusID", "hap", "read_name", "read_length"]),
        (summary_path, ["LocusID", "cell_type", "hap", "n_reads", "mean_length", "median_length"]),
        (clustered_path, ["LocusID", "cell_type", "hap", "is_two_tight_groups", "has_cross_hap_mixing"]),
        (changed_alleles_path, _CHANGED_ALLELE_COLS),
    ):
        _write_tsv(path, [], fields)
    summary = {
        "status": status,
        "reason": reason,
        "sample_id": args.sample_id,
        "group_a": args.group_a,
        "group_b": args.group_b,
        "sample_a_label": args.sample_a_label,
        "sample_b_label": args.sample_b_label,
        "split_dir": str(args.split_dir),
        "merged_tdb": str(args.merged_tdb),
        "group_a_fasta": str(args.group_a_fasta),
        "group_b_fasta": str(args.group_b_fasta),
        "n_targets": 0,
        "n_read_rows": 0,
        "n_summary_rows": 0,
        "n_changed_allele_rows": 0,
        "targets_tsv": str(targets_path),
        "read_lengths_tsv": str(reads_path),
        "summary_tsv": str(summary_path),
        "summary_clustered_tsv": str(clustered_path),
        "changed_alleles_tsv": str(changed_alleles_path),
        "plots": [],
    }
    _write_json(args.output_dir / "summary.json", summary)
    return summary


def _import_analysis_modules() -> dict[str, Any]:
    import importlib

    modules: dict[str, Any] = {}
    for name in ("numpy", "pandas", "tdb"):
        modules[name] = importlib.import_module(name)
    try:
        modules["matplotlib.pyplot"] = importlib.import_module("matplotlib.pyplot")
        modules["seaborn"] = importlib.import_module("seaborn")
    except Exception:
        modules["matplotlib.pyplot"] = None
        modules["seaborn"] = None
    return modules


def _find_pairing_targets(data, sample_a_label: str, sample_b_label: str, min_expansion_bp: int, *, np, pd):
    sample_tables = data.get("sample", {})
    if sample_a_label not in sample_tables or sample_b_label not in sample_tables:
        return pd.DataFrame(
            columns=[
                "LocusID",
                "pairing",
                "pairing_confidence",
                "group_a_hap1_length",
                "group_a_hap2_length",
                "group_b_hap1_length",
                "group_b_hap2_length",
                "hap1_delta_bp",
                "hap2_delta_bp",
                "max_abs_delta_bp",
                "chrom",
                "start",
                "end",
                "igv_locus",
                "region",
            ]
        )

    calls = pd.concat(
        [
            sample_tables[sample_a_label][["LocusID", "allele_number", "haplotype"]].assign(sample=sample_a_label),
            sample_tables[sample_b_label][["LocusID", "allele_number", "haplotype"]].assign(sample=sample_b_label),
        ],
        ignore_index=True,
    )
    allele_info = data["allele"][["LocusID", "allele_number", "allele_length"]].copy()
    calls = calls.merge(allele_info, on=["LocusID", "allele_number"], how="left")
    calls = calls[calls["haplotype"].isin([0, 1])].copy()

    valid_haps = (
        calls.groupby(["sample", "LocusID", "haplotype"])
        .size()
        .reset_index(name="n")
        .query("n == 1")[["sample", "LocusID", "haplotype"]]
    )
    calls = calls.merge(valid_haps, on=["sample", "LocusID", "haplotype"], how="inner")

    valid_loci = (
        calls.groupby(["sample", "LocusID"])["haplotype"]
        .nunique()
        .reset_index(name="n_haps")
        .query("n_haps == 2")[["sample", "LocusID"]]
    )
    calls = calls.merge(valid_loci, on=["sample", "LocusID"], how="inner")

    shared_loci = calls.groupby("LocusID")["sample"].nunique().loc[lambda x: x == 2].index
    calls = calls[calls["LocusID"].isin(shared_loci)].copy()
    if calls.empty:
        return pd.DataFrame(
            columns=[
                "LocusID",
                "pairing",
                "pairing_confidence",
                "group_a_hap1_length",
                "group_a_hap2_length",
                "group_b_hap1_length",
                "group_b_hap2_length",
                "hap1_delta_bp",
                "hap2_delta_bp",
                "max_abs_delta_bp",
                "chrom",
                "start",
                "end",
                "igv_locus",
                "region",
            ]
        )

    wide = calls.pivot_table(
        index="LocusID",
        columns=["sample", "haplotype"],
        values="allele_length",
        aggfunc="first",
    )
    wide.columns = [f"{sample_name}_h{hap}" for sample_name, hap in wide.columns]
    wide = wide.reset_index()
    required_cols = [
        f"{sample_a_label}_h0",
        f"{sample_a_label}_h1",
        f"{sample_b_label}_h0",
        f"{sample_b_label}_h1",
    ]
    wide = wide.dropna(subset=required_cols).copy()
    if wide.empty:
        return wide

    a0 = wide[f"{sample_a_label}_h0"].astype("int32")
    a1 = wide[f"{sample_a_label}_h1"].astype("int32")
    b0 = wide[f"{sample_b_label}_h0"].astype("int32")
    b1 = wide[f"{sample_b_label}_h1"].astype("int32")

    direct_cost = (a0 - b0).abs() + (a1 - b1).abs()
    swap_cost = (a0 - b1).abs() + (a1 - b0).abs()
    use_swap = swap_cost < direct_cost
    total_cost = direct_cost + swap_cost

    wide["pairing_confidence"] = np.where(
        total_cost == 0,
        1.0,
        (total_cost - 2 * np.minimum(direct_cost, swap_cost)) / total_cost,
    )
    wide["pairing"] = np.where(use_swap, "swapped", "direct")
    wide["group_a_hap1_length"] = a0
    wide["group_a_hap2_length"] = a1
    wide["group_b_hap1_length"] = np.where(use_swap, b1, b0).astype("int32")
    wide["group_b_hap2_length"] = np.where(use_swap, b0, b1).astype("int32")
    wide["hap1_delta_bp"] = wide["group_b_hap1_length"] - wide["group_a_hap1_length"]
    wide["hap2_delta_bp"] = wide["group_b_hap2_length"] - wide["group_a_hap2_length"]
    wide["max_abs_delta_bp"] = wide[["hap1_delta_bp", "hap2_delta_bp"]].abs().max(axis=1)

    targets = wide[wide["max_abs_delta_bp"] >= min_expansion_bp].copy()
    if targets.empty:
        return targets

    targets = targets.merge(
        data["locus"][["LocusID", "chrom", "start", "end"]],
        on="LocusID",
        how="left",
    )
    targets["igv_locus"] = (
        targets["chrom"].astype(str)
        + ":"
        + targets["start"].astype(str)
        + "-"
        + targets["end"].astype(str)
    )
    targets["region"] = (
        targets["chrom"].astype(str)
        + "_"
        + targets["start"].astype(str)
        + "_"
        + targets["end"].astype(str)
    )
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    targets["chrom"] = pd.Categorical(targets["chrom"], categories=chrom_order, ordered=True)
    return targets.sort_values(["max_abs_delta_bp", "chrom", "start"], ascending=[False, True, True]).reset_index(drop=True)


def _parse_target_reads(fasta_path: Path, cell_type: str, target_set: set[str], region_to_locus: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    keep = False
    region = None
    hap = "unknown"
    read_name = ""
    seq_len = 0
    with fasta_path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if keep and region in region_to_locus:
                    rows.append(
                        {
                            "cell_type": cell_type,
                            "region": region,
                            "LocusID": region_to_locus[region],
                            "hap": hap,
                            "read_name": read_name,
                            "read_length": seq_len,
                        }
                    )
                keep = False
                region = None
                hap = "unknown"
                read_name = line[1:].split()[0] if len(line) > 1 else ""
                seq_len = 0
                region_match = _REGION_RE.search(line)
                if region_match:
                    region = region_match.group(1)
                    if region in target_set:
                        keep = True
                        hap_match = _HAP_RE.search(line)
                        if hap_match:
                            hap = hap_match.group(1)
            elif keep:
                seq_len += len(line)
    if keep and region in region_to_locus:
        rows.append(
            {
                "cell_type": cell_type,
                "region": region,
                "LocusID": region_to_locus[region],
                "hap": hap,
                "read_name": read_name,
                "read_length": seq_len,
            }
        )
    return rows


def _collect_target_reads(args: TrPostArgs, targets, *, pd):
    target_set = set(targets["region"])
    region_to_locus = targets.set_index("region")["LocusID"].to_dict()
    rows = _parse_target_reads(args.group_a_fasta, args.sample_a_label, target_set, region_to_locus)
    rows.extend(_parse_target_reads(args.group_b_fasta, args.sample_b_label, target_set, region_to_locus))
    if not rows:
        return pd.DataFrame(columns=["cell_type", "region", "LocusID", "hap", "read_name", "read_length"])
    df_reads = pd.DataFrame(rows)
    return df_reads.merge(
        targets[["LocusID", "chrom", "start", "end", "igv_locus"]],
        on="LocusID",
        how="left",
    )


def _summarize_reads(df_reads, *, np, pd):
    if df_reads.empty:
        return pd.DataFrame(
            columns=[
                "LocusID",
                "cell_type",
                "hap",
                "chrom",
                "start",
                "end",
                "n_reads",
                "mean_length",
                "median_length",
                "std_length",
                "mad_length",
                "iqr_length",
                "read_lengths",
            ]
        )
    return (
        df_reads.groupby(["LocusID", "cell_type", "hap"], as_index=False)
        .agg(
            chrom=("chrom", "first"),
            start=("start", "first"),
            end=("end", "first"),
            n_reads=("read_length", "size"),
            mean_length=("read_length", "mean"),
            median_length=("read_length", "median"),
            std_length=("read_length", "std"),
            mad_length=("read_length", lambda x: float(np.median(np.abs(x - np.median(x))))),
            iqr_length=("read_length", lambda x: float(np.percentile(x, 75) - np.percentile(x, 25))),
            read_lengths=("read_length", list),
        )
    )


def _df_reads_fix_hp(df_reads, *, min_reads_per_hap: int, np):
    def _mad(x) -> float:
        return float(np.median(np.abs(x - np.median(x))))

    def _center_scale(x, min_scale: float = 1.0) -> tuple[float, float]:
        center = float(np.median(x))
        scale = max(1.4826 * _mad(x), min_scale)
        return center, scale

    out = df_reads.copy().reset_index(drop=True)
    out["hap_fixed"] = out["hap"]
    out["hp_changed"] = False
    out["own_center"] = np.nan
    out["other_center"] = np.nan
    out["own_scale"] = np.nan
    out["score_margin"] = np.nan

    for (_, _), group in out.groupby(["LocusID", "cell_type"], sort=False):
        hap1 = group.loc[group["hap"] == "hap1", "read_length"].to_numpy(dtype=float)
        hap2 = group.loc[group["hap"] == "hap2", "read_length"].to_numpy(dtype=float)
        if len(hap1) < min_reads_per_hap or len(hap2) < min_reads_per_hap:
            continue
        c1, s1 = _center_scale(hap1)
        c2, s2 = _center_scale(hap2)
        center_sep = abs(c1 - c2)
        pooled_scale = max((s1 + s2) / 2.0, 1.0)
        if center_sep < 25.0 or center_sep / pooled_scale < 2.0:
            continue
        idx = group.index
        read_lengths = out.loc[idx, "read_length"].to_numpy(dtype=float)
        hap = out.loc[idx, "hap"].to_numpy()
        is_h1 = hap == "hap1"
        is_h2 = hap == "hap2"
        is_phased = is_h1 | is_h2
        own_center = np.where(is_h1, c1, np.where(is_h2, c2, np.nan))
        other_center = np.where(is_h1, c2, np.where(is_h2, c1, np.nan))
        own_scale = np.where(is_h1, s1, np.where(is_h2, s2, np.nan))
        dist_own = np.abs(read_lengths - own_center)
        dist_other = np.abs(read_lengths - other_center)
        margin = dist_own - dist_other
        should_flip = (
            (dist_own >= 1.5 * own_scale)
            & (margin >= 20.0)
            & (margin >= 1.0 * pooled_scale)
            & is_phased
        )
        if should_flip.any():
            new_hap = np.where(is_h1[should_flip], "hap2", "hap1")
            out.loc[idx[should_flip], "hap_fixed"] = new_hap
            out.loc[idx[should_flip], "hp_changed"] = True
        out.loc[idx, "own_center"] = own_center
        out.loc[idx, "other_center"] = other_center
        out.loc[idx, "own_scale"] = own_scale
        out.loc[idx, "score_margin"] = margin
    return out


def _detect_supported_expansion(lengths, *, np) -> dict[str, Any]:
    x = np.asarray(lengths, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 8:
        return {
            "is_two_tight_groups": False,
            "largest_gap": None,
            "left_n": None,
            "right_n": None,
            "expanded_n": None,
            "expanded_fraction": None,
            "baseline_median": None,
            "expanded_median": None,
            "expanded_span": None,
        }
    x = np.sort(x)
    gaps = np.diff(x)
    split_idx = int(np.argmax(gaps))
    largest_gap = float(gaps[split_idx])
    left = x[: split_idx + 1]
    right = x[split_idx + 1 :]
    if len(left) == 0 or len(right) == 0:
        return {
            "is_two_tight_groups": False,
            "largest_gap": largest_gap,
            "left_n": len(left),
            "right_n": len(right),
            "expanded_n": 0,
            "expanded_fraction": 0.0,
            "baseline_median": None,
            "expanded_median": None,
            "expanded_span": None,
        }
    left_med = float(np.median(left))
    right_med = float(np.median(right))
    left_mad = float(np.median(np.abs(left - left_med)))
    right_mad = float(np.median(np.abs(right - right_med)))
    pooled_mad = max((left_mad + right_mad) / 2.0, 1.0)
    if right_med > left_med:
        baseline, expanded = left, right
        baseline_med, expanded_med = left_med, right_med
        baseline_mad, expanded_mad = left_mad, right_mad
    else:
        baseline, expanded = right, left
        baseline_med, expanded_med = right_med, left_med
        baseline_mad, expanded_mad = right_mad, left_mad
    expanded_n = len(expanded)
    expanded_fraction = expanded_n / len(x)
    expanded_span = float(np.max(expanded) - np.min(expanded)) if expanded_n > 1 else 0.0
    is_supported = (
        largest_gap >= 100
        and largest_gap >= 4.0 * pooled_mad
        and len(baseline) >= 2
        and expanded_n >= 2
        and expanded_fraction >= 0.20
        and baseline_mad <= 80.0
        and expanded_mad <= 80.0
    )
    return {
        "is_two_tight_groups": bool(is_supported),
        "largest_gap": largest_gap,
        "left_n": len(left),
        "right_n": len(right),
        "expanded_n": expanded_n,
        "expanded_fraction": expanded_fraction,
        "baseline_median": baseline_med,
        "expanded_median": expanded_med,
        "expanded_span": expanded_span,
    }


def _cross_haplotype_mixing_check(group, *, np) -> dict[str, Any]:
    hap1 = group.loc[group["hap"] == "hap1", "read_length"].to_numpy(dtype=float)
    hap2 = group.loc[group["hap"] == "hap2", "read_length"].to_numpy(dtype=float)
    default = {
        "cross_hap_checkable": False,
        "hap1_median": None,
        "hap2_median": None,
        "hap_center_gap": None,
        "n_cross_assigned": None,
        "frac_cross_assigned": None,
        "has_cross_hap_mixing": False,
    }
    if len(hap1) == 0 or len(hap2) == 0:
        return default
    med1 = float(np.median(hap1))
    med2 = float(np.median(hap2))
    center_gap = abs(med1 - med2)
    if center_gap < 100:
        default["hap1_median"] = med1
        default["hap2_median"] = med2
        default["hap_center_gap"] = center_gap
        return default
    cross = 0
    total = 0
    for value in hap1:
        total += 1
        if abs(value - med2) < abs(value - med1):
            cross += 1
    for value in hap2:
        total += 1
        if abs(value - med1) < abs(value - med2):
            cross += 1
    frac = cross / total if total else None
    return {
        "cross_hap_checkable": True,
        "hap1_median": med1,
        "hap2_median": med2,
        "hap_center_gap": center_gap,
        "n_cross_assigned": cross,
        "frac_cross_assigned": frac,
        "has_cross_hap_mixing": bool(frac is not None and frac >= 0.10),
    }


def _build_clustered_summary(final_table, df_reads_clean, *, np, pd):
    if final_table.empty:
        return final_table.copy()
    cluster_rows = [
        _detect_supported_expansion(lengths, np=np)
        for lengths in final_table["read_lengths"].tolist()
    ]
    cluster_flags = pd.DataFrame(cluster_rows)
    mixing_rows: list[dict[str, Any]] = []
    phased_reads = df_reads_clean[df_reads_clean["hap"].isin(["hap1", "hap2"])]
    for (locus_id, cell_type), group in phased_reads.groupby(["LocusID", "cell_type"], sort=False):
        row = _cross_haplotype_mixing_check(group, np=np)
        row["LocusID"] = locus_id
        row["cell_type"] = cell_type
        mixing_rows.append(row)
    if mixing_rows:
        cross_hap_flags = pd.DataFrame(mixing_rows)
    else:
        cross_hap_flags = pd.DataFrame(columns=["LocusID", "cell_type"])
    return (
        final_table.join(cluster_flags)
        .merge(cross_hap_flags, on=["LocusID", "cell_type"], how="left")
        .sort_values(["std_length", "n_reads"], ascending=[False, False])
        .reset_index(drop=True)
    )


def _build_changed_allele_table(
    targets,
    df_reads,
    *,
    sample_a_label: str,
    sample_b_label: str,
    min_expansion_bp: int,
    np,
    pd,
) -> Any:
    """One row per (locus, allele) where |TDB delta| >= min_expansion_bp.

    Identifies which group carries the expanded/contracted allele and attaches
    supporting read names and lengths from that group alongside baseline reads
    for direct comparison.
    """
    if targets.empty:
        return pd.DataFrame(columns=_CHANGED_ALLELE_COLS)

    # Build lookup: (LocusID, cell_type, hap) -> (read_names, read_lengths) sorted longest-first
    read_lookup: dict[tuple[Any, str, str], tuple[list[str], list[int]]] = {}
    if not df_reads.empty:
        for (locus_id, cell_type, hap), grp in df_reads.groupby(
            ["LocusID", "cell_type", "hap"], sort=False
        ):
            grp_sorted = grp.sort_values("read_length", ascending=False)
            read_lookup[(locus_id, cell_type, hap)] = (
                grp_sorted["read_name"].tolist(),
                grp_sorted["read_length"].tolist(),
            )

    rows: list[dict[str, Any]] = []
    for _, tgt in targets.iterrows():
        locus_id = tgt["LocusID"]
        for allele_key, delta, a_len, b_len in [
            ("hap1", int(tgt["hap1_delta_bp"]), int(tgt["group_a_hap1_length"]), int(tgt["group_b_hap1_length"])),
            ("hap2", int(tgt["hap2_delta_bp"]), int(tgt["group_a_hap2_length"]), int(tgt["group_b_hap2_length"])),
        ]:
            if abs(delta) < min_expansion_bp:
                continue
            # change_celltype = the group with the longer (expanded) allele
            if delta > 0:
                change_ct, baseline_ct = sample_b_label, sample_a_label
                change_tdb, baseline_tdb = b_len, a_len
            else:
                change_ct, baseline_ct = sample_a_label, sample_b_label
                change_tdb, baseline_tdb = a_len, b_len

            change_names, change_lengths = read_lookup.get((locus_id, change_ct, allele_key), ([], []))
            baseline_names, baseline_lengths = read_lookup.get((locus_id, baseline_ct, allele_key), ([], []))

            rows.append({
                "chrom": tgt["chrom"],
                "start": int(tgt["start"]),
                "end": int(tgt["end"]),
                "LocusID": locus_id,
                "change_allele": allele_key,
                "change_celltype": change_ct,
                "baseline_celltype": baseline_ct,
                "tdb_change_length": change_tdb,
                "tdb_baseline_length": baseline_tdb,
                "change_length_bp": abs(delta),
                "read_mean_change_length": float(np.mean(change_lengths)) if change_lengths else None,
                "read_mean_baseline_length": float(np.mean(baseline_lengths)) if baseline_lengths else None,
                "n_change_reads": len(change_names),
                "n_baseline_reads": len(baseline_names),
                "change_read_names": change_names,
                "change_read_lengths": change_lengths,
                "baseline_read_names": baseline_names,
                "baseline_read_lengths": baseline_lengths,
                "pairing": tgt.get("pairing"),
                "pairing_confidence": tgt.get("pairing_confidence"),
                "max_abs_delta_bp": int(tgt["max_abs_delta_bp"]),
            })

    if not rows:
        return pd.DataFrame(columns=_CHANGED_ALLELE_COLS)

    df = pd.DataFrame(rows)
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    df["chrom"] = pd.Categorical(df["chrom"], categories=chrom_order, ordered=True)
    return df.sort_values(
        ["max_abs_delta_bp", "chrom", "start"], ascending=[False, True, True]
    ).reset_index(drop=True)


def _plot_overview(
    df_reads,
    targets,
    *,
    changed_allele_table,
    sample_a_label: str,
    sample_b_label: str,
    group_a: str,
    group_b: str,
    max_fold: float,
    outpath: Path,
    np,
    pd,
    plt,
    sns,
):
    if df_reads.empty or targets.empty:
        return
    palette = {"hap1": "#4C78A8", "hap2": "#F58518", "unknown": "#888888"}
    changed_color = "#D62728"
    cell_offset = {sample_a_label: -0.18, sample_b_label: 0.18}
    rng = np.random.default_rng(42)
    locus_order = targets["LocusID"].tolist()
    locus_to_center = {lid: idx for idx, lid in enumerate(locus_order)}
    locus_to_label = targets.set_index("LocusID")["igv_locus"].to_dict()

    # (LocusID, cell_type, hap) triples that carry the expansion/contraction
    changed_keys: set[tuple] = set()
    if changed_allele_table is not None and not changed_allele_table.empty:
        changed_keys = set(zip(
            changed_allele_table["LocusID"],
            changed_allele_table["change_celltype"],
            changed_allele_table["change_allele"],
        ))

    def _fold_clean(group):
        values = group["read_length"]
        if values.empty or max_fold == float("inf"):
            return group
        return group[values <= max_fold * values.median()]

    cleaned_groups = [_fold_clean(group) for _, group in df_reads.groupby(["LocusID", "cell_type", "hap"], sort=False)]
    plot_df = pd.concat(cleaned_groups, ignore_index=False).copy() if cleaned_groups else df_reads.iloc[0:0].copy()
    plot_df["locus_center"] = plot_df["LocusID"].map(locus_to_center)
    plot_df = plot_df.dropna(subset=["locus_center"]).copy()
    plot_df["x_plot"] = (
        plot_df["locus_center"].astype(float)
        + plot_df["cell_type"].map(cell_offset).astype(float)
        + rng.normal(0, 0.03, size=len(plot_df))
    )
    plot_df["is_changed"] = [
        (row.LocusID, row.cell_type, row.hap) in changed_keys
        for row in plot_df.itertuples()
    ]
    fig_width = max(14, len(locus_order) * 0.55)
    fig, ax = plt.subplots(figsize=(fig_width, 7))
    # Baseline reads (not expanded/contracted) — hap-coloured, lower alpha
    for hap in ["hap1", "hap2", "unknown"]:
        sub = plot_df[(plot_df["hap"] == hap) & ~plot_df["is_changed"]]
        if sub.empty:
            continue
        ax.scatter(
            sub["x_plot"],
            sub["read_length"],
            s=30,
            alpha=0.50,
            c=palette[hap],
            edgecolor="black",
            linewidth=0.15,
            label=hap,
            zorder=3,
        )
    # Changed (expanded/contracted) reads — red, on top
    sub_changed = plot_df[plot_df["is_changed"]]
    if not sub_changed.empty:
        ax.scatter(
            sub_changed["x_plot"],
            sub_changed["read_length"],
            s=38,
            alpha=0.88,
            c=changed_color,
            edgecolor="black",
            linewidth=0.3,
            label="expanded",
            zorder=4,
        )
    for (locus_id, cell_type, hap), group in plot_df.groupby(["LocusID", "cell_type", "hap"], sort=False):
        if locus_id not in locus_to_center:
            continue
        center = locus_to_center[locus_id]
        offset = cell_offset.get(cell_type, 0.0)
        mean_value = group["read_length"].mean()
        is_changed_group = (locus_id, cell_type, hap) in changed_keys
        ax.hlines(
            mean_value,
            center + offset - 0.11,
            center + offset + 0.11,
            colors=changed_color if is_changed_group else palette.get(hap, "#888888"),
            linewidth=2.5,
            zorder=5,
        )
    for idx in range(len(locus_order) - 1):
        ax.axvline(idx + 0.5, color="lightgray", linewidth=0.5, alpha=0.5, zorder=0)
    xticks: list[float] = []
    xticklabels: list[str] = []
    for locus_id in locus_order:
        center = locus_to_center[locus_id]
        xticks.extend([center + cell_offset[sample_a_label], center + cell_offset[sample_b_label]])
        xticklabels.extend([group_a, group_b])
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=90, fontsize=8)
    ax.set_xlabel(f"Per locus: {group_a} vs {group_b}")
    ax.set_ylabel("Trimmed read length (bp)")
    ymin, ymax = ax.get_ylim()
    label_y = ymax + (ymax - ymin) * 0.03
    for locus_id in locus_order:
        ax.text(locus_to_center[locus_id], label_y, locus_to_label[locus_id], rotation=90, ha="center", va="bottom", fontsize=7)
    ax.set_ylim(ymin, label_y + (ymax - ymin) * 0.20)
    sns.despine(ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    unique_handles = []
    unique_labels = []
    seen = set()
    for handle, label in zip(handles, labels):
        if label in seen:
            continue
        unique_handles.append(handle)
        unique_labels.append(label)
        seen.add(label)
    ax.legend(unique_handles, unique_labels, title="Haplotype / status", bbox_to_anchor=(1.01, 1), loc="upper left", frameon=True)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight", dpi=150)
    plt.close(fig)



def _write_dataframe_tsv(df, path: Path, *, json_cols: set[str] | None = None) -> None:
    json_cols = json_cols or set()
    out = df.copy()
    for col in json_cols:
        if col in out.columns:
            out[col] = out[col].apply(lambda value: json.dumps(value) if isinstance(value, list) else value)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t", index=False)


def tr_post_processing_main(cli_args=None) -> dict[str, Any]:
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(message)s")
    parser = _build_arg_parser()
    args = _resolve_args(parser.parse_args(cli_args))

    missing_inputs = []
    if not args.split_dir.exists():
        missing_inputs.append(str(args.split_dir))
    if not args.merged_tdb.exists():
        missing_inputs.append(str(args.merged_tdb))
    if not args.group_a_fasta.exists():
        missing_inputs.append(str(args.group_a_fasta))
    if not args.group_b_fasta.exists():
        missing_inputs.append(str(args.group_b_fasta))
    if missing_inputs:
        return _write_empty_outputs(args, status="skipped", reason="Missing required input(s): " + ", ".join(missing_inputs))
    if not args.merged_tdb.is_dir():
        return _write_empty_outputs(args, status="skipped", reason=f"Merged TDB is not a directory: {args.merged_tdb}")

    try:
        modules = _import_analysis_modules()
    except Exception as exc:
        return _write_empty_outputs(args, status="skipped", reason=f"Required analysis dependency missing: {exc}")

    np = modules["numpy"]
    pd = modules["pandas"]
    tdb = modules["tdb"]
    plt = modules["matplotlib.pyplot"]
    sns = modules["seaborn"]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    targets_path = args.output_dir / "targets.tsv"
    reads_path = args.output_dir / "read_lengths.tsv"
    summary_path = args.output_dir / "summary.tsv"
    clustered_path = args.output_dir / "summary_clustered.tsv"
    changed_alleles_path = args.output_dir / "changed_alleles.tsv"
    plot_dir = args.output_dir / "plots"
    overview_plot = plot_dir / "all_loci_overview.png"

    data = tdb.load_tdb(str(args.merged_tdb))
    targets = _find_pairing_targets(
        data,
        args.sample_a_label,
        args.sample_b_label,
        args.min_expansion_bp,
        np=np,
        pd=pd,
    )
    _write_dataframe_tsv(targets, targets_path)
    if targets.empty:
        summary = {
            "status": "completed_empty",
            "reason": "No loci passed the paired-haplotype delta filter",
            "sample_id": args.sample_id,
            "group_a": args.group_a,
            "group_b": args.group_b,
            "sample_a_label": args.sample_a_label,
            "sample_b_label": args.sample_b_label,
            "split_dir": str(args.split_dir),
            "merged_tdb": str(args.merged_tdb),
            "group_a_fasta": str(args.group_a_fasta),
            "group_b_fasta": str(args.group_b_fasta),
            "n_targets": 0,
            "n_read_rows": 0,
            "n_summary_rows": 0,
            "n_changed_allele_rows": 0,
            "targets_tsv": str(targets_path),
            "read_lengths_tsv": str(reads_path),
            "summary_tsv": str(summary_path),
            "summary_clustered_tsv": str(clustered_path),
            "changed_alleles_tsv": str(changed_alleles_path),
            "plots": [],
        }
        _write_tsv(reads_path, [], ["cell_type", "region", "LocusID", "hap", "read_name", "read_length"])
        _write_tsv(summary_path, [], ["LocusID", "cell_type", "hap", "n_reads", "mean_length", "median_length"])
        _write_tsv(clustered_path, [], ["LocusID", "cell_type", "hap", "is_two_tight_groups", "has_cross_hap_mixing"])
        _write_tsv(changed_alleles_path, [], _CHANGED_ALLELE_COLS)
        _write_json(args.output_dir / "summary.json", summary)
        return summary

    df_reads = _collect_target_reads(args, targets, pd=pd)
    if df_reads.empty:
        return _write_empty_outputs(args, status="completed_empty", reason="No trimmed reads matched the retained target loci")
    df_reads = df_reads[df_reads["read_length"] > 1].copy()
    df_reads = _df_reads_fix_hp(df_reads, min_reads_per_hap=args.min_reads_per_hap, np=np)
    df_reads["hap"] = df_reads["hap_fixed"]
    final_table = _summarize_reads(df_reads, np=np, pd=pd)
    clustered_summary = _build_clustered_summary(final_table, df_reads, np=np, pd=pd)
    changed_allele_table = _build_changed_allele_table(
        targets,
        df_reads,
        sample_a_label=args.sample_a_label,
        sample_b_label=args.sample_b_label,
        min_expansion_bp=args.min_expansion_bp,
        np=np,
        pd=pd,
    )

    _READS_FRONT = ["LocusID", "chrom", "start", "end", "cell_type", "hap", "read_name", "read_length"]
    _reads_rest = [c for c in df_reads.columns if c not in _READS_FRONT]
    _write_dataframe_tsv(df_reads[_READS_FRONT + _reads_rest], reads_path)
    _write_dataframe_tsv(final_table, summary_path, json_cols={"read_lengths"})
    _write_dataframe_tsv(clustered_summary, clustered_path, json_cols={"read_lengths"})
    _write_dataframe_tsv(changed_allele_table, changed_alleles_path, json_cols=_CHANGED_ALLELE_JSON_COLS)

    plots: list[str] = []
    if args.make_plots and plt is not None and sns is not None:
        try:
            sns.set_theme(style="whitegrid", context="talk")
            _plot_overview(
                df_reads,
                targets,
                changed_allele_table=changed_allele_table,
                sample_a_label=args.sample_a_label,
                sample_b_label=args.sample_b_label,
                group_a=args.group_a,
                group_b=args.group_b,
                max_fold=args.max_fold,
                outpath=overview_plot,
                np=np,
                pd=pd,
                plt=plt,
                sns=sns,
            )
            if overview_plot.exists():
                plots.append(str(overview_plot))
        except Exception as exc:
            logging.warning("TR post-processing plots were skipped: %s", exc)

    summary = {
        "status": "completed",
        "sample_id": args.sample_id,
        "group_a": args.group_a,
        "group_b": args.group_b,
        "sample_a_label": args.sample_a_label,
        "sample_b_label": args.sample_b_label,
        "split_dir": str(args.split_dir),
        "merged_tdb": str(args.merged_tdb),
        "group_a_fasta": str(args.group_a_fasta),
        "group_b_fasta": str(args.group_b_fasta),
        "params": {
            "min_expansion_bp": args.min_expansion_bp,
            "min_reads_per_hap": args.min_reads_per_hap,
            "max_fold": args.max_fold,
            "make_plots": args.make_plots,
        },
        "n_targets": int(len(targets)),
        "n_read_rows": int(len(df_reads)),
        "n_summary_rows": int(len(clustered_summary)),
        "n_changed_allele_rows": int(len(changed_allele_table)),
        "targets_tsv": str(targets_path),
        "read_lengths_tsv": str(reads_path),
        "summary_tsv": str(summary_path),
        "summary_clustered_tsv": str(clustered_path),
        "changed_alleles_tsv": str(changed_alleles_path),
        "plots": plots,
    }
    _write_json(args.output_dir / "summary.json", summary)
    return summary


def main() -> int:
    tr_post_processing_main()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
