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

from sniffcell.discover.discover import (
    _expand_path,
    _infer_sample_id,
    _sanitize_token,
    _write_json,
)


_HAP_RE = re.compile(r"_(hap\d+)_phased-set")
_REGION_RE = re.compile(r"(chr[^_]+_\d+_\d+)")

_TARGET_COLS: list[str] = [
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

_CHANGED_ALLELE_COLS: list[str] = [
    "chrom", "start", "end", "LocusID",
    "change_allele", "change_celltype", "baseline_celltype",
    "tdb_change_length", "tdb_baseline_length", "change_length_bp",
    "read_mean_change_length", "read_mean_baseline_length",
    "n_change_reads", "n_baseline_reads",
    "n_change_support_reads", "n_baseline_support_reads",
    "change_read_names", "change_read_lengths",
    "baseline_read_names", "baseline_read_lengths",
    "change_support_read_names", "change_support_read_lengths",
    "baseline_support_read_names", "baseline_support_read_lengths",
    "pairing", "pairing_confidence", "max_abs_delta_bp",
    "signal_class",
    "tail_read_count", "tail_far_read_count",
    "tail_baseline_same_hap_max_bp", "tail_change_other_hap_upper_bp",
    "tail_anchor_bp", "tail_max_excess_bp",
    "sample_change_lower_bp", "sample_change_upper_bp",
    "sample_baseline_lower_bp", "sample_baseline_upper_bp",
    "sample_lower_delta_bp", "sample_upper_delta_bp",
    "sample_range_supports_tail",
]

_TR_BED_COLS: list[str] = [
    "chrom", "start", "end", "trid",
    "change_allele", "change_type",
    "change_group", "baseline_group",
    "n_change_reads", "n_baseline_reads",
    "n_change_support_reads", "n_baseline_support_reads",
    "tdb_change_length", "tdb_baseline_length", "change_length_bp",
    "change_read_mean", "baseline_read_mean",
    "change_read_range", "baseline_read_range",
    "change_read_names", "baseline_read_names",
    "change_support_read_names", "baseline_support_read_names",
    "pairing", "pairing_confidence",
    "signal_class",
    "tail_read_count", "tail_far_read_count",
    "tail_baseline_same_hap_max_bp", "tail_change_other_hap_upper_bp",
    "tail_anchor_bp", "tail_max_excess_bp",
    "sample_change_lower_bp", "sample_change_upper_bp",
    "sample_baseline_lower_bp", "sample_baseline_upper_bp",
    "sample_lower_delta_bp", "sample_upper_delta_bp",
    "sample_range_supports_tail",
    "tr_tier", "tr_pass_for_harmonized",
    "change_median_bp", "baseline_median_bp",
    "median_shift_bp", "median_shift_ratio",
    "other_hap_median_delta_bp",
    "n_hp_changed", "hp_changed_fraction",
    "change_cross_hap_mixing", "baseline_cross_hap_mixing",
]

_TR_TIER_ORDER: dict[str, int] = {
    "strong": 0,
    "supportive": 1,
    "weak": 2,
}

_TAIL_SIGNAL_CLASS = "tail_expansion"
_TDB_SIGNAL_CLASS = "tdb_delta"
_TAIL_MIN_READS_PER_GROUP = 10
_TAIL_MARGIN_BP = 15
_TAIL_FAR_MARGIN_BP = 40
_TAIL_MAX_OTHER_HAP_MEDIAN_DELTA_BP = 10
_TAIL_MIN_SAMPLE_UPPER_DELTA_BP = 5


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
    tail_expansion_rescue: bool
    tail_require_sample_range_support: bool


@dataclass(frozen=True)
class TrPostOutputs:
    targets_tsv: Path
    read_lengths_tsv: Path
    summary_tsv: Path
    summary_clustered_tsv: Path
    tr_bed_tsv: Path
    plots_dir: Path
    overview_plot: Path
    summary_json: Path


def _build_arg_parser(
    *,
    prog: str = "python -m sniffcell.discover.tr_post_processing",
    add_help: bool = True,
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Summarize pairwise tandem-repeat length differences between two split groups "
            "using merged TDB output and Medaka trimmed reads."
        ),
        add_help=add_help,
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
    parser.add_argument(
        "--tail-expansion-rescue",
        action="store_true",
        default=False,
        help=(
            "Experimental: rescue loci with little/no TDB hap delta when one group shows a clean, "
            "hap-specific expansion tail in the trimmed reads."
        ),
    )
    parser.add_argument(
        "--tail-expansion-require-sample-range-support",
        dest="tail_require_sample_range_support",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Require tail-expansion rescue rows to also show same-haplotype TDB sample-range support. "
            "Use --no-tail-expansion-require-sample-range-support to allow trimmed-read-only tail rescue."
        ),
    )
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
        tail_expansion_rescue=bool(raw_args.tail_expansion_rescue),
        tail_require_sample_range_support=bool(raw_args.tail_require_sample_range_support),
    )


def _write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _build_output_paths(output_dir: Path) -> TrPostOutputs:
    plots_dir = output_dir / "plots"
    return TrPostOutputs(
        targets_tsv=output_dir / "targets.tsv",
        read_lengths_tsv=output_dir / "read_lengths.tsv",
        summary_tsv=output_dir / "summary.tsv",
        summary_clustered_tsv=output_dir / "summary_clustered.tsv",
        tr_bed_tsv=output_dir / "tr_changes.bed.tsv",
        plots_dir=plots_dir,
        overview_plot=plots_dir / "all_loci_overview.png",
        summary_json=output_dir / "summary.json",
    )


def _build_summary_payload(
    args: TrPostArgs,
    outputs: TrPostOutputs,
    *,
    status: str,
    reason: str | None = None,
    n_targets: int = 0,
    n_read_rows: int = 0,
    n_summary_rows: int = 0,
    n_tr_bed_rows: int = 0,
    n_tr_strong_rows: int = 0,
    n_tr_supportive_rows: int = 0,
    n_tr_weak_rows: int = 0,
    n_tail_rescue_rows: int = 0,
    plots: list[str] | None = None,
) -> dict[str, Any]:
    summary = {
        "status": status,
        "sample_id": args.sample_id,
        "group_a": args.group_a,
        "group_b": args.group_b,
        "sample_a_label": args.sample_a_label,
        "sample_b_label": args.sample_b_label,
        "split_dir": str(args.split_dir),
        "merged_tdb": str(args.merged_tdb),
        "group_a_fasta": str(args.group_a_fasta),
        "group_b_fasta": str(args.group_b_fasta),
        "n_targets": int(n_targets),
        "n_read_rows": int(n_read_rows),
        "n_summary_rows": int(n_summary_rows),
        "n_tr_bed_rows": int(n_tr_bed_rows),
        "n_tr_strong_rows": int(n_tr_strong_rows),
        "n_tr_supportive_rows": int(n_tr_supportive_rows),
        "n_tr_weak_rows": int(n_tr_weak_rows),
        "n_tail_rescue_rows": int(n_tail_rescue_rows),
        "targets_tsv": str(outputs.targets_tsv),
        "read_lengths_tsv": str(outputs.read_lengths_tsv),
        "summary_tsv": str(outputs.summary_tsv),
        "summary_clustered_tsv": str(outputs.summary_clustered_tsv),
        "tr_bed_tsv": str(outputs.tr_bed_tsv),
        "plots": plots or [],
    }
    if reason is not None:
        summary["reason"] = reason
    return summary


def _write_empty_outputs(args: TrPostArgs, *, status: str, reason: str) -> dict[str, str]:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = _build_output_paths(args.output_dir)
    for path, fields in (
        (outputs.targets_tsv, ["LocusID", "chrom", "start", "end"]),
        (outputs.read_lengths_tsv, ["cell_type", "region", "LocusID", "hap", "read_name", "read_length"]),
        (outputs.summary_tsv, ["LocusID", "cell_type", "hap", "n_reads", "mean_length", "median_length"]),
        (outputs.summary_clustered_tsv, ["LocusID", "cell_type", "hap", "is_two_tight_groups", "has_cross_hap_mixing"]),
        (outputs.tr_bed_tsv, _TR_BED_COLS),
    ):
        _write_tsv(path, [], fields)
    summary = _build_summary_payload(args, outputs, status=status, reason=reason)
    _write_json(outputs.summary_json, summary)
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


def _empty_targets_frame(*, pd):
    return pd.DataFrame(columns=_TARGET_COLS)


def _collect_valid_haplotype_calls(data, sample_a_label: str, sample_b_label: str, *, pd):
    sample_tables = data.get("sample", {})
    if sample_a_label not in sample_tables or sample_b_label not in sample_tables:
        return pd.DataFrame(columns=["LocusID", "allele_number", "haplotype", "sample", "allele_length"])

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
    return calls[calls["LocusID"].isin(shared_loci)].copy()


def _pair_haplotype_lengths(calls, sample_a_label: str, sample_b_label: str, *, np, pd):
    if calls.empty:
        return _empty_targets_frame(pd=pd)

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
        return _empty_targets_frame(pd=pd)

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
    return wide


def _annotate_pairing_candidates(wide, data, *, pd):
    if wide.empty:
        return _empty_targets_frame(pd=pd)

    targets = wide.merge(
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
    return targets[_TARGET_COLS].sort_values(
        ["max_abs_delta_bp", "chrom", "start"],
        ascending=[False, True, True],
    ).reset_index(drop=True)


def _filter_pairing_targets(candidates, min_expansion_bp: int, *, pd):
    if candidates.empty:
        return _empty_targets_frame(pd=pd)
    targets = candidates[candidates["max_abs_delta_bp"] >= min_expansion_bp].copy()
    if targets.empty:
        return _empty_targets_frame(pd=pd)
    return targets.reset_index(drop=True)


def _find_pairing_targets(data, sample_a_label: str, sample_b_label: str, min_expansion_bp: int, *, np, pd):
    calls = _collect_valid_haplotype_calls(data, sample_a_label, sample_b_label, pd=pd)
    wide = _pair_haplotype_lengths(calls, sample_a_label, sample_b_label, np=np, pd=pd)
    candidates = _annotate_pairing_candidates(wide, data, pd=pd)
    return _filter_pairing_targets(candidates, min_expansion_bp, pd=pd)


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


def _collect_reads_for_loci(args: TrPostArgs, loci, *, pd):
    target_set = set(loci["region"])
    region_to_locus = loci.set_index("region")["LocusID"].to_dict()
    rows = _parse_target_reads(args.group_a_fasta, args.sample_a_label, target_set, region_to_locus)
    rows.extend(_parse_target_reads(args.group_b_fasta, args.sample_b_label, target_set, region_to_locus))
    if not rows:
        return pd.DataFrame(columns=["cell_type", "region", "LocusID", "hap", "read_name", "read_length"])
    df_reads = pd.DataFrame(rows)
    return df_reads.merge(
        loci[["LocusID", "chrom", "start", "end", "igv_locus"]],
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


def _union_targets(primary_targets, rescue_targets, *, pd):
    if primary_targets.empty:
        return rescue_targets.reset_index(drop=True)
    if rescue_targets.empty:
        return primary_targets.reset_index(drop=True)
    combined = pd.concat([primary_targets, rescue_targets], ignore_index=True)
    combined = combined.drop_duplicates(subset=["LocusID"], keep="first")
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    combined["chrom"] = pd.Categorical(combined["chrom"], categories=chrom_order, ordered=True)
    return combined.sort_values(
        ["max_abs_delta_bp", "chrom", "start"],
        ascending=[False, True, True],
    ).reset_index(drop=True)


def _empty_tail_metrics() -> dict[str, Any]:
    return {
        "tail_read_count": 0,
        "tail_far_read_count": 0,
        "tail_baseline_same_hap_max_bp": None,
        "tail_change_other_hap_upper_bp": None,
        "tail_anchor_bp": None,
        "tail_max_excess_bp": None,
        "sample_change_lower_bp": None,
        "sample_change_upper_bp": None,
        "sample_baseline_lower_bp": None,
        "sample_baseline_upper_bp": None,
        "sample_lower_delta_bp": None,
        "sample_upper_delta_bp": None,
        "sample_range_supports_tail": None,
    }


def _compute_tail_expansion_metrics(
    change_lengths,
    baseline_lengths,
    change_other_lengths,
    *,
    np,
) -> dict[str, Any]:
    if not change_lengths or not baseline_lengths or not change_other_lengths:
        return _empty_tail_metrics()

    baseline_same_hap_max = int(max(baseline_lengths))
    change_other_hap_upper = int(np.ceil(np.percentile(change_other_lengths, 95)))
    anchor = max(baseline_same_hap_max, change_other_hap_upper)
    tail_excesses = [
        int(length) - anchor
        for length in change_lengths
        if int(length) > anchor + _TAIL_MARGIN_BP
    ]
    far_excesses = [excess for excess in tail_excesses if excess >= _TAIL_FAR_MARGIN_BP]
    return {
        "tail_read_count": len(tail_excesses),
        "tail_far_read_count": len(far_excesses),
        "tail_baseline_same_hap_max_bp": baseline_same_hap_max,
        "tail_change_other_hap_upper_bp": change_other_hap_upper,
        "tail_anchor_bp": int(anchor),
        "tail_max_excess_bp": int(max(tail_excesses)) if tail_excesses else None,
    }


def _tail_support_score(metrics: dict[str, Any]) -> tuple[int, int, int]:
    return (
        int(metrics.get("tail_far_read_count") or 0),
        int(metrics.get("tail_read_count") or 0),
        int(metrics.get("tail_max_excess_bp") or 0),
    )


def _build_sample_range_lookup(data, sample_labels: tuple[str, str], *, pd):
    columns = [
        "LocusID",
        "sample",
        "hap",
        "length_range_lower",
        "length_range_upper",
    ]
    sample_tables = data.get("sample", {})
    rows: list[Any] = []
    for sample_label in sample_labels:
        table = sample_tables.get(sample_label)
        if table is None or table.empty:
            continue
        required = {"LocusID", "haplotype", "length_range_lower", "length_range_upper"}
        if not required.issubset(set(table.columns)):
            continue
        subset = table[list(required)].copy()
        subset = subset[subset["haplotype"].isin([0, 1])].copy()
        subset["sample"] = sample_label
        subset["hap"] = subset["haplotype"].map({0: "hap1", 1: "hap2"})
        valid = (
            subset.groupby(["sample", "LocusID", "hap"], as_index=False)
            .size()
            .query("size == 1")[["sample", "LocusID", "hap"]]
        )
        subset = subset.merge(valid, on=["sample", "LocusID", "hap"], how="inner")
        rows.append(subset[columns])
    if not rows:
        return {}
    combined = pd.concat(rows, ignore_index=True)
    return {
        (row["LocusID"], row["sample"], row["hap"]): (
            int(row["length_range_lower"]),
            int(row["length_range_upper"]),
        )
        for _, row in combined.iterrows()
    }


def _evaluate_tail_direction(
    *,
    change_lengths,
    baseline_lengths,
    change_other_lengths,
    baseline_other_lengths,
    change_sample_range,
    baseline_sample_range,
    require_sample_range_support: bool,
    np,
) -> dict[str, Any]:
    metrics = _compute_tail_expansion_metrics(
        change_lengths,
        baseline_lengths,
        change_other_lengths,
        np=np,
    )
    metrics["n_change_reads"] = len(change_lengths)
    metrics["n_baseline_reads"] = len(baseline_lengths)
    if change_other_lengths and baseline_other_lengths:
        metrics["other_hap_median_delta_bp"] = abs(
            float(np.median(change_other_lengths)) - float(np.median(baseline_other_lengths))
        )
    else:
        metrics["other_hap_median_delta_bp"] = None
    if change_sample_range is not None and baseline_sample_range is not None:
        change_lower, change_upper = change_sample_range
        baseline_lower, baseline_upper = baseline_sample_range
        sample_lower_delta_bp = int(change_lower) - int(baseline_lower)
        sample_upper_delta_bp = int(change_upper) - int(baseline_upper)
        metrics["sample_change_lower_bp"] = int(change_lower)
        metrics["sample_change_upper_bp"] = int(change_upper)
        metrics["sample_baseline_lower_bp"] = int(baseline_lower)
        metrics["sample_baseline_upper_bp"] = int(baseline_upper)
        metrics["sample_lower_delta_bp"] = sample_lower_delta_bp
        metrics["sample_upper_delta_bp"] = sample_upper_delta_bp
        metrics["sample_range_supports_tail"] = bool(sample_upper_delta_bp >= _TAIL_MIN_SAMPLE_UPPER_DELTA_BP)
    else:
        metrics["sample_change_lower_bp"] = None
        metrics["sample_change_upper_bp"] = None
        metrics["sample_baseline_lower_bp"] = None
        metrics["sample_baseline_upper_bp"] = None
        metrics["sample_lower_delta_bp"] = None
        metrics["sample_upper_delta_bp"] = None
        metrics["sample_range_supports_tail"] = False
    metrics["eligible"] = bool(
        metrics["n_change_reads"] >= _TAIL_MIN_READS_PER_GROUP
        and metrics["n_baseline_reads"] >= _TAIL_MIN_READS_PER_GROUP
        and metrics["tail_read_count"] >= 2
        and metrics["tail_far_read_count"] >= 1
        and metrics["tail_max_excess_bp"] is not None
        and metrics["tail_max_excess_bp"] >= _TAIL_FAR_MARGIN_BP
        and metrics["other_hap_median_delta_bp"] is not None
        and metrics["other_hap_median_delta_bp"] <= _TAIL_MAX_OTHER_HAP_MEDIAN_DELTA_BP
        and (
            bool(metrics["sample_range_supports_tail"])
            if require_sample_range_support
            else True
        )
    )
    return metrics


def _build_tail_rescue_specs(
    data,
    pairing_candidates,
    df_reads,
    *,
    sample_a_label: str,
    sample_b_label: str,
    min_expansion_bp: int,
    require_sample_range_support: bool,
    np,
    pd,
):
    columns = [
        "LocusID",
        "change_allele",
        "change_celltype",
        "baseline_celltype",
        "signal_class",
        "sample_change_lower_bp",
        "sample_change_upper_bp",
        "sample_baseline_lower_bp",
        "sample_baseline_upper_bp",
        "sample_lower_delta_bp",
        "sample_upper_delta_bp",
        "sample_range_supports_tail",
    ]
    if pairing_candidates.empty or df_reads.empty:
        return pd.DataFrame(columns=columns)

    length_lookup: dict[tuple[Any, str, str], list[int]] = {}
    for (locus_id, cell_type, hap), group in df_reads.groupby(["LocusID", "cell_type", "hap"], sort=False):
        length_lookup[(locus_id, cell_type, hap)] = group["read_length"].astype(int).tolist()
    sample_range_lookup = _build_sample_range_lookup(
        data,
        (sample_a_label, sample_b_label),
        pd=pd,
    )

    rows: list[dict[str, Any]] = []
    for _, locus in pairing_candidates.iterrows():
        if int(locus["max_abs_delta_bp"]) >= min_expansion_bp:
            continue
        locus_id = locus["LocusID"]
        for allele_key, other_hap in (("hap1", "hap2"), ("hap2", "hap1")):
            a_same = length_lookup.get((locus_id, sample_a_label, allele_key), [])
            b_same = length_lookup.get((locus_id, sample_b_label, allele_key), [])
            a_other = length_lookup.get((locus_id, sample_a_label, other_hap), [])
            b_other = length_lookup.get((locus_id, sample_b_label, other_hap), [])
            a_same_range = sample_range_lookup.get((locus_id, sample_a_label, allele_key))
            b_same_range = sample_range_lookup.get((locus_id, sample_b_label, allele_key))
            if not a_same or not b_same or not a_other or not b_other:
                continue

            candidate_a = _evaluate_tail_direction(
                change_lengths=a_same,
                baseline_lengths=b_same,
                change_other_lengths=a_other,
                baseline_other_lengths=b_other,
                change_sample_range=a_same_range,
                baseline_sample_range=b_same_range,
                require_sample_range_support=require_sample_range_support,
                np=np,
            )
            candidate_b = _evaluate_tail_direction(
                change_lengths=b_same,
                baseline_lengths=a_same,
                change_other_lengths=b_other,
                baseline_other_lengths=a_other,
                change_sample_range=b_same_range,
                baseline_sample_range=a_same_range,
                require_sample_range_support=require_sample_range_support,
                np=np,
            )
            if not candidate_a["eligible"] and not candidate_b["eligible"]:
                continue

            score_a = _tail_support_score(candidate_a)
            score_b = _tail_support_score(candidate_b)
            if candidate_a["eligible"] and candidate_b["eligible"] and score_a == score_b:
                continue

            if candidate_a["eligible"] and (not candidate_b["eligible"] or score_a > score_b):
                change_celltype, baseline_celltype = sample_a_label, sample_b_label
                chosen_metrics = candidate_a
            elif candidate_b["eligible"] and (not candidate_a["eligible"] or score_b > score_a):
                change_celltype, baseline_celltype = sample_b_label, sample_a_label
                chosen_metrics = candidate_b
            else:
                continue

            rows.append(
                {
                    "LocusID": locus_id,
                    "change_allele": allele_key,
                    "change_celltype": change_celltype,
                    "baseline_celltype": baseline_celltype,
                    "signal_class": _TAIL_SIGNAL_CLASS,
                    "sample_change_lower_bp": chosen_metrics["sample_change_lower_bp"],
                    "sample_change_upper_bp": chosen_metrics["sample_change_upper_bp"],
                    "sample_baseline_lower_bp": chosen_metrics["sample_baseline_lower_bp"],
                    "sample_baseline_upper_bp": chosen_metrics["sample_baseline_upper_bp"],
                    "sample_lower_delta_bp": chosen_metrics["sample_lower_delta_bp"],
                    "sample_upper_delta_bp": chosen_metrics["sample_upper_delta_bp"],
                    "sample_range_supports_tail": chosen_metrics["sample_range_supports_tail"],
                }
            )

    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows).drop_duplicates(
        subset=["LocusID", "change_allele", "change_celltype", "baseline_celltype"],
        keep="first",
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
    tail_rescue_specs=None,
    tail_rescue_reads=None,
    np,
    pd,
) -> Any:
    """One row per retained signal-bearing (locus, allele).

    Primary rows come from the TDB paired-haplotype delta filter. Optional
    rescue rows can be added for clean read-tail expansions even when the merged
    TDB allele lengths are unchanged between groups.
    """
    if targets.empty and (tail_rescue_specs is None or tail_rescue_specs.empty):
        return pd.DataFrame(columns=_CHANGED_ALLELE_COLS)

    def _build_read_lookup(frame) -> dict[tuple[Any, str, str], tuple[list[str], list[int]]]:
        lookup: dict[tuple[Any, str, str], tuple[list[str], list[int]]] = {}
        if frame is None or frame.empty:
            return lookup
        for (locus_id, cell_type, hap), grp in frame.groupby(
            ["LocusID", "cell_type", "hap"], sort=False
        ):
            grp_sorted = grp.sort_values("read_length", ascending=False)
            lookup[(locus_id, cell_type, hap)] = (
                grp_sorted["read_name"].tolist(),
                grp_sorted["read_length"].tolist(),
            )
        return lookup

    def _select_support_reads(
        *,
        change_names: list[str],
        change_lengths: list[int],
        baseline_names: list[str],
        baseline_lengths: list[int],
        change_tdb: int,
        baseline_tdb: int,
        signal_class: str,
        tail_metrics: dict[str, Any] | None,
    ) -> tuple[list[str], list[int], list[str], list[int]]:
        if not change_lengths:
            return [], [], [], []

        signal = str(signal_class).strip().lower()
        change_pairs = list(zip(change_names, change_lengths, strict=False))
        baseline_pairs = list(zip(baseline_names, baseline_lengths, strict=False))
        is_expansion = signal == _TAIL_SIGNAL_CLASS or int(change_tdb) >= int(baseline_tdb)
        delta_bp = abs(float(change_tdb) - float(baseline_tdb))
        min_shift = max(25.0, delta_bp * 0.25)

        if signal == _TAIL_SIGNAL_CLASS:
            anchor = tail_metrics.get("tail_anchor_bp") if tail_metrics else None
            if anchor is None:
                threshold = (
                    min(float(change_tdb), float(baseline_tdb) + min_shift)
                    if is_expansion
                    else max(float(change_tdb), float(baseline_tdb) - min_shift)
                )
            else:
                threshold = float(anchor)
        else:
            threshold = (
                min(float(change_tdb), float(baseline_tdb) + min_shift)
                if is_expansion
                else max(float(change_tdb), float(baseline_tdb) - min_shift)
            )

        if is_expansion:
            change_support = [(name, length) for name, length in change_pairs if float(length) >= threshold]
            baseline_support = [(name, length) for name, length in baseline_pairs if float(length) <= threshold]
        else:
            change_support = [(name, length) for name, length in change_pairs if float(length) <= threshold]
            baseline_support = [(name, length) for name, length in baseline_pairs if float(length) >= threshold]

        if not change_support and change_pairs:
            if is_expansion:
                best_name, best_length = max(change_pairs, key=lambda item: float(item[1]))
            else:
                best_name, best_length = min(change_pairs, key=lambda item: float(item[1]))
            change_support = [(best_name, best_length)]

        return (
            [name for name, _ in change_support],
            [int(length) for _, length in change_support],
            [name for name, _ in baseline_support],
            [int(length) for _, length in baseline_support],
        )

    # Cleaned hap assignments are the default support table. Tail rescues keep
    # their raw-haplotype support reads so the evidence is not erased by hp-fix.
    read_lookup = _build_read_lookup(df_reads)
    tail_read_lookup = _build_read_lookup(tail_rescue_reads) if tail_rescue_reads is not None else read_lookup

    target_lookup = {
        row["LocusID"]: row
        for _, row in targets.iterrows()
    }
    rows: list[dict[str, Any]] = []
    seen_keys: set[tuple[Any, str, str, str]] = set()

    def _append_row(
        *,
        tgt,
        allele_key: str,
        change_ct: str,
        baseline_ct: str,
        change_tdb: int,
        baseline_tdb: int,
        change_length_bp: int,
        signal_class: str,
        support_lookup=None,
        tail_metrics: dict[str, Any] | None = None,
    ) -> None:
        key = (tgt["LocusID"], allele_key, change_ct, baseline_ct)
        if key in seen_keys:
            return

        active_lookup = support_lookup or read_lookup
        change_names, change_lengths = active_lookup.get((tgt["LocusID"], change_ct, allele_key), ([], []))
        baseline_names, baseline_lengths = active_lookup.get((tgt["LocusID"], baseline_ct, allele_key), ([], []))
        metrics = tail_metrics or _empty_tail_metrics()
        (
            change_support_names,
            change_support_lengths,
            baseline_support_names,
            baseline_support_lengths,
        ) = _select_support_reads(
            change_names=change_names,
            change_lengths=change_lengths,
            baseline_names=baseline_names,
            baseline_lengths=baseline_lengths,
            change_tdb=change_tdb,
            baseline_tdb=baseline_tdb,
            signal_class=signal_class,
            tail_metrics=metrics,
        )
        rows.append({
            "chrom": tgt["chrom"],
            "start": int(tgt["start"]),
            "end": int(tgt["end"]),
            "LocusID": tgt["LocusID"],
            "change_allele": allele_key,
            "change_celltype": change_ct,
            "baseline_celltype": baseline_ct,
            "tdb_change_length": int(change_tdb),
            "tdb_baseline_length": int(baseline_tdb),
            "change_length_bp": int(change_length_bp),
            "read_mean_change_length": float(np.mean(change_lengths)) if change_lengths else None,
            "read_mean_baseline_length": float(np.mean(baseline_lengths)) if baseline_lengths else None,
            "n_change_reads": len(change_names),
            "n_baseline_reads": len(baseline_names),
            "n_change_support_reads": len(change_support_names),
            "n_baseline_support_reads": len(baseline_support_names),
            "change_read_names": change_names,
            "change_read_lengths": change_lengths,
            "baseline_read_names": baseline_names,
            "baseline_read_lengths": baseline_lengths,
            "change_support_read_names": change_support_names,
            "change_support_read_lengths": change_support_lengths,
            "baseline_support_read_names": baseline_support_names,
            "baseline_support_read_lengths": baseline_support_lengths,
            "pairing": tgt.get("pairing"),
            "pairing_confidence": tgt.get("pairing_confidence"),
            "max_abs_delta_bp": int(tgt["max_abs_delta_bp"]),
            "signal_class": signal_class,
            "tail_read_count": metrics.get("tail_read_count"),
            "tail_far_read_count": metrics.get("tail_far_read_count"),
            "tail_baseline_same_hap_max_bp": metrics.get("tail_baseline_same_hap_max_bp"),
            "tail_change_other_hap_upper_bp": metrics.get("tail_change_other_hap_upper_bp"),
            "tail_anchor_bp": metrics.get("tail_anchor_bp"),
            "tail_max_excess_bp": metrics.get("tail_max_excess_bp"),
            "sample_change_lower_bp": metrics.get("sample_change_lower_bp"),
            "sample_change_upper_bp": metrics.get("sample_change_upper_bp"),
            "sample_baseline_lower_bp": metrics.get("sample_baseline_lower_bp"),
            "sample_baseline_upper_bp": metrics.get("sample_baseline_upper_bp"),
            "sample_lower_delta_bp": metrics.get("sample_lower_delta_bp"),
            "sample_upper_delta_bp": metrics.get("sample_upper_delta_bp"),
            "sample_range_supports_tail": metrics.get("sample_range_supports_tail"),
        })
        seen_keys.add(key)

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

            _append_row(
                tgt=tgt,
                allele_key=allele_key,
                change_ct=change_ct,
                baseline_ct=baseline_ct,
                change_tdb=change_tdb,
                baseline_tdb=baseline_tdb,
                change_length_bp=abs(delta),
                signal_class=_TDB_SIGNAL_CLASS,
                support_lookup=read_lookup,
            )

    if tail_rescue_specs is not None and not tail_rescue_specs.empty:
        opposite_hap = {"hap1": "hap2", "hap2": "hap1"}
        for rescue in tail_rescue_specs.itertuples(index=False):
            tgt = target_lookup.get(rescue.LocusID)
            if tgt is None:
                continue
            allele_key = rescue.change_allele
            other_hap = opposite_hap.get(allele_key)
            if other_hap is None:
                continue

            change_names, change_lengths = tail_read_lookup.get((rescue.LocusID, rescue.change_celltype, allele_key), ([], []))
            baseline_names, baseline_lengths = tail_read_lookup.get((rescue.LocusID, rescue.baseline_celltype, allele_key), ([], []))
            _, change_other_lengths = tail_read_lookup.get((rescue.LocusID, rescue.change_celltype, other_hap), ([], []))
            if not change_lengths or not baseline_lengths or not change_other_lengths:
                continue

            tail_metrics = _compute_tail_expansion_metrics(
                change_lengths,
                baseline_lengths,
                change_other_lengths,
                np=np,
            )
            tail_metrics.update(
                {
                    "sample_change_lower_bp": getattr(rescue, "sample_change_lower_bp", None),
                    "sample_change_upper_bp": getattr(rescue, "sample_change_upper_bp", None),
                    "sample_baseline_lower_bp": getattr(rescue, "sample_baseline_lower_bp", None),
                    "sample_baseline_upper_bp": getattr(rescue, "sample_baseline_upper_bp", None),
                    "sample_lower_delta_bp": getattr(rescue, "sample_lower_delta_bp", None),
                    "sample_upper_delta_bp": getattr(rescue, "sample_upper_delta_bp", None),
                    "sample_range_supports_tail": getattr(rescue, "sample_range_supports_tail", None),
                }
            )
            if tail_metrics["tail_read_count"] < 2 or tail_metrics["tail_far_read_count"] < 1:
                continue

            group_prefix = "group_a" if rescue.change_celltype == sample_a_label else "group_b"
            base_prefix = "group_a" if rescue.baseline_celltype == sample_a_label else "group_b"
            if allele_key == "hap1":
                change_tdb = int(tgt[f"{group_prefix}_hap1_length"])
                baseline_tdb = int(tgt[f"{base_prefix}_hap1_length"])
            else:
                change_tdb = int(tgt[f"{group_prefix}_hap2_length"])
                baseline_tdb = int(tgt[f"{base_prefix}_hap2_length"])

            _append_row(
                tgt=tgt,
                allele_key=allele_key,
                change_ct=rescue.change_celltype,
                baseline_ct=rescue.baseline_celltype,
                change_tdb=change_tdb,
                baseline_tdb=baseline_tdb,
                change_length_bp=int(tail_metrics["tail_max_excess_bp"] or 0),
                signal_class=_TAIL_SIGNAL_CLASS,
                support_lookup=tail_read_lookup,
                tail_metrics=tail_metrics,
            )

    if not rows:
        return pd.DataFrame(columns=_CHANGED_ALLELE_COLS)

    df = pd.DataFrame(rows)
    chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    df["chrom"] = pd.Categorical(df["chrom"], categories=chrom_order, ordered=True)
    return df.sort_values(
        ["max_abs_delta_bp", "chrom", "start"], ascending=[False, True, True]
    ).reset_index(drop=True)


def _build_plot_dataframe(
    df_reads,
    targets,
    *,
    changed_allele_table,
    sample_a_label: str,
    sample_b_label: str,
    max_fold: float,
    np,
    pd,
):
    if df_reads.empty or targets.empty:
        return None

    cell_offset = {sample_a_label: -0.18, sample_b_label: 0.18}
    rng = np.random.default_rng(42)
    locus_order = targets["LocusID"].tolist()
    locus_to_center = {lid: idx for idx, lid in enumerate(locus_order)}
    locus_to_label = targets.set_index("LocusID")["igv_locus"].to_dict()

    changed_keys: set[tuple[Any, str, str]] = set()
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

    cleaned_groups = [
        _fold_clean(group)
        for _, group in df_reads.groupby(["LocusID", "cell_type", "hap"], sort=False)
    ]
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
    return {
        "plot_df": plot_df,
        "changed_keys": changed_keys,
        "cell_offset": cell_offset,
        "locus_order": locus_order,
        "locus_to_center": locus_to_center,
        "locus_to_label": locus_to_label,
    }


def _render_overview_plot(
    plot_context: dict[str, Any],
    *,
    group_a: str,
    group_b: str,
    sample_a_label: str,
    sample_b_label: str,
    outpath: Path,
    plt,
    sns,
):
    plot_df = plot_context["plot_df"]
    changed_keys = plot_context["changed_keys"]
    cell_offset = plot_context["cell_offset"]
    locus_order = plot_context["locus_order"]
    locus_to_center = plot_context["locus_to_center"]
    locus_to_label = plot_context["locus_to_label"]

    palette = {"hap1": "#4C78A8", "hap2": "#F58518", "unknown": "#888888"}
    changed_color = "#D62728"
    fig_width = max(14, len(locus_order) * 0.55)
    fig, ax = plt.subplots(figsize=(fig_width, 7))

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
        ax.text(
            locus_to_center[locus_id],
            label_y,
            locus_to_label[locus_id],
            rotation=90,
            ha="center",
            va="bottom",
            fontsize=7,
        )
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
    ax.legend(
        unique_handles,
        unique_labels,
        title="Haplotype / status",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=True,
    )
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight", dpi=150)
    plt.close(fig)


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
    plot_context = _build_plot_dataframe(
        df_reads,
        targets,
        changed_allele_table=changed_allele_table,
        sample_a_label=sample_a_label,
        sample_b_label=sample_b_label,
        max_fold=max_fold,
        np=np,
        pd=pd,
    )
    if plot_context is None:
        return
    _render_overview_plot(
        plot_context,
        group_a=group_a,
        group_b=group_b,
        sample_a_label=sample_a_label,
        sample_b_label=sample_b_label,
        outpath=outpath,
        plt=plt,
        sns=sns,
    )



def _write_dataframe_tsv(df, path: Path, *, json_cols: set[str] | None = None) -> None:
    json_cols = json_cols or set()
    out = df.copy()
    for col in json_cols:
        if col in out.columns:
            out[col] = out[col].apply(lambda value: json.dumps(value) if isinstance(value, list) else value)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, sep="\t", index=False)


def _metric_or_none(value) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _bool_or_none(value) -> bool | None:
    if value is None:
        return None
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    return None


def _build_haplotype_qc_table(clustered_summary, df_reads, *, pd):
    base_cols = [
        "LocusID",
        "cell_type",
        "hap",
        "median_length",
        "has_cross_hap_mixing",
    ]
    if clustered_summary.empty:
        return pd.DataFrame(
            columns=base_cols + ["n_hp_changed", "hp_changed_fraction"]
        )

    hp_changed = (
        df_reads.groupby(["LocusID", "cell_type", "hap"], as_index=False)
        .agg(
            n_hp_changed=("hp_changed", lambda s: int(s.fillna(False).sum())),
            hp_changed_fraction=("hp_changed", lambda s: float(s.fillna(False).mean()) if len(s) else 0.0),
        )
    )
    out = clustered_summary[base_cols].merge(
        hp_changed,
        on=["LocusID", "cell_type", "hap"],
        how="left",
    )
    out["n_hp_changed"] = out["n_hp_changed"].fillna(0).astype(int)
    out["hp_changed_fraction"] = out["hp_changed_fraction"].fillna(0.0).astype(float)
    return out


def _classify_tr_row(row, *, require_sample_range_support: bool = True) -> tuple[str, bool]:
    signal_class = str(row.get("signal_class", _TDB_SIGNAL_CLASS)).strip().lower()
    change_type = str(row.get("change_type", "")).strip().lower()
    change_length_bp = _metric_or_none(row.get("change_length_bp"))
    n_change_reads = _metric_or_none(row.get("n_change_reads"))
    n_baseline_reads = _metric_or_none(row.get("n_baseline_reads"))
    median_shift_bp = _metric_or_none(row.get("median_shift_bp"))
    median_shift_ratio = _metric_or_none(row.get("median_shift_ratio"))
    other_hap_median_delta_bp = _metric_or_none(row.get("other_hap_median_delta_bp"))
    hp_changed_fraction = _metric_or_none(row.get("hp_changed_fraction"))
    pairing_confidence = _metric_or_none(row.get("pairing_confidence"))
    change_cross = _bool_or_none(row.get("change_cross_hap_mixing"))
    baseline_cross = _bool_or_none(row.get("baseline_cross_hap_mixing"))
    sample_range_supports_tail = _bool_or_none(row.get("sample_range_supports_tail"))
    tail_read_count = _metric_or_none(row.get("tail_read_count"))
    tail_far_read_count = _metric_or_none(row.get("tail_far_read_count"))
    tail_max_excess_bp = _metric_or_none(row.get("tail_max_excess_bp"))

    if signal_class == _TAIL_SIGNAL_CLASS:
        supportive = (
            change_type == "expansion"
            and n_change_reads is not None and n_change_reads >= _TAIL_MIN_READS_PER_GROUP
            and n_baseline_reads is not None and n_baseline_reads >= _TAIL_MIN_READS_PER_GROUP
            and tail_read_count is not None and tail_read_count >= 2.0
            and tail_far_read_count is not None and tail_far_read_count >= 1.0
            and tail_max_excess_bp is not None and tail_max_excess_bp >= _TAIL_FAR_MARGIN_BP
            and other_hap_median_delta_bp is not None and other_hap_median_delta_bp <= _TAIL_MAX_OTHER_HAP_MEDIAN_DELTA_BP
            and (
                sample_range_supports_tail is True
                if require_sample_range_support
                else True
            )
            and hp_changed_fraction is not None and hp_changed_fraction <= 0.15
            and change_cross is False
            and baseline_cross is False
        )
        if not supportive:
            return "weak", False
        strong = (
            tail_read_count >= 3.0
            and tail_far_read_count >= 2.0
            and tail_max_excess_bp >= 100.0
            and hp_changed_fraction == 0.0
        )
        return ("strong", False) if strong else ("supportive", False)

    _no_cross_mix = change_cross is False and baseline_cross is False
    supportive = (
        change_type == "expansion"
        and change_length_bp is not None and change_length_bp >= 100.0
        and n_change_reads is not None and n_change_reads >= 8.0
        and n_baseline_reads is not None and n_baseline_reads >= 8.0
        and median_shift_bp is not None and median_shift_bp >= 50.0
        and other_hap_median_delta_bp is not None and other_hap_median_delta_bp <= 50.0
        and hp_changed_fraction is not None and hp_changed_fraction <= 0.15
        and _no_cross_mix
    )
    # Path A: moderate expansion with fewer reads/smaller shift, gated by pairing confidence.
    # pairing_confidence == 0.0 means both haplotype orientations are equally plausible (no
    # directional signal); every noise/artifact row in benchmarking had pairing_confidence=0.0
    # while true positives had pairing_confidence >= 0.735.
    if not supportive:
        supportive = (
            change_type == "expansion"
            and change_length_bp is not None and change_length_bp >= 50.0
            and n_change_reads is not None and n_change_reads >= 5.0
            and n_baseline_reads is not None and n_baseline_reads >= 10.0
            and median_shift_bp is not None and median_shift_bp >= 30.0
            and other_hap_median_delta_bp is not None and other_hap_median_delta_bp <= 30.0
            and hp_changed_fraction is not None and hp_changed_fraction == 0.0
            and _no_cross_mix
            and pairing_confidence is not None and pairing_confidence >= 0.05
        )
    # Path B: small confirmed expansion, very high pairing confidence required.
    if not supportive:
        supportive = (
            change_type == "expansion"
            and change_length_bp is not None and change_length_bp >= 15.0
            and n_change_reads is not None and n_change_reads >= 8.0
            and n_baseline_reads is not None and n_baseline_reads >= 10.0
            and median_shift_bp is not None and median_shift_bp >= 15.0
            and other_hap_median_delta_bp is not None and other_hap_median_delta_bp <= 20.0
            and hp_changed_fraction is not None and hp_changed_fraction == 0.0
            and _no_cross_mix
            and pairing_confidence is not None and pairing_confidence >= 0.75
        )
    if not supportive:
        return "weak", False

    strong = (
        n_change_reads >= 10.0
        and n_baseline_reads >= 10.0
        and median_shift_bp >= 100.0
        and other_hap_median_delta_bp <= 25.0
        and hp_changed_fraction == 0.0
        and (
            (pairing_confidence is not None and pairing_confidence >= 0.05)
            or (median_shift_ratio is not None and median_shift_ratio >= 1.15)
        )
    )
    return ("strong", True) if strong else ("supportive", True)


def _build_tr_bed_table(changed_allele_table, clustered_summary, df_reads, *, pd, require_sample_range_support: bool = True):
    if changed_allele_table.empty:
        return pd.DataFrame(columns=_TR_BED_COLS)

    n_change_support_reads = changed_allele_table.get("n_change_support_reads", changed_allele_table["n_change_reads"])
    n_baseline_support_reads = changed_allele_table.get("n_baseline_support_reads", changed_allele_table["n_baseline_reads"])
    change_support_read_names = changed_allele_table.get("change_support_read_names", changed_allele_table["change_read_names"])
    baseline_support_read_names = changed_allele_table.get("baseline_support_read_names", changed_allele_table["baseline_read_names"])

    hap_qc = _build_haplotype_qc_table(clustered_summary, df_reads, pd=pd)
    changed_metrics = hap_qc.rename(
        columns={
            "LocusID": "trid",
            "cell_type": "change_group",
            "hap": "change_allele",
            "median_length": "change_median_bp",
            "has_cross_hap_mixing": "change_cross_hap_mixing",
        }
    )
    baseline_metrics = hap_qc.rename(
        columns={
            "LocusID": "trid",
            "cell_type": "baseline_group",
            "hap": "change_allele",
            "median_length": "baseline_median_bp",
            "has_cross_hap_mixing": "baseline_cross_hap_mixing",
        }
    )
    other_metrics = hap_qc.rename(
        columns={
            "LocusID": "trid",
            "cell_type": "metric_group",
            "hap": "other_hap",
            "median_length": "other_hap_median_bp",
        }
    )

    bed = pd.DataFrame(
        {
            "chrom": changed_allele_table["chrom"].astype(str),
            "start": changed_allele_table["start"].astype(int),
            "end": changed_allele_table["end"].astype(int),
            "trid": changed_allele_table["LocusID"].astype(int),
            "change_allele": changed_allele_table["change_allele"],
            "change_type": [
                "expansion"
                if str(signal_class).strip().lower() == _TAIL_SIGNAL_CLASS
                or int(change_len) > int(base_len)
                else "contraction"
                for signal_class, change_len, base_len in zip(
                    changed_allele_table["signal_class"],
                    changed_allele_table["tdb_change_length"],
                    changed_allele_table["tdb_baseline_length"],
                )
            ],
            "change_group": changed_allele_table["change_celltype"],
            "baseline_group": changed_allele_table["baseline_celltype"],
            "n_change_reads": changed_allele_table["n_change_reads"].astype(int),
            "n_baseline_reads": changed_allele_table["n_baseline_reads"].astype(int),
            "n_change_support_reads": pd.to_numeric(n_change_support_reads, errors="coerce").fillna(0).astype(int),
            "n_baseline_support_reads": pd.to_numeric(n_baseline_support_reads, errors="coerce").fillna(0).astype(int),
            "tdb_change_length": changed_allele_table["tdb_change_length"].astype(int),
            "tdb_baseline_length": changed_allele_table["tdb_baseline_length"].astype(int),
            "change_length_bp": changed_allele_table["change_length_bp"].astype(int),
            "change_read_mean": [
                round(float(value), 1) if value is not None else "."
                for value in changed_allele_table["read_mean_change_length"]
            ],
            "baseline_read_mean": [
                round(float(value), 1) if value is not None else "."
                for value in changed_allele_table["read_mean_baseline_length"]
            ],
            "change_read_range": [
                f"{min(lengths)}-{max(lengths)}" if lengths else "."
                for lengths in changed_allele_table["change_read_lengths"]
            ],
            "baseline_read_range": [
                f"{min(lengths)}-{max(lengths)}" if lengths else "."
                for lengths in changed_allele_table["baseline_read_lengths"]
            ],
            "change_read_names": [value or [] for value in changed_allele_table["change_read_names"]],
            "baseline_read_names": [value or [] for value in changed_allele_table["baseline_read_names"]],
            "change_support_read_names": [value or [] for value in change_support_read_names],
            "baseline_support_read_names": [value or [] for value in baseline_support_read_names],
            "pairing": changed_allele_table["pairing"],
            "pairing_confidence": changed_allele_table["pairing_confidence"],
            "signal_class": changed_allele_table["signal_class"].fillna(_TDB_SIGNAL_CLASS),
            "tail_read_count": changed_allele_table["tail_read_count"],
            "tail_far_read_count": changed_allele_table["tail_far_read_count"],
            "tail_baseline_same_hap_max_bp": changed_allele_table["tail_baseline_same_hap_max_bp"],
            "tail_change_other_hap_upper_bp": changed_allele_table["tail_change_other_hap_upper_bp"],
            "tail_anchor_bp": changed_allele_table["tail_anchor_bp"],
            "tail_max_excess_bp": changed_allele_table["tail_max_excess_bp"],
            "sample_change_lower_bp": changed_allele_table["sample_change_lower_bp"],
            "sample_change_upper_bp": changed_allele_table["sample_change_upper_bp"],
            "sample_baseline_lower_bp": changed_allele_table["sample_baseline_lower_bp"],
            "sample_baseline_upper_bp": changed_allele_table["sample_baseline_upper_bp"],
            "sample_lower_delta_bp": changed_allele_table["sample_lower_delta_bp"],
            "sample_upper_delta_bp": changed_allele_table["sample_upper_delta_bp"],
            "sample_range_supports_tail": changed_allele_table["sample_range_supports_tail"],
        }
    )

    bed = bed.merge(
        changed_metrics[
            [
                "trid",
                "change_group",
                "change_allele",
                "change_median_bp",
                "n_hp_changed",
                "hp_changed_fraction",
                "change_cross_hap_mixing",
            ]
        ],
        on=["trid", "change_group", "change_allele"],
        how="left",
    )
    bed = bed.merge(
        baseline_metrics[
            [
                "trid",
                "baseline_group",
                "change_allele",
                "baseline_median_bp",
                "baseline_cross_hap_mixing",
            ]
        ],
        on=["trid", "baseline_group", "change_allele"],
        how="left",
    )

    opposite_hap = {"hap1": "hap2", "hap2": "hap1"}
    bed["other_hap"] = bed["change_allele"].map(opposite_hap).fillna("unknown")
    bed = bed.merge(
        other_metrics[["trid", "metric_group", "other_hap", "other_hap_median_bp"]].rename(
            columns={
                "metric_group": "change_group",
                "other_hap_median_bp": "change_other_hap_median_bp",
            }
        ),
        on=["trid", "change_group", "other_hap"],
        how="left",
    )
    bed = bed.merge(
        other_metrics[["trid", "metric_group", "other_hap", "other_hap_median_bp"]].rename(
            columns={
                "metric_group": "baseline_group",
                "other_hap_median_bp": "baseline_other_hap_median_bp",
            }
        ),
        on=["trid", "baseline_group", "other_hap"],
        how="left",
    )

    bed["median_shift_bp"] = bed["change_median_bp"] - bed["baseline_median_bp"]
    bed["median_shift_ratio"] = bed["change_median_bp"] / bed["baseline_median_bp"]
    bed["other_hap_median_delta_bp"] = (
        bed["change_other_hap_median_bp"] - bed["baseline_other_hap_median_bp"]
    ).abs()

    tiers = bed.apply(
        lambda row: _classify_tr_row(
            row,
            require_sample_range_support=require_sample_range_support,
        ),
        axis=1,
        result_type="expand",
    )
    bed[["tr_tier", "tr_pass_for_harmonized"]] = tiers
    bed["tr_tier"] = bed["tr_tier"].fillna("weak")
    bed["tr_pass_for_harmonized"] = bed["tr_pass_for_harmonized"].fillna(False).astype(bool)
    bed["_tier_rank"] = bed["tr_tier"].map(_TR_TIER_ORDER).fillna(99).astype(int)
    bed["_pairing_conf_sort"] = pd.to_numeric(bed["pairing_confidence"], errors="coerce").fillna(-1.0)
    bed["_tail_excess_sort"] = pd.to_numeric(bed["tail_max_excess_bp"], errors="coerce").fillna(-1.0)

    bed = bed.sort_values(
        [
            "tr_pass_for_harmonized",
            "_tier_rank",
            "_tail_excess_sort",
            "change_length_bp",
            "n_change_reads",
            "_pairing_conf_sort",
        ],
        ascending=[False, True, False, False, False, False],
    ).reset_index(drop=True)
    return bed[_TR_BED_COLS]


def _write_tr_bed_tsv(tr_bed_table, path: Path) -> None:
    _write_dataframe_tsv(
        tr_bed_table,
        path,
        json_cols={
            "change_read_names",
            "baseline_read_names",
            "change_support_read_names",
            "baseline_support_read_names",
        },
    )


def tr_post_processing_main(cli_args=None) -> dict[str, Any]:
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(message)s")
    parser = _build_arg_parser()
    args = _resolve_args(parser.parse_args(cli_args))
    outputs = _build_output_paths(args.output_dir)

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
    data = tdb.load_tdb(str(args.merged_tdb))
    calls = _collect_valid_haplotype_calls(data, args.sample_a_label, args.sample_b_label, pd=pd)
    wide = _pair_haplotype_lengths(calls, args.sample_a_label, args.sample_b_label, np=np, pd=pd)
    pairing_candidates = _annotate_pairing_candidates(wide, data, pd=pd)
    primary_targets = _filter_pairing_targets(pairing_candidates, args.min_expansion_bp, pd=pd)

    all_reads = None
    tail_rescue_specs = pd.DataFrame(
        columns=["LocusID", "change_allele", "change_celltype", "baseline_celltype", "signal_class"]
    )
    targets = primary_targets
    if args.tail_expansion_rescue and not pairing_candidates.empty:
        all_reads = _collect_reads_for_loci(args, pairing_candidates, pd=pd)
        if not all_reads.empty:
            all_reads = all_reads[all_reads["read_length"] > 1].copy()
            tail_rescue_specs = _build_tail_rescue_specs(
                data,
                pairing_candidates,
                all_reads,
                sample_a_label=args.sample_a_label,
                sample_b_label=args.sample_b_label,
                min_expansion_bp=args.min_expansion_bp,
                require_sample_range_support=args.tail_require_sample_range_support,
                np=np,
                pd=pd,
            )
            if not tail_rescue_specs.empty:
                rescue_targets = pairing_candidates[
                    pairing_candidates["LocusID"].isin(tail_rescue_specs["LocusID"])
                ].copy()
                targets = _union_targets(primary_targets, rescue_targets, pd=pd)
    _write_dataframe_tsv(targets, outputs.targets_tsv)
    if targets.empty:
        summary = _build_summary_payload(
            args,
            outputs,
            status="completed_empty",
            reason=(
                "No loci passed the paired-haplotype delta filter"
                if not args.tail_expansion_rescue
                else "No loci passed the paired-haplotype delta filter or tail-expansion rescue"
            ),
        )
        _write_tsv(outputs.read_lengths_tsv, [], ["cell_type", "region", "LocusID", "hap", "read_name", "read_length"])
        _write_tsv(outputs.summary_tsv, [], ["LocusID", "cell_type", "hap", "n_reads", "mean_length", "median_length"])
        _write_tsv(outputs.summary_clustered_tsv, [], ["LocusID", "cell_type", "hap", "is_two_tight_groups", "has_cross_hap_mixing"])
        _write_tsv(outputs.tr_bed_tsv, [], _TR_BED_COLS)
        _write_json(outputs.summary_json, summary)
        return summary

    if all_reads is None:
        df_reads = _collect_reads_for_loci(args, targets, pd=pd)
    else:
        df_reads = all_reads[all_reads["LocusID"].isin(set(targets["LocusID"]))].copy()
    if df_reads.empty:
        return _write_empty_outputs(args, status="completed_empty", reason="No trimmed reads matched the retained target loci")
    df_reads = df_reads[df_reads["read_length"] > 1].copy()
    df_reads_raw = df_reads.copy()
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
        tail_rescue_specs=tail_rescue_specs,
        tail_rescue_reads=df_reads_raw,
        np=np,
        pd=pd,
    )
    tr_bed_table = _build_tr_bed_table(
        changed_allele_table,
        clustered_summary,
        df_reads,
        pd=pd,
        require_sample_range_support=args.tail_require_sample_range_support,
    )

    _READS_FRONT = ["LocusID", "chrom", "start", "end", "cell_type", "hap", "read_name", "read_length"]
    _reads_rest = [c for c in df_reads.columns if c not in _READS_FRONT]
    _write_dataframe_tsv(df_reads[_READS_FRONT + _reads_rest], outputs.read_lengths_tsv)
    _write_dataframe_tsv(final_table, outputs.summary_tsv, json_cols={"read_lengths"})
    _write_dataframe_tsv(clustered_summary, outputs.summary_clustered_tsv, json_cols={"read_lengths"})
    _write_tr_bed_tsv(tr_bed_table, outputs.tr_bed_tsv)

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
                outpath=outputs.overview_plot,
                np=np,
                pd=pd,
                plt=plt,
                sns=sns,
            )
            if outputs.overview_plot.exists():
                plots.append(str(outputs.overview_plot))
        except Exception as exc:
            logging.warning("TR post-processing plots were skipped: %s", exc)

    summary = _build_summary_payload(
        args,
        outputs,
        status="completed",
        n_targets=len(targets),
        n_read_rows=len(df_reads),
        n_summary_rows=len(clustered_summary),
        n_tr_bed_rows=len(tr_bed_table),
        n_tr_strong_rows=int((tr_bed_table["tr_tier"] == "strong").sum()) if not tr_bed_table.empty else 0,
        n_tr_supportive_rows=int((tr_bed_table["tr_tier"] == "supportive").sum()) if not tr_bed_table.empty else 0,
        n_tr_weak_rows=int((tr_bed_table["tr_tier"] == "weak").sum()) if not tr_bed_table.empty else 0,
        n_tail_rescue_rows=int((tr_bed_table["signal_class"] == _TAIL_SIGNAL_CLASS).sum()) if not tr_bed_table.empty else 0,
        plots=plots,
    )
    summary["params"] = {
        "min_expansion_bp": args.min_expansion_bp,
        "min_reads_per_hap": args.min_reads_per_hap,
        "max_fold": args.max_fold,
        "make_plots": args.make_plots,
        "tail_expansion_rescue": args.tail_expansion_rescue,
        "tail_require_sample_range_support": args.tail_require_sample_range_support,
    }
    _write_json(outputs.summary_json, summary)
    return summary


def main(argv: list[str] | None = None) -> int:
    tr_post_processing_main(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
