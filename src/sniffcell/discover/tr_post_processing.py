"""Read-length paired-comparison TR caller for the discover pipeline.

This module implements the "somatic scan" TR-discovery method: instead of
modelling tandem-repeat allele lengths through a merged TDB, it works directly
off the per-group spanning reads that ``medaka tandem`` already trims out
(``trimmed_reads.fasta``) and uses each read's length as a proxy for the TR
allele length at its locus.

For every locus that has spanning reads in either of the two split groups we
pool the read lengths per group (haplotype-agnostic) and apply a paired,
directional comparison:

    A locus is a confident TR change in favour of one group ("change group")
    when its top ``min_supporting_reads`` reads are EACH longer than the other
    group's longest read by at least ``margin_bp``.

The longest read of the comparator group is the "anchor"; the comparator must
contribute at least one read (loci where the baseline group has no spanning
reads are skipped to avoid coverage-driven false positives). The method only
calls expansions of the change group relative to the baseline group, which is
exactly the signal the genome-wide somatic scan was built around.

``margin_bp`` and ``min_supporting_reads`` are the two tunable parameters.

The output schema (``tr_changes.bed.tsv`` + ``summary.json``) is kept
compatible with :mod:`sniffcell.discover.harmonize_variants` so the rest of the
discover pipeline is unaffected.
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import statistics
from collections import defaultdict
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


DEFAULT_MARGIN_BP = 100
DEFAULT_MIN_SUPPORTING_READS = 2

_READLEN_SIGNAL_CLASS = "readlen_paired"
_CHANGE_ALLELE = "all"  # the scan pools reads across haplotypes
_CHANGE_TYPE = "expansion"

# First ``chrN_start_end`` triple anywhere in a trimmed-read FASTA header.
_REGION_RE = re.compile(r"(chr[0-9A-Za-z]+)_(\d+)_(\d+)")

_TR_TIER_ORDER: dict[str, int] = {
    "strong": 0,
    "supportive": 1,
    "weak": 2,
}

_CHROM_ORDER: dict[str, int] = {f"chr{i}": i for i in range(1, 23)}
_CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25, "chrMT": 25})

# Columns written to tr_changes.bed.tsv. The first block is what
# harmonize_variants consumes; the rest are evidence/QC fields.
_TR_BED_COLS: list[str] = [
    "chrom", "start", "end", "trid",
    "change_allele", "change_type",
    "change_group", "baseline_group",
    "n_change_reads", "n_baseline_reads",
    "n_change_support_reads", "n_baseline_support_reads",
    "change_max_bp", "baseline_max_bp", "baseline_anchor_bp",
    "margin_bp", "min_supporting_reads",
    "change_length_bp", "change_support_min_excess_bp",
    "change_read_mean", "baseline_read_mean",
    "change_read_range", "baseline_read_range",
    "change_read_names", "baseline_read_names",
    "change_support_read_names", "baseline_support_read_names",
    "signal_class",
    "tr_tier", "tr_pass_for_harmonized",
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
    group_a_fasta: Path
    group_b_fasta: Path
    margin_bp: int
    min_supporting_reads: int
    make_plots: bool


@dataclass(frozen=True)
class TrPostOutputs:
    read_lengths_tsv: Path
    targets_tsv: Path
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
            "Call tandem-repeat changes between two split groups using a paired "
            "read-length comparison of medaka trimmed reads. A locus is called when "
            "one group's top --min-supporting-reads reads each exceed the other "
            "group's longest read by at least --margin-bp."
        ),
        add_help=add_help,
    )
    parser.add_argument("--split-dir", required=True, help="deconv_requested_group_splits directory")
    parser.add_argument("--groups", required=True, help="Exactly two group names, comma-separated")
    parser.add_argument("--output-dir", default=None, help="Output directory for the TR post-processing report")
    parser.add_argument("--sample-id", default=None, help="Optional sample ID override")
    parser.add_argument("--sample-a-label", default=None, help="Read-group label for the first group (default <sample>.<group_a>)")
    parser.add_argument("--sample-b-label", default=None, help="Read-group label for the second group (default <sample>.<group_b>)")
    parser.add_argument(
        "--group-a-fasta",
        default=None,
        help="trimmed_reads.fasta for the first group (default <split-dir>/medaka_tandem/<group_a>.medaka/trimmed_reads.fasta)",
    )
    parser.add_argument(
        "--group-b-fasta",
        default=None,
        help="trimmed_reads.fasta for the second group (default <split-dir>/medaka_tandem/<group_b>.medaka/trimmed_reads.fasta)",
    )
    parser.add_argument(
        "--margin-bp",
        type=int,
        default=DEFAULT_MARGIN_BP,
        help=(
            "Minimum length (bp) by which each supporting read of the change group must exceed "
            f"the baseline group's longest read. Default={DEFAULT_MARGIN_BP}."
        ),
    )
    parser.add_argument(
        "--min-supporting-reads",
        type=int,
        default=DEFAULT_MIN_SUPPORTING_READS,
        help=(
            "Number of the change group's longest reads that must each clear the baseline-max + margin "
            f"threshold for a locus to be called. Default={DEFAULT_MIN_SUPPORTING_READS}."
        ),
    )
    parser.add_argument(
        "--merged-tdb",
        default=None,
        help="Deprecated and ignored. The read-length method does not use a TDB; kept for backward compatibility.",
    )
    parser.add_argument("--skip-plots", action="store_true", default=False, help="Skip PNG generation")
    return parser


def _resolve_args(raw_args) -> TrPostArgs:
    split_dir = _expand_path(raw_args.split_dir)
    tokens = [x.strip() for x in str(raw_args.groups).split(",") if x.strip()]
    if len(tokens) != 2:
        raise ValueError("--groups must contain exactly two group names")
    group_a, group_b = tokens
    sample_id = raw_args.sample_id or _infer_sample_id(split_dir.parent)
    output_dir = (
        _expand_path(raw_args.output_dir)
        if raw_args.output_dir
        else split_dir / "postprocess" / f"tr_post_processing_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
    sample_a_label = raw_args.sample_a_label or f"{sample_id}.{group_a}"
    sample_b_label = raw_args.sample_b_label or f"{sample_id}.{group_b}"
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
    margin_bp = int(raw_args.margin_bp)
    min_supporting_reads = int(raw_args.min_supporting_reads)
    if margin_bp < 0:
        raise ValueError("--margin-bp must be >= 0")
    if min_supporting_reads < 1:
        raise ValueError("--min-supporting-reads must be >= 1")
    return TrPostArgs(
        split_dir=split_dir,
        output_dir=output_dir,
        sample_id=sample_id,
        group_a=group_a,
        group_b=group_b,
        sample_a_label=sample_a_label,
        sample_b_label=sample_b_label,
        group_a_fasta=group_a_fasta,
        group_b_fasta=group_b_fasta,
        margin_bp=margin_bp,
        min_supporting_reads=min_supporting_reads,
        make_plots=not bool(raw_args.skip_plots),
    )


def _write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _build_output_paths(output_dir: Path) -> TrPostOutputs:
    plots_dir = output_dir / "plots"
    return TrPostOutputs(
        read_lengths_tsv=output_dir / "read_lengths.tsv",
        targets_tsv=output_dir / "targets.tsv",
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
    n_loci_scanned: int = 0,
    n_targets: int = 0,
    n_tr_strong_rows: int = 0,
    n_tr_supportive_rows: int = 0,
    n_tr_weak_rows: int = 0,
    plots: list[str] | None = None,
) -> dict[str, Any]:
    summary = {
        "status": status,
        "method": "readlen_paired_scan",
        "sample_id": args.sample_id,
        "group_a": args.group_a,
        "group_b": args.group_b,
        "sample_a_label": args.sample_a_label,
        "sample_b_label": args.sample_b_label,
        "split_dir": str(args.split_dir),
        "group_a_fasta": str(args.group_a_fasta),
        "group_b_fasta": str(args.group_b_fasta),
        "n_loci_scanned": int(n_loci_scanned),
        "n_targets": int(n_targets),
        "n_tr_strong_rows": int(n_tr_strong_rows),
        "n_tr_supportive_rows": int(n_tr_supportive_rows),
        "n_tr_weak_rows": int(n_tr_weak_rows),
        "read_lengths_tsv": str(outputs.read_lengths_tsv),
        "targets_tsv": str(outputs.targets_tsv),
        "tr_bed_tsv": str(outputs.tr_bed_tsv),
        "plots": plots or [],
        "params": {
            "margin_bp": args.margin_bp,
            "min_supporting_reads": args.min_supporting_reads,
            "make_plots": args.make_plots,
        },
    }
    if reason is not None:
        summary["reason"] = reason
    return summary


def _write_empty_outputs(args: TrPostArgs, *, status: str, reason: str) -> dict[str, Any]:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = _build_output_paths(args.output_dir)
    _write_tsv(
        outputs.read_lengths_tsv,
        [],
        ["trid", "chrom", "start", "end", "cell_type", "read_name", "read_length"],
    )
    _write_tsv(outputs.targets_tsv, [], _TR_BED_COLS)
    _write_tsv(outputs.tr_bed_tsv, [], _TR_BED_COLS)
    summary = _build_summary_payload(args, outputs, status=status, reason=reason)
    _write_json(outputs.summary_json, summary)
    return summary


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

LocusKey = tuple[str, int, int]


def _parse_fasta_lengths(path: Path) -> dict[LocusKey, list[tuple[str, int]]]:
    """Map each TR locus to its ``(read_name, read_length)`` records.

    Read length is summed across wrapped sequence lines. The locus is taken
    from the first ``chrN_start_end`` triple in the read header (matching the
    medaka trimmed-reads naming convention).
    """
    loci: dict[LocusKey, list[tuple[str, int]]] = defaultdict(list)
    read_name: str | None = None
    key: LocusKey | None = None
    length = 0

    def _flush() -> None:
        nonlocal read_name, key, length
        if read_name is not None and key is not None and length > 0:
            loci[key].append((read_name, length))

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                _flush()
                header = line[1:]
                read_name = header.split()[0] if header.split() else header
                length = 0
                key = None
                match = _REGION_RE.search(header)
                if match:
                    key = (match.group(1), int(match.group(2)), int(match.group(3)))
            elif key is not None:
                length += len(line.strip())
        _flush()
    return loci


# ---------------------------------------------------------------------------
# Paired read-length scan
# ---------------------------------------------------------------------------

def _direction_excess(
    change_lengths_desc: list[int],
    baseline_lengths_desc: list[int],
    *,
    margin_bp: int,
    min_supporting_reads: int,
) -> int | None:
    """Excess (bp) of the change group's longest read over the baseline anchor,
    or ``None`` if this direction does not pass the paired-comparison rule.

    Passes when the baseline contributes >=1 read, the change group has at least
    ``min_supporting_reads`` reads, and each of its top ``min_supporting_reads``
    reads exceeds ``baseline_max + margin_bp``.
    """
    if not baseline_lengths_desc:
        return None
    if len(change_lengths_desc) < min_supporting_reads:
        return None
    anchor = baseline_lengths_desc[0]
    threshold = anchor + margin_bp
    top = change_lengths_desc[:min_supporting_reads]
    if all(length > threshold for length in top):
        return change_lengths_desc[0] - anchor
    return None


def _assign_tier(*, n_change_support_reads: int, n_baseline_reads: int, min_supporting_reads: int) -> str:
    """Confident calls are 'supportive'; promote to 'strong' when there is an
    extra corroborating change read beyond the minimum and a well-anchored
    baseline."""
    if n_change_support_reads >= min_supporting_reads + 1 and n_baseline_reads >= min_supporting_reads:
        return "strong"
    return "supportive"


def _chrom_rank(chrom: str) -> int:
    return _CHROM_ORDER.get(str(chrom), 99)


def _scan_loci(
    a_loci: dict[LocusKey, list[tuple[str, int]]],
    b_loci: dict[LocusKey, list[tuple[str, int]]],
    *,
    sample_a_label: str,
    sample_b_label: str,
    margin_bp: int,
    min_supporting_reads: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for key in set(a_loci) | set(b_loci):
        a_pairs = sorted(a_loci.get(key, []), key=lambda item: item[1], reverse=True)
        b_pairs = sorted(b_loci.get(key, []), key=lambda item: item[1], reverse=True)
        a_lengths = [length for _, length in a_pairs]
        b_lengths = [length for _, length in b_pairs]

        excess_a = _direction_excess(
            a_lengths, b_lengths, margin_bp=margin_bp, min_supporting_reads=min_supporting_reads
        )
        excess_b = _direction_excess(
            b_lengths, a_lengths, margin_bp=margin_bp, min_supporting_reads=min_supporting_reads
        )
        if excess_a is None and excess_b is None:
            continue
        if excess_a is not None and excess_b is not None:
            # Both groups carry reads longer than the other by the margin — an
            # ambiguous, bidirectional signal. Skip unless one clearly dominates.
            if excess_a == excess_b:
                continue
            change_is_a = excess_a > excess_b
        else:
            change_is_a = excess_a is not None

        if change_is_a:
            change_pairs, baseline_pairs = a_pairs, b_pairs
            change_group, baseline_group = sample_a_label, sample_b_label
        else:
            change_pairs, baseline_pairs = b_pairs, a_pairs
            change_group, baseline_group = sample_b_label, sample_a_label

        chrom, start, end = key
        anchor = baseline_pairs[0][1]
        threshold = anchor + margin_bp
        change_support = [(name, length) for name, length in change_pairs if length > threshold]
        change_lengths = [length for _, length in change_pairs]
        baseline_lengths = [length for _, length in baseline_pairs]
        support_lengths = [length for _, length in change_support]
        change_max = change_lengths[0]

        tier = _assign_tier(
            n_change_support_reads=len(change_support),
            n_baseline_reads=len(baseline_pairs),
            min_supporting_reads=min_supporting_reads,
        )
        rows.append(
            {
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "trid": f"{chrom}_{start}_{end}",
                "change_allele": _CHANGE_ALLELE,
                "change_type": _CHANGE_TYPE,
                "change_group": change_group,
                "baseline_group": baseline_group,
                "n_change_reads": len(change_pairs),
                "n_baseline_reads": len(baseline_pairs),
                "n_change_support_reads": len(change_support),
                "n_baseline_support_reads": len(baseline_pairs),
                "change_max_bp": int(change_max),
                "baseline_max_bp": int(anchor),
                "baseline_anchor_bp": int(threshold),
                "margin_bp": int(margin_bp),
                "min_supporting_reads": int(min_supporting_reads),
                "change_length_bp": int(change_max - anchor),
                "change_support_min_excess_bp": int(min(support_lengths) - anchor) if support_lengths else 0,
                "change_read_mean": round(statistics.fmean(change_lengths), 1) if change_lengths else ".",
                "baseline_read_mean": round(statistics.fmean(baseline_lengths), 1) if baseline_lengths else ".",
                "change_read_range": f"{min(change_lengths)}-{max(change_lengths)}" if change_lengths else ".",
                "baseline_read_range": f"{min(baseline_lengths)}-{max(baseline_lengths)}" if baseline_lengths else ".",
                "change_read_names": [name for name, _ in change_pairs],
                "baseline_read_names": [name for name, _ in baseline_pairs],
                "change_support_read_names": [name for name, _ in change_support],
                "baseline_support_read_names": [name for name, _ in baseline_pairs],
                "signal_class": _READLEN_SIGNAL_CLASS,
                "tr_tier": tier,
                "tr_pass_for_harmonized": True,
            }
        )

    rows.sort(
        key=lambda row: (
            _TR_TIER_ORDER.get(row["tr_tier"], 99),
            -int(row["change_length_bp"]),
            _chrom_rank(row["chrom"]),
            int(row["start"]),
        )
    )
    return rows


def _read_length_rows(
    rows: list[dict[str, Any]],
    a_loci: dict[LocusKey, list[tuple[str, int]]],
    b_loci: dict[LocusKey, list[tuple[str, int]]],
    *,
    sample_a_label: str,
    sample_b_label: str,
) -> list[dict[str, Any]]:
    """Per-read length records for every called locus (both groups)."""
    out: list[dict[str, Any]] = []
    for row in rows:
        key: LocusKey = (row["chrom"], int(row["start"]), int(row["end"]))
        for label, loci in ((sample_a_label, a_loci), (sample_b_label, b_loci)):
            for name, length in loci.get(key, []):
                out.append(
                    {
                        "trid": row["trid"],
                        "chrom": row["chrom"],
                        "start": int(row["start"]),
                        "end": int(row["end"]),
                        "cell_type": label,
                        "read_name": name,
                        "read_length": int(length),
                    }
                )
    return out


def _serialize_bed_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Render list-valued read-name columns as JSON-ish strings for the TSV."""
    import json

    json_cols = {
        "change_read_names",
        "baseline_read_names",
        "change_support_read_names",
        "baseline_support_read_names",
    }
    serialized: list[dict[str, Any]] = []
    for row in rows:
        new_row = dict(row)
        for col in json_cols:
            new_row[col] = json.dumps(new_row.get(col, []))
        serialized.append(new_row)
    return serialized


def _plot_overview(rows: list[dict[str, Any]], read_rows: list[dict[str, Any]], outpath: Path) -> bool:
    if not rows or not read_rows:
        return False
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # pragma: no cover - matplotlib optional
        return False

    order = {row["trid"]: idx for idx, row in enumerate(rows)}
    by_group: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for rec in read_rows:
        idx = order.get(rec["trid"])
        if idx is None:
            continue
        by_group[rec["cell_type"]].append((idx, int(rec["read_length"])))

    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig_width = max(6.0, 0.4 * len(rows) + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, 5.0))
    for group_label, points in sorted(by_group.items()):
        xs = [x for x, _ in points]
        ys = [y for _, y in points]
        ax.scatter(xs, ys, s=14, alpha=0.6, label=group_label)
    ax.set_xticks(range(len(rows)))
    ax.set_xticklabels([row["trid"] for row in rows], rotation=90, fontsize=6)
    ax.set_ylabel("trimmed read length (bp)")
    ax.set_xlabel("called TR locus")
    ax.set_title("Read-length paired scan: called TR loci")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    try:
        fig.savefig(outpath, dpi=120)
    finally:
        plt.close(fig)
    return outpath.exists()


def tr_post_processing_main(cli_args=None) -> dict[str, Any]:
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(message)s")
    parser = _build_arg_parser()
    args = _resolve_args(parser.parse_args(cli_args))
    outputs = _build_output_paths(args.output_dir)

    missing_inputs = [
        str(path)
        for path in (args.split_dir, args.group_a_fasta, args.group_b_fasta)
        if not path.exists()
    ]
    if missing_inputs:
        return _write_empty_outputs(
            args, status="skipped", reason="Missing required input(s): " + ", ".join(missing_inputs)
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    a_loci = _parse_fasta_lengths(args.group_a_fasta)
    b_loci = _parse_fasta_lengths(args.group_b_fasta)
    n_loci_scanned = len(set(a_loci) | set(b_loci))

    rows = _scan_loci(
        a_loci,
        b_loci,
        sample_a_label=args.sample_a_label,
        sample_b_label=args.sample_b_label,
        margin_bp=args.margin_bp,
        min_supporting_reads=args.min_supporting_reads,
    )

    read_rows = _read_length_rows(
        rows,
        a_loci,
        b_loci,
        sample_a_label=args.sample_a_label,
        sample_b_label=args.sample_b_label,
    )
    _write_tsv(
        outputs.read_lengths_tsv,
        read_rows,
        ["trid", "chrom", "start", "end", "cell_type", "read_name", "read_length"],
    )

    bed_rows = _serialize_bed_rows(rows)
    _write_tsv(outputs.targets_tsv, bed_rows, _TR_BED_COLS)
    _write_tsv(outputs.tr_bed_tsv, bed_rows, _TR_BED_COLS)

    plots: list[str] = []
    if args.make_plots and _plot_overview(rows, read_rows, outputs.overview_plot):
        plots.append(str(outputs.overview_plot))

    n_strong = sum(1 for row in rows if row["tr_tier"] == "strong")
    n_supportive = sum(1 for row in rows if row["tr_tier"] == "supportive")
    n_weak = sum(1 for row in rows if row["tr_tier"] == "weak")
    summary = _build_summary_payload(
        args,
        outputs,
        status="completed" if rows else "completed_empty",
        reason=None if rows else "No loci passed the paired read-length comparison",
        n_loci_scanned=n_loci_scanned,
        n_targets=len(rows),
        n_tr_strong_rows=n_strong,
        n_tr_supportive_rows=n_supportive,
        n_tr_weak_rows=n_weak,
        plots=plots,
    )
    _write_json(outputs.summary_json, summary)
    return summary


def main(argv: list[str] | None = None) -> int:
    tr_post_processing_main(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
