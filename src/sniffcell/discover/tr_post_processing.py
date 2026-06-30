"""Read-length paired-comparison TR caller for the discover pipeline.

This module implements the "somatic scan" TR-discovery method: instead of
modelling tandem-repeat allele lengths through a merged TDB, it works directly
off the per-group spanning reads emitted by the TR genotyper. For medaka tandem
this is ``trimmed_reads.fasta``. For TRGT this is the spanning BAM; the read
length proxy is flank-adjusted using TRGT's ``FL`` tag when available.

For every locus that has spanning reads in either of the two split groups we
pool the read lengths per group (haplotype-agnostic) and apply a paired,
directional comparison:

    A locus is a confident TR change in favour of one group ("change group")
    when both groups have at least ``min_total_reads`` spanning reads and the
    change group's top ``min_supporting_reads`` reads are EACH longer than the
    other group's longest read by at least ``margin_bp``.

The longest read of the comparator group is the "anchor"; the comparator must
contribute at least one read (loci where the baseline group has no spanning
reads are skipped to avoid coverage-driven false positives). The method only
calls expansions of the change group relative to the baseline group, which is
exactly the signal the genome-wide somatic scan was built around.

``margin_bp``, ``min_supporting_reads``, and ``min_total_reads`` are the tunable
parameters.

A locus is additionally dropped when the repeat unit driving the call is shorter
than ``min_motif_size`` bp. The motif is read directly off the change group's
longest spanning read, so homopolymer (1 bp) tracts -- where ONT length
estimation slips badly and inflates a handful of reads into spurious multi-kb
"expansions" -- are excluded by default (``min_motif_size=2``). Set
``--min-motif-size 3`` to also drop dinucleotide (e.g. AT/TA) tracts, or
``--min-motif-size 1`` to recover the previous, unfiltered behaviour.

Called loci are additionally flagged ``hap_dropout_low_conf`` (and downgraded to
the ``weak`` tier with ``tr_pass_for_harmonized=False``) when the expansion is
confined to a haplotype that the baseline group has no spanning reads from -- for
example every supporting read is ``hap2`` but the baseline group contributes no
``hap2`` reads at the locus. The haplotype is read straight off the medaka
trimmed-read name (``..._hap1_phased-set...``); such calls cannot be told apart
from haplotype-coverage dropout in the baseline group, so they are marked low
confidence rather than dropped. Read names without a phased ``hap1``/``hap2``
token (e.g. ``hap0``) never trigger the flag, and TRGT spanning-BAM inputs --
whose read names lack the token -- are left unflagged.

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
    _discover_groups,
    _expand_path,
    _infer_sample_id,
    _resolve_two_group_names,
    _sanitize_token,
    _write_json,
)


DEFAULT_MARGIN_BP = 50
DEFAULT_MIN_SUPPORTING_READS = 3
DEFAULT_MIN_TOTAL_READS = 5
DEFAULT_TRGT_FALLBACK_FLANK_BP = 50
DEFAULT_MIN_MOTIF_SIZE = 2

# Motif detection on the change group's longest read: the smallest repeat unit
# (1..``_MOTIF_MAX_PERIOD`` bp) is accepted only when at least
# ``_MOTIF_MIN_MATCH_FRAC`` of bases match the base ``period`` positions ahead,
# evaluated over the most repetitive ``_MOTIF_PROBE_WINDOW_BP`` window so flanks
# do not dilute the signal. Below the match floor the motif is "undetermined"
# and the locus is kept (the filter only ever drops *confident* short motifs).
_MOTIF_MAX_PERIOD = 6
_MOTIF_MIN_MATCH_FRAC = 0.7
_MOTIF_PROBE_WINDOW_BP = 200

_READLEN_SIGNAL_CLASS = "readlen_paired"
_CHANGE_ALLELE = "all"  # the scan pools reads across haplotypes
_CHANGE_TYPE = "expansion"

# First ``chrN_start_end`` triple anywhere in a trimmed-read FASTA header.
_REGION_RE = re.compile(r"(chr[0-9A-Za-z]+)_(\d+)_(\d+)")

# Phased haplotype token in a medaka trimmed-read name, e.g. ``..._hap2_phased-set...``.
# Only ``hap1``/``hap2`` are treated as phased; ``hap0`` (unphased) is ignored.
_HAP_RE = re.compile(r"_hap([12])_")

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
    "margin_bp", "min_supporting_reads", "min_total_reads",
    "change_length_bp", "change_support_min_excess_bp",
    "change_read_mean", "baseline_read_mean",
    "change_read_range", "baseline_read_range",
    "change_read_names", "baseline_read_names",
    "change_support_read_names", "baseline_support_read_names",
    "change_support_haps", "baseline_haps", "hap_dropout_low_conf",
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
    group_a_spanning_bam: Path | None
    group_b_spanning_bam: Path | None
    tr_bed: Path | None
    trgt_fallback_flank_bp: int
    margin_bp: int
    min_supporting_reads: int
    min_total_reads: int
    min_motif_size: int
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
            "both groups have at least --min-total-reads reads and one group's top "
            "--min-supporting-reads reads each exceed the other group's longest read "
            "by at least --margin-bp."
        ),
        add_help=add_help,
    )
    parser.add_argument("--split-dir", required=True, help="deconv_requested_group_splits directory")
    parser.add_argument(
        "--groups",
        default=None,
        help="Exactly two group names, comma-separated. Default: infer from a two-group split manifest.",
    )
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
        "--discover-run-id",
        default=None,
        help="Optional discover run ID to search for medaka_tandem inputs when group FASTAs are omitted.",
    )
    parser.add_argument(
        "--group-a-spanning-bam",
        default=None,
        help="trgt spanning BAM for the first group. Used when --group-a-fasta is unavailable.",
    )
    parser.add_argument(
        "--group-b-spanning-bam",
        default=None,
        help="trgt spanning BAM for the second group. Used when --group-b-fasta is unavailable.",
    )
    parser.add_argument(
        "--tr-bed",
        default=None,
        help="TR BED used by trgt. Required when using --group-a/--group-b-spanning-bam.",
    )
    parser.add_argument(
        "--trgt-fallback-flank-bp",
        type=int,
        default=DEFAULT_TRGT_FALLBACK_FLANK_BP,
        help=(
            "Fallback per-side TRGT flank length to subtract when a spanning BAM read lacks "
            f"the FL tag. Default={DEFAULT_TRGT_FALLBACK_FLANK_BP}."
        ),
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
        "--min-total-reads",
        type=int,
        default=DEFAULT_MIN_TOTAL_READS,
        help=(
            "Minimum total spanning reads required in each group before testing a locus. "
            f"Default={DEFAULT_MIN_TOTAL_READS}."
        ),
    )
    parser.add_argument(
        "--min-motif-size",
        type=int,
        default=DEFAULT_MIN_MOTIF_SIZE,
        help=(
            "Drop a called locus when the repeat unit of the change group's longest read is "
            "shorter than this many bp. Excludes homopolymer/dinucleotide tracts where ONT "
            "length estimation slips and inflates a few reads into spurious multi-kb expansions. "
            f"Use 1 to disable. Default={DEFAULT_MIN_MOTIF_SIZE}."
        ),
    )
    parser.add_argument(
        "--merged-tdb",
        default=None,
        help="Deprecated and ignored. The read-length method does not use a TDB; kept for backward compatibility.",
    )
    parser.add_argument("--skip-plots", action="store_true", default=False, help="Skip PNG generation")
    return parser


def _group_path_tokens(split_dir: Path, group_name: str) -> list[str]:
    tokens = [group_name, _sanitize_token(group_name)]
    try:
        for group in _discover_groups(split_dir):
            if group.name == group_name:
                tokens.append(Path(group.bam_path).stem)
                break
    except (FileNotFoundError, ValueError):
        pass
    return list(dict.fromkeys(token for token in tokens if token))


def _trimmed_reads_path(base_dir: Path, group_token: str) -> Path:
    return base_dir / "medaka_tandem" / f"{group_token}.medaka" / "trimmed_reads.fasta"


def _discover_run_dirs(split_dir: Path, discover_run_id: str | None) -> list[Path]:
    discover_root = split_dir / "discover"
    if discover_run_id:
        return [discover_root / discover_run_id]
    if not discover_root.exists():
        return []
    run_dirs = [path for path in discover_root.iterdir() if path.is_dir()]
    return sorted(run_dirs, key=lambda path: path.stat().st_mtime, reverse=True)


def _trimmed_read_pairs(
    base_dir: Path,
    group_a_tokens: list[str],
    group_b_tokens: list[str],
) -> list[tuple[Path, Path]]:
    pairs: list[tuple[Path, Path]] = []
    for group_a_token in group_a_tokens:
        for group_b_token in group_b_tokens:
            pairs.append(
                (
                    _trimmed_reads_path(base_dir, group_a_token),
                    _trimmed_reads_path(base_dir, group_b_token),
                )
            )
    return pairs


def _resolve_default_trimmed_read_pair(
    *,
    split_dir: Path,
    group_a: str,
    group_b: str,
    discover_run_id: str | None = None,
) -> tuple[Path, Path]:
    group_a_tokens = _group_path_tokens(split_dir, group_a)
    group_b_tokens = _group_path_tokens(split_dir, group_b)
    direct_pairs = _trimmed_read_pairs(split_dir, group_a_tokens, group_b_tokens)
    for group_a_fasta, group_b_fasta in direct_pairs:
        if group_a_fasta.exists() and group_b_fasta.exists():
            return group_a_fasta, group_b_fasta

    discover_pairs: list[tuple[Path, Path]] = []
    for run_dir in _discover_run_dirs(split_dir, discover_run_id):
        run_pairs = _trimmed_read_pairs(run_dir, group_a_tokens, group_b_tokens)
        for group_a_fasta, group_b_fasta in run_pairs:
            if group_a_fasta.exists() and group_b_fasta.exists():
                return group_a_fasta, group_b_fasta
        discover_pairs.extend(run_pairs)

    if discover_run_id and discover_pairs:
        return discover_pairs[0]
    return direct_pairs[0]


def _resolve_args(raw_args) -> TrPostArgs:
    split_dir = _expand_path(raw_args.split_dir)
    group_a, group_b = _resolve_two_group_names(split_dir, raw_args.groups)
    sample_id = raw_args.sample_id or _infer_sample_id(split_dir.parent)
    output_dir = (
        _expand_path(raw_args.output_dir)
        if raw_args.output_dir
        else split_dir / "postprocess" / f"tr_post_processing_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
    sample_a_label = raw_args.sample_a_label or f"{sample_id}.{group_a}"
    sample_b_label = raw_args.sample_b_label or f"{sample_id}.{group_b}"
    if raw_args.group_a_fasta or raw_args.group_b_fasta:
        group_a_fasta = (
            _expand_path(raw_args.group_a_fasta)
            if raw_args.group_a_fasta
            else _trimmed_reads_path(split_dir, group_a)
        )
        group_b_fasta = (
            _expand_path(raw_args.group_b_fasta)
            if raw_args.group_b_fasta
            else _trimmed_reads_path(split_dir, group_b)
        )
    else:
        group_a_fasta, group_b_fasta = _resolve_default_trimmed_read_pair(
            split_dir=split_dir,
            group_a=group_a,
            group_b=group_b,
            discover_run_id=getattr(raw_args, "discover_run_id", None),
        )
    group_a_spanning_bam = _expand_path(raw_args.group_a_spanning_bam) if raw_args.group_a_spanning_bam else None
    group_b_spanning_bam = _expand_path(raw_args.group_b_spanning_bam) if raw_args.group_b_spanning_bam else None
    tr_bed = _expand_path(raw_args.tr_bed) if raw_args.tr_bed else None
    trgt_fallback_flank_bp = int(raw_args.trgt_fallback_flank_bp)
    margin_bp = int(raw_args.margin_bp)
    min_supporting_reads = int(raw_args.min_supporting_reads)
    min_total_reads = int(raw_args.min_total_reads)
    min_motif_size = int(raw_args.min_motif_size)
    if trgt_fallback_flank_bp < 0:
        raise ValueError("--trgt-fallback-flank-bp must be >= 0")
    if margin_bp < 0:
        raise ValueError("--margin-bp must be >= 0")
    if min_supporting_reads < 1:
        raise ValueError("--min-supporting-reads must be >= 1")
    if min_total_reads < min_supporting_reads:
        raise ValueError("--min-total-reads must be >= --min-supporting-reads")
    if min_motif_size < 1:
        raise ValueError("--min-motif-size must be >= 1")
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
        group_a_spanning_bam=group_a_spanning_bam,
        group_b_spanning_bam=group_b_spanning_bam,
        tr_bed=tr_bed,
        trgt_fallback_flank_bp=trgt_fallback_flank_bp,
        margin_bp=margin_bp,
        min_supporting_reads=min_supporting_reads,
        min_total_reads=min_total_reads,
        min_motif_size=min_motif_size,
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
    n_tr_motif_filtered: int = 0,
    n_tr_hap_dropout: int = 0,
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
        "group_a_spanning_bam": str(args.group_a_spanning_bam) if args.group_a_spanning_bam else "",
        "group_b_spanning_bam": str(args.group_b_spanning_bam) if args.group_b_spanning_bam else "",
        "tr_bed": str(args.tr_bed) if args.tr_bed else "",
        "n_loci_scanned": int(n_loci_scanned),
        "n_targets": int(n_targets),
        "n_tr_strong_rows": int(n_tr_strong_rows),
        "n_tr_supportive_rows": int(n_tr_supportive_rows),
        "n_tr_weak_rows": int(n_tr_weak_rows),
        "n_tr_motif_filtered": int(n_tr_motif_filtered),
        "n_tr_hap_dropout": int(n_tr_hap_dropout),
        "read_lengths_tsv": str(outputs.read_lengths_tsv),
        "targets_tsv": str(outputs.targets_tsv),
        "tr_bed_tsv": str(outputs.tr_bed_tsv),
        "plots": plots or [],
        "params": {
            "margin_bp": args.margin_bp,
            "min_supporting_reads": args.min_supporting_reads,
            "min_total_reads": args.min_total_reads,
            "min_motif_size": args.min_motif_size,
            "trgt_fallback_flank_bp": args.trgt_fallback_flank_bp,
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


def _dominant_motif_period(seq: str) -> tuple[int, float]:
    """Smallest repeat-unit length (bp) of ``seq`` and its support fraction.

    For each candidate period ``p`` in ``1.._MOTIF_MAX_PERIOD`` we measure the
    fraction of bases equal to the base ``p`` positions ahead, evaluated over the
    most periodic ``_MOTIF_PROBE_WINDOW_BP`` window so long flanks cannot mask a
    short repeat. Returns ``(period, match_fraction)`` for the best period, or
    ``(0, 0.0)`` when the sequence is too short to assess.
    """
    seq = seq.upper()
    n = len(seq)
    if n < 2:
        return (0, 0.0)
    window = min(_MOTIF_PROBE_WINDOW_BP, n)
    step = max(1, window // 2)
    starts = list(range(0, n - window + 1, step)) or [0]
    best_period = 0
    best_frac = 0.0
    for start in starts:
        sub = seq[start:start + window]
        for period in range(1, _MOTIF_MAX_PERIOD + 1):
            denom = len(sub) - period
            if denom <= 0:
                continue
            matches = sum(1 for i in range(denom) if sub[i] == sub[i + period])
            frac = matches / denom
            if frac > best_frac:
                best_frac = frac
                best_period = period
    return (best_period, best_frac)


def _parse_fasta_longest_seqs(path: Path, wanted: set[LocusKey]) -> dict[LocusKey, str]:
    """Longest spanning-read sequence per locus, restricted to ``wanted`` keys.

    Used after the scan to read the motif off the change group's longest read
    only for loci that were actually called, so the extra FASTA pass stays cheap.
    """
    if not wanted:
        return {}
    longest: dict[LocusKey, tuple[int, str]] = {}
    key: LocusKey | None = None
    chunks: list[str] = []

    def _flush() -> None:
        nonlocal key, chunks
        if key is not None and key in wanted and chunks:
            seq = "".join(chunks)
            previous = longest.get(key)
            if previous is None or len(seq) > previous[0]:
                longest[key] = (len(seq), seq)

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                _flush()
                header = line[1:]
                chunks = []
                key = None
                match = _REGION_RE.search(header)
                if match:
                    candidate = (match.group(1), int(match.group(2)), int(match.group(3)))
                    if candidate in wanted:
                        key = candidate
            elif key is not None:
                chunks.append(line.strip())
        _flush()
    return {key: value[1] for key, value in longest.items()}


def _apply_motif_size_filter(
    rows: list[dict[str, Any]],
    *,
    sample_a_label: str,
    sample_b_label: str,
    group_a_fasta: Path,
    group_b_fasta: Path,
    group_a_source: str,
    group_b_source: str,
    min_motif_size: int,
) -> tuple[list[dict[str, Any]], int]:
    """Drop called loci whose change-group repeat unit is < ``min_motif_size`` bp.

    The motif is taken from the change group's longest spanning read. Loci whose
    motif cannot be confidently determined (``match_fraction`` below
    ``_MOTIF_MIN_MATCH_FRAC``) are kept -- the filter only removes confident short
    motifs. Returns ``(kept_rows, n_dropped)``. Sequences are only available for
    FASTA (medaka) inputs; TRGT spanning-BAM groups are left unfiltered.
    """
    if min_motif_size <= 1 or not rows:
        return rows, 0

    wanted_a = {
        (row["chrom"], int(row["start"]), int(row["end"]))
        for row in rows
        if row["change_group"] == sample_a_label
    }
    wanted_b = {
        (row["chrom"], int(row["start"]), int(row["end"]))
        for row in rows
        if row["change_group"] == sample_b_label
    }
    seqs_a = (
        _parse_fasta_longest_seqs(group_a_fasta, wanted_a)
        if group_a_source == "fasta"
        else {}
    )
    seqs_b = (
        _parse_fasta_longest_seqs(group_b_fasta, wanted_b)
        if group_b_source == "fasta"
        else {}
    )

    kept: list[dict[str, Any]] = []
    dropped = 0
    for row in rows:
        key: LocusKey = (row["chrom"], int(row["start"]), int(row["end"]))
        seq = seqs_a.get(key) if row["change_group"] == sample_a_label else seqs_b.get(key)
        if seq:
            period, frac = _dominant_motif_period(seq)
            if period > 0 and frac >= _MOTIF_MIN_MATCH_FRAC and period < min_motif_size:
                dropped += 1
                continue
        kept.append(row)
    return kept, dropped


def _phased_haps_from_read_names(names: list[str]) -> set[str]:
    """Phased haplotype labels ({"1", "2"}) carried by ``names``.

    The label is read off the medaka trimmed-read name (``..._hap1_phased-set...``).
    Unphased reads (``hap0``) and names without the token contribute nothing.
    """
    haps: set[str] = set()
    for name in names:
        match = _HAP_RE.search(name)
        if match:
            haps.add(match.group(1))
    return haps


def _apply_haplotype_dropout_filter(rows: list[dict[str, Any]]) -> int:
    """Flag calls whose expansion sits on a haplotype the baseline group lacks.

    For each called locus the change group's *supporting* (expanded) read names
    and the baseline group's full read names are parsed for phased ``hap1``/
    ``hap2`` labels. When every supporting read is phased and none of those
    haplotypes appear among the baseline reads (e.g. all supporting reads are
    ``hap2`` but the baseline group has no ``hap2`` reads at the locus), the
    apparent expansion cannot be distinguished from haplotype-coverage dropout in
    the baseline group, so the row is marked low confidence: ``tr_tier`` is set to
    ``weak`` and ``tr_pass_for_harmonized`` to ``False``.

    Records ``change_support_haps``, ``baseline_haps`` and ``hap_dropout_low_conf``
    on every row (set in-place). Returns the number of flagged rows. Loci whose
    supporting reads are unphased -- including all TRGT spanning-BAM inputs, whose
    read names lack the token -- are never flagged.
    """
    flagged = 0
    for row in rows:
        support_haps = _phased_haps_from_read_names(row.get("change_support_read_names", []))
        baseline_haps = _phased_haps_from_read_names(row.get("baseline_read_names", []))
        dropout = bool(support_haps) and support_haps.isdisjoint(baseline_haps)
        row["change_support_haps"] = ",".join(sorted(support_haps)) if support_haps else "."
        row["baseline_haps"] = ",".join(sorted(baseline_haps)) if baseline_haps else "."
        row["hap_dropout_low_conf"] = dropout
        if dropout:
            row["tr_tier"] = "weak"
            row["tr_pass_for_harmonized"] = False
            flagged += 1
    return flagged


def _parse_tr_bed_loci(path: Path) -> dict[str, LocusKey]:
    """Map TRGT BED IDs to ``(chrom, start, end)`` locus keys."""
    loci: dict[str, LocusKey] = {}
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            trid = ""
            for token in fields[3].split(";"):
                if token.startswith("ID="):
                    trid = token[3:]
                    break
            if not trid:
                trid = fields[3]
            loci[trid] = (chrom, start, end)
    return loci


def _parse_trgt_flanks(fl_tag_value: Any, fallback_flank_bp: int) -> tuple[int, int]:
    """Return TRGT left/right flank lengths from an FL tag.

    TRGT writes ``FL:B:I,left,right`` in spanning BAMs. The string branch keeps
    this compatible with SAM text values and mirrors mTRplotter's parser.
    """
    if fl_tag_value is not None:
        if isinstance(fl_tag_value, str):
            parts = fl_tag_value.split(",")
            value_parts = parts[1:] if parts and parts[0].isalpha() else parts
        else:
            try:
                value_parts = list(fl_tag_value)
            except TypeError:
                value_parts = []
        if len(value_parts) >= 2:
            try:
                return int(value_parts[0]), int(value_parts[1])
            except (TypeError, ValueError):
                pass
    return fallback_flank_bp, fallback_flank_bp


def _parse_trgt_spanning_bam_lengths(
    path: Path,
    tr_bed: Path,
    *,
    fallback_flank_bp: int,
) -> dict[LocusKey, list[tuple[str, int]]]:
    """Map each TRGT locus to unique ``(read_name, read_length)`` records.

    TRGT's spanning BAM stores the repeat ID in the ``TR`` tag and the
    repeat-spanning sequence in the BAM query sequence. Following mTRplotter,
    the per-read allele-length proxy is ``len(query_sequence) - FL_left -
    FL_right`` where TRGT's ``FL`` tag is present, with a configurable fallback
    flank length for older or malformed records.
    """
    import pysam

    tr_loci = _parse_tr_bed_loci(tr_bed)
    by_locus: dict[LocusKey, dict[str, int]] = defaultdict(dict)
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            seq = read.query_sequence
            if not seq:
                continue
            try:
                trid = read.get_tag("TR")
            except KeyError:
                continue
            locus = tr_loci.get(str(trid))
            if locus is None:
                continue
            read_name = str(read.query_name)
            try:
                fl_tag = read.get_tag("FL")
            except KeyError:
                fl_tag = None
            flank_left, flank_right = _parse_trgt_flanks(fl_tag, fallback_flank_bp)
            read_length = len(seq) - flank_left - flank_right
            previous = by_locus[locus].get(read_name)
            if previous is None or read_length > previous:
                by_locus[locus][read_name] = read_length
    return {locus: sorted(reads.items()) for locus, reads in by_locus.items()}


def _load_group_loci(
    *,
    fasta: Path,
    spanning_bam: Path | None,
    tr_bed: Path | None,
    trgt_fallback_flank_bp: int,
) -> tuple[dict[LocusKey, list[tuple[str, int]]] | None, str, str | None]:
    if fasta.exists():
        return _parse_fasta_lengths(fasta), "fasta", None
    if spanning_bam is not None and spanning_bam.exists():
        if tr_bed is None or not tr_bed.exists():
            return None, "trgt_spanning_bam", f"Missing TR BED for spanning BAM input: {tr_bed or ''}"
        return (
            _parse_trgt_spanning_bam_lengths(
                spanning_bam,
                tr_bed,
                fallback_flank_bp=trgt_fallback_flank_bp,
            ),
            "trgt_spanning_bam",
            None,
        )
    candidates = [str(fasta)]
    if spanning_bam is not None:
        candidates.append(str(spanning_bam))
    return None, "missing", "Missing required input(s): " + ", ".join(candidates)


# ---------------------------------------------------------------------------
# Paired read-length scan
# ---------------------------------------------------------------------------

def _direction_excess(
    change_lengths_desc: list[int],
    baseline_lengths_desc: list[int],
    *,
    margin_bp: int,
    min_supporting_reads: int,
    min_total_reads: int,
) -> int | None:
    """Excess (bp) of the change group's longest read over the baseline anchor,
    or ``None`` if this direction does not pass the paired-comparison rule.

    Passes when both groups contribute at least ``min_total_reads`` reads and
    each of the change group's top ``min_supporting_reads`` reads exceeds
    ``baseline_max + margin_bp``.
    """
    if len(baseline_lengths_desc) < min_total_reads:
        return None
    if len(change_lengths_desc) < min_total_reads:
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
    min_total_reads: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for key in set(a_loci) | set(b_loci):
        a_pairs = sorted(a_loci.get(key, []), key=lambda item: item[1], reverse=True)
        b_pairs = sorted(b_loci.get(key, []), key=lambda item: item[1], reverse=True)
        a_lengths = [length for _, length in a_pairs]
        b_lengths = [length for _, length in b_pairs]

        excess_a = _direction_excess(
            a_lengths,
            b_lengths,
            margin_bp=margin_bp,
            min_supporting_reads=min_supporting_reads,
            min_total_reads=min_total_reads,
        )
        excess_b = _direction_excess(
            b_lengths,
            a_lengths,
            margin_bp=margin_bp,
            min_supporting_reads=min_supporting_reads,
            min_total_reads=min_total_reads,
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
                "min_total_reads": int(min_total_reads),
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

    missing_inputs = [str(args.split_dir)] if not args.split_dir.exists() else []
    a_loci: dict[LocusKey, list[tuple[str, int]]] | None = None
    b_loci: dict[LocusKey, list[tuple[str, int]]] | None = None
    group_a_source = "missing"
    group_b_source = "missing"
    if not missing_inputs:
        a_loci, group_a_source, a_error = _load_group_loci(
            fasta=args.group_a_fasta,
            spanning_bam=args.group_a_spanning_bam,
            tr_bed=args.tr_bed,
            trgt_fallback_flank_bp=args.trgt_fallback_flank_bp,
        )
        b_loci, group_b_source, b_error = _load_group_loci(
            fasta=args.group_b_fasta,
            spanning_bam=args.group_b_spanning_bam,
            tr_bed=args.tr_bed,
            trgt_fallback_flank_bp=args.trgt_fallback_flank_bp,
        )
        missing_inputs = [err for err in (a_error, b_error) if err]
    if missing_inputs or a_loci is None or b_loci is None:
        return _write_empty_outputs(
            args, status="skipped", reason="; ".join(missing_inputs)
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    n_loci_scanned = len(set(a_loci) | set(b_loci))

    rows = _scan_loci(
        a_loci,
        b_loci,
        sample_a_label=args.sample_a_label,
        sample_b_label=args.sample_b_label,
        margin_bp=args.margin_bp,
        min_supporting_reads=args.min_supporting_reads,
        min_total_reads=args.min_total_reads,
    )

    rows, n_motif_filtered = _apply_motif_size_filter(
        rows,
        sample_a_label=args.sample_a_label,
        sample_b_label=args.sample_b_label,
        group_a_fasta=args.group_a_fasta,
        group_b_fasta=args.group_b_fasta,
        group_a_source=group_a_source,
        group_b_source=group_b_source,
        min_motif_size=args.min_motif_size,
    )
    if n_motif_filtered:
        logging.info(
            "Dropped %d locus call(s) with repeat unit < %d bp (--min-motif-size)",
            n_motif_filtered,
            args.min_motif_size,
        )

    n_hap_dropout = _apply_haplotype_dropout_filter(rows)
    if n_hap_dropout:
        logging.info(
            "Flagged %d locus call(s) as low confidence: expansion on a haplotype "
            "absent from the baseline group (hap_dropout_low_conf)",
            n_hap_dropout,
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
        n_tr_motif_filtered=n_motif_filtered,
        n_tr_hap_dropout=n_hap_dropout,
        plots=plots,
    )
    summary["group_a_source"] = group_a_source
    summary["group_b_source"] = group_b_source
    _write_json(outputs.summary_json, summary)
    return summary


def main(argv: list[str] | None = None) -> int:
    tr_post_processing_main(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
