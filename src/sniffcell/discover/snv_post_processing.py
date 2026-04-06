from __future__ import annotations

import argparse
import csv
import gzip
import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterator, TextIO

from sniffcell.discover.discover import (
    _expand_path,
    _infer_sample_id,
    _sanitize_token,
    _write_json,
)


@dataclass(frozen=True)
class TwoSampleSnvArgs:
    split_dir: Path
    output_dir: Path
    group_a: str
    group_b: str
    group_a_gvcf: Path
    group_b_gvcf: Path
    sample_id: str
    min_dp: int
    max_dp: int
    min_dp_absence: int
    min_gq: int
    min_other_af: float


@dataclass(frozen=True)
class CallRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    filt: str
    gt: str
    gq: int
    dp: int
    ad: tuple[int, ...]
    af: str


@dataclass(frozen=True)
class PositionSummary:
    chrom: str
    pos: int
    records: tuple[CallRecord, ...]
    nonref_call: CallRecord | None
    nonref_snp_call: CallRecord | None


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _parse_int(value: str) -> int:
    if value in {"", "."}:
        return 0
    return int(float(value))


def _parse_float_str(value: str) -> str:
    if value in {"", "."}:
        return "."
    return value


def _parse_ad(value: str) -> tuple[int, ...]:
    if value in {"", "."}:
        return ()
    return tuple(_parse_int(part) for part in value.split(","))


def _is_ref_gt(gt: str) -> bool:
    return gt in {"0/0", "0|0"}


def _is_nonref_gt(gt: str) -> bool:
    return gt not in {"0/0", "0|0", "./.", ".|.", ".", "ref", "mis"}


def _is_snp(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and alt != "." and "," not in alt


def _record_rank(rec: CallRecord) -> tuple[int, int]:
    return (rec.gq, rec.dp)


def _alt_read_count(rec: CallRecord, alt: str) -> int:
    if rec.alt in {"", "."}:
        return 0
    alts = rec.alt.split(",")
    try:
        alt_index = alts.index(alt)
    except ValueError:
        return 0
    ad_index = alt_index + 1
    if len(rec.ad) <= ad_index:
        return 0
    return rec.ad[ad_index]


def _better_record(current: CallRecord | None, candidate: CallRecord) -> CallRecord:
    if current is None or _record_rank(candidate) > _record_rank(current):
        return candidate
    return current


def _parse_record(line: str) -> CallRecord | None:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 10:
        return None
    chrom, pos, _variant_id, ref, alt, _qual, filt, _info, fmt, sample = fields[:10]
    fmt_keys = fmt.split(":")
    sample_vals = sample.split(":")
    fmt_map = dict(zip(fmt_keys, sample_vals))
    return CallRecord(
        chrom=chrom,
        pos=int(pos),
        ref=ref,
        alt=alt,
        filt=filt,
        gt=fmt_map.get("GT", "."),
        gq=_parse_int(fmt_map.get("GQ", ".")),
        dp=_parse_int(fmt_map.get("DP", ".")),
        ad=_parse_ad(fmt_map.get("AD", ".")),
        af=_parse_float_str(fmt_map.get("AF", ".")),
    )


def _summarize_position(records: list[CallRecord], *, min_dp: int, max_dp: int, min_gq: int) -> PositionSummary:
    chrom = records[0].chrom
    pos = records[0].pos
    nonref_call: CallRecord | None = None
    nonref_snp_call: CallRecord | None = None
    for rec in records:
        if (
            _is_nonref_gt(rec.gt)
            and rec.filt == "PASS"
            and rec.dp >= min_dp
            and rec.dp <= max_dp
            and rec.gq >= min_gq
            and _alt_read_count(rec, rec.alt) > 0
        ):
            nonref_call = _better_record(nonref_call, rec)
            if _is_snp(rec.ref, rec.alt):
                nonref_snp_call = _better_record(nonref_snp_call, rec)
    return PositionSummary(
        chrom=chrom,
        pos=pos,
        records=tuple(records),
        nonref_call=nonref_call,
        nonref_snp_call=nonref_snp_call,
    )


def _iter_position_summaries(path: Path, *, min_dp: int, max_dp: int, min_gq: int) -> Iterator[PositionSummary]:
    with _open_text(path) as handle:
        current: list[CallRecord] = []
        current_key: tuple[str, int] | None = None
        for line in handle:
            if not line or line.startswith("#"):
                continue
            rec = _parse_record(line)
            if rec is None:
                continue
            key = (rec.chrom, rec.pos)
            if current_key is None:
                current_key = key
                current = [rec]
                continue
            if key == current_key:
                current.append(rec)
                continue
            yield _summarize_position(current, min_dp=min_dp, max_dp=max_dp, min_gq=min_gq)
            current_key = key
            current = [rec]
        if current:
            yield _summarize_position(current, min_dp=min_dp, max_dp=max_dp, min_gq=min_gq)


def _natural_contig_key(chrom: str) -> tuple[int, int, str]:
    if chrom.startswith("chr"):
        tail = chrom[3:]
        if tail.isdigit():
            return (0, int(tail), chrom)
        if tail == "X":
            return (0, 23, chrom)
        if tail == "Y":
            return (0, 24, chrom)
        if tail in {"M", "MT"}:
            return (0, 25, chrom)
    return (1, 0, chrom)


def _compare_positions(left: PositionSummary, right: PositionSummary) -> int:
    left_key = (_natural_contig_key(left.chrom), left.pos)
    right_key = (_natural_contig_key(right.chrom), right.pos)
    if left_key < right_key:
        return -1
    if left_key > right_key:
        return 1
    return 0


def _merged_headers() -> list[str]:
    return [
        "direction",
        "chrom",
        "pos",
        "ref",
        "alt",
        "change_group",
        "other_group",
        "target_gt",
        "target_gq",
        "target_dp",
        "target_alt_ad",
        "target_af",
        "other_gt",
        "other_gq",
        "other_dp",
        "other_alt_ad",
        "other_af",
    ]


def _row_from_match(
    *,
    direction: str,
    target_group: str,
    other_group: str,
    target_call: CallRecord,
    other_call: CallRecord,
) -> dict[str, Any]:
    return {
        "direction": direction,
        "chrom": target_call.chrom,
        "pos": target_call.pos,
        "ref": target_call.ref,
        "alt": target_call.alt,
        "change_group": target_group,
        "other_group": other_group,
        "target_gt": target_call.gt,
        "target_gq": target_call.gq,
        "target_dp": target_call.dp,
        "target_alt_ad": _alt_read_count(target_call, target_call.alt),
        "target_af": target_call.af,
        "other_gt": other_call.gt,
        "other_gq": other_call.gq,
        "other_dp": other_call.dp,
        "other_alt_ad": _alt_read_count(other_call, target_call.alt),
        "other_af": other_call.af,
    }


def _open_writer(path: Path) -> tuple[Any, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = path.open("w", newline="", encoding="utf-8")
    writer = csv.DictWriter(handle, fieldnames=_merged_headers(), delimiter="\t")
    writer.writeheader()
    return handle, writer


def _absence_supporting_record(summary: PositionSummary, alt: str, *, min_dp_absence: int, max_dp: int, min_other_af: float) -> CallRecord | None:
    best: CallRecord | None = None
    for rec in summary.records:
        if rec.dp < min_dp_absence or rec.dp > max_dp:
            continue
        if _alt_read_count(rec, alt) != 0:
            continue
        if rec.af != "." and float(rec.af) < min_other_af:
            continue
        best = _better_record(best, rec)
    return best


def _has_any_alt_support(summary: PositionSummary, alt: str) -> bool:
    for rec in summary.records:
        if _alt_read_count(rec, alt) > 0:
            return True
    return False


def compare_group_specific_snvs(
    *,
    group_a_label: str,
    group_b_label: str,
    group_a_gvcf: Path,
    group_b_gvcf: Path,
    output_dir: Path,
    min_dp: int,
    max_dp: int,
    min_dp_absence: int,
    min_gq: int,
    min_other_af: float,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    sets_dir = output_dir / "sets"
    sets_dir.mkdir(parents=True, exist_ok=True)
    merged_path = output_dir / "snv_changes.tsv"
    group_a_only_path = sets_dir / f"{_sanitize_token(group_a_label)}.only.tsv"
    group_b_only_path = sets_dir / f"{_sanitize_token(group_b_label)}.only.tsv"

    iter_a = iter(_iter_position_summaries(group_a_gvcf, min_dp=min_dp, max_dp=max_dp, min_gq=min_gq))
    iter_b = iter(_iter_position_summaries(group_b_gvcf, min_dp=min_dp, max_dp=max_dp, min_gq=min_gq))
    current_a = next(iter_a, None)
    current_b = next(iter_b, None)

    group_a_rows: list[dict[str, Any]] = []
    group_b_rows: list[dict[str, Any]] = []

    while current_a is not None and current_b is not None:
        cmp = _compare_positions(current_a, current_b)
        if cmp < 0:
            current_a = next(iter_a, None)
            continue
        if cmp > 0:
            current_b = next(iter_b, None)
            continue

        other_for_a = None
        if current_a.nonref_snp_call is not None and not _has_any_alt_support(current_b, current_a.nonref_snp_call.alt):
            other_for_a = _absence_supporting_record(current_b, current_a.nonref_snp_call.alt, min_dp_absence=min_dp_absence, max_dp=max_dp, min_other_af=min_other_af)
        if current_a.nonref_snp_call is not None and other_for_a is not None:
            group_a_rows.append(_row_from_match(
                direction="group_a_only",
                target_group=group_a_label,
                other_group=group_b_label,
                target_call=current_a.nonref_snp_call,
                other_call=other_for_a,
            ))

        other_for_b = None
        if current_b.nonref_snp_call is not None and not _has_any_alt_support(current_a, current_b.nonref_snp_call.alt):
            other_for_b = _absence_supporting_record(current_a, current_b.nonref_snp_call.alt, min_dp_absence=min_dp_absence, max_dp=max_dp, min_other_af=min_other_af)
        if current_b.nonref_snp_call is not None and other_for_b is not None:
            group_b_rows.append(_row_from_match(
                direction="group_b_only",
                target_group=group_b_label,
                other_group=group_a_label,
                target_call=current_b.nonref_snp_call,
                other_call=other_for_b,
            ))

        current_a = next(iter_a, None)
        current_b = next(iter_b, None)

    def _sort_key(row: dict[str, Any]) -> float:
        af = row.get("other_af", ".")
        return float(af) if af != "." else 0.0

    group_a_rows.sort(key=_sort_key, reverse=True)
    group_b_rows.sort(key=_sort_key, reverse=True)
    all_rows = sorted(group_a_rows + group_b_rows, key=_sort_key, reverse=True)

    merged_handle, merged_writer = _open_writer(merged_path)
    group_a_handle, group_a_writer = _open_writer(group_a_only_path)
    group_b_handle, group_b_writer = _open_writer(group_b_only_path)
    try:
        for row in all_rows:
            merged_writer.writerow(row)
        for row in group_a_rows:
            group_a_writer.writerow(row)
        for row in group_b_rows:
            group_b_writer.writerow(row)
    finally:
        merged_handle.close()
        group_a_handle.close()
        group_b_handle.close()

    summary: dict[str, Any] = {
        "group_a_only_count": len(group_a_rows),
        "group_b_only_count": len(group_b_rows),
        "group_a_only_top": f"{group_a_rows[0]['chrom']}:{group_a_rows[0]['pos']}:{group_a_rows[0]['ref']}:{group_a_rows[0]['alt']}" if group_a_rows else "",
        "group_b_only_top": f"{group_b_rows[0]['chrom']}:{group_b_rows[0]['pos']}:{group_b_rows[0]['ref']}:{group_b_rows[0]['alt']}" if group_b_rows else "",
        "merged_count": len(all_rows),
    }

    summary["group_a_only_tsv"] = str(group_a_only_path)
    summary["group_b_only_tsv"] = str(group_b_only_path)
    summary["merged_tsv"] = str(merged_path)
    return summary


def _build_arg_parser(
    *,
    prog: str = "python -m sniffcell.discover.snv_post_processing",
    add_help: bool = True,
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Compare two Clair3 gVCF outputs and emit group-specific SNPs supported by explicit reference calls in the other group.",
        add_help=add_help,
    )
    parser.add_argument("--split-dir", required=True, help="deconv_requested_group_splits directory")
    parser.add_argument("--groups", required=True, help="Exactly two group names, comma-separated")
    parser.add_argument("--group-a-gvcf", required=True, help="gVCF or pileup gVCF for the first group")
    parser.add_argument("--group-b-gvcf", required=True, help="gVCF or pileup gVCF for the second group")
    parser.add_argument("--output-dir", default=None, help="Output directory. Default: <split-dir>/postprocess/snv_post_processing_<timestamp>")
    parser.add_argument("--sample-id", default=None, help="Optional sample ID override")
    parser.add_argument("--min-dp", type=int, default=5, help="Minimum DP for target existence calls. Default=5")
    parser.add_argument("--max-dp", type=int, default=50, help="Maximum DP for target existence calls. Sites above this are likely repetitive/low-complexity regions. Default=50")
    parser.add_argument("--min-dp-absence", type=int, default=15, help="Minimum DP in the other group to confidently assert absence of an ALT allele. Default=15")
    parser.add_argument("--min-gq", type=int, default=21, help="Minimum GQ for target non-reference SNV calls. Default=21")
    parser.add_argument("--min-other-af", type=float, default=0.9, help="Minimum reference allele fraction in the other group's absence record. For Clair3 RefCall records, AF=ref_reads/DP, so 0.9 requires 90%% of reads to support the reference. Default=0.9")
    return parser


def _resolve_args(raw_args: Any) -> TwoSampleSnvArgs:
    split_dir = _expand_path(raw_args.split_dir)
    tokens = [x.strip() for x in str(raw_args.groups).split(",") if x.strip()]
    if len(tokens) != 2:
        raise ValueError("--groups must contain exactly two group names")
    output_dir = (
        _expand_path(raw_args.output_dir)
        if raw_args.output_dir
        else split_dir / "postprocess" / f"snv_post_processing_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
    return TwoSampleSnvArgs(
        split_dir=split_dir,
        output_dir=output_dir,
        group_a=tokens[0],
        group_b=tokens[1],
        group_a_gvcf=_expand_path(raw_args.group_a_gvcf),
        group_b_gvcf=_expand_path(raw_args.group_b_gvcf),
        sample_id=raw_args.sample_id or _infer_sample_id(split_dir.parent),
        min_dp=int(raw_args.min_dp),
        max_dp=int(raw_args.max_dp),
        min_dp_absence=int(raw_args.min_dp_absence),
        min_gq=int(raw_args.min_gq),
        min_other_af=float(raw_args.min_other_af),
    )


def snv_post_processing_main(cli_args: list[str] | None = None) -> dict[str, Any]:
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(message)s")
    parser = _build_arg_parser()
    args = _resolve_args(parser.parse_args(cli_args))

    for path in (args.split_dir, args.group_a_gvcf, args.group_b_gvcf):
        if not path.exists():
            raise FileNotFoundError(path)

    group_a_label = f"{args.sample_id}.{args.group_a}"
    group_b_label = f"{args.sample_id}.{args.group_b}"
    summary = compare_group_specific_snvs(
        group_a_label=group_a_label,
        group_b_label=group_b_label,
        group_a_gvcf=args.group_a_gvcf,
        group_b_gvcf=args.group_b_gvcf,
        output_dir=args.output_dir,
        min_dp=args.min_dp,
        max_dp=args.max_dp,
        min_dp_absence=args.min_dp_absence,
        min_gq=args.min_gq,
        min_other_af=args.min_other_af,
    )
    summary.update(
        {
            "split_dir": str(args.split_dir),
            "sample_id": args.sample_id,
            "group_a": args.group_a,
            "group_b": args.group_b,
            "group_a_gvcf": str(args.group_a_gvcf),
            "group_b_gvcf": str(args.group_b_gvcf),
            "params": {
                "min_dp": args.min_dp,
                "max_dp": args.max_dp,
                "min_dp_absence": args.min_dp_absence,
                "min_gq": args.min_gq,
                "min_other_af": args.min_other_af,
            },
        }
    )

    readme_path = args.output_dir / "README.txt"
    readme_path.write_text(
        "\n".join(
            [
                "Two-group SNP post-processing from Clair3 gVCF outputs.",
                "Logic: keep PASS SNP non-reference calls in one group only when the target call has AD-supported alt reads, GQ above threshold, and the other group has zero AD support for that same ALT at the same coordinate.",
                f"Minimum DP for target existence calls: {args.min_dp}",
                f"Maximum DP for target existence calls (above = likely repetitive region): {args.max_dp}",
                f"Minimum DP for absence support in other group: {args.min_dp_absence}",
                f"Minimum GQ for target non-reference SNP calls: {args.min_gq}",
                f"Minimum reference allele fraction (AF) in other group's absence record: {args.min_other_af}",
                f"Group A label: {group_a_label}",
                f"Group B label: {group_b_label}",
                f"Group A gVCF: {args.group_a_gvcf}",
                f"Group B gVCF: {args.group_b_gvcf}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    summary["readme"] = str(readme_path)
    _write_json(args.output_dir / "summary.json", summary)
    return summary


def main(argv: list[str] | None = None) -> int:
    snv_post_processing_main(argv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
