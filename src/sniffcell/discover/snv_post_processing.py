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
    min_gq: int


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
    af: str


@dataclass(frozen=True)
class PositionSummary:
    chrom: str
    pos: int
    nonref_call: CallRecord | None
    nonref_snp_call: CallRecord | None
    ref_call: CallRecord | None


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


def _is_ref_gt(gt: str) -> bool:
    return gt in {"0/0", "0|0"}


def _is_nonref_gt(gt: str) -> bool:
    return gt not in {"0/0", "0|0", "./.", ".|.", ".", "ref", "mis"}


def _is_snp(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and alt != "." and "," not in alt


def _record_rank(rec: CallRecord) -> tuple[int, int]:
    return (rec.gq, rec.dp)


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
        af=_parse_float_str(fmt_map.get("AF", ".")),
    )


def _summarize_position(records: list[CallRecord], *, min_dp: int, min_gq: int) -> PositionSummary:
    chrom = records[0].chrom
    pos = records[0].pos
    nonref_call: CallRecord | None = None
    nonref_snp_call: CallRecord | None = None
    ref_call: CallRecord | None = None
    for rec in records:
        if _is_nonref_gt(rec.gt) and rec.filt == "PASS" and rec.dp >= min_dp and rec.gq >= min_gq:
            nonref_call = _better_record(nonref_call, rec)
            if _is_snp(rec.ref, rec.alt):
                nonref_snp_call = _better_record(nonref_snp_call, rec)
        elif (
            rec.filt == "RefCall"
            and rec.alt == "."
            and _is_ref_gt(rec.gt)
            and rec.dp >= min_dp
            and rec.gq >= min_gq
        ):
            ref_call = _better_record(ref_call, rec)
    return PositionSummary(
        chrom=chrom,
        pos=pos,
        nonref_call=nonref_call,
        nonref_snp_call=nonref_snp_call,
        ref_call=ref_call,
    )


def _iter_position_summaries(path: Path, *, min_dp: int, min_gq: int) -> Iterator[PositionSummary]:
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
            yield _summarize_position(current, min_dp=min_dp, min_gq=min_gq)
            current_key = key
            current = [rec]
        if current:
            yield _summarize_position(current, min_dp=min_dp, min_gq=min_gq)


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
        "target_af",
        "other_gt",
        "other_gq",
        "other_dp",
        "other_af",
    ]


def _row_from_match(
    *,
    direction: str,
    target_group: str,
    other_group: str,
    target_call: CallRecord,
    other_ref: CallRecord,
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
        "target_af": target_call.af,
        "other_gt": other_ref.gt,
        "other_gq": other_ref.gq,
        "other_dp": other_ref.dp,
        "other_af": other_ref.af,
    }


def _open_writer(path: Path) -> tuple[Any, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = path.open("w", newline="", encoding="utf-8")
    writer = csv.DictWriter(handle, fieldnames=_merged_headers(), delimiter="\t")
    writer.writeheader()
    return handle, writer


def compare_group_specific_snvs(
    *,
    group_a_label: str,
    group_b_label: str,
    group_a_gvcf: Path,
    group_b_gvcf: Path,
    output_dir: Path,
    min_dp: int,
    min_gq: int,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    sets_dir = output_dir / "sets"
    sets_dir.mkdir(parents=True, exist_ok=True)
    merged_path = output_dir / "snv_changes.tsv"
    group_a_only_path = sets_dir / f"{_sanitize_token(group_a_label)}.only.tsv"
    group_b_only_path = sets_dir / f"{_sanitize_token(group_b_label)}.only.tsv"

    merged_handle, merged_writer = _open_writer(merged_path)
    group_a_handle, group_a_writer = _open_writer(group_a_only_path)
    group_b_handle, group_b_writer = _open_writer(group_b_only_path)
    try:
        iter_a = iter(_iter_position_summaries(group_a_gvcf, min_dp=min_dp, min_gq=min_gq))
        iter_b = iter(_iter_position_summaries(group_b_gvcf, min_dp=min_dp, min_gq=min_gq))
        current_a = next(iter_a, None)
        current_b = next(iter_b, None)

        summary: dict[str, Any] = {
            "group_a_only_count": 0,
            "group_b_only_count": 0,
            "group_a_only_top": "",
            "group_b_only_top": "",
            "merged_count": 0,
        }

        while current_a is not None and current_b is not None:
            cmp = _compare_positions(current_a, current_b)
            if cmp < 0:
                current_a = next(iter_a, None)
                continue
            if cmp > 0:
                current_b = next(iter_b, None)
                continue

            if current_a.nonref_snp_call is not None and current_b.nonref_call is None and current_b.ref_call is not None:
                row = _row_from_match(
                    direction="group_a_only",
                    target_group=group_a_label,
                    other_group=group_b_label,
                    target_call=current_a.nonref_snp_call,
                    other_ref=current_b.ref_call,
                )
                group_a_writer.writerow(row)
                merged_writer.writerow(row)
                summary["group_a_only_count"] += 1
                summary["merged_count"] += 1
                if not summary["group_a_only_top"]:
                    summary["group_a_only_top"] = f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"

            if current_b.nonref_snp_call is not None and current_a.nonref_call is None and current_a.ref_call is not None:
                row = _row_from_match(
                    direction="group_b_only",
                    target_group=group_b_label,
                    other_group=group_a_label,
                    target_call=current_b.nonref_snp_call,
                    other_ref=current_a.ref_call,
                )
                group_b_writer.writerow(row)
                merged_writer.writerow(row)
                summary["group_b_only_count"] += 1
                summary["merged_count"] += 1
                if not summary["group_b_only_top"]:
                    summary["group_b_only_top"] = f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"

            current_a = next(iter_a, None)
            current_b = next(iter_b, None)
    finally:
        merged_handle.close()
        group_a_handle.close()
        group_b_handle.close()

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
    parser.add_argument("--min-dp", type=int, default=5, help="Minimum DP for both target existence and explicit absence. Default=5")
    parser.add_argument("--min-gq", type=int, default=0, help="Minimum GQ for both target existence and explicit absence. Default=0")
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
        min_gq=int(raw_args.min_gq),
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
        min_gq=args.min_gq,
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
                "min_gq": args.min_gq,
            },
        }
    )

    readme_path = args.output_dir / "README.txt"
    readme_path.write_text(
        "\n".join(
            [
                "Two-group SNP post-processing from Clair3 gVCF outputs.",
                "Logic: keep PASS SNP non-reference calls in one group only when the other group has no passing non-reference call at that coordinate and does have an explicit RefCall 0/0.",
                f"Minimum DP for both existence and absence: {args.min_dp}",
                f"Minimum GQ for both existence and absence: {args.min_gq}",
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
