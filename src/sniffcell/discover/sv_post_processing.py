from __future__ import annotations

import argparse
import csv
import gzip
import logging
import shutil
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from sniffcell.postprocess.postprocess import (
    DEFAULT_MOSAIC_FILTER_EXPR,
    REPO_ROOT,
    SplitGroup,
    _ensure_tool,
    _expand_path,
    _infer_sample_id,
    _read_json,
    _sanitize_token,
    _write_json,
)


@dataclass(frozen=True)
class TwoSampleSvArgs:
    split_dir: Path
    reference: Path
    output_dir: Path
    group_a: str
    group_b: str
    mosaic_filter_expression: str
    min_total_ad: int
    min_target_alt_ad: int
    other_max_alt_ad: int
    bcftools_bin: str
    truvari_bin: str
    kanpig_bin: str
    threads: int
    kanpig_seqsim: float
    kanpig_sizesim: float
    sample_id: str


def _discover_groups(split_dir: Path) -> list[SplitGroup]:
    manifest_path = split_dir / "requested_group_splits.tsv"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing requested split manifest: {manifest_path}")
    with manifest_path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    groups: list[SplitGroup] = []
    for row in rows:
        group_name = (row.get("requested_group") or "").strip()
        bam_path = (row.get("bam_path") or "").strip()
        if not group_name or not bam_path:
            raise ValueError(f"Malformed row in {manifest_path}: {row}")
        bam = _expand_path(bam_path)
        groups.append(
            SplitGroup(
                name=group_name,
                bam_path=str(bam),
                bai_path=str(_expand_path(str(bam) + ".bai")),
                read_summary_path=row.get("read_summary_path") or None,
            )
        )
    return groups


def _group_lookup(groups: list[SplitGroup], name: str) -> SplitGroup:
    for group in groups:
        if group.name == name:
            return group
    raise KeyError(name)


def _sniffles_vcf_path(split_dir: Path, group_name: str) -> Path:
    return split_dir / f"{group_name}.sniffles.vcf.gz"


def _run(cmd: list[str], *, stdout_path: Path | None = None, stderr_path: Path | None = None) -> None:
    if stdout_path is not None:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
    if stderr_path is not None:
        stderr_path.parent.mkdir(parents=True, exist_ok=True)
    stdout_handle = stdout_path.open("a") if stdout_path is not None else None
    stderr_handle = stderr_path.open("a") if stderr_path is not None else None
    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=stdout_handle if stdout_handle is not None else subprocess.DEVNULL,
            stderr=stderr_handle if stderr_handle is not None else subprocess.DEVNULL,
        )
    finally:
        if stdout_handle is not None:
            stdout_handle.close()
        if stderr_handle is not None:
            stderr_handle.close()


def _run_shell(command: str, *, stdout_path: Path, stderr_path: Path) -> None:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("a") as stdout_handle, stderr_path.open("a") as stderr_handle:
        subprocess.run(command, shell=True, check=True, stdout=stdout_handle, stderr=stderr_handle)


def _filter_sniffles_vcf(
    *,
    bcftools_bin: str,
    input_vcf: Path,
    output_vcf: Path,
    expr: str,
    stdout_path: Path,
    stderr_path: Path,
) -> None:
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    _run(
        [
            bcftools_bin,
            "view",
            "-f",
            "PASS",
            "-i",
            expr,
            "-Oz",
            "-o",
            str(output_vcf),
            str(input_vcf),
        ],
        stdout_path=stdout_path,
        stderr_path=stderr_path,
    )
    _run(
        [bcftools_bin, "index", "-t", "-f", str(output_vcf)],
        stdout_path=stdout_path,
        stderr_path=stderr_path,
    )


def _collapse_two_vcfs(
    *,
    bcftools_bin: str,
    truvari_bin: str,
    reference: Path,
    group_a_vcf: Path,
    group_b_vcf: Path,
    stage_dir: Path,
    stdout_path: Path,
    stderr_path: Path,
) -> tuple[Path, Path]:
    stage_dir.mkdir(parents=True, exist_ok=True)
    sampleless_a = stage_dir / f"{group_a_vcf.stem.replace('.vcf', '')}.sites.vcf.gz"
    sampleless_b = stage_dir / f"{group_b_vcf.stem.replace('.vcf', '')}.sites.vcf.gz"
    merged_input = stage_dir / "collapse.inputs.vcf.gz"
    raw_output = stage_dir / "collapsed.vcf"
    removed_vcf = stage_dir / "removed.vcf"
    sorted_output = stage_dir / "collapsed.sorted.vcf.gz"

    merge_cmd = "\n".join(
        [
            "set -euo pipefail",
            " ".join([bcftools_bin, "view", "-G", "-Oz", "-o", str(sampleless_a), str(group_a_vcf)]),
            " ".join([bcftools_bin, "index", "-t", "-f", str(sampleless_a)]),
            " ".join([bcftools_bin, "view", "-G", "-Oz", "-o", str(sampleless_b), str(group_b_vcf)]),
            " ".join([bcftools_bin, "index", "-t", "-f", str(sampleless_b)]),
            " ".join([bcftools_bin, "concat", "-a", "-Oz", "-o", str(merged_input), str(sampleless_a), str(sampleless_b)]),
            " ".join([bcftools_bin, "index", "-t", "-f", str(merged_input)]),
        ]
    )
    _run_shell(merge_cmd, stdout_path=stdout_path, stderr_path=stderr_path)
    _run(
        [
            truvari_bin,
            "collapse",
            "-i",
            str(merged_input),
            "-o",
            str(raw_output),
            "-c",
            str(removed_vcf),
            "-f",
            str(reference),
            "-r",
            "500",
            "-p",
            "0.95",
            "-P",
            "0.95",
            "--passonly",
        ],
        stdout_path=stdout_path,
        stderr_path=stderr_path,
    )
    _run(
        [bcftools_bin, "sort", "-Oz", "-o", str(sorted_output), str(raw_output)],
        stdout_path=stdout_path,
        stderr_path=stderr_path,
    )
    _run(
        [bcftools_bin, "index", "-t", "-f", str(sorted_output)],
        stdout_path=stdout_path,
        stderr_path=stderr_path,
    )
    return sorted_output, removed_vcf


def _run_kanpig_on_merged_sites(
    *,
    bcftools_bin: str,
    kanpig_bin: str,
    merged_vcf: Path,
    group_a: SplitGroup,
    group_b: SplitGroup,
    reference: Path,
    output_vcf: Path,
    rnames_tsv: Path,
    sample_id: str,
    threads: int,
    seqsim: float,
    sizesim: float,
    stdout_path: Path,
    stderr_path: Path,
) -> Path:
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    raw_output = output_vcf.with_suffix("")
    _run(
        [
            kanpig_bin,
            "mosaic",
            "--input",
            str(merged_vcf),
            "--reads",
            group_a.bam_path,
            "--reads",
            group_b.bam_path,
            "--sample",
            f"{sample_id}_{group_a.name}",
            "--sample",
            f"{sample_id}_{group_b.name}",
            "--reference",
            str(reference),
            "--threads",
            str(threads),
            "--passonly",
            "--seqsim",
            str(seqsim),
            "--sizesim",
            str(sizesim),
            "--rnames",
            str(rnames_tsv),
            "--out",
            str(raw_output),
        ],
        stdout_path=stdout_path,
        stderr_path=stderr_path,
    )
    _run(
        [bcftools_bin, "sort", "-Oz", "-o", str(output_vcf), str(raw_output)],
        stdout_path=stdout_path,
        stderr_path=stderr_path,
    )
    _run(
        [bcftools_bin, "index", "-t", "-f", str(output_vcf)],
        stdout_path=stdout_path,
        stderr_path=stderr_path,
    )
    return output_vcf


def _alt_present(gt: str) -> bool:
    return gt in {"0/1", "1/0", "1/1", "0|1", "1|0", "1|1"}


def _parse_ad(sample_field: str) -> tuple[str, int, int] | None:
    parts = sample_field.rstrip().split(":")
    if len(parts) < 7:
        return None
    gt = parts[0]
    ad = parts[6]
    if ad in {".", "./."}:
        return gt, 0, 0
    vals = ad.split(",")
    if len(vals) < 2:
        return gt, 0, 0
    ref = int(vals[0]) if vals[0] != "." else 0
    alt = int(vals[1]) if vals[1] != "." else 0
    return gt, ref, alt


def _filter_sample_specific_by_ad(
    *,
    kanpig_vcf_gz: Path,
    sample_a_label: str,
    sample_b_label: str,
    min_total_ad: int,
    min_target_alt_ad: int,
    other_max_alt_ad: int,
    output_dir: Path,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    out_a = output_dir / f"{_sanitize_token(sample_a_label)}.only.AD_logic.vcf"
    out_b = output_dir / f"{_sanitize_token(sample_b_label)}.only.AD_logic.vcf"
    out_shared = output_dir / "shared.both_samples.AD_logic.vcf"

    opener = gzip.open
    with kanpig_vcf_gz.open("rb") as probe:
        if probe.read(2) != b"\x1f\x8b":
            opener = open

    with opener(kanpig_vcf_gz, "rt", encoding="utf-8") as inp, \
        out_a.open("w", encoding="utf-8") as a_handle, \
        out_b.open("w", encoding="utf-8") as b_handle, \
        out_shared.open("w", encoding="utf-8") as shared_handle:
        for line in inp:
            if line.startswith("#"):
                a_handle.write(line)
                b_handle.write(line)
                shared_handle.write(line)
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) < 11:
                continue
            parsed_a = _parse_ad(row[9])
            parsed_b = _parse_ad(row[10])
            if parsed_a is None or parsed_b is None:
                continue
            gt_a, ref_a, alt_a = parsed_a
            gt_b, ref_b, alt_b = parsed_b
            if (ref_a + alt_a) <= min_total_ad or (ref_b + alt_b) <= min_total_ad:
                continue
            if alt_a >= min_target_alt_ad and alt_b <= other_max_alt_ad:
                a_handle.write(line)
            elif alt_b >= min_target_alt_ad and alt_a <= other_max_alt_ad:
                b_handle.write(line)
            elif _alt_present(gt_a) and _alt_present(gt_b):
                shared_handle.write(line)

    for plain in (out_a, out_b, out_shared):
        subprocess.run(["bgzip", "-f", str(plain)], check=True)
        subprocess.run(["tabix", "-f", "-p", "vcf", str(plain) + ".gz"], check=True)
    return {
        "sample_a_only": Path(str(out_a) + ".gz"),
        "sample_b_only": Path(str(out_b) + ".gz"),
        "shared": Path(str(out_shared) + ".gz"),
    }


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m sniffcell.postprocess.sv_post_processing",
        description="Compare two per-BAM Sniffles VCFs, genotype merged sites with Kanpig, and emit AD-based sample-specific SV sets.",
    )
    parser.add_argument("--split-dir", required=True, help="deconv_requested_group_splits directory")
    parser.add_argument("--reference", required=True, help="Reference FASTA")
    parser.add_argument("--groups", required=True, help="Exactly two group names, comma-separated")
    parser.add_argument("--output-dir", default=None, help="Output directory. Default: <split-dir>/postprocess/sv_post_processing_<timestamp>")
    parser.add_argument("--sample-id", default=None, help="Optional sample ID override")
    parser.add_argument("--mosaic-filter-expression", default=DEFAULT_MOSAIC_FILTER_EXPR)
    parser.add_argument("--min-total-ad", type=int, default=5, help="Require AD_ref + AD_alt > this threshold in both samples. Default=5")
    parser.add_argument("--min-target-alt-ad", type=int, default=1, help="Minimum AD_alt for the kept sample. Default=1")
    parser.add_argument("--other-max-alt-ad", type=int, default=0, help="Maximum AD_alt allowed in the absent sample. Default=0")
    parser.add_argument("--bcftools-bin", default=None)
    parser.add_argument("--truvari-bin", default=None)
    parser.add_argument("--kanpig-bin", default=None)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--kanpig-seqsim", type=float, default=0.8)
    parser.add_argument("--kanpig-sizesim", type=float, default=0.85)
    return parser


def _resolve_args(raw_args) -> TwoSampleSvArgs:
    split_dir = _expand_path(raw_args.split_dir)
    reference = _expand_path(raw_args.reference)
    tokens = [x.strip() for x in str(raw_args.groups).split(",") if x.strip()]
    if len(tokens) != 2:
        raise ValueError("--groups must contain exactly two group names")
    sample_id = raw_args.sample_id or _infer_sample_id(split_dir.parent)
    output_dir = _expand_path(raw_args.output_dir) if raw_args.output_dir else split_dir / "postprocess" / f"sv_post_processing_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    return TwoSampleSvArgs(
        split_dir=split_dir,
        reference=reference,
        output_dir=output_dir,
        group_a=tokens[0],
        group_b=tokens[1],
        mosaic_filter_expression=str(raw_args.mosaic_filter_expression),
        min_total_ad=int(raw_args.min_total_ad),
        min_target_alt_ad=int(raw_args.min_target_alt_ad),
        other_max_alt_ad=int(raw_args.other_max_alt_ad),
        bcftools_bin=_ensure_tool("bcftools", raw_args.bcftools_bin),
        truvari_bin=_ensure_tool("truvari", raw_args.truvari_bin),
        kanpig_bin=_ensure_tool("kanpig", raw_args.kanpig_bin),
        threads=int(raw_args.threads),
        kanpig_seqsim=float(raw_args.kanpig_seqsim),
        kanpig_sizesim=float(raw_args.kanpig_sizesim),
        sample_id=sample_id,
    )


def sv_post_processing_main(cli_args=None) -> dict[str, str]:
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s %(message)s")
    parser = _build_arg_parser()
    args = _resolve_args(parser.parse_args(cli_args))

    groups = _discover_groups(args.split_dir)
    group_a = _group_lookup(groups, args.group_a)
    group_b = _group_lookup(groups, args.group_b)
    for path in (
        Path(group_a.bam_path),
        Path(group_a.bai_path),
        Path(group_b.bam_path),
        Path(group_b.bai_path),
        _sniffles_vcf_path(args.split_dir, group_a.name),
        _sniffles_vcf_path(args.split_dir, group_b.name),
        args.reference,
    ):
        if not path.exists():
            raise FileNotFoundError(path)

    logs_dir = args.output_dir / "logs"
    filtered_dir = args.output_dir / "sv" / "sniffles"
    collapse_dir = args.output_dir / "sv" / "truvari_collapse" / f"{_sanitize_token(group_a.name)}_vs_{_sanitize_token(group_b.name)}"
    kanpig_dir = args.output_dir / "kanpig_merged"
    ad_dir = args.output_dir / "ad_filtered"
    for path in (args.output_dir, logs_dir, filtered_dir, collapse_dir, kanpig_dir, ad_dir):
        path.mkdir(parents=True, exist_ok=True)

    filtered_a = filtered_dir / _sanitize_token(group_a.name) / "sniffles.mosaic_only.vcf.gz"
    filtered_b = filtered_dir / _sanitize_token(group_b.name) / "sniffles.mosaic_only.vcf.gz"
    _filter_sniffles_vcf(
        bcftools_bin=args.bcftools_bin,
        input_vcf=_sniffles_vcf_path(args.split_dir, group_a.name),
        output_vcf=filtered_a,
        expr=args.mosaic_filter_expression,
        stdout_path=logs_dir / f"sniffles_filter.{_sanitize_token(group_a.name)}.out",
        stderr_path=logs_dir / f"sniffles_filter.{_sanitize_token(group_a.name)}.err",
    )
    _filter_sniffles_vcf(
        bcftools_bin=args.bcftools_bin,
        input_vcf=_sniffles_vcf_path(args.split_dir, group_b.name),
        output_vcf=filtered_b,
        expr=args.mosaic_filter_expression,
        stdout_path=logs_dir / f"sniffles_filter.{_sanitize_token(group_b.name)}.out",
        stderr_path=logs_dir / f"sniffles_filter.{_sanitize_token(group_b.name)}.err",
    )

    collapsed_sorted, removed_vcf = _collapse_two_vcfs(
        bcftools_bin=args.bcftools_bin,
        truvari_bin=args.truvari_bin,
        reference=args.reference,
        group_a_vcf=filtered_a,
        group_b_vcf=filtered_b,
        stage_dir=collapse_dir,
        stdout_path=logs_dir / "collapse.out",
        stderr_path=logs_dir / "collapse.err",
    )

    kanpig_merged = _run_kanpig_on_merged_sites(
        bcftools_bin=args.bcftools_bin,
        kanpig_bin=args.kanpig_bin,
        merged_vcf=collapsed_sorted,
        group_a=group_a,
        group_b=group_b,
        reference=args.reference,
        output_vcf=kanpig_dir / "kanpig_merged.sorted.vcf.gz",
        rnames_tsv=kanpig_dir / "kanpig_merged.rnames.tsv",
        sample_id=args.sample_id,
        threads=args.threads,
        seqsim=args.kanpig_seqsim,
        sizesim=args.kanpig_sizesim,
        stdout_path=logs_dir / "kanpig_merged.out",
        stderr_path=logs_dir / "kanpig_merged.err",
    )

    filtered_outputs = _filter_sample_specific_by_ad(
        kanpig_vcf_gz=kanpig_merged,
        sample_a_label=group_a.name,
        sample_b_label=group_b.name,
        min_total_ad=args.min_total_ad,
        min_target_alt_ad=args.min_target_alt_ad,
        other_max_alt_ad=args.other_max_alt_ad,
        output_dir=ad_dir,
    )

    summary = {
        "split_dir": str(args.split_dir),
        "reference": str(args.reference),
        "group_a": group_a.name,
        "group_b": group_b.name,
        "filtered_group_a_vcf": str(filtered_a),
        "filtered_group_b_vcf": str(filtered_b),
        "collapsed_sorted_vcf": str(collapsed_sorted),
        "removed_vcf": str(removed_vcf),
        "kanpig_merged_vcf": str(kanpig_merged),
        "sample_a_only_vcf": str(filtered_outputs["sample_a_only"]),
        "sample_b_only_vcf": str(filtered_outputs["sample_b_only"]),
        "shared_vcf": str(filtered_outputs["shared"]),
        "params": {
            "min_total_ad": args.min_total_ad,
            "min_target_alt_ad": args.min_target_alt_ad,
            "other_max_alt_ad": args.other_max_alt_ad,
            "mosaic_filter_expression": args.mosaic_filter_expression,
            "threads": args.threads,
            "kanpig_seqsim": args.kanpig_seqsim,
            "kanpig_sizesim": args.kanpig_sizesim,
        },
    }
    _write_json(args.output_dir / "summary.json", summary)
    return summary


def main() -> int:
    sv_post_processing_main()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
