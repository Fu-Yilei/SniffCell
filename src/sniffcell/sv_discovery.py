#!/usr/bin/env python3
"""
Two-pass Sniffles discovery wrapper for confident SV inputs to SniffCell.

Design:
1. Run Sniffles in TR regions (high sensitivity defaults).
2. Run Sniffles in non-TR regions (mosaic defaults with configurable AF floor).
3. Apply confidence filters (SVTYPE, PASS, optional Q100 BED).
4. Optionally merge TR/non-TR confident calls into one VCF.
"""

from __future__ import annotations

import argparse
import json
import logging
import shlex
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TR_BED = str(REPO_ROOT / "atlas" / "adotto.v2.trgt.bed")
DEFAULT_Q100_BED = str(REPO_ROOT / "atlas" / "GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed")


@dataclass
class CommandRunner:
    dry_run: bool = False
    commands: list[str] | None = None

    def __post_init__(self) -> None:
        if self.commands is None:
            self.commands = []

    def run(self, cmd: Sequence[str] | str, shell: bool = False) -> None:
        if shell:
            if not isinstance(cmd, str):
                raise TypeError("shell=True requires a string command")
            cmd_str = cmd
        else:
            if isinstance(cmd, str):
                raise TypeError("shell=False requires a sequence command")
            cmd_str = shlex.join(cmd)

        self.commands.append(cmd_str)
        logging.info("CMD: %s", cmd_str)
        if self.dry_run:
            return

        subprocess.run(cmd, shell=shell, check=True)


def _expand_path(value: str | Path) -> Path:
    return Path(value).expanduser().resolve()


def _ensure_executable(tool: str) -> str:
    found = shutil.which(tool)
    if not found:
        raise FileNotFoundError(f"Required executable not found in PATH: {tool}")
    return found


def _parse_svtypes(text: str) -> tuple[str, ...]:
    cleaned = [x.strip().upper() for x in text.split(",") if x.strip()]
    return tuple(cleaned)


def _build_svtype_expr(svtypes: Iterable[str]) -> str | None:
    tokens = [t for t in svtypes if t]
    if not tokens:
        return None
    terms = [f'INFO/SVTYPE="{svtype}"' for svtype in tokens]
    return " || ".join(terms)


def _build_additional_filter_expr(
    *,
    max_vaf: float | None,
    max_del_len: int | None,
    min_support: int | None,
    max_support: int | None,
) -> str | None:
    terms: list[str] = []
    if max_vaf is not None:
        # Prefer VAF when present; fall back to AF when VAF is missing.
        terms.append(f'((INFO/VAF<={max_vaf:g}) || (INFO/VAF="." && INFO/AF<={max_vaf:g}))')
    if max_del_len is not None:
        terms.append(f'(INFO/SVTYPE!="DEL" || ABS(INFO/SVLEN)<={max_del_len})')
    if min_support is not None:
        terms.append(f"(INFO/SUPPORT>={min_support})")
    if max_support is not None:
        terms.append(f"(INFO/SUPPORT<={max_support})")
    if not terms:
        return None
    return " && ".join(terms)


def _index_vcf(vcf_path: Path, bcftools_bin: str, runner: CommandRunner) -> None:
    runner.run([bcftools_bin, "index", "-t", "-f", str(vcf_path)])


def _ensure_fasta_index(reference_fa: Path, samtools_bin: str, runner: CommandRunner) -> Path:
    fai = Path(str(reference_fa) + ".fai")
    if fai.exists():
        return fai
    logging.info("Reference index missing, creating: %s", fai)
    runner.run([samtools_bin, "faidx", str(reference_fa)])
    return fai


def _write_complement_bed(
    tr_bed: Path,
    genome_fai: Path,
    output_bed: Path,
    bedtools_bin: str,
    runner: CommandRunner,
    force_recompute: bool = False,
) -> None:
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    if output_bed.exists() and not force_recompute:
        logging.info("Using existing non-TR complement BED: %s", output_bed)
        return

    cmd = (
        "set -euo pipefail; "
        f"{shlex.quote(bedtools_bin)} sort -faidx {shlex.quote(str(genome_fai))} -i {shlex.quote(str(tr_bed))} | "
        f"{shlex.quote(bedtools_bin)} merge -i - | "
        f"{shlex.quote(bedtools_bin)} complement -i - -g {shlex.quote(str(genome_fai))} "
        f"> {shlex.quote(str(output_bed))}"
    )
    runner.run(cmd, shell=True)


def _build_sniffles_command(
    *,
    sniffles_bin: str,
    input_bam: Path,
    reference_fa: Path,
    output_vcf: Path,
    output_snf: Path,
    threads: int,
    mosaic_af_min: float,
    include_germline: bool,
    no_qc: bool,
    regions_bed: Path | None = None,
    cluster_merge_len: float | None = None,
    extra_args: Sequence[str] = (),
) -> list[str]:
    cmd = [
        sniffles_bin,
        "--input",
        str(input_bam),
        "--reference",
        str(reference_fa),
        "--vcf",
        str(output_vcf),
        "--snf",
        str(output_snf),
        "--threads",
        str(threads),
        "--mosaic",
        "--mosaic-af-min",
        str(mosaic_af_min),
        "--output-rnames",
        "--allow-overwrite",
    ]
    if include_germline:
        cmd.append("--mosaic-include-germline")
    if no_qc:
        cmd.append("--no-qc")
    if regions_bed is not None:
        cmd.extend(["--regions", str(regions_bed)])
    if cluster_merge_len is not None:
        cmd.extend(["--cluster-merge-len", str(cluster_merge_len)])
    cmd.extend(extra_args)
    return cmd


def _split_onepass_vcf_with_bedtools(
    *,
    bedtools_bin: str,
    bcftools_bin: str,
    input_vcf: Path,
    tr_bed: Path,
    tr_output_vcf: Path,
    nontr_output_vcf: Path,
    genome_fai: Path,
    runner: CommandRunner,
) -> None:
    tr_cmd = (
        "set -euo pipefail; "
        f"{shlex.quote(bedtools_bin)} intersect -header -sorted -g {shlex.quote(str(genome_fai))} "
        f"-a {shlex.quote(str(input_vcf))} -b {shlex.quote(str(tr_bed))} | "
        f"{shlex.quote(bcftools_bin)} view -Oz -o {shlex.quote(str(tr_output_vcf))} -"
    )
    runner.run(tr_cmd, shell=True)
    _index_vcf(tr_output_vcf, bcftools_bin, runner)

    nontr_cmd = (
        "set -euo pipefail; "
        f"{shlex.quote(bedtools_bin)} intersect -header -sorted -v -g {shlex.quote(str(genome_fai))} "
        f"-a {shlex.quote(str(input_vcf))} -b {shlex.quote(str(tr_bed))} | "
        f"{shlex.quote(bcftools_bin)} view -Oz -o {shlex.quote(str(nontr_output_vcf))} -"
    )
    runner.run(nontr_cmd, shell=True)
    _index_vcf(nontr_output_vcf, bcftools_bin, runner)


def _bcftools_view(
    *,
    bcftools_bin: str,
    input_vcf: Path,
    output_vcf: Path,
    runner: CommandRunner,
    include_bed: Path | None = None,
    exclude_bed: Path | None = None,
    pass_only: bool = False,
    svtypes: tuple[str, ...] = (),
    include_expr: str | None = None,
) -> None:
    if include_bed is not None and exclude_bed is not None:
        raise ValueError("Use at most one of include_bed/exclude_bed per bcftools view call")

    cmd: list[str] = [bcftools_bin, "view"]
    if include_bed is not None:
        cmd.extend(["-R", str(include_bed)])
    if exclude_bed is not None:
        cmd.extend(["-T", f"^{exclude_bed}"])
    if pass_only:
        cmd.extend(["-f", "PASS"])

    exprs: list[str] = []
    svtype_expr = _build_svtype_expr(svtypes)
    if svtype_expr:
        exprs.append(f"({svtype_expr})")
    if include_expr:
        exprs.append(f"({include_expr})")
    if exprs:
        cmd.extend(["-i", " && ".join(exprs)])

    cmd.extend(["-Oz", "-o", str(output_vcf), str(input_vcf)])
    runner.run(cmd)
    _index_vcf(output_vcf, bcftools_bin, runner)


def _apply_confidence_filters(
    *,
    bcftools_bin: str,
    input_vcf: Path,
    output_vcf: Path,
    runner: CommandRunner,
    svtypes: tuple[str, ...],
    pass_only: bool,
    include_bed: Path | None,
    exclude_bed: Path | None,
    include_expr: str | None,
) -> None:
    with tempfile.TemporaryDirectory(prefix="sv_filter_", dir=str(output_vcf.parent)) as td:
        tmpdir = Path(td)
        current = input_vcf

        if exclude_bed is not None:
            next_vcf = tmpdir / "step1.exclude.vcf.gz"
            _bcftools_view(
                bcftools_bin=bcftools_bin,
                input_vcf=current,
                output_vcf=next_vcf,
                runner=runner,
                exclude_bed=exclude_bed,
            )
            current = next_vcf

        if include_bed is not None:
            next_vcf = tmpdir / "step2.include.vcf.gz"
            _bcftools_view(
                bcftools_bin=bcftools_bin,
                input_vcf=current,
                output_vcf=next_vcf,
                runner=runner,
                include_bed=include_bed,
            )
            current = next_vcf

        _bcftools_view(
            bcftools_bin=bcftools_bin,
            input_vcf=current,
            output_vcf=output_vcf,
            runner=runner,
            pass_only=pass_only,
            svtypes=svtypes,
            include_expr=include_expr,
        )


def _merge_vcfs(
    *,
    bcftools_bin: str,
    tr_vcf: Path,
    nontr_vcf: Path,
    merged_vcf: Path,
    runner: CommandRunner,
) -> None:
    merged_vcf.parent.mkdir(parents=True, exist_ok=True)
    tmp_concat = merged_vcf.parent / "merged.concat.unsorted.vcf.gz"
    runner.run(
        [
            bcftools_bin,
            "concat",
            "-a",
            "-Oz",
            "-o",
            str(tmp_concat),
            str(tr_vcf),
            str(nontr_vcf),
        ]
    )
    _index_vcf(tmp_concat, bcftools_bin, runner)

    runner.run(
        [
            bcftools_bin,
            "sort",
            "-Oz",
            "-o",
            str(merged_vcf),
            str(tmp_concat),
        ]
    )
    _index_vcf(merged_vcf, bcftools_bin, runner)

    if not runner.dry_run:
        tmp_concat.unlink(missing_ok=True)
        tmp_index = Path(str(tmp_concat) + ".tbi")
        tmp_index.unlink(missing_ok=True)


def build_parser(
    *,
    prog: str = "sniffcell-discover-sv",
    add_help: bool = True,
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Push-button dual-pass Sniffles SV discovery before SniffCell: "
            "TR-focused pass + non-TR pass + optional Q100 confidence filter."
        ),
        add_help=add_help,
    )
    parser.add_argument("-i", "--input", required=True, help="Input BAM/CRAM.")
    parser.add_argument("-r", "--reference", required=True, help="Reference FASTA used for mapping.")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory.")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Threads passed to Sniffles.")

    parser.add_argument("--tr-bed", default=DEFAULT_TR_BED, help="Tandem repeat BED for TR pass.")
    parser.add_argument(
        "--q100-bed",
        default=DEFAULT_Q100_BED,
        help="Q100 benchmark BED used as final confidence include filter.",
    )
    parser.add_argument(
        "--disable-q100-filter",
        action="store_true",
        help="Disable final Q100 BED include filtering.",
    )

    parser.add_argument(
        "--run-mode",
        choices=["two-pass", "one-pass"],
        default="two-pass",
        help=(
            "SV discovery strategy: 'two-pass' runs Sniffles separately on TR and non-TR; "
            "'one-pass' runs Sniffles once genome-wide then splits TR/non-TR."
        ),
    )
    parser.add_argument(
        "--onepass-splitter",
        choices=["bedtools", "bcftools"],
        default="bedtools",
        help="Tool used to split one-pass raw VCF into TR/non-TR sets.",
    )
    parser.add_argument(
        "--no-onepass-reuse-raw",
        action="store_true",
        help="Force rerunning one-pass Sniffles even if raw VCF already exists.",
    )
    parser.add_argument(
        "--nontr-region-mode",
        choices=["auto", "complement", "postfilter"],
        default="auto",
        help=(
            "How to restrict non-TR pass: "
            "'complement' builds genome\\TR BED and runs Sniffles only there; "
            "'postfilter' runs genome-wide then excludes TR with bcftools."
        ),
    )
    parser.add_argument(
        "--force-recompute-complement",
        action="store_true",
        help="Regenerate non-TR complement BED even if cached.",
    )

    parser.add_argument("--svtypes", default="INS,DEL", help="Comma-separated SVTYPEs to keep. Empty keeps all.")
    parser.add_argument(
        "--keep-nonpass",
        action="store_true",
        help="Keep non-PASS calls (default keeps PASS only for confidence).",
    )
    parser.add_argument(
        "--max-vaf",
        type=float,
        default=None,
        help="Optional maximum allowed VAF/AF in confident calls (e.g., 0.2 keeps VAF<=0.2).",
    )
    parser.add_argument(
        "--max-del-len",
        type=int,
        default=None,
        help="Optional maximum absolute DEL size in bp; larger DELs are excluded.",
    )
    parser.add_argument(
        "--min-support",
        type=int,
        default=None,
        help="Optional minimum INFO/SUPPORT to keep.",
    )
    parser.add_argument(
        "--max-support",
        type=int,
        default=None,
        help="Optional maximum INFO/SUPPORT to keep.",
    )

    parser.add_argument("--tr-mosaic-af-min", type=float, default=0.001, help="Mosaic AF minimum for TR pass.")
    parser.add_argument(
        "--nontr-mosaic-af-min",
        type=float,
        default=0.01,
        help="Mosaic AF minimum for non-TR pass.",
    )

    parser.add_argument(
        "--tr-with-qc",
        action="store_true",
        help="Enable Sniffles QC in TR pass (default TR pass uses --no-qc).",
    )
    parser.add_argument(
        "--nontr-no-qc",
        action="store_true",
        help="Disable Sniffles QC in non-TR pass (default keeps QC enabled).",
    )

    parser.add_argument(
        "--tr-exclude-germline",
        action="store_true",
        help="Do not include germline calls in TR pass.",
    )
    parser.add_argument(
        "--nontr-exclude-germline",
        action="store_true",
        help="Do not include germline calls in non-TR pass.",
    )

    parser.add_argument(
        "--tr-cluster-merge-len",
        type=float,
        default=0.01,
        help="TR pass value for Sniffles --cluster-merge-len.",
    )
    parser.add_argument(
        "--extra-tr-args",
        default="",
        help="Extra raw arguments appended to TR Sniffles command.",
    )
    parser.add_argument(
        "--extra-nontr-args",
        default="",
        help="Extra raw arguments appended to non-TR Sniffles command.",
    )

    parser.add_argument("--no-merge", action="store_true", help="Skip merged VCF generation.")

    parser.add_argument("--sniffles-bin", default="sniffles")
    parser.add_argument("--bcftools-bin", default="bcftools")
    parser.add_argument("--bedtools-bin", default="bedtools")
    parser.add_argument("--samtools-bin", default="samtools")

    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing.")
    parser.add_argument("--log-level", default="INFO", help="DEBUG/INFO/WARNING/ERROR")
    return parser


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    return build_parser().parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    sniffles_bin = _ensure_executable(args.sniffles_bin)
    bcftools_bin = _ensure_executable(args.bcftools_bin)
    samtools_bin = _ensure_executable(args.samtools_bin)

    nontr_mode = args.nontr_region_mode
    bedtools_present = shutil.which(args.bedtools_bin) is not None
    if args.run_mode == "two-pass":
        if nontr_mode == "auto":
            nontr_mode = "complement" if bedtools_present else "postfilter"
        if nontr_mode == "complement":
            if not bedtools_present:
                raise FileNotFoundError(
                    f"Requested non-TR mode '{nontr_mode}' requires bedtools, but it is not available: {args.bedtools_bin}"
                )
            bedtools_bin = _ensure_executable(args.bedtools_bin)
        else:
            bedtools_bin = args.bedtools_bin
    else:
        # one-pass mode does not build non-TR complement BEDs.
        nontr_mode = f"one-pass-split-{args.onepass_splitter}"
        if args.onepass_splitter == "bedtools":
            if not bedtools_present:
                raise FileNotFoundError(
                    f"Requested one-pass splitter '{args.onepass_splitter}' requires bedtools, "
                    f"but it is not available: {args.bedtools_bin}"
                )
            bedtools_bin = _ensure_executable(args.bedtools_bin)
        else:
            bedtools_bin = args.bedtools_bin

    input_bam = _expand_path(args.input)
    reference = _expand_path(args.reference)
    tr_bed = _expand_path(args.tr_bed)
    q100_bed = _expand_path(args.q100_bed)

    output_dir = _expand_path(args.output_dir)
    tr_dir = output_dir / "tr"
    nontr_dir = output_dir / "non_tr"
    merged_dir = output_dir / "merged"
    regions_dir = output_dir / "regions"
    for path in [output_dir, tr_dir, nontr_dir, merged_dir, regions_dir]:
        path.mkdir(parents=True, exist_ok=True)

    if not input_bam.exists():
        raise FileNotFoundError(f"Input BAM/CRAM not found: {input_bam}")
    if not reference.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {reference}")
    if not tr_bed.exists():
        raise FileNotFoundError(f"TR BED not found: {tr_bed}")
    if not args.disable_q100_filter and not q100_bed.exists():
        raise FileNotFoundError(f"Q100 BED not found: {q100_bed}")

    svtypes = _parse_svtypes(args.svtypes)
    pass_only = not args.keep_nonpass
    if args.max_vaf is not None and (args.max_vaf < 0 or args.max_vaf > 1):
        raise ValueError(f"--max-vaf must be in [0,1], got: {args.max_vaf}")
    if args.max_del_len is not None and args.max_del_len <= 0:
        raise ValueError(f"--max-del-len must be > 0, got: {args.max_del_len}")
    if args.min_support is not None and args.min_support <= 0:
        raise ValueError(f"--min-support must be > 0, got: {args.min_support}")
    if args.max_support is not None and args.max_support <= 0:
        raise ValueError(f"--max-support must be > 0, got: {args.max_support}")
    if (
        args.min_support is not None
        and args.max_support is not None
        and args.min_support > args.max_support
    ):
        raise ValueError(
            f"--min-support ({args.min_support}) cannot exceed --max-support ({args.max_support})"
        )
    additional_filter_expr = _build_additional_filter_expr(
        max_vaf=args.max_vaf,
        max_del_len=args.max_del_len,
        min_support=args.min_support,
        max_support=args.max_support,
    )
    tr_no_qc = not args.tr_with_qc
    tr_include_germline = not args.tr_exclude_germline
    nontr_include_germline = not args.nontr_exclude_germline
    runner = CommandRunner(dry_run=args.dry_run)

    logging.info("Run mode resolved to: %s", args.run_mode)
    if args.run_mode == "one-pass":
        logging.info("One-pass splitter resolved to: %s", args.onepass_splitter)
    logging.info("Non-TR mode resolved to: %s", nontr_mode)
    if args.disable_q100_filter:
        logging.info("Q100 filter disabled")
    else:
        logging.info("Q100 include BED enabled: %s", q100_bed)
    if additional_filter_expr:
        logging.info("Additional confidence expression enabled: %s", additional_filter_expr)

    tr_raw_vcf = tr_dir / "sniffles.tr.raw.vcf.gz"
    tr_raw_snf = tr_dir / "sniffles.tr.raw.snf"
    tr_conf_vcf = tr_dir / "sniffles.tr.confident.vcf.gz"

    nontr_raw_vcf = nontr_dir / "sniffles.nontr.raw.vcf.gz"
    nontr_raw_snf = nontr_dir / "sniffles.nontr.raw.snf"
    nontr_conf_vcf = nontr_dir / "sniffles.nontr.confident.vcf.gz"

    nontr_regions_bed: Path | None = None
    if args.run_mode == "two-pass":
        tr_cmd = _build_sniffles_command(
            sniffles_bin=sniffles_bin,
            input_bam=input_bam,
            reference_fa=reference,
            output_vcf=tr_raw_vcf,
            output_snf=tr_raw_snf,
            threads=args.threads,
            mosaic_af_min=args.tr_mosaic_af_min,
            include_germline=tr_include_germline,
            no_qc=tr_no_qc,
            regions_bed=tr_bed,
            cluster_merge_len=args.tr_cluster_merge_len,
            extra_args=tuple(shlex.split(args.extra_tr_args)),
        )
        runner.run(tr_cmd)
        _index_vcf(tr_raw_vcf, bcftools_bin, runner)

        _apply_confidence_filters(
            bcftools_bin=bcftools_bin,
            input_vcf=tr_raw_vcf,
            output_vcf=tr_conf_vcf,
            runner=runner,
            svtypes=svtypes,
            pass_only=pass_only,
            include_bed=None if args.disable_q100_filter else q100_bed,
            exclude_bed=None,
            include_expr=additional_filter_expr,
        )

        nontr_exclude_tr_bed: Path | None = None
        if nontr_mode == "complement":
            genome_fai = _ensure_fasta_index(reference, samtools_bin, runner)
            nontr_regions_bed = regions_dir / "non_tr.complement.bed"
            _write_complement_bed(
                tr_bed=tr_bed,
                genome_fai=genome_fai,
                output_bed=nontr_regions_bed,
                bedtools_bin=bedtools_bin,
                runner=runner,
                force_recompute=args.force_recompute_complement,
            )
        else:
            nontr_exclude_tr_bed = tr_bed

        nontr_cmd = _build_sniffles_command(
            sniffles_bin=sniffles_bin,
            input_bam=input_bam,
            reference_fa=reference,
            output_vcf=nontr_raw_vcf,
            output_snf=nontr_raw_snf,
            threads=args.threads,
            mosaic_af_min=args.nontr_mosaic_af_min,
            include_germline=nontr_include_germline,
            no_qc=args.nontr_no_qc,
            regions_bed=nontr_regions_bed,
            extra_args=tuple(shlex.split(args.extra_nontr_args)),
        )
        runner.run(nontr_cmd)
        _index_vcf(nontr_raw_vcf, bcftools_bin, runner)

        _apply_confidence_filters(
            bcftools_bin=bcftools_bin,
            input_vcf=nontr_raw_vcf,
            output_vcf=nontr_conf_vcf,
            runner=runner,
            svtypes=svtypes,
            pass_only=pass_only,
            include_bed=None if args.disable_q100_filter else q100_bed,
            exclude_bed=nontr_exclude_tr_bed,
            include_expr=additional_filter_expr,
        )
    else:
        raw_dir = output_dir / "raw"
        raw_dir.mkdir(parents=True, exist_ok=True)
        raw_vcf = raw_dir / "sniffles.raw.vcf.gz"
        raw_snf = raw_dir / "sniffles.raw.snf"
        reuse_onepass_raw = not args.no_onepass_reuse_raw

        # one-pass mode uses the most sensitive AF floor and TR cluster merge setting.
        onepass_af_min = min(args.tr_mosaic_af_min, args.nontr_mosaic_af_min)
        # one-pass cannot differ by TR/non-TR germline flag; include if either pass would include.
        onepass_include_germline = tr_include_germline or nontr_include_germline
        # one-pass cannot differ by TR/non-TR QC flag; default to TR behavior for sensitivity.
        onepass_no_qc = tr_no_qc

        if reuse_onepass_raw and raw_vcf.exists():
            logging.info("Reusing existing one-pass raw VCF: %s", raw_vcf)
            _index_vcf(raw_vcf, bcftools_bin, runner)
        else:
            onepass_cmd = _build_sniffles_command(
                sniffles_bin=sniffles_bin,
                input_bam=input_bam,
                reference_fa=reference,
                output_vcf=raw_vcf,
                output_snf=raw_snf,
                threads=args.threads,
                mosaic_af_min=onepass_af_min,
                include_germline=onepass_include_germline,
                no_qc=onepass_no_qc,
                regions_bed=None,
                cluster_merge_len=args.tr_cluster_merge_len,
                extra_args=tuple(shlex.split(args.extra_tr_args)),
            )
            runner.run(onepass_cmd)
            _index_vcf(raw_vcf, bcftools_bin, runner)

        if args.onepass_splitter == "bedtools":
            genome_fai = _ensure_fasta_index(reference, samtools_bin, runner)
            _split_onepass_vcf_with_bedtools(
                bedtools_bin=bedtools_bin,
                bcftools_bin=bcftools_bin,
                input_vcf=raw_vcf,
                tr_bed=tr_bed,
                tr_output_vcf=tr_raw_vcf,
                nontr_output_vcf=nontr_raw_vcf,
                genome_fai=genome_fai,
                runner=runner,
            )
        else:
            _bcftools_view(
                bcftools_bin=bcftools_bin,
                input_vcf=raw_vcf,
                output_vcf=tr_raw_vcf,
                runner=runner,
                include_bed=tr_bed,
            )
            _bcftools_view(
                bcftools_bin=bcftools_bin,
                input_vcf=raw_vcf,
                output_vcf=nontr_raw_vcf,
                runner=runner,
                exclude_bed=tr_bed,
            )

        _apply_confidence_filters(
            bcftools_bin=bcftools_bin,
            input_vcf=tr_raw_vcf,
            output_vcf=tr_conf_vcf,
            runner=runner,
            svtypes=svtypes,
            pass_only=pass_only,
            include_bed=None if args.disable_q100_filter else q100_bed,
            exclude_bed=None,
            include_expr=additional_filter_expr,
        )

        _apply_confidence_filters(
            bcftools_bin=bcftools_bin,
            input_vcf=nontr_raw_vcf,
            output_vcf=nontr_conf_vcf,
            runner=runner,
            svtypes=svtypes,
            pass_only=pass_only,
            include_bed=None if args.disable_q100_filter else q100_bed,
            exclude_bed=None,
            include_expr=additional_filter_expr,
        )

    merged_vcf: Path | None = None
    if not args.no_merge:
        merged_vcf = merged_dir / "sniffles.confident.merged.vcf.gz"
        _merge_vcfs(
            bcftools_bin=bcftools_bin,
            tr_vcf=tr_conf_vcf,
            nontr_vcf=nontr_conf_vcf,
            merged_vcf=merged_vcf,
            runner=runner,
        )

    manifest = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "bam": str(input_bam),
            "reference": str(reference),
            "tr_bed": str(tr_bed),
            "q100_bed": None if args.disable_q100_filter else str(q100_bed),
        },
        "settings": {
            "run_mode": args.run_mode,
            "onepass_splitter": args.onepass_splitter,
            "onepass_reuse_raw": reuse_onepass_raw if args.run_mode == "one-pass" else None,
            "threads": args.threads,
            "svtypes": list(svtypes),
            "pass_only": pass_only,
            "nontr_region_mode": nontr_mode,
            "tr_mosaic_af_min": args.tr_mosaic_af_min,
            "nontr_mosaic_af_min": args.nontr_mosaic_af_min,
            "max_vaf": args.max_vaf,
            "max_del_len": args.max_del_len,
            "min_support": args.min_support,
            "max_support": args.max_support,
            "additional_filter_expr": additional_filter_expr,
            "tr_no_qc": tr_no_qc,
            "nontr_no_qc": args.nontr_no_qc,
            "tr_include_germline": tr_include_germline,
            "nontr_include_germline": nontr_include_germline,
            "tr_cluster_merge_len": args.tr_cluster_merge_len,
        },
        "outputs": {
            "tr_raw_vcf": str(tr_raw_vcf),
            "tr_confident_vcf": str(tr_conf_vcf),
            "nontr_raw_vcf": str(nontr_raw_vcf),
            "nontr_confident_vcf": str(nontr_conf_vcf),
            "nontr_regions_bed": str(nontr_regions_bed) if nontr_regions_bed else None,
            "merged_confident_vcf": str(merged_vcf) if merged_vcf else None,
        },
        "dry_run": args.dry_run,
        "executed_commands": runner.commands,
    }
    manifest_path = output_dir / "sv_discovery_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"TR confident VCF: {tr_conf_vcf}")
    print(f"non-TR confident VCF: {nontr_conf_vcf}")
    if merged_vcf:
        print(f"Merged confident VCF: {merged_vcf}")
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
