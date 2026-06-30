import dataclasses
import json
import os
import shlex
import stat
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

from sniffcell.parse_args import parse_args
from sniffcell.discover.discover import (
    DEFAULT_STAGE_ORDER,
    GROUP_SCOPED_STAGES,
    RunContext,
    VALID_STAGES,
    _apply_platform_tr_substitution,
    _build_context,
    _build_recursive_cli,
    _clear_force_rerun_state,
    _discover_groups,
    _parse_stages,
    _select_groups,
    _render_slurm,
    _render_submit_script,
    _sanitize_token,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_exec(path: Path) -> None:
    path.write_text("#!/bin/bash\nexit 0\n")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def _make_python_exec(path: Path, body: str) -> None:
    path.write_text(f"#!/usr/bin/env python3\n{body}", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def _build_fake_runtime_env(root: Path) -> tuple[Path, Path, Path, dict[str, str]]:
    deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
    tool_dir = Path(next(iter(tool_paths.values()))).parent

    _make_python_exec(
        tool_dir / "sniffles",
        """
import sys
from pathlib import Path
args = sys.argv[1:]
vcf = Path(args[args.index("--vcf") + 1])
snf = Path(args[args.index("--snf") + 1])
vcf.parent.mkdir(parents=True, exist_ok=True)
vcf.write_text("##fileformat=VCFv4.2\\n", encoding="utf-8")
snf.parent.mkdir(parents=True, exist_ok=True)
snf.write_text("snf\\n", encoding="utf-8")
""",
    )

    _make_python_exec(
        tool_dir / "bcftools",
        """
import sys
from pathlib import Path
args = sys.argv[1:]
sub = args[0]
if sub == "index":
    target = Path(args[-1])
    target.parent.mkdir(parents=True, exist_ok=True)
    Path(str(target) + ".tbi").write_text("tbi\\n", encoding="utf-8")
elif sub in {"view", "sort", "concat"}:
    out = Path(args[args.index("-o") + 1])
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("##fileformat=VCFv4.2\\n", encoding="utf-8")
else:
    raise SystemExit(f"unsupported bcftools subcommand: {sub}")
""",
    )

    _make_python_exec(
        tool_dir / "kanpig",
        """
import sys
from pathlib import Path
args = sys.argv[1:]
raw_vcf = Path(args[args.index("-o") + 1] if "-o" in args else args[args.index("--out") + 1])
rnames = Path(args[args.index("--rnames") + 1])
samples = []
for idx, arg in enumerate(args):
    if arg in {"--sample", "-s"} and idx + 1 < len(args):
        samples.append(args[idx + 1])
raw_vcf.parent.mkdir(parents=True, exist_ok=True)
if len(samples) >= 2:
    raw_vcf.write_text(
        "##fileformat=VCFv4.2\\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"GT\\">\\n"
        "##FORMAT=<ID=FT,Number=1,Type=Integer,Description=\\"FT\\">\\n"
        "##FORMAT=<ID=SQ,Number=1,Type=Integer,Description=\\"SQ\\">\\n"
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\\"GQ\\">\\n"
        "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\\"PS\\">\\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\\"DP\\">\\n"
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\\"AD\\">\\n"
        "##FORMAT=<ID=KS,Number=.,Type=Integer,Description=\\"KS\\">\\n"
        f"#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t{samples[0]}\\t{samples[1]}\\n"
        "chr1\\t10\\tsv1\\tA\\t<DEL>\\t60\\tPASS\\t.\\tGT:FT:SQ:GQ:PS:DP:AD:KS\\t0|1:0:50:50:1:10:0,10:99\\t0|0:0:10:10:1:12:12,0:99\\n",
        encoding="utf-8",
    )
else:
    raw_vcf.write_text("##fileformat=VCFv4.2\\n", encoding="utf-8")
rnames.parent.mkdir(parents=True, exist_ok=True)
rnames.write_text("sv_id\\tread_name\\n", encoding="utf-8")
""",
    )

    _make_python_exec(
        tool_dir / "truvari",
        """
import sys
from pathlib import Path
args = sys.argv[1:]
if args[0] != "collapse":
    raise SystemExit("expected truvari collapse")
raw_vcf = Path(args[args.index("-o") + 1])
removed_vcf = Path(args[args.index("-c") + 1])
raw_vcf.parent.mkdir(parents=True, exist_ok=True)
raw_vcf.write_text("##fileformat=VCFv4.2\\n", encoding="utf-8")
removed_vcf.parent.mkdir(parents=True, exist_ok=True)
removed_vcf.write_text("##fileformat=VCFv4.2\\n", encoding="utf-8")
""",
    )

    _make_python_exec(
        tool_dir / "medaka",
        """
import sys
from pathlib import Path
args = sys.argv[1:]
if args[0] != "tandem":
    raise SystemExit("expected medaka tandem")
stage_dir = Path(args[-1])
output_vcf = stage_dir / "medaka_to_ref.TR.vcf"
output_vcf.parent.mkdir(parents=True, exist_ok=True)
output_vcf.write_text("##fileformat=VCFv4.2\\n", encoding="utf-8")
# Trimmed spanning reads consumed by the read-length TR scan. Make Neuron's
# reads clearly longer than Oligodendrocyte's so the scan emits one call.
group = stage_dir.name.replace(".medaka", "")
lengths = [400, 410, 420] if group == "Neuron" else [100, 105]
fasta = stage_dir / "trimmed_reads.fasta"
with fasta.open("w", encoding="utf-8") as handle:
    for idx, length in enumerate(lengths):
        handle.write(f">r{idx}_chr1_0_10_pad_0_0_fwd_hap1_phased-set1_ploidy2\\n")
        handle.write("A" * length + "\\n")
""",
    )

    _make_python_exec(
        tool_dir / "tdb",
        """
import sys
from pathlib import Path
args = sys.argv[1:]
sub = args[0]
output = Path(args[args.index("-o") + 1])
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(f"{sub}\\n", encoding="utf-8")
""",
    )

    _make_python_exec(
        tool_dir / "modkit",
        """
import sys
from pathlib import Path
args = sys.argv[1:]
if args[0] != "pileup":
    raise SystemExit("expected modkit pileup")
output = Path(args[2])
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text("bedmethyl\\n", encoding="utf-8")
""",
    )

    _make_python_exec(
        tool_dir / "tabix",
        """
import sys
from pathlib import Path
target = Path(sys.argv[-1])
target.parent.mkdir(parents=True, exist_ok=True)
Path(str(target) + ".tbi").write_text("tbi\\n", encoding="utf-8")
""",
    )

    _make_python_exec(
        tool_dir / "run_clair3.sh",
        """
import sys
import gzip
from pathlib import Path
args = sys.argv[1:]
output_dir = None
for arg in args:
    if arg.startswith("--output="):
        output_dir = Path(arg.split("=", 1)[1])
        break
if output_dir is None:
    raise SystemExit("missing --output")
output_dir.mkdir(parents=True, exist_ok=True)
merge_output = output_dir / "merge_output.vcf.gz"
pileup_output = output_dir / "pileup.vcf.gz"
group_name = output_dir.name

if group_name == "Neuron":
    records = [
        "chr1\\t10\\t.\\tA\\tG\\t60\\tPASS\\tP\\tGT:GQ:DP:AD:AF\\t0/1:30:8:4,4:0.5",
        "chr1\\t20\\t.\\tC\\t.\\t20\\tRefCall\\tP\\tGT:GQ:DP:AD:AF\\t0/0:20:18:18:1.0",
    ]
else:
    records = [
        "chr1\\t10\\t.\\tA\\t.\\t20\\tRefCall\\tP\\tGT:GQ:DP:AD:AF\\t0/0:20:19:19:1.0",
        "chr1\\t20\\t.\\tC\\tT\\t60\\tPASS\\tP\\tGT:GQ:DP:AD:AF\\t0/1:30:9:4,5:0.5",
    ]

for path in (merge_output, pileup_output):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\\n")
        handle.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\\"GT\\">\\n")
        handle.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\\"GQ\\">\\n")
        handle.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\\"DP\\">\\n")
        handle.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\\"AD\\">\\n")
        handle.write("##FORMAT=<ID=AF,Number=1,Type=Float,Description=\\"AF\\">\\n")
        handle.write("#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tSAMPLE\\n")
        for record in records:
            handle.write(record + "\\n")
""",
    )

    return deconv_dir, ref, tr_bed, tool_paths


def _build_minimal_env(root: Path) -> tuple[Path, Path, Path, dict[str, str]]:
    """
    Create a minimal on-disk environment for _build_context tests.
    Returns (deconv_dir, ref, tr_bed, tool_paths).
    """
    sample_dir = root / "sample1"
    deconv_dir = sample_dir / "deconv"
    split_dir = deconv_dir / "deconv_requested_group_splits"
    split_dir.mkdir(parents=True)
    for group_name in ("Neuron", "Oligodendrocyte"):
        (split_dir / f"{group_name}.bam").write_text("")
        (split_dir / f"{group_name}.bam.bai").write_text("")
    (split_dir / "requested_group_splits.tsv").write_text(
        "requested_group\tbam_path\tread_summary_path\n"
        f"Neuron\t{split_dir / 'Neuron.bam'}\t{split_dir / 'Neuron.read_summary.tsv'}\n"
        f"Oligodendrocyte\t{split_dir / 'Oligodendrocyte.bam'}\t"
        f"{split_dir / 'Oligodendrocyte.read_summary.tsv'}\n"
    )
    ref = root / "ref.fa"
    tr_bed = root / "tr.bed"
    ref.write_text(">chr1\nA\n")
    tr_bed.write_text("chr1\t0\t10\n")
    tool_dir = root / "bin"
    tool_dir.mkdir()
    tool_paths: dict[str, str] = {}
    for tool_name in (
        "sniffles", "bcftools", "kanpig", "truvari", "medaka",
        "tdb", "modkit", "tabix", "trgt", "samtools", "run_clair3.sh",
    ):
        p = tool_dir / tool_name
        _make_exec(p)
        tool_paths[tool_name] = str(p)
    return deconv_dir, ref, tr_bed, tool_paths


def _base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths) -> list[str]:
    return [
        "discover",
        "--deconv-dir", str(deconv_dir),
        "--reference", str(ref),
        "--tr-bed", str(tr_bed),
        "--sex", "male",
        "--scheduler", "slurm",
        "--run-id", "testrun",
        "--sniffles-bin", tool_paths["sniffles"],
        "--bcftools-bin", tool_paths["bcftools"],
        "--kanpig-bin", tool_paths["kanpig"],
        "--truvari-bin", tool_paths["truvari"],
        "--medaka-bin", tool_paths["medaka"],
        "--tdb-bin", tool_paths["tdb"],
        "--modkit-bin", tool_paths["modkit"],
        "--tabix-bin", tool_paths["tabix"],
        "--trgt-bin", tool_paths["trgt"],
        "--samtools-bin", tool_paths["samtools"],
        "--clair3-bin", tool_paths["run_clair3.sh"],
        "--clair3-model-path", "/tmp/clair3_model",
    ]


def _base_local_argv(deconv_dir, ref, tr_bed, tool_paths) -> list[str]:
    argv = _base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths)
    scheduler_idx = argv.index("--scheduler")
    argv[scheduler_idx + 1] = "local"
    return argv


# ---------------------------------------------------------------------------
# parse_args tests
# ---------------------------------------------------------------------------

class TestDiscoverParseArgs(unittest.TestCase):

    def test_discover_accepts_core_arguments(self):
        args = parse_args(
            [
                "discover",
                "--deconv-dir", "/tmp/sample/deconv",
                "--reference", "/tmp/ref.fa",
                "--tr-bed", "/tmp/tr.bed",
                "--sex", "male",
            ]
        )
        self.assertEqual(args.command, "discover")
        self.assertEqual(args.scheduler, "local")
        self.assertIsNone(args.slurm_account)
        self.assertEqual(args.threads, 16)
        self.assertEqual(args.mods_mode, "separate")

    def test_threads_default_is_16(self):
        args = parse_args([
            "discover",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "female",
        ])
        self.assertEqual(args.threads, 16)

    def test_threads_override(self):
        args = parse_args([
            "discover",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "female",
            "--threads", "32",
        ])
        self.assertEqual(args.threads, 32)

    def test_slurm_account_optional_none_by_default(self):
        args = parse_args([
            "discover",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "male",
        ])
        self.assertIsNone(args.slurm_account)

    def test_slurm_account_can_be_set(self):
        args = parse_args([
            "discover",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "male",
            "--slurm-account", "mylab",
        ])
        self.assertEqual(args.slurm_account, "mylab")

    def test_hidden_args_still_parseable(self):
        """--run-id, --stages, and --groups are suppressed but functional."""
        args = parse_args([
            "discover",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "male",
            "--run-id", "myrun",
            "--stages", "sv",
            "--groups", "A,B",
        ])
        self.assertEqual(args.run_id, "myrun")
        self.assertEqual(args.stages, "sv")
        self.assertEqual(args.groups, "A,B")

    def test_removed_args_are_absent(self):
        """Removed arguments should not be attributes on args."""
        args = parse_args([
            "discover",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "male",
        ])
        self.assertFalse(hasattr(args, "keep_going"))
        self.assertFalse(hasattr(args, "post_tdb_script"))
        self.assertFalse(hasattr(args, "post_tdb_args"))
        self.assertFalse(hasattr(args, "slurm_partition"))
        self.assertFalse(hasattr(args, "submit"))
        # Per-tool thread args should not exist
        for attr in (
            "sniffles_threads", "kanpig_threads", "medaka_workers",
            "tdb_merge_threads", "modkit_threads", "clair3_threads",
        ):
            self.assertFalse(hasattr(args, attr), f"unexpected attr: {attr}")

    def test_kanpig_passonly_default_true(self):
        args = parse_args([
            "discover",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "male",
        ])
        self.assertTrue(args.kanpig_passonly)

    def test_mods_mode_default_separate(self):
        args = parse_args([
            "discover",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "male",
        ])
        self.assertEqual(args.mods_mode, "separate")

    def test_platform_default(self):
        args = parse_args([
            "discover",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "male",
        ])
        self.assertEqual(args.platform, "ont")

    def test_platform_flag_and_deprecated_alias(self):
        base = [
            "discover",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "male",
        ]
        args = parse_args(base + ["--platform", "hifi"])
        self.assertEqual(args.platform, "hifi")
        # --clair3-platform remains a deprecated alias for the same dest.
        alias = parse_args(base + ["--clair3-platform", "hifi"])
        self.assertEqual(alias.platform, "hifi")

    def test_medaka_phasing_accepts_abpoa(self):
        args = parse_args([
            "discover",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "male",
            "--medaka-phasing", "abpoa",
        ])
        self.assertEqual(args.medaka_phasing, "abpoa")

    def test_discover_tools_run_accepts_discover_arguments(self):
        args = parse_args([
            "discover",
            "tools",
            "run",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "female",
            "--threads", "32",
        ])
        self.assertEqual(args.command, "discover")
        self.assertEqual(args.discover_section, "tools")
        self.assertEqual(args.discover_tools_command, "run")

    def test_discover_ctprocessing_snv_accepts_arguments(self):
        args = parse_args([
            "discover",
            "ctprocessing",
            "snv",
            "--split-dir", "/tmp/split",
            "--groups", "A,B",
            "--group-a-gvcf", "/tmp/a.vcf.gz",
            "--group-b-gvcf", "/tmp/b.vcf.gz",
        ])
        self.assertEqual(args.command, "discover")
        self.assertEqual(args.discover_section, "ctprocessing")
        self.assertEqual(args.discover_ctprocessing_command, "snv")
        self.assertEqual(args.min_dp, 5)

    def test_discover_tools_check_parses_envcheck_arguments(self):
        args = parse_args([
            "discover",
            "tools",
            "check",
            "--stages", "sv,mods",
            "--json",
        ])
        self.assertEqual(args.command, "discover")
        self.assertEqual(args.discover_section, "tools")
        self.assertEqual(args.discover_tools_command, "check")
        self.assertEqual(args.stages, "sv,mods")
        self.assertTrue(args.json)

    def test_discover_tools_sv_parses(self):
        args = parse_args([
            "discover",
            "tools",
            "sv",
            "-i", "/tmp/in.bam",
            "-r", "/tmp/ref.fa",
            "-o", "/tmp/out",
        ])
        self.assertEqual(args.command, "discover")
        self.assertEqual(args.discover_section, "tools")
        self.assertEqual(args.discover_tools_command, "sv")
        self.assertEqual(args.input, "/tmp/in.bam")

    def test_discover_ctprocessing_sv_parses(self):
        args = parse_args([
            "discover",
            "ctprocessing",
            "sv",
            "--split-dir", "/tmp/splits",
            "--reference", "/tmp/ref.fa",
            "--groups", "A,B",
        ])
        self.assertEqual(args.command, "discover")
        self.assertEqual(args.discover_section, "ctprocessing")
        self.assertEqual(args.discover_ctprocessing_command, "sv")
        self.assertEqual(args.groups, "A,B")

    def test_discover_ctprocessing_tr_parses(self):
        args = parse_args([
            "discover",
            "ctprocessing",
            "tr",
            "--split-dir", "/tmp/splits",
            "--groups", "A,B",
            "--discover-run-id", "run1",
        ])
        self.assertEqual(args.command, "discover")
        self.assertEqual(args.discover_section, "ctprocessing")
        self.assertEqual(args.discover_ctprocessing_command, "tr")
        self.assertEqual(args.groups, "A,B")
        self.assertEqual(args.discover_run_id, "run1")

    def test_discover_ctprocessing_groups_are_optional(self):
        cases = [
            (
                ["snv", "--split-dir", "/tmp/splits", "--group-a-gvcf", "/tmp/a.vcf.gz", "--group-b-gvcf", "/tmp/b.vcf.gz"],
                "snv",
            ),
            (["sv", "--split-dir", "/tmp/splits", "--reference", "/tmp/ref.fa"], "sv"),
            (["tr", "--split-dir", "/tmp/splits"], "tr"),
        ]
        for argv_tail, command in cases:
            with self.subTest(command=command):
                args = parse_args(["discover", "ctprocessing", *argv_tail])
                self.assertEqual(args.discover_ctprocessing_command, command)
                self.assertIsNone(args.groups)

    def test_discover_ctprocessing_harmonize_parses(self):
        args = parse_args([
            "discover",
            "ctprocessing",
            "harmonize",
            "--tr-bed", "/tmp/tr.tsv",
            "--sv-bed", "/tmp/sv.tsv",
            "--group-a-label", "sample.A",
            "--group-b-label", "sample.B",
            "--output", "/tmp/harmonized.tsv",
        ])
        self.assertEqual(args.command, "discover")
        self.assertEqual(args.discover_section, "ctprocessing")
        self.assertEqual(args.discover_ctprocessing_command, "harmonize")
        self.assertEqual(args.output, "/tmp/harmonized.tsv")


# ---------------------------------------------------------------------------
# _parse_stages tests
# ---------------------------------------------------------------------------

class TestParseStages(unittest.TestCase):

    def test_none_returns_default_order(self):
        self.assertEqual(_parse_stages(None), DEFAULT_STAGE_ORDER)
        self.assertNotIn("clair3", DEFAULT_STAGE_ORDER)

    def test_empty_string_returns_default_order(self):
        self.assertEqual(_parse_stages(""), DEFAULT_STAGE_ORDER)

    def test_sv_alias(self):
        stages = _parse_stages("sv")
        self.assertEqual(stages, ("sniffles", "sniffles_filter", "kanpig", "collapse"))

    def test_tdb_alias(self):
        stages = _parse_stages("tdb")
        self.assertEqual(stages, ("tdb_create", "tdb_merge"))

    def test_snv_alias(self):
        stages = _parse_stages("snv")
        self.assertEqual(stages, ("clair3",))

    def test_mods_alias(self):
        stages = _parse_stages("mods")
        self.assertEqual(stages, ("modkit",))

    def test_all_alias_returns_full_order(self):
        stages = _parse_stages("all")
        self.assertEqual(stages, DEFAULT_STAGE_ORDER)
        self.assertNotIn("clair3", stages)

    def test_sv_tdb_modkit_combo(self):
        stages = _parse_stages("sv,tdb,modkit")
        self.assertEqual(
            stages,
            (
                "sniffles",
                "sniffles_filter",
                "kanpig",
                "collapse",
                "tdb_create",
                "tdb_merge",
                "modkit",
            ),
        )

    def test_individual_stage(self):
        self.assertEqual(_parse_stages("sniffles"), ("sniffles",))
        self.assertEqual(_parse_stages("modkit"), ("modkit",))
        self.assertEqual(_parse_stages("clair3"), ("clair3",))

    def test_deduplication(self):
        stages = _parse_stages("sniffles,sniffles")
        self.assertEqual(stages, ("sniffles",))

    def test_invalid_stage_raises(self):
        with self.assertRaises(ValueError):
            _parse_stages("foobar")

    def test_post_tdb_not_a_valid_stage(self):
        with self.assertRaises(ValueError):
            _parse_stages("post_tdb")

    def test_all_stages_in_default_order_are_valid(self):
        for stage in DEFAULT_STAGE_ORDER:
            result = _parse_stages(stage)
            self.assertIn(stage, result)

    def test_post_tdb_not_in_default_order(self):
        self.assertNotIn("post_tdb", DEFAULT_STAGE_ORDER)

    def test_trgt_is_valid_stage(self):
        self.assertIn("trgt", VALID_STAGES)
        self.assertIn("trgt", DEFAULT_STAGE_ORDER)
        self.assertEqual(_parse_stages("trgt"), ("trgt",))

    def test_medaka_is_valid_but_not_default(self):
        self.assertIn("medaka", VALID_STAGES)
        self.assertNotIn("medaka", DEFAULT_STAGE_ORDER)
        self.assertEqual(_parse_stages("medaka"), ("medaka",))


class TestPlatformTrSubstitution(unittest.TestCase):

    def test_ont_platform_keeps_user_requested_medaka(self):
        stages = ("sniffles", "medaka", "tdb_create")
        self.assertEqual(_apply_platform_tr_substitution(stages, "ont"), stages)

    def test_hifi_platform_keeps_user_requested_medaka(self):
        stages = ("sniffles", "medaka", "tdb_create")
        self.assertEqual(
            _apply_platform_tr_substitution(stages, "hifi"),
            stages,
        )

    def test_hifi_case_insensitive_noop(self):
        self.assertEqual(
            _apply_platform_tr_substitution(("medaka",), "HiFi"),
            ("medaka",),
        )

    def test_no_medaka_means_no_substitution(self):
        stages = ("sniffles", "kanpig")
        self.assertEqual(_apply_platform_tr_substitution(stages, "hifi"), stages)

    def test_explicit_trgt_already_present_is_left_alone(self):
        stages = ("sniffles", "medaka", "trgt")
        self.assertEqual(_apply_platform_tr_substitution(stages, "hifi"), stages)

    def test_none_platform_is_noop(self):
        stages = ("medaka",)
        self.assertEqual(_apply_platform_tr_substitution(stages, None), stages)


# ---------------------------------------------------------------------------
# _sanitize_token tests
# ---------------------------------------------------------------------------

class TestSanitizeToken(unittest.TestCase):

    def test_plain_name(self):
        from sniffcell.discover.discover import _sanitize_token
        self.assertEqual(_sanitize_token("Neuron"), "Neuron")
        self.assertEqual(_sanitize_token("Oligodendrocyte"), "Oligodendrocyte")

    def test_spaces_become_underscore(self):
        from sniffcell.discover.discover import _sanitize_token
        self.assertEqual(_sanitize_token("T Cell"), "T_Cell")

    def test_special_chars_become_underscore(self):
        from sniffcell.discover.discover import _sanitize_token
        self.assertEqual(_sanitize_token("B/NK Cell"), "B_NK_Cell")

    def test_multiple_underscores_collapsed(self):
        from sniffcell.discover.discover import _sanitize_token
        self.assertEqual(_sanitize_token("A__B"), "A_B")

    def test_leading_trailing_stripped(self):
        from sniffcell.discover.discover import _sanitize_token
        result = _sanitize_token("!Neuron!")
        self.assertFalse(result.startswith("_"))
        self.assertFalse(result.endswith("_"))

    def test_empty_falls_back_to_group(self):
        from sniffcell.discover.discover import _sanitize_token
        self.assertEqual(_sanitize_token(""), "group")
        self.assertEqual(_sanitize_token("!!!"), "group")

    def test_alphanumeric_passthrough(self):
        from sniffcell.discover.discover import _sanitize_token
        self.assertEqual(_sanitize_token("abc123"), "abc123")


# ---------------------------------------------------------------------------
# _select_groups tests
# ---------------------------------------------------------------------------

class TestSelectGroups(unittest.TestCase):

    def _make_groups(self, names):
        from sniffcell.discover.discover import SplitGroup
        return [SplitGroup(name=n, bam_path=f"/tmp/{n}.bam", bai_path=f"/tmp/{n}.bam.bai") for n in names]

    def test_no_text_two_groups_returns_both(self):
        groups = self._make_groups(["A", "B"])
        result = _select_groups(groups, None)
        self.assertEqual(result, ["A", "B"])

    def test_no_text_non_two_raises(self):
        groups = self._make_groups(["A", "B", "C"])
        with self.assertRaises(ValueError):
            _select_groups(groups, None)

    def test_explicit_groups_selection(self):
        groups = self._make_groups(["A", "B", "C"])
        result = _select_groups(groups, "A,C")
        self.assertEqual(result, ["A", "C"])

    def test_unknown_group_raises(self):
        groups = self._make_groups(["A", "B"])
        with self.assertRaises(ValueError):
            _select_groups(groups, "A,X")

    def test_empty_text_raises(self):
        groups = self._make_groups(["A", "B"])
        with self.assertRaises(ValueError):
            _select_groups(groups, "  ,  ")


# ---------------------------------------------------------------------------
# RunContext shape tests
# ---------------------------------------------------------------------------

class TestRunContextShape(unittest.TestCase):

    def test_no_keep_going_field(self):
        fields = {f.name for f in dataclasses.fields(RunContext)}
        self.assertNotIn("keep_going", fields)

    def test_has_params_and_tool_paths(self):
        fields = {f.name for f in dataclasses.fields(RunContext)}
        self.assertIn("params", fields)
        self.assertIn("tool_paths", fields)

    def test_has_scheduler_dry_run_force(self):
        fields = {f.name for f in dataclasses.fields(RunContext)}
        for f in ("scheduler", "dry_run", "force", "rerun_failed"):
            self.assertIn(f, fields)


# ---------------------------------------------------------------------------
# _build_context integration tests
# ---------------------------------------------------------------------------

class TestBuildContext(unittest.TestCase):

    def test_build_context_sample_id_inferred(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            self.assertEqual(ctx.sample_id, "sample1")

    def test_build_context_selected_groups_two(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            self.assertEqual(ctx.selected_groups, ["Neuron", "Oligodendrocyte"])

    def test_build_context_threads_in_params(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            argv = _base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths) + ["--threads", "24"]
            args = parse_args(argv)
            ctx = _build_context(args)
            self.assertEqual(ctx.params["threads"], 24)

    def test_build_context_no_keep_going_in_params(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            self.assertNotIn("keep_going", ctx.params)

    def test_build_context_slurm_account_in_params(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            argv = _base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths) + ["--slurm-account", "labacct"]
            args = parse_args(argv)
            ctx = _build_context(args)
            self.assertEqual(ctx.params["slurm_account"], "labacct")

    def test_build_context_mods_mode_propagated(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            argv = _base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths) + ["--mods-mode", "combined"]
            args = parse_args(argv)
            ctx = _build_context(args)
            self.assertEqual(ctx.params["mods_mode"], "combined")

    def test_build_context_medaka_phasing_propagated(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            argv = _base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths) + ["--medaka-phasing", "abpoa"]
            args = parse_args(argv)
            ctx = _build_context(args)
            self.assertEqual(ctx.params["medaka_phasing"], "abpoa")


# ---------------------------------------------------------------------------
# _build_recursive_cli tests
# ---------------------------------------------------------------------------

class TestBuildRecursiveCli(unittest.TestCase):

    def _get_ctx(self, root, extra_argv=None):
        deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
        argv = _base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths)
        if extra_argv:
            argv += extra_argv
        args = parse_args(argv)
        return _build_context(args)

    def test_contains_scheduler_local(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td))
            cli = _build_recursive_cli(ctx, "sniffles")
            self.assertIn("--scheduler", cli)
            idx = cli.index("--scheduler")
            self.assertEqual(cli[idx + 1], "local")

    def test_contains_threads(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td), ["--threads", "8"])
            cli = _build_recursive_cli(ctx, "sniffles")
            self.assertIn("--threads", cli)
            idx = cli.index("--threads")
            self.assertEqual(cli[idx + 1], "8")

    def test_no_per_tool_thread_flags(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td))
            cli = _build_recursive_cli(ctx, "sniffles")
            for flag in (
                "--sniffles-threads", "--kanpig-threads", "--medaka-workers",
                "--tdb-merge-threads", "--modkit-threads",
            ):
                self.assertNotIn(flag, cli)

    def test_contains_run_id(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td))
            cli = _build_recursive_cli(ctx, "sniffles")
            self.assertIn("--run-id", cli)
            idx = cli.index("--run-id")
            self.assertEqual(cli[idx + 1], ctx.run_id)

    def test_groups_included_when_provided(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td))
            cli = _build_recursive_cli(ctx, "sniffles", groups="Neuron,Oligodendrocyte")
            self.assertIn("--groups", cli)
            idx = cli.index("--groups")
            self.assertEqual(cli[idx + 1], "Neuron,Oligodendrocyte")

    def test_mods_mode_forwarded(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td), ["--mods-mode", "combined"])
            cli = _build_recursive_cli(ctx, "modkit")
            self.assertIn("--mods-mode", cli)
            idx = cli.index("--mods-mode")
            self.assertEqual(cli[idx + 1], "combined")

    def test_kanpig_sample_name_template_forwarded(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td), [
                "--kanpig-sample-name-template", "{sample_id}__{group}"
            ])
            cli = _build_recursive_cli(ctx, "kanpig")
            self.assertIn("--kanpig-sample-name-template", cli)
            idx = cli.index("--kanpig-sample-name-template")
            self.assertEqual(cli[idx + 1], "{sample_id}__{group}")

    def test_medaka_sample_name_template_forwarded(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td), [
                "--medaka-sample-name-template", "{sample_id}_{group}"
            ])
            cli = _build_recursive_cli(ctx, "medaka")
            self.assertIn("--medaka-sample-name-template", cli)
            idx = cli.index("--medaka-sample-name-template")
            self.assertEqual(cli[idx + 1], "{sample_id}_{group}")

    def test_medaka_phasing_forwarded_when_set(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td), ["--medaka-phasing", "abpoa"])
            cli = _build_recursive_cli(ctx, "medaka")
            self.assertIn("--medaka-phasing", cli)
            idx = cli.index("--medaka-phasing")
            self.assertEqual(cli[idx + 1], "abpoa")

    def test_medaka_phasing_absent_when_not_set(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td))
            cli = _build_recursive_cli(ctx, "medaka")
            self.assertNotIn("--medaka-phasing", cli)

    def test_tdb_create_force_forwarded_when_set(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td), ["--tdb-create-force"])
            cli = _build_recursive_cli(ctx, "tdb_create")
            self.assertIn("--tdb-create-force", cli)

    def test_tdb_create_force_absent_when_not_set(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td))
            cli = _build_recursive_cli(ctx, "tdb_create")
            self.assertNotIn("--tdb-create-force", cli)

    def test_collapse_use_default_kanpig_not_forwarded(self):
        """collapse_use=kanpig is the default so not forwarded (saves CLI clutter)."""
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td))
            cli = _build_recursive_cli(ctx, "collapse")
            self.assertNotIn("--collapse-use", cli)

    def test_collapse_use_non_default_forwarded(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td), ["--collapse-use", "sniffles"])
            cli = _build_recursive_cli(ctx, "collapse")
            self.assertIn("--collapse-use", cli)
            idx = cli.index("--collapse-use")
            self.assertEqual(cli[idx + 1], "sniffles")

    def test_force_forwarded_when_set(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td), ["--force"])
            cli = _build_recursive_cli(ctx, "sniffles")
            self.assertIn("--force", cli)

    def test_force_absent_when_not_set(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td))
            cli = _build_recursive_cli(ctx, "sniffles")
            self.assertNotIn("--force", cli)

    def test_force_rerun_clears_selected_task_state(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(
                _base_local_argv(deconv_dir, ref, tr_bed, tool_paths)
                + ["--stages", "collapse,medaka", "--force"]
            )
            ctx = _build_context(args)
            ctx.run_root.mkdir(parents=True, exist_ok=True)
            ctx.status_dir.mkdir(parents=True, exist_ok=True)
            (ctx.status_dir / "collapse.done.json").write_text("{}", encoding="utf-8")
            (ctx.status_dir / "collapse_inputs.done.json").write_text("{}", encoding="utf-8")
            (ctx.status_dir / "sv_post_processing.done.json").write_text("{}", encoding="utf-8")
            (ctx.status_dir / "medaka.Neuron.done.json").write_text("{}", encoding="utf-8")
            (ctx.status_dir / "medaka.Oligodendrocyte.failed.json").write_text("{}", encoding="utf-8")
            (ctx.status_dir / "sniffles.Neuron.done.json").write_text("{}", encoding="utf-8")
            (ctx.status_dir / "discover_status.json").write_text(
                json.dumps(
                    {
                        "collapse": {"state": "completed"},
                        "collapse_inputs": {"state": "completed"},
                        "sv_post_processing": {"state": "completed"},
                        "medaka.Neuron": {"state": "running"},
                        "medaka.Oligodendrocyte": {"state": "failed"},
                        "sniffles.Neuron": {"state": "completed"},
                    }
                ),
                encoding="utf-8",
            )

            _clear_force_rerun_state(ctx)

            status_payload = json.loads((ctx.status_dir / "discover_status.json").read_text(encoding="utf-8"))
            self.assertNotIn("collapse", status_payload)
            self.assertNotIn("collapse_inputs", status_payload)
            self.assertNotIn("sv_post_processing", status_payload)
            self.assertNotIn("medaka.Neuron", status_payload)
            self.assertNotIn("medaka.Oligodendrocyte", status_payload)
            self.assertIn("sniffles.Neuron", status_payload)
            self.assertFalse((ctx.status_dir / "collapse.done.json").exists())
            self.assertFalse((ctx.status_dir / "collapse_inputs.done.json").exists())
            self.assertFalse((ctx.status_dir / "sv_post_processing.done.json").exists())
            self.assertFalse((ctx.status_dir / "medaka.Neuron.done.json").exists())
            self.assertFalse((ctx.status_dir / "medaka.Oligodendrocyte.failed.json").exists())
            self.assertTrue((ctx.status_dir / "sniffles.Neuron.done.json").exists())

    def test_clair3_model_path_forwarded_when_set(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td))
            cli = _build_recursive_cli(ctx, "clair3")
            self.assertIn("--clair3-model-path", cli)
            idx = cli.index("--clair3-model-path")
            self.assertEqual(cli[idx + 1], "/tmp/clair3_model")

    def test_clair3_command_enables_gvcf_by_default(self):
        from sniffcell.discover.discover import _run_clair3

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            argv = [
                "discover",
                "--deconv-dir", str(deconv_dir),
                "--reference", str(ref),
                "--tr-bed", str(tr_bed),
                "--sex", "male",
                "--scheduler", "local",
                "--dry-run",
                "--sniffles-bin", tool_paths["sniffles"],
                "--bcftools-bin", tool_paths["bcftools"],
                "--kanpig-bin", tool_paths["kanpig"],
                "--truvari-bin", tool_paths["truvari"],
                "--medaka-bin", tool_paths["medaka"],
                "--tdb-bin", tool_paths["tdb"],
                "--modkit-bin", tool_paths["modkit"],
                "--tabix-bin", tool_paths["tabix"],
                "--clair3-bin", tool_paths["run_clair3.sh"],
                "--clair3-model-path", "/tmp/clair3_model",
            ]
            args = parse_args(argv)
            ctx = _build_context(args)
            _run_clair3(ctx, "Neuron")
            cmd_text = (ctx.commands_dir / "clair3.Neuron.command.txt").read_text()
            self.assertIn("--gvcf", cmd_text)

    def test_force_rerun_clears_snv_post_processing_when_clair3_selected(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(
                _base_local_argv(deconv_dir, ref, tr_bed, tool_paths)
                + ["--stages", "clair3", "--force"]
            )
            ctx = _build_context(args)
            ctx.run_root.mkdir(parents=True, exist_ok=True)
            ctx.status_dir.mkdir(parents=True, exist_ok=True)
            (ctx.status_dir / "clair3.Neuron.done.json").write_text("{}", encoding="utf-8")
            (ctx.status_dir / "clair3.Oligodendrocyte.done.json").write_text("{}", encoding="utf-8")
            (ctx.status_dir / "snv_post_processing.done.json").write_text("{}", encoding="utf-8")
            (ctx.status_dir / "discover_status.json").write_text(
                json.dumps(
                    {
                        "clair3.Neuron": {"state": "completed"},
                        "clair3.Oligodendrocyte": {"state": "completed"},
                        "snv_post_processing": {"state": "completed"},
                    }
                ),
                encoding="utf-8",
            )

            _clear_force_rerun_state(ctx)

            status_payload = json.loads((ctx.status_dir / "discover_status.json").read_text(encoding="utf-8"))
            self.assertNotIn("clair3.Neuron", status_payload)
            self.assertNotIn("clair3.Oligodendrocyte", status_payload)
            self.assertNotIn("snv_post_processing", status_payload)
            self.assertFalse((ctx.status_dir / "clair3.Neuron.done.json").exists())
            self.assertFalse((ctx.status_dir / "clair3.Oligodendrocyte.done.json").exists())
            self.assertFalse((ctx.status_dir / "snv_post_processing.done.json").exists())


class TestPostprocessContextAndSlurm(unittest.TestCase):

    def test_build_context_discovers_two_groups_and_renders_slurm(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))

            ctx = _build_context(args)
            self.assertEqual(ctx.sample_id, "sample1")
            self.assertEqual(ctx.selected_groups, ["Neuron", "Oligodendrocyte"])

            _render_slurm(ctx)
            sniffles_script = ctx.slurm_dir / "sniffles.array.sbatch.sh"
            self.assertTrue(sniffles_script.exists())
            script_text = sniffles_script.read_text()
            self.assertNotIn("#SBATCH --partition=", script_text)
            self.assertNotIn("#SBATCH --account=", script_text)
            self.assertIn('${GROUP_NAME}', script_text)
            self.assertIn("--scheduler local", script_text)

            submit_script = ctx.slurm_dir / "submit_pipeline.sh"
            self.assertTrue(submit_script.exists())
            submit_text = submit_script.read_text()
            self.assertIn('PARTITION=', submit_text)
            self.assertIn('_acct()', submit_text)

    def test_all_group_scoped_scripts_are_generated(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            _render_slurm(ctx)
            # Only the group-scoped stages that were actually requested for this
            # run should have sbatch scripts emitted.
            for stage in GROUP_SCOPED_STAGES & set(ctx.stages):
                script = ctx.slurm_dir / f"{stage}.array.sbatch.sh"
                self.assertTrue(script.exists(), f"Missing SLURM script for {stage}")

    def test_sample_scoped_scripts_are_generated(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            _render_slurm(ctx)
            for script_name in ("collapse.sbatch.sh", "tr_post_processing.sbatch.sh"):
                self.assertTrue(
                    (ctx.slurm_dir / script_name).exists(),
                    f"Missing SLURM script: {script_name}",
                )
            self.assertFalse((ctx.slurm_dir / "tdb_merge.sbatch.sh").exists())
            self.assertFalse((ctx.slurm_dir / "snv_post_processing.sbatch.sh").exists())

    def test_submit_script_has_concurrent_first_tier(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            _render_slurm(ctx)
            submit_text = (ctx.slurm_dir / "submit_pipeline.sh").read_text()
            # First-tier independent jobs
            self.assertIn("SNIFFLES_JID", submit_text)
            self.assertIn("TRGT_JID", submit_text)
            self.assertIn("MODKIT_JID", submit_text)
            self.assertNotIn("MEDAKA_JID", submit_text)
            self.assertNotIn("CLAIR3_JID", submit_text)
            self.assertNotIn("SNV_POST_PROCESSING_JID", submit_text)

    def test_submit_script_has_dependent_chain(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            _render_slurm(ctx)
            submit_text = (ctx.slurm_dir / "submit_pipeline.sh").read_text()
            # sniffles_filter depends on sniffles
            self.assertIn("SNIFFLES_FILTER_JID", submit_text)
            self.assertIn("afterok:${SNIFFLES_JID}", submit_text)
            # kanpig depends on sniffles_filter
            self.assertIn("afterok:${SNIFFLES_FILTER_JID}", submit_text)
            # tr_post_processing depends on trgt
            self.assertIn("afterok:${TRGT_JID}", submit_text)
            # Medaka/TDB are opt-in, not part of default submission
            self.assertNotIn("afterok:${MEDAKA_JID}", submit_text)
            self.assertNotIn("afterok:${TDB_CREATE_JID}", submit_text)

    def test_submit_script_account_prefilled_when_given(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            argv = _base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths) + ["--slurm-account", "mylab"]
            args = parse_args(argv)
            ctx = _build_context(args)
            _render_slurm(ctx)
            submit_text = (ctx.slurm_dir / "submit_pipeline.sh").read_text()
            self.assertIn('ACCOUNT="mylab"', submit_text)

    def test_submit_script_account_empty_when_not_given(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            _render_slurm(ctx)
            submit_text = (ctx.slurm_dir / "submit_pipeline.sh").read_text()
            self.assertIn('ACCOUNT=""', submit_text)

    def test_submit_script_is_executable(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            _render_slurm(ctx)
            submit_script = ctx.slurm_dir / "submit_pipeline.sh"
            self.assertTrue(os.access(submit_script, os.X_OK))

    def test_sbatch_scripts_are_executable(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            _render_slurm(ctx)
            for script in ctx.slurm_dir.glob("*.sbatch.sh"):
                self.assertTrue(os.access(script, os.X_OK), f"Not executable: {script.name}")

    def test_groups_tsv_written(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            _render_slurm(ctx)
            groups_tsv = ctx.slurm_dir / "groups.tsv"
            self.assertTrue(groups_tsv.exists())
            content = groups_tsv.read_text()
            self.assertIn("Neuron", content)
            self.assertIn("Oligodendrocyte", content)

    def test_sniffles_command_omits_sample_name_when_unsupported(self):
        """_run_sniffles should not pass --sample-name when the sniffles binary lacks that option."""
        from sniffcell.discover.discover import _run_sniffles
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            # Local dry-run so _run_sniffles writes the command without executing
            argv = [
                "discover",
                "--deconv-dir", str(deconv_dir),
                "--reference", str(ref),
                "--tr-bed", str(tr_bed),
                "--sex", "male",
                "--scheduler", "local",
                "--dry-run",
                "--sniffles-bin", tool_paths["sniffles"],
                "--bcftools-bin", tool_paths["bcftools"],
                "--kanpig-bin", tool_paths["kanpig"],
                "--truvari-bin", tool_paths["truvari"],
                "--medaka-bin", tool_paths["medaka"],
                "--tdb-bin", tool_paths["tdb"],
                "--modkit-bin", tool_paths["modkit"],
                "--tabix-bin", tool_paths["tabix"],
                "--clair3-bin", tool_paths["run_clair3.sh"],
            ]
            args = parse_args(argv)
            ctx = _build_context(args)
            _run_sniffles(ctx, "Neuron")
            cmd_file = ctx.commands_dir / "sniffles.Neuron.command.txt"
            cmd_text = cmd_file.read_text()
            self.assertNotIn("--sample-name", cmd_text)
            self.assertNotIn("sample1_Neuron", cmd_text)

    def test_sniffles_command_contains_sample_name_when_supported(self):
        """_run_sniffles should pass --sample-name when the sniffles binary supports it."""
        from sniffcell.discover.discover import _run_sniffles
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            argv = [
                "discover",
                "--deconv-dir", str(deconv_dir),
                "--reference", str(ref),
                "--tr-bed", str(tr_bed),
                "--sex", "male",
                "--scheduler", "local",
                "--dry-run",
                "--sniffles-bin", tool_paths["sniffles"],
                "--bcftools-bin", tool_paths["bcftools"],
                "--kanpig-bin", tool_paths["kanpig"],
                "--truvari-bin", tool_paths["truvari"],
                "--medaka-bin", tool_paths["medaka"],
                "--tdb-bin", tool_paths["tdb"],
                "--modkit-bin", tool_paths["modkit"],
                "--tabix-bin", tool_paths["tabix"],
                "--clair3-bin", tool_paths["run_clair3.sh"],
            ]
            args = parse_args(argv)
            ctx = _build_context(args)
            with patch("sniffcell.discover.discover._sniffles_supports_sample_name", return_value=True):
                _run_sniffles(ctx, "Neuron")
            cmd_file = ctx.commands_dir / "sniffles.Neuron.command.txt"
            cmd_text = cmd_file.read_text()
            self.assertIn("--sample-name", cmd_text)
            self.assertIn("sample1_Neuron", cmd_text)

    def test_medaka_command_contains_phasing_when_requested(self):
        from sniffcell.discover.discover import _run_medaka
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            argv = [
                "discover",
                "--deconv-dir", str(deconv_dir),
                "--reference", str(ref),
                "--tr-bed", str(tr_bed),
                "--sex", "male",
                "--scheduler", "local",
                "--dry-run",
                "--medaka-phasing", "abpoa",
                "--sniffles-bin", tool_paths["sniffles"],
                "--bcftools-bin", tool_paths["bcftools"],
                "--kanpig-bin", tool_paths["kanpig"],
                "--truvari-bin", tool_paths["truvari"],
                "--medaka-bin", tool_paths["medaka"],
                "--tdb-bin", tool_paths["tdb"],
                "--modkit-bin", tool_paths["modkit"],
                "--tabix-bin", tool_paths["tabix"],
                "--clair3-bin", tool_paths["run_clair3.sh"],
            ]
            args = parse_args(argv)
            ctx = _build_context(args)
            _run_medaka(ctx, "Neuron")
            cmd_file = ctx.commands_dir / "medaka.Neuron.command.txt"
            cmd_text = cmd_file.read_text()
            self.assertIn("--phasing", cmd_text)
            self.assertIn("abpoa", cmd_text)

    def test_sniffles_script_contains_threads_param(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            argv = _base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths) + ["--threads", "8"]
            args = parse_args(argv)
            ctx = _build_context(args)
            _render_slurm(ctx)
            script_text = (ctx.slurm_dir / "sniffles.array.sbatch.sh").read_text()
            self.assertIn("--threads 8", script_text)
            self.assertIn("#SBATCH --cpus-per-task=8", script_text)

    def test_mods_mode_forwarded_in_modkit_script(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            argv = _base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths) + ["--mods-mode", "combined"]
            args = parse_args(argv)
            ctx = _build_context(args)
            _render_slurm(ctx)
            script_text = (ctx.slurm_dir / "modkit.array.sbatch.sh").read_text()
            self.assertIn("--mods-mode combined", script_text)


class TestTrPostProcessingScan(unittest.TestCase):

    def test_direction_excess_passes_when_topk_clear_margin(self):
        from sniffcell.discover.tr_post_processing import _direction_excess

        excess = _direction_excess(
            [1300, 1280, 1260, 900, 850],
            [1000, 990, 980, 970, 960],
            margin_bp=100,
            min_supporting_reads=3,
            min_total_reads=5,
        )
        self.assertEqual(excess, 300)  # 1300 - 1000

    def test_direction_excess_requires_min_supporting_reads(self):
        from sniffcell.discover.tr_post_processing import _direction_excess

        # only one read clears baseline_max(1000) + margin(100)
        self.assertIsNone(
            _direction_excess(
                [1300, 1050, 1040, 1030, 1020],
                [1000, 990, 980, 970, 960],
                margin_bp=100,
                min_supporting_reads=3,
                min_total_reads=5,
            )
        )

    def test_direction_excess_skips_when_baseline_empty(self):
        from sniffcell.discover.tr_post_processing import _direction_excess

        self.assertIsNone(
            _direction_excess([2000, 1900, 1800, 1700, 1600], [], margin_bp=100, min_supporting_reads=3, min_total_reads=5)
        )

    def test_direction_excess_requires_total_reads_in_both_groups(self):
        from sniffcell.discover.tr_post_processing import _direction_excess

        self.assertIsNone(
            _direction_excess(
                [1300, 1280, 1260, 1240],
                [1000, 990, 980, 970, 960],
                margin_bp=100,
                min_supporting_reads=3,
                min_total_reads=5,
            )
        )
        self.assertIsNone(
            _direction_excess(
                [1300, 1280, 1260, 1240, 1220],
                [1000, 990, 980, 970],
                margin_bp=100,
                min_supporting_reads=3,
                min_total_reads=5,
            )
        )

    def test_assign_tier_strong_vs_supportive(self):
        from sniffcell.discover.tr_post_processing import _assign_tier

        self.assertEqual(_assign_tier(n_change_support_reads=3, n_baseline_reads=2, min_supporting_reads=2), "strong")
        self.assertEqual(_assign_tier(n_change_support_reads=2, n_baseline_reads=2, min_supporting_reads=2), "supportive")
        self.assertEqual(_assign_tier(n_change_support_reads=5, n_baseline_reads=1, min_supporting_reads=2), "supportive")

    def test_scan_loci_calls_expansion_and_picks_change_group(self):
        from sniffcell.discover.tr_post_processing import _scan_loci

        a_loci = {("chr1", 100, 200): [("a0", 1300), ("a1", 1280), ("a2", 1260), ("a3", 1240), ("a4", 850)]}
        b_loci = {("chr1", 100, 200): [("b0", 1000), ("b1", 1010), ("b2", 990), ("b3", 980), ("b4", 970)]}
        rows = _scan_loci(
            a_loci, b_loci,
            sample_a_label="s.Neuron", sample_b_label="s.Oligo",
            margin_bp=100, min_supporting_reads=3, min_total_reads=5,
        )
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["change_group"], "s.Neuron")
        self.assertEqual(row["baseline_group"], "s.Oligo")
        self.assertEqual(row["change_type"], "expansion")
        self.assertEqual(row["baseline_max_bp"], 1010)
        self.assertEqual(row["change_max_bp"], 1300)
        self.assertEqual(row["change_length_bp"], 290)
        self.assertEqual(row["tr_tier"], "strong")
        self.assertTrue(row["tr_pass_for_harmonized"])

    def test_scan_loci_skips_within_margin_and_empty_baseline(self):
        from sniffcell.discover.tr_post_processing import _scan_loci

        a_loci = {
            ("chr3", 0, 10): [("a0", 1050), ("a1", 1040), ("a2", 1030), ("a3", 1020), ("a4", 1010)],
            ("chr4", 0, 10): [("a0", 2000), ("a1", 1900), ("a2", 1800), ("a3", 1700), ("a4", 1600)],
        }
        b_loci = {("chr3", 0, 10): [("b0", 1000), ("b1", 1005), ("b2", 990), ("b3", 980), ("b4", 970)]}
        rows = _scan_loci(
            a_loci, b_loci,
            sample_a_label="s.Neuron", sample_b_label="s.Oligo",
            margin_bp=100, min_supporting_reads=3, min_total_reads=5,
        )
        self.assertEqual(rows, [])

    def test_parse_fasta_lengths_sums_wrapped_sequence(self):
        from sniffcell.discover.tr_post_processing import _parse_fasta_lengths

        with tempfile.TemporaryDirectory() as td:
            fasta = Path(td) / "trimmed_reads.fasta"
            fasta.write_text(
                ">r0_chr1_100_200_pad_0_0_fwd_hap1_phased-set1_ploidy2\n"
                "AAAA\nAAA\n"  # 7 bp across two lines
                ">r1_chr1_100_200_pad_0_0_fwd_hap2_phased-set1_ploidy2\n"
                "AAAAA\n",
                encoding="utf-8",
            )
            loci = _parse_fasta_lengths(fasta)
            lengths = sorted(length for _, length in loci[("chr1", 100, 200)])
            self.assertEqual(lengths, [5, 7])


class TestTrPostProcessingMain(unittest.TestCase):

    def _write_fasta(self, path: Path, records: list[tuple[str, int]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as handle:
            for header, length in records:
                handle.write(f">{header}\n")
                handle.write("A" * int(length) + "\n")

    def _records(self, tag: str, region: str, lengths: list[int]) -> list[tuple[str, int]]:
        return [
            (f"{tag}_{idx}_{region}_pad_0_0_fwd_hap1_phased-set1_ploidy2", length)
            for idx, length in enumerate(lengths)
        ]

    def test_tr_resolve_args_infers_groups_from_manifest(self):
        from sniffcell.discover.tr_post_processing import _resolve_args

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, _ref, _tr_bed, _tool_paths = _build_minimal_env(root)
            split_dir = deconv_dir / "deconv_requested_group_splits"
            args = _resolve_args(SimpleNamespace(
                split_dir=str(split_dir),
                groups=None,
                output_dir=str(root / "out"),
                sample_id=None,
                sample_a_label=None,
                sample_b_label=None,
                group_a_fasta=None,
                group_b_fasta=None,
                group_a_spanning_bam=None,
                group_b_spanning_bam=None,
                tr_bed=None,
                trgt_fallback_flank_bp=50,
                margin_bp=100,
                min_supporting_reads=3,
                min_total_reads=5,
                min_motif_size=1,
                skip_plots=True,
            ))
            self.assertEqual((args.group_a, args.group_b), ("Neuron", "Oligodendrocyte"))
            self.assertEqual(args.sample_id, "sample1")
            self.assertEqual(
                args.group_a_fasta,
                split_dir / "medaka_tandem" / "Neuron.medaka" / "trimmed_reads.fasta",
            )

    def test_tr_resolve_args_finds_discover_run_medaka_inputs(self):
        from sniffcell.discover.tr_post_processing import _resolve_args

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, _ref, _tr_bed, _tool_paths = _build_minimal_env(root)
            split_dir = deconv_dir / "deconv_requested_group_splits"
            stale_run = split_dir / "discover" / "newer_without_medaka"
            medaka_dir = split_dir / "discover" / "run_with_medaka" / "medaka_tandem"
            stale_run.mkdir(parents=True)
            group_a_fasta = medaka_dir / "Neuron.medaka" / "trimmed_reads.fasta"
            group_b_fasta = medaka_dir / "Oligodendrocyte.medaka" / "trimmed_reads.fasta"
            self._write_fasta(group_a_fasta, self._records("a", "chr1_100_200", [100]))
            self._write_fasta(group_b_fasta, self._records("b", "chr1_100_200", [100]))

            args = _resolve_args(SimpleNamespace(
                split_dir=str(split_dir),
                groups=None,
                output_dir=str(root / "out"),
                sample_id=None,
                sample_a_label=None,
                sample_b_label=None,
                group_a_fasta=None,
                group_b_fasta=None,
                discover_run_id=None,
                group_a_spanning_bam=None,
                group_b_spanning_bam=None,
                tr_bed=None,
                trgt_fallback_flank_bp=50,
                margin_bp=100,
                min_supporting_reads=3,
                min_total_reads=5,
                min_motif_size=1,
                skip_plots=True,
            ))
            self.assertEqual(args.group_a_fasta, group_a_fasta)
            self.assertEqual(args.group_b_fasta, group_b_fasta)

    def test_tr_post_processing_main_writes_tiered_tr_changes(self):
        from sniffcell.discover.tr_post_processing import tr_post_processing_main

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            split_dir = root / "sample1" / "deconv" / "deconv_requested_group_splits"
            medaka_dir = split_dir / "medaka_tandem"
            group_a_fasta = medaka_dir / "Neuron.medaka" / "trimmed_reads.fasta"
            group_b_fasta = medaka_dir / "Oligodendrocyte.medaka" / "trimmed_reads.fasta"

            a_records: list[tuple[str, int]] = []
            b_records: list[tuple[str, int]] = []
            # chr1: Neuron expanded, many supporting reads -> strong
            a_records += self._records("a_strong", "chr1_100_200", [1300, 1280, 1260, 1240, 1220])
            b_records += self._records("b_strong", "chr1_100_200", [1000, 1005, 1010, 1008, 995])
            # chr2: Oligo expanded, exactly min supporting reads -> supportive
            a_records += self._records("a_supp", "chr2_200_300", [1000, 990, 980, 970, 960])
            b_records += self._records("b_supp", "chr2_200_300", [1250, 1240, 1230, 900, 880])
            # chr3: within margin -> not called
            a_records += self._records("a_neg", "chr3_300_400", [1050, 1040, 1030])
            b_records += self._records("b_neg", "chr3_300_400", [1000, 1005])
            # chr4: baseline has no reads -> skipped
            a_records += self._records("a_skip", "chr4_400_500", [2000, 1900])

            self._write_fasta(group_a_fasta, a_records)
            self._write_fasta(group_b_fasta, b_records)

            summary = tr_post_processing_main(
                [
                    "--split-dir", str(split_dir),
                    "--groups", "Neuron,Oligodendrocyte",
                    "--output-dir", str(root / "out"),
                    "--sample-id", "sample1",
                    "--sample-a-label", "sample1.Neuron",
                    "--sample-b-label", "sample1.Oligodendrocyte",
                    "--group-a-fasta", str(group_a_fasta),
                    "--group-b-fasta", str(group_b_fasta),
                    # Homopolymer placeholder reads; this test covers the length
                    # tiering logic, not the motif filter, so disable it here.
                    "--min-motif-size", "1",
                    "--skip-plots",
                ]
            )

            self.assertEqual(summary["status"], "completed")
            self.assertEqual(summary["n_targets"], 2)
            self.assertEqual(summary["n_tr_strong_rows"], 1)
            self.assertEqual(summary["n_tr_supportive_rows"], 1)
            self.assertEqual(summary["n_tr_weak_rows"], 0)
            self.assertEqual(summary["params"]["margin_bp"], 50)
            self.assertEqual(summary["params"]["min_supporting_reads"], 3)
            self.assertEqual(summary["params"]["min_total_reads"], 5)

            tr_bed = pd.read_csv(root / "out" / "tr_changes.bed.tsv", sep="\t")
            for col in (
                "tr_tier", "tr_pass_for_harmonized", "change_group", "baseline_group",
                "change_type", "change_length_bp", "n_change_support_reads",
                "change_support_read_names",
            ):
                self.assertIn(col, tr_bed.columns)
            self.assertEqual(tr_bed["tr_tier"].tolist(), ["strong", "supportive"])
            strong_row = tr_bed.iloc[0]
            self.assertEqual(strong_row["change_group"], "sample1.Neuron")
            self.assertEqual(strong_row["change_type"], "expansion")
            self.assertEqual(int(strong_row["change_length_bp"]), 290)
            supportive_row = tr_bed.iloc[1]
            self.assertEqual(supportive_row["change_group"], "sample1.Oligodendrocyte")
            self.assertEqual(int(supportive_row["n_change_support_reads"]), 3)

    def test_tr_post_processing_main_margin_and_min_reads_are_parameters(self):
        from sniffcell.discover.tr_post_processing import tr_post_processing_main

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            split_dir = root / "sample1" / "deconv" / "deconv_requested_group_splits"
            medaka_dir = split_dir / "medaka_tandem"
            group_a_fasta = medaka_dir / "Neuron.medaka" / "trimmed_reads.fasta"
            group_b_fasta = medaka_dir / "Oligodendrocyte.medaka" / "trimmed_reads.fasta"

            # Neuron tops out 60 bp over Oligo's longest read: called only when margin <= 50.
            self._write_fasta(group_a_fasta, self._records("a", "chr1_100_200", [1060, 1055, 1052, 900, 890]))
            self._write_fasta(group_b_fasta, self._records("b", "chr1_100_200", [1000, 990, 980, 970, 960]))

            base_argv = [
                "--split-dir", str(split_dir),
                "--groups", "Neuron,Oligodendrocyte",
                "--sample-id", "sample1",
                "--sample-a-label", "sample1.Neuron",
                "--sample-b-label", "sample1.Oligodendrocyte",
                "--group-a-fasta", str(group_a_fasta),
                "--group-b-fasta", str(group_b_fasta),
                # Homopolymer placeholder reads; this test covers the margin /
                # min-reads parameters, not the motif filter, so disable it here.
                "--min-motif-size", "1",
                "--skip-plots",
            ]

            strict = tr_post_processing_main(base_argv + ["--output-dir", str(root / "strict"), "--margin-bp", "100"])
            self.assertEqual(strict["n_targets"], 0)

            loose = tr_post_processing_main(base_argv + ["--output-dir", str(root / "loose"), "--margin-bp", "50"])
            self.assertEqual(loose["n_targets"], 1)

            # Requiring 4 supporting reads drops the call (only 3 reads clear the threshold).
            need3 = tr_post_processing_main(
                base_argv + ["--output-dir", str(root / "need4"), "--margin-bp", "50", "--min-supporting-reads", "4"]
            )
            self.assertEqual(need3["n_targets"], 0)

    def test_tr_post_processing_main_min_motif_size_filters_dinucleotide(self):
        from sniffcell.discover.tr_post_processing import tr_post_processing_main

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            split_dir = root / "sample1" / "deconv" / "deconv_requested_group_splits"
            medaka_dir = split_dir / "medaka_tandem"
            group_a_fasta = medaka_dir / "Neuron.medaka" / "trimmed_reads.fasta"
            group_b_fasta = medaka_dir / "Oligodendrocyte.medaka" / "trimmed_reads.fasta"

            def write(path: Path, records: list[tuple[str, str]]) -> None:
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("w", encoding="utf-8") as handle:
                    for header, seq in records:
                        handle.write(f">{header}\n{seq}\n")

            def at(length: int) -> str:
                return ("AT" * (length // 2 + 1))[:length]

            # chr1: Neuron carries an AT-dinucleotide expansion vs the short Oligo
            # baseline -> a valid length call whose motif is 2 bp.
            a_records = [
                (f"a_{i}_chr1_100_200_pad_0_0_fwd_hap1_phased-set1_ploidy2", at(n))
                for i, n in enumerate([1300, 1280, 1260, 1240, 1220])
            ]
            b_records = [
                (f"b_{i}_chr1_100_200_pad_0_0_fwd_hap1_phased-set1_ploidy2", at(n))
                for i, n in enumerate([1000, 1005, 1010, 1008, 995])
            ]
            write(group_a_fasta, a_records)
            write(group_b_fasta, b_records)

            base = [
                "--split-dir", str(split_dir),
                "--groups", "Neuron,Oligodendrocyte",
                "--sample-id", "sample1",
                "--sample-a-label", "sample1.Neuron",
                "--sample-b-label", "sample1.Oligodendrocyte",
                "--group-a-fasta", str(group_a_fasta),
                "--group-b-fasta", str(group_b_fasta),
                "--skip-plots",
            ]

            # Default (min-motif-size=2) keeps the dinucleotide locus: only
            # homopolymers are dropped by default.
            default = tr_post_processing_main(base + ["--output-dir", str(root / "default")])
            self.assertEqual(default["n_targets"], 1)
            self.assertEqual(default["n_tr_motif_filtered"], 0)
            self.assertEqual(default["params"]["min_motif_size"], 2)

            # Raising the floor to 3 drops the dinucleotide locus.
            strict = tr_post_processing_main(
                base + ["--output-dir", str(root / "strict"), "--min-motif-size", "3"]
            )
            self.assertEqual(strict["n_targets"], 0)
            self.assertEqual(strict["n_tr_motif_filtered"], 1)

            # Disabling the filter (=1) also keeps the call.
            disabled = tr_post_processing_main(
                base + ["--output-dir", str(root / "off"), "--min-motif-size", "1"]
            )
            self.assertEqual(disabled["n_targets"], 1)
            self.assertEqual(disabled["n_tr_motif_filtered"], 0)

    def test_tr_post_processing_main_flags_haplotype_dropout(self):
        from sniffcell.discover.tr_post_processing import tr_post_processing_main

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            split_dir = root / "sample1" / "deconv" / "deconv_requested_group_splits"
            medaka_dir = split_dir / "medaka_tandem"
            group_a_fasta = medaka_dir / "Neuron.medaka" / "trimmed_reads.fasta"
            group_b_fasta = medaka_dir / "Oligodendrocyte.medaka" / "trimmed_reads.fasta"

            def recs(tag: str, region: str, lengths: list[int], hap: int):
                return [
                    (f"{tag}_{i}_{region}_pad_0_0_fwd_hap{hap}_phased-set1_ploidy2", n)
                    for i, n in enumerate(lengths)
                ]

            # chr1: Neuron expansion supported only by hap2 reads while the Oligo
            # baseline has hap1 reads only -> haplotype dropout -> low confidence.
            a = recs("a_drop", "chr1_100_200", [1300, 1280, 1260, 1240, 1220], 2)
            b = recs("b_drop", "chr1_100_200", [1000, 1005, 1010, 1008, 995], 1)
            # chr2: same expansion on hap2 but the baseline also carries hap2 reads
            # -> not a dropout.
            a += recs("a_ok", "chr2_200_300", [1300, 1280, 1260, 1240, 1220], 2)
            b += recs("b_ok", "chr2_200_300", [1000, 1005, 1010, 1008, 995], 2)

            self._write_fasta(group_a_fasta, a)
            self._write_fasta(group_b_fasta, b)

            summary = tr_post_processing_main([
                "--split-dir", str(split_dir),
                "--groups", "Neuron,Oligodendrocyte",
                "--output-dir", str(root / "out"),
                "--sample-id", "sample1",
                "--sample-a-label", "sample1.Neuron",
                "--sample-b-label", "sample1.Oligodendrocyte",
                "--group-a-fasta", str(group_a_fasta),
                "--group-b-fasta", str(group_b_fasta),
                "--min-motif-size", "1",
                "--skip-plots",
            ])

            self.assertEqual(summary["n_targets"], 2)
            self.assertEqual(summary["n_tr_hap_dropout"], 1)

            tr_bed = pd.read_csv(root / "out" / "tr_changes.bed.tsv", sep="\t")
            for col in ("change_support_haps", "baseline_haps", "hap_dropout_low_conf"):
                self.assertIn(col, tr_bed.columns)
            rows = {r["trid"]: r for _, r in tr_bed.iterrows()}

            drop = rows["chr1_100_200"]
            self.assertTrue(bool(drop["hap_dropout_low_conf"]))
            self.assertEqual(str(drop["change_support_haps"]), "2")
            self.assertEqual(str(drop["baseline_haps"]), "1")
            self.assertEqual(drop["tr_tier"], "weak")
            self.assertFalse(bool(drop["tr_pass_for_harmonized"]))

            ok = rows["chr2_200_300"]
            self.assertFalse(bool(ok["hap_dropout_low_conf"]))
            self.assertTrue(bool(ok["tr_pass_for_harmonized"]))

    def test_tr_post_processing_main_skips_when_fasta_missing(self):
        from sniffcell.discover.tr_post_processing import tr_post_processing_main

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            split_dir = root / "sample1" / "deconv" / "deconv_requested_group_splits"
            split_dir.mkdir(parents=True)
            summary = tr_post_processing_main(
                [
                    "--split-dir", str(split_dir),
                    "--groups", "Neuron,Oligodendrocyte",
                    "--output-dir", str(root / "out"),
                    "--sample-id", "sample1",
                    "--group-a-fasta", str(root / "missing_a.fasta"),
                    "--group-b-fasta", str(root / "missing_b.fasta"),
                    "--skip-plots",
                ]
            )
            self.assertEqual(summary["status"], "skipped")
            tr_bed = pd.read_csv(root / "out" / "tr_changes.bed.tsv", sep="\t")
            self.assertIn("tr_pass_for_harmonized", tr_bed.columns)
            self.assertEqual(len(tr_bed), 0)

    def test_tr_post_processing_main_accepts_trgt_spanning_bam(self):
        from array import array

        import pysam

        from sniffcell.discover.tr_post_processing import tr_post_processing_main

        def write_bam(path: Path, lengths: list[int]) -> None:
            path.parent.mkdir(parents=True, exist_ok=True)
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
            with pysam.AlignmentFile(str(path), "wb", header=header) as bam_out:
                for idx, length in enumerate(lengths):
                    read = pysam.AlignedSegment()
                    read.query_name = f"read{idx}"
                    read.query_sequence = "A" * length
                    read.flag = 0
                    read.reference_id = 0
                    read.reference_start = 50
                    read.mapping_quality = 60
                    read.cigar = ((0, length),)
                    read.query_qualities = pysam.qualitystring_to_array("I" * length)
                    read.set_tag("TR", "tr1")
                    read.set_tag("FL", array("I", [50, 50]))
                    bam_out.write(read)

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            split_dir = root / "sample1" / "deconv" / "deconv_requested_group_splits"
            split_dir.mkdir(parents=True)
            tr_bed = root / "repeats.bed"
            tr_bed.write_text("chr1\t100\t200\tID=tr1;MOTIFS=A;STRUC=(A)n\n", encoding="utf-8")
            group_a_bam = root / "a.spanning.bam"
            group_b_bam = root / "b.spanning.bam"
            write_bam(group_a_bam, [330, 325, 320, 210, 205])
            write_bam(group_b_bam, [200, 198, 196, 194, 192])

            summary = tr_post_processing_main(
                [
                    "--split-dir", str(split_dir),
                    "--groups", "Neuron,Oligodendrocyte",
                    "--output-dir", str(root / "out"),
                    "--sample-id", "sample1",
                    "--sample-a-label", "sample1.Neuron",
                    "--sample-b-label", "sample1.Oligodendrocyte",
                    "--group-a-spanning-bam", str(group_a_bam),
                    "--group-b-spanning-bam", str(group_b_bam),
                    "--tr-bed", str(tr_bed),
                    "--skip-plots",
                ]
            )

            self.assertEqual(summary["status"], "completed")
            self.assertEqual(summary["group_a_source"], "trgt_spanning_bam")
            self.assertEqual(summary["n_targets"], 1)
            tr_rows = pd.read_csv(root / "out" / "tr_changes.bed.tsv", sep="\t")
            self.assertEqual(tr_rows["change_group"].tolist(), ["sample1.Neuron"])
            self.assertEqual(int(tr_rows.iloc[0]["change_length_bp"]), 130)
            read_lengths = sorted(pd.read_csv(root / "out" / "read_lengths.tsv", sep="\t")["read_length"].tolist())
            self.assertEqual(read_lengths, [92, 94, 96, 98, 100, 105, 110, 220, 225, 230])


class TestPostprocessLocalIntegration(unittest.TestCase):

    def test_discover_main_local_second_run_skips_completed_stage(self):
        from sniffcell.discover.discover import discover_main

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_fake_runtime_env(root)
            argv = _base_local_argv(deconv_dir, ref, tr_bed, tool_paths) + ["--stages", "sniffles"]

            discover_main(parse_args(argv))
            discover_main(parse_args(argv))

            run_root = deconv_dir / "deconv_requested_group_splits" / "discover" / "testrun"
            status_text = (run_root / "status" / "discover_status.json").read_text(encoding="utf-8")
            self.assertIn('"state": "skipped"', status_text)
            self.assertTrue((run_root / "status" / "sniffles.Neuron.done.json").exists())
            self.assertTrue((run_root / "status" / "sniffles.Oligodendrocyte.done.json").exists())


if __name__ == "__main__":
    unittest.main()
