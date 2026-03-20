import dataclasses
import os
import shlex
import stat
import sys
import tempfile
import unittest
from pathlib import Path

from sniffcell.parse_args import parse_args
from sniffcell.discover.discover import (
    DEFAULT_STAGE_ORDER,
    GROUP_SCOPED_STAGES,
    RunContext,
    _build_context,
    _build_recursive_cli,
    _parse_stages,
    _render_slurm,
    _render_submit_script,
    _sanitize_token,
    _select_groups,
    _discover_groups,
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
from pathlib import Path
args = sys.argv[1:]
output_dir = None
for arg in args:
    if arg.startswith("--output="):
        output_dir = Path(arg.split("=", 1)[1])
        break
if output_dir is None:
    raise SystemExit("missing --output")
output = output_dir / "merge_output.vcf.gz"
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text("##fileformat=VCFv4.2\\n", encoding="utf-8")
""",
    )

    _make_python_exec(
        tool_dir / "run_clairs",
        """
import sys
from pathlib import Path
args = sys.argv[1:]
output_dir = Path(args[args.index("-o") + 1])
output = output_dir / "output.vcf.gz"
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text("##fileformat=VCFv4.2\\n", encoding="utf-8")
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
        "tdb", "modkit", "tabix", "run_clair3.sh", "run_clairs",
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
        "--clair3-bin", tool_paths["run_clair3.sh"],
        "--clair3-model-path", "/tmp/clair3_model",
        "--clairs-bin", tool_paths["run_clairs"],
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
        """--run-id, --stages, --groups, --clairs-tumor-group are suppressed but functional."""
        args = parse_args([
            "discover",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "male",
            "--run-id", "myrun",
            "--stages", "sv",
            "--groups", "A,B",
            "--clairs-tumor-group", "A",
        ])
        self.assertEqual(args.run_id, "myrun")
        self.assertEqual(args.stages, "sv")
        self.assertEqual(args.groups, "A,B")
        self.assertEqual(args.clairs_tumor_group, "A")

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
            "tdb_merge_threads", "modkit_threads", "clair3_threads", "clairs_threads",
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

    def test_clairs_platform_default(self):
        args = parse_args([
            "discover",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "male",
        ])
        self.assertEqual(args.clairs_platform, "ont_r10_dorado_sup_5khz")

    def test_clair3_platform_default(self):
        args = parse_args([
            "discover",
            "--deconv-dir", "/tmp/d",
            "--reference", "/tmp/r",
            "--tr-bed", "/tmp/t",
            "--sex", "male",
        ])
        self.assertEqual(args.clair3_platform, "ont")


# ---------------------------------------------------------------------------
# _parse_stages tests
# ---------------------------------------------------------------------------

class TestParseStages(unittest.TestCase):

    def test_none_returns_default_order(self):
        self.assertEqual(_parse_stages(None), DEFAULT_STAGE_ORDER)

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
        self.assertEqual(stages, ("clair3", "clairs"))

    def test_mods_alias(self):
        stages = _parse_stages("mods")
        self.assertEqual(stages, ("modkit",))

    def test_all_alias_returns_full_order(self):
        stages = _parse_stages("all")
        self.assertEqual(stages, DEFAULT_STAGE_ORDER)

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
        self.assertEqual(_parse_stages("clairs"), ("clairs",))

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

    def test_build_context_clairs_tumor_group_in_params(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            argv = _base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths) + [
                "--clairs-tumor-group", "Neuron",
                "--groups", "Neuron,Oligodendrocyte",
                "--scheduler", "local",
            ]
            args = parse_args(argv)
            ctx = _build_context(args)
            self.assertEqual(ctx.params["clairs_tumor_group"], "Neuron")

    def test_build_context_mods_mode_propagated(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            argv = _base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths) + ["--mods-mode", "combined"]
            args = parse_args(argv)
            ctx = _build_context(args)
            self.assertEqual(ctx.params["mods_mode"], "combined")


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

    def test_clairs_tumor_group_included_when_specified(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td))
            cli = _build_recursive_cli(ctx, "clairs", clairs_tumor_group="Neuron")
            self.assertIn("--clairs-tumor-group", cli)
            idx = cli.index("--clairs-tumor-group")
            self.assertEqual(cli[idx + 1], "Neuron")

    def test_clairs_tumor_group_absent_when_not_specified(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td))
            cli = _build_recursive_cli(ctx, "clairs")
            self.assertNotIn("--clairs-tumor-group", cli)

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

    def test_clair3_model_path_forwarded_when_set(self):
        with tempfile.TemporaryDirectory() as td:
            ctx = self._get_ctx(Path(td))
            cli = _build_recursive_cli(ctx, "clair3")
            self.assertIn("--clair3-model-path", cli)
            idx = cli.index("--clair3-model-path")
            self.assertEqual(cli[idx + 1], "/tmp/clair3_model")


# ---------------------------------------------------------------------------
# SLURM script generation tests
# ---------------------------------------------------------------------------

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
            for stage in GROUP_SCOPED_STAGES:
                script = ctx.slurm_dir / f"{stage}.array.sbatch.sh"
                self.assertTrue(script.exists(), f"Missing SLURM script for {stage}")

    def test_sample_scoped_scripts_are_generated(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            _render_slurm(ctx)
            for script_name in ("collapse.sbatch.sh", "tdb_merge.sbatch.sh"):
                self.assertTrue(
                    (ctx.slurm_dir / script_name).exists(),
                    f"Missing SLURM script: {script_name}",
                )

    def test_clairs_both_direction_scripts_generated(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            _render_slurm(ctx)
            # Both directions must exist
            fwd = ctx.slurm_dir / "clairs_Neuron_vs_Oligodendrocyte.sbatch.sh"
            rev = ctx.slurm_dir / "clairs_Oligodendrocyte_vs_Neuron.sbatch.sh"
            self.assertTrue(fwd.exists(), "Missing forward ClairS SLURM script")
            self.assertTrue(rev.exists(), "Missing reverse ClairS SLURM script")

    def test_clairs_scripts_contain_clairs_tumor_group(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            _render_slurm(ctx)
            fwd = (ctx.slurm_dir / "clairs_Neuron_vs_Oligodendrocyte.sbatch.sh").read_text()
            rev = (ctx.slurm_dir / "clairs_Oligodendrocyte_vs_Neuron.sbatch.sh").read_text()
            self.assertIn("--clairs-tumor-group Neuron", fwd)
            self.assertIn("--clairs-tumor-group Oligodendrocyte", rev)

    def test_clairs_scripts_have_no_partition_or_account(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            _render_slurm(ctx)
            for script in (
                "clairs_Neuron_vs_Oligodendrocyte.sbatch.sh",
                "clairs_Oligodendrocyte_vs_Neuron.sbatch.sh",
            ):
                text = (ctx.slurm_dir / script).read_text()
                self.assertNotIn("#SBATCH --partition=", text)
                self.assertNotIn("#SBATCH --account=", text)

    def test_submit_script_has_both_clairs_jids(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_minimal_env(root)
            args = parse_args(_base_slurm_argv(deconv_dir, ref, tr_bed, tool_paths))
            ctx = _build_context(args)
            _render_slurm(ctx)
            submit_text = (ctx.slurm_dir / "submit_pipeline.sh").read_text()
            self.assertIn("CLAIRS_NEURON_VS_OLIGODENDROCYTE_JID", submit_text)
            self.assertIn("CLAIRS_OLIGODENDROCYTE_VS_NEURON_JID", submit_text)

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
            self.assertIn("CLAIR3_JID", submit_text)
            self.assertIn("MEDAKA_JID", submit_text)
            self.assertIn("MODKIT_JID", submit_text)

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
            # tdb_create depends on medaka
            self.assertIn("afterok:${MEDAKA_JID}", submit_text)
            # tdb_merge depends on tdb_create
            self.assertIn("afterok:${TDB_CREATE_JID}", submit_text)

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

    def test_sniffles_command_contains_sample_name(self):
        """_run_sniffles must pass --sample-name <sampleid>_<group> to sniffles."""
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
                "--clairs-bin", tool_paths["run_clairs"],
            ]
            args = parse_args(argv)
            ctx = _build_context(args)
            _run_sniffles(ctx, "Neuron")
            cmd_file = ctx.commands_dir / "sniffles.Neuron.command.txt"
            cmd_text = cmd_file.read_text()
            self.assertIn("--sample-name", cmd_text)
            self.assertIn("sample1_Neuron", cmd_text)

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


class TestPostprocessLocalIntegration(unittest.TestCase):

    def test_discover_main_local_full_pipeline_simulated(self):
        from sniffcell.discover.discover import discover_main

        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            deconv_dir, ref, tr_bed, tool_paths = _build_fake_runtime_env(root)
            argv = _base_local_argv(deconv_dir, ref, tr_bed, tool_paths) + [
                "--stages",
                "sniffles,sniffles_filter,kanpig,collapse,medaka,tdb_create,tdb_merge,clair3,clairs,modkit",
            ]
            args = parse_args(argv)
            discover_main(args)

            run_root = deconv_dir / "deconv_requested_group_splits" / "discover" / "testrun"
            self.assertTrue((run_root / "run_summary.json").exists())
            self.assertTrue((run_root / "manifest" / "discover_run_manifest.json").exists())
            self.assertTrue((run_root / "manifest" / "discover_task_manifest.tsv").exists())

            for group in ("Neuron", "Oligodendrocyte"):
                self.assertTrue((run_root / "sv" / "sniffles" / group / "sniffles.raw.vcf.gz").exists())
                self.assertTrue((run_root / "sv" / "sniffles" / group / "sniffles.mosaic_only.vcf.gz").exists())
                self.assertTrue((run_root / "sv" / "kanpig" / group / "kanpig.mosaic.vcf.gz").exists())
                self.assertTrue((run_root / "medaka_tandem" / f"{group}.medaka" / "medaka_to_ref.TR.vcf").exists())
                self.assertTrue((run_root / "medaka_tandem" / "tdb" / f"{group}.tdb").exists())
                self.assertTrue((run_root / "modkit" / group / f"{group}.cpg.bedmethyl.gz").exists())
                self.assertTrue((run_root / "snv" / "clair3" / group / "merge_output.vcf.gz").exists())

            self.assertTrue(
                (run_root / "sv" / "truvari_collapse" / "Neuron_vs_Oligodendrocyte" / "collapsed.sorted.vcf.gz").exists()
            )
            self.assertTrue((run_root / "medaka_tandem" / "tdb" / "sample1.merged.tdb").exists())
            self.assertTrue(
                (run_root / "medaka_tandem" / "tr_post_processing" / "Neuron_vs_Oligodendrocyte" / "summary.json").exists()
            )
            self.assertTrue((run_root / "snv" / "clairs" / "Neuron_vs_Oligodendrocyte" / "output.vcf.gz").exists())
            self.assertTrue((run_root / "snv" / "clairs" / "Oligodendrocyte_vs_Neuron" / "output.vcf.gz").exists())
            self.assertTrue(
                (run_root / "sv" / "sv_post_processing" / "Neuron_vs_Oligodendrocyte" / "summary.json").exists()
            )

            status_text = (run_root / "status" / "discover_status.json").read_text(encoding="utf-8")
            self.assertIn('"state": "completed"', status_text)

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
