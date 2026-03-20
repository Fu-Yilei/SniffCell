import os
import stat
import tempfile
import unittest
from pathlib import Path

from sniffcell.parse_args import parse_args
from sniffcell.postprocess.postprocess import _build_context, _parse_stages, _render_slurm


class TestPostprocessParseArgs(unittest.TestCase):
    def test_postprocess_accepts_core_arguments(self):
        args = parse_args(
            [
                "postprocess",
                "--deconv-dir",
                "/tmp/sample/deconv",
                "--reference",
                "/tmp/ref.fa",
                "--tr-bed",
                "/tmp/tr.bed",
                "--sex",
                "male",
            ]
        )

        self.assertEqual(args.command, "postprocess")
        self.assertEqual(args.scheduler, "local")
        self.assertEqual(args.slurm_partition, "medium")
        self.assertEqual(args.slurm_account, "proj-fs0006")
        self.assertEqual(args.sniffles_threads, 24)
        self.assertEqual(args.medaka_workers, 8)
        self.assertEqual(args.mods_mode, "separate")

    def test_parse_stages_expands_aliases(self):
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


class TestPostprocessContextAndSlurm(unittest.TestCase):
    def _make_exec(self, path: Path) -> None:
        path.write_text("#!/bin/bash\nexit 0\n")
        path.chmod(path.stat().st_mode | stat.S_IXUSR)

    def test_build_context_discovers_two_groups_and_renders_slurm(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            sample_dir = root / "sample1"
            deconv_dir = sample_dir / "deconv"
            split_dir = deconv_dir / "deconv_requested_group_splits"
            split_dir.mkdir(parents=True)
            for group_name in ("Neuron", "Oligodendrocyte"):
                bam = split_dir / f"{group_name}.bam"
                bai = split_dir / f"{group_name}.bam.bai"
                bam.write_text("")
                bai.write_text("")
            (split_dir / "requested_group_splits.tsv").write_text(
                "requested_group\tbam_path\tread_summary_path\n"
                f"Neuron\t{split_dir / 'Neuron.bam'}\t{split_dir / 'Neuron.read_summary.tsv'}\n"
                f"Oligodendrocyte\t{split_dir / 'Oligodendrocyte.bam'}\t{split_dir / 'Oligodendrocyte.read_summary.tsv'}\n"
            )
            ref = root / "ref.fa"
            tr_bed = root / "tr.bed"
            ref.write_text(">chr1\nA\n")
            tr_bed.write_text("chr1\t0\t10\n")
            tool_dir = root / "bin"
            tool_dir.mkdir()
            tool_paths = {}
            for tool_name in ("sniffles", "bcftools", "kanpig", "truvari", "medaka", "tdb", "modkit", "tabix"):
                tool_path = tool_dir / tool_name
                self._make_exec(tool_path)
                tool_paths[tool_name] = str(tool_path)

            args = parse_args(
                [
                    "postprocess",
                    "--deconv-dir",
                    str(deconv_dir),
                    "--reference",
                    str(ref),
                    "--tr-bed",
                    str(tr_bed),
                    "--sex",
                    "male",
                    "--scheduler",
                    "slurm",
                    "--run-id",
                    "testrun",
                    "--sniffles-bin",
                    tool_paths["sniffles"],
                    "--bcftools-bin",
                    tool_paths["bcftools"],
                    "--kanpig-bin",
                    tool_paths["kanpig"],
                    "--truvari-bin",
                    tool_paths["truvari"],
                    "--medaka-bin",
                    tool_paths["medaka"],
                    "--tdb-bin",
                    tool_paths["tdb"],
                    "--modkit-bin",
                    tool_paths["modkit"],
                    "--tabix-bin",
                    tool_paths["tabix"],
                ]
            )

            ctx = _build_context(args)
            self.assertEqual(ctx.sample_id, "sample1")
            self.assertEqual(ctx.selected_groups, ["Neuron", "Oligodendrocyte"])

            _render_slurm(ctx, submit=False)
            sniffles_script = ctx.slurm_dir / "sniffles.array.sbatch.sh"
            self.assertTrue(sniffles_script.exists())
            script_text = sniffles_script.read_text()
            self.assertIn("#SBATCH --partition=medium", script_text)
            self.assertIn("#SBATCH --account=proj-fs0006", script_text)
            self.assertIn('${GROUP_NAME}', script_text)
            self.assertIn("--scheduler local", script_text)


if __name__ == "__main__":
    unittest.main()
