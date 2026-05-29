import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from sniffcell.deconv.regions import (
    TargetRegion,
    load_target_regions,
    parse_region_spec,
    regions_main,
    resolve_region_plan,
    select_ctdmrs_for_target,
)
from sniffcell.parse_args import parse_args


class TestDeconvRegionParseArgs(unittest.TestCase):
    def test_regions_command_accepts_bed_target(self):
        args = parse_args(
            [
                "regions",
                "-b",
                "ctdmr.tsv",
                "--regions",
                "targets.bed",
                "-o",
                "out_dir",
                "--regions-ctdmrs",
                "3",
                "--regions-left-ctdmrs",
                "1",
            ]
        )

        self.assertEqual(args.command, "regions")
        self.assertEqual(args.bed, "ctdmr.tsv")
        self.assertEqual(args.regions, "targets.bed")
        self.assertEqual(args.output, "out_dir")
        self.assertEqual(args.regions_ctdmrs, 3)
        self.assertEqual(args.regions_left_ctdmrs, 1)
        self.assertIsNone(args.regions_right_ctdmrs)

    def test_deconv_accepts_regions_string(self):
        args = parse_args(
            [
                "deconv",
                "-i",
                "input.bam",
                "-r",
                "ref.fa",
                "-b",
                "ctdmr.tsv",
                "-o",
                "out_dir",
                "--regions",
                "chr13:100-200",
                "--regions-left-ctdmrs",
                "4",
                "--regions-right-ctdmrs",
                "5",
            ]
        )

        self.assertEqual(args.command, "deconv")
        self.assertEqual(args.regions, "chr13:100-200")
        self.assertEqual(args.regions_left_ctdmrs, 4)
        self.assertEqual(args.regions_right_ctdmrs, 5)


class TestDeconvRegionParsing(unittest.TestCase):
    def test_parse_region_spec(self):
        target = parse_region_spec("chr13:102161550-102161589")

        self.assertEqual(
            target,
            TargetRegion(chrom="chr13", start=102161550, end=102161589, name=None),
        )

    def test_load_target_regions_from_bed(self):
        with tempfile.TemporaryDirectory() as td:
            bed_path = Path(td) / "targets.bed"
            bed_path.write_text("chr1\t10\t20\talpha\nchr2\t30\t40\n", encoding="utf-8")

            targets = load_target_regions(str(bed_path))

        self.assertEqual(
            targets,
            [
                TargetRegion(chrom="chr1", start=10, end=20, name="alpha"),
                TargetRegion(chrom="chr2", start=30, end=40, name=None),
            ],
        )


class TestDeconvRegionSelection(unittest.TestCase):
    def setUp(self):
        self.ctdmr_df = pd.DataFrame(
            [
                {"chr": "chr1", "start": 10, "end": 20, "best_group": "A", "best_group_leaves": "A", "other_group_leaves": "B", "code_order": "A|B", "best_dir": "hyper"},
                {"chr": "chr1", "start": 30, "end": 40, "best_group": "A", "best_group_leaves": "A", "other_group_leaves": "B", "code_order": "A|B", "best_dir": "hyper"},
                {"chr": "chr1", "start": 48, "end": 62, "best_group": "A", "best_group_leaves": "A", "other_group_leaves": "B", "code_order": "A|B", "best_dir": "hyper"},
                {"chr": "chr1", "start": 70, "end": 75, "best_group": "B", "best_group_leaves": "B", "other_group_leaves": "A", "code_order": "A|B", "best_dir": "hyper"},
                {"chr": "chr1", "start": 90, "end": 100, "best_group": "B", "best_group_leaves": "B", "other_group_leaves": "A", "code_order": "A|B", "best_dir": "hyper"},
                {"chr": "chr1", "start": 120, "end": 130, "best_group": "B", "best_group_leaves": "B", "other_group_leaves": "A", "code_order": "A|B", "best_dir": "hyper"},
            ]
        )

    def test_select_nearest_ctdmrs_includes_overlap_and_flanks(self):
        selected = select_ctdmrs_for_target(
            self.ctdmr_df,
            target=TargetRegion(chrom="chr1", start=50, end=60, name=None),
            left_ctdmrs=2,
            right_ctdmrs=2,
        )

        self.assertEqual(selected["start"].tolist(), [10, 30, 48, 70, 90])

    def test_select_nearest_ctdmrs_works_without_overlap(self):
        selected = select_ctdmrs_for_target(
            self.ctdmr_df,
            target=TargetRegion(chrom="chr1", start=80, end=85, name=None),
            left_ctdmrs=2,
            right_ctdmrs=1,
        )

        self.assertEqual(selected["start"].tolist(), [48, 70, 90])


class TestRegionPlanOutputs(unittest.TestCase):
    def setUp(self):
        self.ctdmr_df = pd.DataFrame(
            [
                {"chr": "chr1", "start": 10, "end": 20, "name": "dmr1", "best_group": "A", "other_group": "B", "best_group_leaves": "A", "other_group_leaves": "B", "code_order": "A|B", "best_dir": "hyper", "mean_A": 0.9, "mean_B": 0.1},
                {"chr": "chr1", "start": 30, "end": 40, "name": "dmr2", "best_group": "A", "other_group": "B", "best_group_leaves": "A", "other_group_leaves": "B", "code_order": "A|B", "best_dir": "hyper", "mean_A": 0.8, "mean_B": 0.2},
                {"chr": "chr1", "start": 48, "end": 62, "name": "dmr3", "best_group": "A", "other_group": "B", "best_group_leaves": "A", "other_group_leaves": "B", "code_order": "A|B", "best_dir": "hypo", "mean_A": 0.2, "mean_B": 0.8},
                {"chr": "chr1", "start": 70, "end": 75, "name": "dmr4", "best_group": "B", "other_group": "A", "best_group_leaves": "B", "other_group_leaves": "A", "code_order": "A|B", "best_dir": "hyper", "mean_A": 0.2, "mean_B": 0.8},
                {"chr": "chr1", "start": 90, "end": 100, "name": "dmr5", "best_group": "B", "other_group": "A", "best_group_leaves": "B", "other_group_leaves": "A", "code_order": "A|B", "best_dir": "hyper", "mean_A": 0.1, "mean_B": 0.9},
            ]
        )

    def test_resolve_region_plan_writes_subset_bed_and_summaries(self):
        with tempfile.TemporaryDirectory() as td:
            target_bed = Path(td) / "targets.bed"
            output_dir = Path(td) / "plan"
            target_bed.write_text("chr1\t50\t60\tlocusA\n", encoding="utf-8")

            plan = resolve_region_plan(
                ctdmr_df=self.ctdmr_df,
                output_dir=str(output_dir),
                regions_arg=str(target_bed),
                left_ctdmrs=2,
                right_ctdmrs=1,
            )

            self.assertEqual(plan.selected_ctdmr_count, 4)
            self.assertTrue((output_dir / "targets.bed").exists())
            self.assertTrue((output_dir / "subset_regions.bed").exists())
            self.assertTrue((output_dir / "ctdmr_subset.tsv").exists())
            self.assertTrue((output_dir / "ctdmr_region_summary.tsv").exists())
            self.assertTrue((output_dir / "ctdmr_selected_summary.tsv").exists())

            self.assertEqual(
                (output_dir / "subset_regions.bed").read_text(encoding="utf-8"),
                "chr1\t10\t75\tregion_1\n",
            )
            subset = pd.read_csv(output_dir / "ctdmr_subset.tsv", sep="\t")
            self.assertEqual(subset["start"].tolist(), [10, 30, 48, 70])

            region_summary = pd.read_csv(output_dir / "ctdmr_region_summary.tsv", sep="\t")
            self.assertEqual(region_summary.loc[0, "target_name"], "locusA")
            self.assertEqual(int(region_summary.loc[0, "selected_ctdmr_count"]), 4)
            self.assertEqual(int(region_summary.loc[0, "left_flank_ctdmr_count"]), 2)
            self.assertEqual(int(region_summary.loc[0, "overlap_ctdmr_count"]), 1)
            self.assertEqual(int(region_summary.loc[0, "right_flank_ctdmr_count"]), 1)
            self.assertEqual(region_summary.loc[0, "best_groups"], "A|B")

            selected_summary = pd.read_csv(output_dir / "ctdmr_selected_summary.tsv", sep="\t")
            self.assertEqual(selected_summary["relation"].tolist(), ["left_flank", "left_flank", "overlap", "right_flank"])
            self.assertEqual(selected_summary["distance_bp"].tolist(), [30, 10, 0, 10])

            manifest = json.loads((output_dir / "region_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["selected_ctdmr_count"], 4)
            self.assertTrue(manifest["outputs"]["subset_regions_bed"].endswith("subset_regions.bed"))

    def test_regions_main_loads_ctdmr_tsv(self):
        with tempfile.TemporaryDirectory() as td:
            ctdmr_path = Path(td) / "ctdmr.tsv"
            target_bed = Path(td) / "targets.bed"
            output_dir = Path(td) / "plan"
            self.ctdmr_df.to_csv(ctdmr_path, sep="\t", index=False)
            target_bed.write_text("chr1\t50\t60\tlocusA\n", encoding="utf-8")
            args = parse_args(
                [
                    "regions",
                    "-b",
                    str(ctdmr_path),
                    "--regions",
                    str(target_bed),
                    "-o",
                    str(output_dir),
                    "--regions-left-ctdmrs",
                    "1",
                    "--regions-right-ctdmrs",
                    "0",
                ]
            )

            plan = regions_main(args)

            subset = pd.read_csv(plan.subset_bed_path, sep="\t")
            self.assertEqual(subset["start"].tolist(), [30, 48])
            self.assertEqual((output_dir / "subset_regions.bed").read_text(encoding="utf-8"), "chr1\t30\t62\tregion_1\n")


if __name__ == "__main__":
    unittest.main()
