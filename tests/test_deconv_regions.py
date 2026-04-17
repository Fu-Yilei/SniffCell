import tempfile
import unittest
from pathlib import Path

import pandas as pd

from sniffcell.deconv.regions import (
    TargetRegion,
    load_target_regions,
    parse_region_spec,
    select_ctdmrs_for_target,
)
from sniffcell.parse_args import parse_args


class TestDeconvRegionParseArgs(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
