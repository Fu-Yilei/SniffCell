import unittest
import tempfile
from pathlib import Path

import pandas as pd
import pysam

from sniffcell.deconv.deconv import (
    _build_deconv_summary,
    _build_read_summary,
    _parse_requested_split_groups,
    _resolve_output_paths,
    _write_requested_split_group_outputs,
    _write_group_split_reads,
)
from sniffcell.parse_args import parse_args


class TestDeconvParseArgs(unittest.TestCase):
    def test_deconv_accepts_threads_and_assignment_mode(self):
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
                "-t",
                "4",
                "--read_assignment_mode",
                "kmeans",
            ]
        )

        self.assertEqual(args.command, "deconv")
        self.assertEqual(args.threads, 4)
        self.assertEqual(args.read_assignment_mode, "kmeans")

    def test_deconv_accepts_requested_bam_split_groups(self):
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
                "--split_bam_groups",
                "t_cell,b_cell,nk_cell;monocyte",
            ]
        )

        self.assertEqual(args.split_bam_groups, "t_cell,b_cell,nk_cell;monocyte")


class TestDeconvSummaries(unittest.TestCase):
    def setUp(self):
        self.read_assignment_df = pd.DataFrame(
            {
                "code_order": ["A|B", "A|B", "A|B", "A|B", "A|B"],
                "code": ["10", "10", "01", "01", "01"],
                "best_group": ["A", "A", "B", "B", "B"],
                "best_group_leaves": ["A", "A", "B", "B", "B"],
                "is_best_group": [True, True, True, True, True],
            },
            index=pd.Index(["r1", "r1", "r2", "r2", "r3"], name="readname"),
        )

    def test_build_read_summary_rolls_up_majority_vote_per_read(self):
        out = _build_read_summary(self.read_assignment_df)

        self.assertEqual(out["readname"].tolist(), ["r1", "r2", "r3"])
        self.assertEqual(out["majority_code"].tolist(), ["10", "01", "01"])
        self.assertEqual(out["primary_celltype"].tolist(), ["A", "B", "B"])
        self.assertEqual(out["linked_celltypes"].tolist(), ["A", "B", "B"])
        self.assertEqual(out["linked_leaf_celltypes"].tolist(), ["A", "B", "B"])
        self.assertEqual(out["n_ctdmrs"].tolist(), [2, 2, 1])
        self.assertAlmostEqual(float(out.loc[out["readname"] == "r1", "majority_pct"].iloc[0]), 1.0)

    def test_build_deconv_summary_reports_all_rows_and_per_read_views(self):
        out = _build_deconv_summary(self.read_assignment_df)

        self.assertEqual(out["summary_mode"].tolist(), ["all_rows", "per_read"])

        all_rows = out.loc[out["summary_mode"] == "all_rows"].iloc[0]
        self.assertEqual(all_rows["primary_celltype"], "B")
        self.assertEqual(all_rows["linked_celltype_counts"], "B:3;A:2")
        self.assertEqual(all_rows["linked_celltype_fractions"], "B:0.600;A:0.400")
        self.assertEqual(int(all_rows["n_evidence_units"]), 5)

        per_read = out.loc[out["summary_mode"] == "per_read"].iloc[0]
        self.assertEqual(per_read["primary_celltype"], "B")
        self.assertEqual(per_read["linked_celltype_counts"], "B:2;A:1")
        self.assertEqual(per_read["linked_celltype_fractions"], "B:0.667;A:0.333")
        self.assertEqual(int(per_read["n_evidence_units"]), 3)


class TestDeconvOutputPaths(unittest.TestCase):
    def test_resolve_output_paths_supports_directory_or_summary_file(self):
        dir_paths = _resolve_output_paths("out_dir")
        self.assertTrue(dir_paths["summary"].endswith("out_dir/deconv_summary.tsv"))
        self.assertTrue(dir_paths["reads"].endswith("out_dir/deconv_reads_classification.tsv"))
        self.assertTrue(dir_paths["group_dir"].endswith("out_dir/deconv_reads_by_group"))

        file_paths = _resolve_output_paths("out_dir/custom_summary.tsv")
        self.assertTrue(file_paths["summary"].endswith("out_dir/custom_summary.tsv"))
        self.assertTrue(file_paths["read_summary"].endswith("out_dir/deconv_read_summary.tsv"))

    def test_write_group_split_reads_creates_one_file_per_best_group(self):
        read_assignment_df = pd.DataFrame(
            {
                "best_group": ["A", "A", "B"],
                "best_group_leaves": ["A", "A", "B"],
                "code_order": ["A|B", "A|B", "A|B"],
                "code": ["10", "10", "01"],
                "is_best_group": [True, True, True],
            },
            index=pd.Index(["r1", "r2", "r3"], name="readname"),
        )

        with tempfile.TemporaryDirectory() as td:
            written = _write_group_split_reads(read_assignment_df, td)

            self.assertEqual(len(written), 2)
            self.assertTrue((Path(td) / "A.tsv").exists())
            self.assertTrue((Path(td) / "B.tsv").exists())

            a_df = pd.read_csv(Path(td) / "A.tsv", sep="\t", index_col=0)
            b_df = pd.read_csv(Path(td) / "B.tsv", sep="\t", index_col=0)
            self.assertEqual(a_df.index.tolist(), ["r1", "r2"])
            self.assertEqual(b_df.index.tolist(), ["r3"])


class TestRequestedBamSplits(unittest.TestCase):
    def _write_test_bam(self, bam_path: Path, read_names: list[str]) -> None:
        header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 100000}]}
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam_out:
            for idx, read_name in enumerate(read_names):
                seg = pysam.AlignedSegment()
                seg.query_name = read_name
                seg.query_sequence = "A" * 100
                seg.flag = 0
                seg.reference_id = 0
                seg.reference_start = 1000 + idx * 10
                seg.mapping_quality = 60
                seg.cigar = ((0, 100),)
                seg.query_qualities = pysam.qualitystring_to_array("I" * 100)
                bam_out.write(seg)
        pysam.index(str(bam_path))

    def _read_bam_names(self, bam_path: Path) -> list[str]:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam_in:
            return [record.query_name for record in bam_in.fetch(until_eof=True)]

    def test_requested_group_splits_expand_parent_labels_to_leaf_matches(self):
        split_specs = _parse_requested_split_groups("t_cell,b_cell,nk_cell;monocyte")
        self.assertEqual([spec["file_stub"] for spec in split_specs], ["t_cell_b_cell_nk_cell", "monocyte"])

        read_assignment_df = pd.DataFrame(
            {
                "best_group": ["lymphocytes", "lymphocytes", "Monocyte"],
                "best_group_leaves": [
                    "T-cell|NK-cell|B-cell",
                    "T-cell|NK-cell|B-cell",
                    "Classical-Monocyte|Nonclassical-Monocyte",
                ],
                "other_group_leaves": [
                    "Classical-Monocyte|Nonclassical-Monocyte",
                    "Classical-Monocyte|Nonclassical-Monocyte",
                    "T-cell|NK-cell|B-cell",
                ],
                "code_order": [
                    "T-cell|NK-cell|B-cell|Classical-Monocyte|Nonclassical-Monocyte",
                    "T-cell|NK-cell|B-cell|Classical-Monocyte|Nonclassical-Monocyte",
                    "T-cell|NK-cell|B-cell|Classical-Monocyte|Nonclassical-Monocyte",
                ],
                "code": ["10100", "00100", "00011"],
                "is_best_group": [True, True, True],
            },
            index=pd.Index(["r1", "r2", "r3"], name="readname"),
        )
        read_summary_df = _build_read_summary(read_assignment_df)

        with tempfile.TemporaryDirectory() as td:
            bam_path = Path(td) / "input.bam"
            self._write_test_bam(bam_path, ["r1", "r2", "r3", "rx"])

            manifest = _write_requested_split_group_outputs(
                bam_path=str(bam_path),
                read_summary_df=read_summary_df,
                read_assignment_df=read_assignment_df,
                split_group_spec="t_cell,b_cell,nk_cell;monocyte",
                output_dir=str(Path(td) / "splits"),
            )

            self.assertEqual(manifest["requested_group"].tolist(), ["t_cell,b_cell,nk_cell", "monocyte"])
            self.assertEqual(manifest["n_reads"].tolist(), [2, 1])

            first_bam = Path(td) / "splits" / "t_cell_b_cell_nk_cell.bam"
            second_bam = Path(td) / "splits" / "monocyte.bam"
            self.assertEqual(self._read_bam_names(first_bam), ["r1", "r2"])
            self.assertEqual(self._read_bam_names(second_bam), ["r3"])

            first_tsv = pd.read_csv(Path(td) / "splits" / "t_cell_b_cell_nk_cell.read_summary.tsv", sep="\t")
            second_tsv = pd.read_csv(Path(td) / "splits" / "monocyte.read_summary.tsv", sep="\t")
            self.assertEqual(first_tsv["readname"].tolist(), ["r1", "r2"])
            self.assertEqual(second_tsv["readname"].tolist(), ["r3"])


if __name__ == "__main__":
    unittest.main()
