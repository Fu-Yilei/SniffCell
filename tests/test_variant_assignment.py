import unittest

import pandas as pd

from sniffcell.anno.anno import _build_sv_readable_reports
from sniffcell.anno.variant_assignment import assign_sv_celltypes
from sniffcell.parse_args import parse_args


class TestVariantAssignment(unittest.TestCase):
    def test_multi_region_votes_per_read_are_aggregated(self):
        sv_df = pd.DataFrame(
            {
                "chr": ["1"],
                "location": [1000],
                "id": ["sv1"],
                "supporting_reads": [["r1", "r2"]],
            }
        )

        read_assignment_df = pd.DataFrame(
            {
                "code_order": ["A|B", "A|B", "A|B", "A|B"],
                "code": ["10", "01", "01", "01"],
                "best_group": ["A", "A", "B", "B"],
                "is_best_group": [True, False, True, True],
                "chr": ["1", "1", "1", "1"],
                "start": [100, 200, 300, 120],
                "end": [150, 250, 350, 170],
            },
            index=pd.Index(["r1", "r1", "r1", "r2"], name="readname"),
        )

        out = assign_sv_celltypes(sv_df, read_assignment_df)
        row = out.iloc[0]

        self.assertEqual(row["n_supporting"], 2)
        self.assertEqual(row["n_overlapped"], 2)
        self.assertEqual(row["majority_code"], "01")
        self.assertEqual(row["assigned_code"], "01")
        self.assertEqual(row["code_counts"], "01:2")
        self.assertEqual(row["majority_schema"], "")
        # With code-based decoding, read-level winning code "01" maps to B.
        self.assertEqual(row["primary_celltype"], "B")
        self.assertEqual(row["majority_linked_celltypes"], "B")
        self.assertEqual(row["majority_excluded_celltypes"], "")
        self.assertEqual(row["linked_celltypes"], "B")
        self.assertEqual(row["linked_celltype_counts"], "B:2")
        self.assertEqual(row["linked_celltype_fractions"], "B:1.000")
        self.assertFalse(bool(row["is_multi_celltype_link"]))
        self.assertAlmostEqual(float(row["majority_pct"]), 1.0)
        self.assertAlmostEqual(float(row["overlap_pct"]), 1.0)

    def test_code_schema_is_qualified_for_mixed_beds(self):
        sv_df = pd.DataFrame(
            {
                "chr": ["1"],
                "location": [1000],
                "id": ["sv1"],
                "supporting_reads": [["r1", "r2"]],
            }
        )

        read_assignment_df = pd.DataFrame(
            {
                "code_order": ["A|B", "C|D"],
                "code": ["10", "10"],
                "best_group": ["A", "C"],
                "is_best_group": [True, True],
                "chr": ["1", "1"],
                "start": [100, 200],
                "end": [150, 250],
            },
            index=pd.Index(["r1", "r2"], name="readname"),
        )

        out = assign_sv_celltypes(sv_df, read_assignment_df)
        row = out.iloc[0]

        self.assertEqual(row["n_supporting"], 2)
        self.assertEqual(row["n_overlapped"], 2)
        self.assertEqual(row["majority_code"], "A|B::10")
        self.assertEqual(row["assigned_code"], "A|B::10")
        self.assertEqual(row["code_counts"], "A|B::10:1;C|D::10:1")
        self.assertEqual(row["majority_schema"], "")
        self.assertEqual(row["primary_celltype"], "A")
        self.assertEqual(row["majority_linked_celltypes"], "A")
        self.assertEqual(row["majority_excluded_celltypes"], "")
        self.assertEqual(row["linked_celltypes"], "A|C")
        self.assertEqual(row["linked_celltype_counts"], "A:1;C:1")
        self.assertEqual(row["linked_celltype_fractions"], "A:0.500;C:0.500")
        self.assertTrue(bool(row["is_multi_celltype_link"]))

    def test_strict_majority_can_still_assign_with_partial_support_overlap(self):
        sv_df = pd.DataFrame(
            {
                "chr": ["1"],
                "location": [1000],
                "id": ["sv1"],
                "supporting_reads": [["r1", "r2", "r3", "r4", "r5"]],
            }
        )

        read_assignment_df = pd.DataFrame(
            {
                "sv_id": ["sv1", "sv1", "sv1", "sv1"],
                "code_order": ["A|B"] * 4,
                "code": ["10", "10", "10", "10"],
                "best_group": ["A", "A", "A", "A"],
                "is_best_group": [True, True, True, True],
                "chr": ["1", "1", "1", "1"],
                "start": [100, 200, 300, 400],
                "end": [150, 250, 350, 450],
            },
            index=pd.Index(["r1", "r2", "r3", "r4"], name="readname"),
        )

        out = assign_sv_celltypes(
            sv_df,
            read_assignment_df,
            unique_reads_for_overlap=True,
            min_overlap_pct=0.8,
            min_agreement_pct=1.0,
        )
        row = out.iloc[0]

        self.assertEqual(row["n_supporting"], 5)
        self.assertEqual(row["n_overlapped"], 4)
        self.assertAlmostEqual(float(row["overlap_pct"]), 0.8, places=6)
        self.assertEqual(row["majority_code"], "10")
        self.assertAlmostEqual(float(row["majority_pct"]), 1.0, places=6)
        self.assertEqual(row["assigned_code"], "10")

    def test_no_usable_evidence_sets_multi_link_na(self):
        sv_df = pd.DataFrame(
            {
                "chr": ["1"],
                "location": [1000],
                "id": ["sv1"],
                "supporting_reads": [["r1"]],
            }
        )

        # read assignment exists but does not overlap supporting read r1
        read_assignment_df = pd.DataFrame(
            {
                "code_order": ["A|B"],
                "code": ["10"],
                "best_group": ["A"],
                "is_best_group": [True],
                "chr": ["1"],
                "start": [100],
                "end": [150],
            },
            index=pd.Index(["rX"], name="readname"),
        )

        out = assign_sv_celltypes(sv_df, read_assignment_df)
        row = out.iloc[0]
        self.assertTrue(pd.isna(row["is_multi_celltype_link"]))
        self.assertTrue(pd.isna(row["has_hard_conflict"]))

    def test_non_supporting_reads_do_not_contribute_to_evidence(self):
        sv_df = pd.DataFrame(
            {
                "chr": ["1"],
                "location": [1000],
                "id": ["sv1"],
                "supporting_reads": [["r1"]],
            }
        )

        read_assignment_df = pd.DataFrame(
            {
                "code_order": ["A|B", "A|B"],
                "code": ["10", "01"],
                "best_group": ["A", "B"],
                "is_best_group": [True, True],
                "chr": ["1", "1"],
                "start": [100, 100],
                "end": [150, 150],
            },
            index=pd.Index(["r1", "rX"], name="readname"),
        )

        out = assign_sv_celltypes(sv_df, read_assignment_df)
        row = out.iloc[0]

        self.assertEqual(row["n_supporting"], 1)
        self.assertEqual(row["n_overlapped"], 1)
        self.assertEqual(row["majority_code"], "10")
        self.assertEqual(row["assigned_code"], "10")
        self.assertEqual(row["code_counts"], "10:1")

    def test_other_cluster_in_hierarchy_is_tracked_by_code(self):
        sv_df = pd.DataFrame(
            {
                "chr": ["1"],
                "location": [1000],
                "id": ["sv1"],
                "supporting_reads": [["r1"]],
            }
        )

        # Region best_group is lymphocytes but this read was assigned to the opposite cluster.
        # With hierarchical code_order/bits, this should decode to Monocyte.
        read_assignment_df = pd.DataFrame(
            {
                "code_order": ["T-cell|NK-cell|B-cell|Monocyte"],
                "code": ["0001"],
                "best_group": ["lymphocytes"],
                "is_best_group": [False],
                "chr": ["1"],
                "start": [100],
                "end": [150],
            },
            index=pd.Index(["r1"], name="readname"),
        )

        out = assign_sv_celltypes(sv_df, read_assignment_df)
        row = out.iloc[0]

        self.assertEqual(row["primary_celltype"], "Monocyte")
        self.assertEqual(row["linked_celltypes"], "Monocyte")
        self.assertEqual(row["linked_celltype_counts"], "Monocyte:1")
        self.assertEqual(row["linked_celltype_fractions"], "Monocyte:1.000")
        self.assertFalse(bool(row["is_multi_celltype_link"]))

    def test_zero_padded_code_is_decoded_for_other_cluster(self):
        sv_df = pd.DataFrame(
            {
                "chr": ["1"],
                "location": [1000],
                "id": ["sv1"],
                "supporting_reads": [["r1"]],
            }
        )

        # Simulate lost leading zero during csv parsing: 0111 -> 111 (int).
        read_assignment_df = pd.DataFrame(
            {
                "code_order": ["A|B|C|D"],
                "code": [111],
                "best_group": ["A"],
                "is_best_group": [False],
                "chr": ["1"],
                "start": [100],
                "end": [150],
            },
            index=pd.Index(["r1"], name="readname"),
        )

        out = assign_sv_celltypes(sv_df, read_assignment_df)
        row = out.iloc[0]

        self.assertEqual(row["majority_code"], "0111")
        self.assertEqual(row["code_counts"], "0111:1")
        self.assertEqual(row["linked_celltypes"], "B|C|D")
        self.assertEqual(row["linked_celltype_counts"], "B:1;C:1;D:1")
        self.assertTrue(bool(row["is_multi_celltype_link"]))

    def test_ambiguous_leaf_subset_is_reported_as_explicit_celltypes(self):
        sv_df = pd.DataFrame(
            {
                "chr": ["1"],
                "location": [1000],
                "id": ["sv1"],
                "supporting_reads": [["r1"]],
            }
        )

        # Ambiguous subset should stay explicit as linked leaf cell types.
        read_assignment_df = pd.DataFrame(
            {
                "code_order": ["T-cell|NK-cell|B-cell|Monocyte"],
                "code": ["1010"],
                "best_group": ["lymphocytes"],
                "best_group_leaves": ["T-cell|NK-cell|B-cell"],
                "is_best_group": [True],
                "chr": ["1"],
                "start": [100],
                "end": [150],
            },
            index=pd.Index(["r1"], name="readname"),
        )

        out = assign_sv_celltypes(sv_df, read_assignment_df)
        row = out.iloc[0]

        self.assertEqual(row["primary_celltype"], "B-cell")
        self.assertEqual(set(row["linked_celltypes"].split("|")), {"T-cell", "B-cell"})
        self.assertEqual(set(row["linked_celltype_counts"].split(";")), {"T-cell:1", "B-cell:1"})
        self.assertEqual(
            set(row["linked_celltype_fractions"].split(";")),
            {"T-cell:1.000", "B-cell:1.000"},
        )
        self.assertTrue(bool(row["is_multi_celltype_link"]))

    def test_all_rows_with_strict_agreement_marks_conflict_unassigned(self):
        sv_df = pd.DataFrame(
            {
                "chr": ["1"],
                "location": [1000],
                "id": ["sv1"],
                "supporting_reads": [["r1", "r2"]],
            }
        )

        # Two reads, each carrying conflicting region-level codes (2x "10", 1x "01").
        read_assignment_df = pd.DataFrame(
            {
                "code_order": ["A|B"] * 6,
                "code": ["10", "10", "01", "10", "10", "01"],
                "best_group": ["A", "A", "B", "A", "A", "B"],
                "is_best_group": [True, True, True, True, True, True],
                "chr": ["1"] * 6,
                "start": [100, 200, 300, 100, 200, 300],
                "end": [150, 250, 350, 150, 250, 350],
            },
            index=pd.Index(["r1", "r1", "r1", "r2", "r2", "r2"], name="readname"),
        )

        out = assign_sv_celltypes(
            sv_df,
            read_assignment_df,
            unique_reads_for_overlap=False,
            min_agreement_pct=1.0,
        )
        row = out.iloc[0]

        self.assertEqual(row["n_supporting"], 2)
        self.assertEqual(row["n_overlapped"], 2)
        self.assertEqual(row["majority_code"], "10")
        self.assertAlmostEqual(float(row["majority_pct"]), 4.0 / 6.0, places=6)
        self.assertTrue(pd.isna(row["assigned_code"]))

    def test_hard_conflict_is_detected_and_blocks_assigned_code(self):
        sv_df = pd.DataFrame(
            {
                "chr": ["1"],
                "location": [1000],
                "id": ["sv1"],
                "supporting_reads": [["r1", "r2", "r3", "r4"]],
            }
        )
        read_assignment_df = pd.DataFrame(
            {
                "code_order": ["A|B|C|D"] * 4,
                "code": ["1110", "0001", "1101", "0010"],
                "best_group": ["A", "D", "A", "B"],
                "is_best_group": [True, True, True, True],
                "chr": ["1"] * 4,
                "start": [100, 200, 300, 400],
                "end": [150, 250, 350, 450],
            },
            index=pd.Index(["r1", "r2", "r3", "r4"], name="readname"),
        )

        out = assign_sv_celltypes(
            sv_df,
            read_assignment_df,
            unique_reads_for_overlap=False,
            min_agreement_pct=0.0,
        )
        row = out.iloc[0]
        self.assertTrue(bool(row["has_hard_conflict"]))
        self.assertTrue(pd.isna(row["assigned_code"]))

    def test_compatible_negative_constraints_are_not_hard_conflict(self):
        sv_df = pd.DataFrame(
            {
                "chr": ["1"],
                "location": [1000],
                "id": ["sv1"],
                "supporting_reads": [["r1", "r2"]],
            }
        )
        read_assignment_df = pd.DataFrame(
            {
                "code_order": ["A|B|C|D", "A|B|C|D"],
                "code": ["1110", "1101"],
                "best_group": ["D", "C"],
                "is_best_group": [False, False],
                "chr": ["1", "1"],
                "start": [100, 200],
                "end": [150, 250],
            },
            index=pd.Index(["r1", "r2"], name="readname"),
        )

        out = assign_sv_celltypes(
            sv_df,
            read_assignment_df,
            unique_reads_for_overlap=False,
            min_agreement_pct=0.0,
        )
        row = out.iloc[0]
        self.assertFalse(bool(row["has_hard_conflict"]))
        self.assertEqual(str(row["intersection_code"]), "1100")


class TestParseArgs(unittest.TestCase):
    def test_anno_accepts_single_bed_path(self):
        args = parse_args(
            [
                "anno",
                "-i",
                "input.bam",
                "-v",
                "input.vcf.gz",
                "-r",
                "ref.fa",
                "-b",
                "a.bed.tsv",
                "-o",
                "out_dir",
            ]
        )
        self.assertEqual(args.bed, "a.bed.tsv")

    def test_viz_subcommand_parses(self):
        args = parse_args(
            [
                "viz",
                "-i",
                "input.bam",
                "-v",
                "input.vcf.gz",
                "-s",
                "sv123",
                "-o",
                "out/sv123",
            ]
        )
        self.assertEqual(args.command, "viz")
        self.assertEqual(args.sv_id, "sv123")
        self.assertEqual(args.format, "png")
        self.assertIsNone(args.read_assignment)
        self.assertIsNone(args.reference)
        self.assertFalse(args.export_tables)

    def test_viz_subcommand_parses_with_anno_output_only(self):
        args = parse_args(
            [
                "viz",
                "--anno_output",
                "anno_out",
                "-s",
                "sv123",
            ]
        )
        self.assertEqual(args.command, "viz")
        self.assertEqual(args.anno_output, "anno_out")
        self.assertEqual(args.sv_id, "sv123")
        self.assertIsNone(args.input)
        self.assertIsNone(args.vcf)
        self.assertIsNone(args.output)

    def test_igvviz_subcommand_parses_with_multiple_bams(self):
        args = parse_args(
            [
                "igvviz",
                "-i",
                "a.bam",
                "b.bam",
                "-v",
                "input.vcf.gz",
                "-s",
                "sv123",
                "-o",
                "out_dir",
            ]
        )
        self.assertEqual(args.command, "igvviz")
        self.assertEqual(args.input, ["a.bam", "b.bam"])
        self.assertEqual(args.vcf, "input.vcf.gz")
        self.assertEqual(args.sv_id, "sv123")
        self.assertEqual(args.output, "out_dir")
        self.assertEqual(args.support_tag, "SC")
        self.assertEqual(args.phase_tag, "HP")
        self.assertEqual(args.snapshot_width, 3600)
        self.assertEqual(args.snapshot_height, 1600)
        self.assertFalse(args.batch_only)

    def test_igvviz_subcommand_parses_with_anno_output_only(self):
        args = parse_args(
            [
                "igvviz",
                "--anno_output",
                "anno_out",
                "-s",
                "sv123",
                "--batch_only",
            ]
        )
        self.assertEqual(args.command, "igvviz")
        self.assertEqual(args.anno_output, "anno_out")
        self.assertEqual(args.sv_id, "sv123")
        self.assertIsNone(args.input)
        self.assertIsNone(args.vcf)
        self.assertTrue(args.batch_only)

    def test_anno_assignment_defaults_are_strict_all_rows(self):
        args = parse_args(
            [
                "anno",
                "-i",
                "input.bam",
                "-v",
                "input.vcf.gz",
                "-r",
                "ref.fa",
                "-b",
                "a.bed.tsv",
                "-o",
                "out_dir",
            ]
        )
        self.assertEqual(args.read_assignment_mode, "closest_reference_mean")
        self.assertEqual(args.evidence_mode, "all_rows")
        self.assertAlmostEqual(args.min_overlap_pct, 0.0)
        self.assertAlmostEqual(args.min_agreement_pct, 1.0)

    def test_anno_assignment_mode_accepts_kmeans(self):
        args = parse_args(
            [
                "anno",
                "-i",
                "input.bam",
                "-v",
                "input.vcf.gz",
                "-r",
                "ref.fa",
                "-b",
                "a.bed.tsv",
                "-o",
                "out_dir",
                "--read_assignment_mode",
                "kmeans",
            ]
        )
        self.assertEqual(args.read_assignment_mode, "kmeans")

    def test_svanno_assignment_defaults_are_strict_all_rows(self):
        args = parse_args(
            [
                "svanno",
                "-v",
                "input.vcf.gz",
                "-i",
                "anno_out/reads_classification.tsv",
                "-o",
                "anno_out",
            ]
        )
        self.assertEqual(args.evidence_mode, "all_rows")
        self.assertAlmostEqual(args.min_overlap_pct, 0.0)
        self.assertAlmostEqual(args.min_agreement_pct, 1.0)

    def test_report_subcommand_parses_with_defaults(self):
        args = parse_args(
            [
                "report",
                "--anno_output",
                "anno_out",
            ]
        )
        self.assertEqual(args.command, "report")
        self.assertEqual(args.anno_output, "anno_out")
        self.assertAlmostEqual(args.min_overlap_pct, 0.8)
        self.assertAlmostEqual(args.min_majority_pct, 1.0)
        self.assertFalse(args.include_unassigned)
        self.assertFalse(args.allow_hard_conflict)
        self.assertEqual(args.max_sv, 0)
        self.assertFalse(args.with_figures)
        self.assertEqual(args.figure_threads, 1)
        self.assertFalse(args.with_igvviz)
        self.assertEqual(args.igv_cmd, "igv.sh")
        self.assertEqual(args.igv_snapshot_format, "png")
        self.assertEqual(args.igv_snapshot_width, 3600)
        self.assertEqual(args.igv_snapshot_height, 1600)
        self.assertEqual(args.format, "png")
        self.assertFalse(hasattr(args, "export_tables"))
        self.assertFalse(hasattr(args, "igv_threads"))
        self.assertFalse(hasattr(args, "igv_visibility_window"))
        self.assertFalse(hasattr(args, "igv_phase_tag"))
        self.assertFalse(hasattr(args, "igv_support_tag"))
        self.assertFalse(hasattr(args, "igv_keep_intermediates"))
        self.assertFalse(hasattr(args, "igv_batch_only"))
        self.assertFalse(hasattr(args, "input"))
        self.assertFalse(hasattr(args, "vcf"))
        self.assertFalse(hasattr(args, "reference"))
        self.assertFalse(hasattr(args, "bed"))
        self.assertFalse(hasattr(args, "read_assignment"))


class TestReadableSvReport(unittest.TestCase):
    def test_build_readable_reports_summarizes_celltypes(self):
        sv_df = pd.DataFrame(
            [
                {
                    "id": "sv1",
                    "sv_chr": "1",
                    "sv_pos": 1000,
                    "sv_len": 200,
                    "vaf": 0.25,
                    "n_supporting": 8,
                    "n_overlapped": 6,
                    "overlap_pct": 0.75,
                    "majority_pct": 0.67,
                    "linked_celltypes": "B-cell|T-cell",
                    "linked_celltype_counts": "B-cell:4;T-cell:2",
                    "linked_celltype_fractions": "B-cell:0.667;T-cell:0.333",
                    "is_multi_celltype_link": True,
                }
            ]
        )

        summary_df, long_df = _build_sv_readable_reports(sv_df)

        self.assertEqual(len(summary_df), 1)
        self.assertEqual(summary_df.iloc[0]["classified_celltypes"], "B-cell|T-cell")
        self.assertEqual(summary_df.iloc[0]["classified_celltype_count"], 2)
        self.assertIn("B-cell (4 reads, 0.667)", summary_df.iloc[0]["classification_summary"])
        self.assertEqual(len(long_df), 2)
        self.assertEqual(long_df.iloc[0]["celltype"], "B-cell")
        self.assertEqual(long_df.iloc[0]["supporting_read_count"], 4)


if __name__ == "__main__":
    unittest.main()
