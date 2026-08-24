import unittest
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd
import pysam

from sniffcell.anno.anno import _run_batch_annotation, _run_compact_annotation
from sniffcell.find.ctdmr import call_ct_combination_dmrs
from sniffcell.report.report import _load_excluded_variant_ids, _read_table, report_main
from sniffcell.tissue_atlas import resolve_tissue_key


class TestFindCombinations(unittest.TestCase):
    def test_call_ct_combination_dmrs_detects_multi_group_patterns(self):
        idx_df = pd.DataFrame(
            {
                "chr": ["1", "1", "1", "1"],
                "start": [100, 200, 300, 400],
                "end": [150, 250, 350, 450],
                "startCpG": [0, 2, 4, 6],
                "endCpG": [2, 4, 6, 8],
            }
        )

        mean_by_group = {
            "A": pd.Series([0.82, 0.80, 0.20, 0.18]),
            "B": pd.Series([0.78, 0.76, 0.22, 0.20]),
            "C": pd.Series([0.20, 0.22, 0.82, 0.80]),
            "D": pd.Series([0.18, 0.20, 0.78, 0.76]),
        }

        dmrs = call_ct_combination_dmrs(
            idx_df=idx_df,
            mean_by_group=mean_by_group,
            diff_threshold=0.40,
            min_rows=2,
            min_cpgs=2,
            min_bp=0,
            direction="both",
            max_gap_bp=100,
            bed_out=None,
        )

        self.assertEqual(len(dmrs), 2)
        first = dmrs.iloc[0]
        second = dmrs.iloc[1]

        self.assertEqual(first["best_group_leaves"], "A|B")
        self.assertEqual(first["other_group_leaves"], "C|D")
        self.assertEqual(first["best_dir"], "hyper")
        self.assertEqual(first["code_order"], "A|B|C|D")

        self.assertEqual(second["best_group_leaves"], "C|D")
        self.assertEqual(second["other_group_leaves"], "A|B")
        self.assertEqual(second["best_dir"], "hyper")
        self.assertEqual(second["code_order"], "A|B|C|D")


class TestTissueAtlasHelpers(unittest.TestCase):
    def test_resolve_tissue_key_by_code_or_name(self):
        atlas = {
            "__tissues__": [
                {
                    "code": "3E",
                    "name": "Colon, Ascending",
                    "key": "3E",
                    "usable": True,
                }
            ],
            "3E": {"A": ["sample_a"], "B": ["sample_b"]},
        }
        by_code, _ = resolve_tissue_key("3E", atlas)
        by_name, _ = resolve_tissue_key("colon ascending", atlas)

        self.assertEqual(by_code, "3E")
        self.assertEqual(by_name, "3E")


class TestCompactAnno(unittest.TestCase):
    def _write_tiny_bam(self, path):
        header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 10000}]}
        with pysam.AlignmentFile(path, "wb", header=header) as bam:
            for i, read_name in enumerate(["read1", "read2", "read_other"]):
                aln = pysam.AlignedSegment()
                aln.query_name = read_name
                aln.query_sequence = "A" * 3000
                aln.flag = 0
                aln.reference_id = 0
                aln.reference_start = 1000 + i * 100
                aln.mapping_quality = 60
                aln.cigar = [(0, 3000)]
                aln.query_qualities = pysam.qualitystring_to_array("I" * 3000)
                bam.write(aln)

    def test_compact_annotation_uses_catalog_and_supporting_reads(self):
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            catalog = tmp / "reads_classification.tsv"
            ctdmr_df = pd.DataFrame(
                {
                    "chr": ["chr1"],
                    "start": [3500],
                    "end": [3600],
                    "best_group": ["A"],
                    "other_group": ["B"],
                    "best_dir": ["hyper"],
                    "mean_best_value": [0.85],
                    "mean_rest_value": [0.15],
                    "best_group_leaves": ["A"],
                    "other_group_leaves": ["B"],
                    "hyper_group_leaves": ["A"],
                    "hypo_group_leaves": ["B"],
                    "code_order": ["A|B"],
                }
            )
            ctdmr_df.to_csv(catalog, sep="\t", index=False)
            bam_path = tmp / "reads.bam"
            self._write_tiny_bam(bam_path)

            outdir = tmp / "anno"
            args = SimpleNamespace(
                command="anno",
                output=str(outdir),
                variant_name="var1",
                variant_location="chr1:1,000-1,100",
                supporting_reads="read1,read2,missing_read",
                catalog=str(catalog),
                input=str(bam_path),
                reference=str(tmp / "ref.fa"),
                kanpig_read_names=None,
                window=100,
                breakpoint_exclusion_frac=0.0,
                evidence_mode="per_read",
                min_overlap_pct=0.0,
                min_agreement_pct=0.0,
                per_read_min_agreement=0.66,
            )

            mm = pd.DataFrame(
                [[0.90], [0.82]],
                index=pd.MultiIndex.from_arrays(
                    [["read1", "read2"], [1, 1]],
                    names=["read_name", "haplotype"],
                ),
                columns=[3550],
            )
            with patch("sniffcell.anno.anno.methyl_matrix_from_bam", return_value=(mm, [850])) as mocked:
                _run_compact_annotation(args)
                mocked.assert_called_once()
            assignment = pd.read_csv(outdir / "variant_assignment.tsv", sep="\t")
            filtered_reads = pd.read_csv(outdir / "reads_classification.tsv", sep="\t")
            mappings = pd.read_csv(outdir / "support_read_mappings.tsv", sep="\t")

        self.assertEqual(len(filtered_reads), 2)
        self.assertEqual(assignment.loc[0, "id"], "var1")
        self.assertEqual(assignment.loc[0, "linked_celltypes"], "A")
        self.assertEqual(int(assignment.loc[0, "n_supporting"]), 3)
        self.assertEqual(int(assignment.loc[0, "n_overlapped"]), 2)
        self.assertEqual(len(mappings), 3)
        self.assertEqual(bool(mappings.loc[mappings["readname"] == "read1", "is_mapped"].iloc[0]), True)

    def test_batch_annotation_uses_bam_and_targeted_ctdmrs(self):
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            catalog = tmp / "ctdmr.tsv"
            pd.DataFrame(
                {
                    "chr": ["chr1"],
                    "start": [3500],
                    "end": [3600],
                    "best_group": ["A"],
                    "other_group": ["B"],
                    "mean_best_value": [0.85],
                    "mean_rest_value": [0.15],
                    "best_group_leaves": ["A"],
                    "other_group_leaves": ["B"],
                    "hyper_group_leaves": ["A"],
                    "hypo_group_leaves": ["B"],
                    "code_order": ["A|B"],
                }
            ).to_csv(catalog, sep="\t", index=False)
            bam_path = tmp / "reads.bam"
            self._write_tiny_bam(bam_path)
            batch = tmp / "batch.tsv"
            pd.DataFrame(
                [
                    {
                        "variant_name": "var1",
                        "variant_location": "chr1:1,000-1,100",
                        "supporting_reads": "read1,read2",
                        "catalog": str(catalog),
                        "bam": str(bam_path),
                    }
                ]
            ).to_csv(batch, sep="\t", index=False)
            outdir = tmp / "batch_out"
            args = SimpleNamespace(
                command="anno",
                output=str(outdir),
                batch=str(batch),
                reference=str(tmp / "ref.fa"),
                window=100,
                breakpoint_exclusion_frac=0.0,
                evidence_mode="per_read",
                min_overlap_pct=0.0,
                min_agreement_pct=0.0,
                per_read_min_agreement=0.66,
            )
            mm = pd.DataFrame(
                [[0.90], [0.82]],
                index=pd.MultiIndex.from_arrays(
                    [["read1", "read2"], [1, 1]],
                    names=["read_name", "haplotype"],
                ),
                columns=[3550],
            )
            with patch("sniffcell.anno.anno.methyl_matrix_from_bam", return_value=(mm, [850])) as mocked:
                _run_batch_annotation(args)
                mocked.assert_called_once()
            assignment = pd.read_csv(outdir / "variant_assignment.tsv", sep="\t")
            mappings = pd.read_csv(outdir / "support_read_mappings.tsv", sep="\t")

        self.assertEqual(assignment.loc[0, "id"], "var1")
        self.assertEqual(assignment.loc[0, "linked_celltypes"], "A")
        self.assertEqual(len(mappings), 2)


class TestLiteReport(unittest.TestCase):
    def test_report_adapter_preserves_leading_zero_assignment_codes(self):
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmpdir:
            table_path = Path(tmpdir) / "assignments.tsv"
            pd.DataFrame(
                {
                    "id": ["var1", "var2"],
                    "majority_code": ["001", "011"],
                    "assigned_code": ["001", "011"],
                    "intersection_code": ["001", "011"],
                }
            ).to_csv(table_path, sep="\t", index=False)

            assignments = _read_table(table_path)

        self.assertEqual(assignments["assigned_code"].tolist(), ["001", "011"])
        self.assertEqual(assignments["majority_code"].tolist(), ["001", "011"])

    def test_exclusion_table_uses_id_column_and_deduplicates(self):
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmpdir:
            exclusion_file = Path(tmpdir) / "excluded.tsv"
            pd.DataFrame(
                {
                    "id": ["var2", "var1", "var2"],
                    "tissue_code": ["3E", "3AD", "3E"],
                }
            ).to_csv(exclusion_file, sep="\t", index=False)

            excluded = _load_excluded_variant_ids(exclusion_file)

        self.assertEqual(excluded, ["var2", "var1"])

    def test_empty_exclusion_table_is_rejected(self):
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmpdir:
            exclusion_file = Path(tmpdir) / "excluded.tsv"
            exclusion_file.write_text("id\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "contains no variant IDs"):
                _load_excluded_variant_ids(exclusion_file)

    def test_report_parser_rejects_include_and_exclude_together(self):
        from sniffcell.parse_args import parse_args

        with self.assertRaises(SystemExit):
            parse_args(
                [
                    "report",
                    "--anno_output",
                    "anno",
                    "--include_variants",
                    "include.tsv",
                    "--exclude_variants",
                    "exclude.tsv",
                ]
            )

    def test_report_writes_html_and_figures_from_lite_batch_output(self):
        import json
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            anno = tmp / "anno"
            anno.mkdir()
            catalog = tmp / "ctdmr.tsv"
            batch = tmp / "batch.tsv"

            pd.DataFrame(
                [
                    {
                        "chr": "chr1",
                        "start": 900,
                        "end": 1200,
                        "best_group": "A",
                        "best_dir": "hyper",
                        "best_group_leaves": "A",
                        "other_group_leaves": "B",
                        "code_order": "A|B",
                        "mean_best_value": 0.85,
                        "mean_rest_value": 0.15,
                    }
                ]
            ).to_csv(catalog, sep="\t", index=False)
            pd.DataFrame(
                [
                    {
                        "variant_name": "var1",
                        "variant_location": "chr1:1001-1001",
                        "supporting_reads": json.dumps(["read1"]),
                        "catalog": str(catalog),
                        "bam": str(tmp / "reads.bam"),
                        "reference": str(tmp / "ref.fa"),
                        "donor": "D1",
                        "tissue_code": "3E",
                        "tissue_name": "Colon, Ascending",
                    }
                ]
            ).to_csv(batch, sep="\t", index=False)
            pd.DataFrame(
                [
                    {
                        "id": "var1",
                        "n_supporting": 1,
                        "n_overlapped": 1,
                        "n_mixed_reads": 0,
                        "overlap_pct": 1.0,
                        "majority_code": "10",
                        "majority_pct": 1.0,
                        "assigned_code": "A|B::10",
                        "linked_celltypes": "A",
                        "linked_celltype_counts": "A:1",
                        "linked_celltype_fractions": "A:1.0",
                        "has_hard_conflict": False,
                        "sv_chr": "chr1",
                        "sv_pos": 1001,
                        "sv_len": 1,
                    }
                ]
            ).to_csv(anno / "variant_assignment.tsv", sep="\t", index=False)
            pd.DataFrame(
                [
                    {
                        "readname": "var1::read1",
                        "chr": "chr1",
                        "start": 900,
                        "end": 1200,
                        "cpgstart": 950,
                        "cpgend": 1150,
                        "best_group": "A",
                        "other_group": "B",
                        "is_best_group": True,
                        "code_order": "A|B",
                        "best_group_leaves": "A",
                        "other_group_leaves": "B",
                        "hyper_group_leaves": "A",
                        "hypo_group_leaves": "B",
                        "code": "10",
                    }
                ]
            ).to_csv(anno / "reads_classification.tsv", sep="\t", index=False)
            pd.DataFrame(
                [
                    {
                        "variant_id": "var1",
                        "readname": "read1",
                        "internal_readname": "var1::read1",
                        "bam": str(tmp / "reads.bam"),
                        "chrom": "chr1",
                        "start": 800,
                        "end": 1300,
                        "mapping_quality": 60,
                        "is_reverse": False,
                        "is_mapped": True,
                    }
                ]
            ).to_csv(anno / "support_read_mappings.tsv", sep="\t", index=False)
            (anno / "anno_batch_manifest.json").write_text(
                json.dumps(
                    {
                        "command": "anno",
                        "mode": "batch",
                        "inputs": {"batch": str(batch), "row_count": 1, "reference": str(tmp / "ref.fa")},
                        "outputs": {},
                    }
                )
                + "\n",
                encoding="utf-8",
            )

            report_dir = tmp / "report"
            args = SimpleNamespace(
                anno_output=str(anno),
                output=str(report_dir),
                include_unassigned=False,
                min_overlap_pct=0.0,
                min_majority_pct=0.0,
                allow_hard_conflict=False,
                max_variants=0,
                with_figures=False,
                window=100,
                figure_format="png",
                figure_dpi=80,
            )
            report_main(args)

            manifest = json.loads((report_dir / "report_manifest.json").read_text(encoding="utf-8"))
            index_exists = (report_dir / "index.html").exists()
            selected_exists = (report_dir / "high_confidence_variants.tsv").exists()
            compat_manifest_exists = (anno / ".sniffcell_full_report_compat" / "anno_run_manifest.json").exists()
            compat_variant_table_exists = (anno / ".sniffcell_full_report_compat" / "lite_variants.tsv").exists()
            figure_count = len(list((report_dir / "figures").glob("*.png")))

        self.assertEqual(manifest["counts"]["variant_selected"], 1)
        self.assertEqual(manifest["counts"]["viz_rendered_or_reused"], 0)
        self.assertTrue(index_exists)
        self.assertTrue(selected_exists)
        self.assertTrue(compat_manifest_exists)
        self.assertTrue(compat_variant_table_exists)
        self.assertEqual(figure_count, 0)


if __name__ == "__main__":
    unittest.main()
