import json
import logging
import tarfile
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd

from sniffcell.report.report import (
    _backfill_sv_fields_from_manifest_vcf,
    _load_sv_assignment,
    _select_high_confidence_svs,
    _viz_supported_for_row,
    report_main,
)


class TestReportSelection(unittest.TestCase):
    def test_load_assignment_preserves_leading_zero_codes(self):
        with tempfile.TemporaryDirectory() as td:
            assignment_path = Path(td) / "variant_assignment.tsv"
            pd.DataFrame(
                {
                    "id": ["var1", "var2"],
                    "majority_code": ["001", "011"],
                    "assigned_code": ["001", "011"],
                }
            ).to_csv(assignment_path, sep="\t", index=False)

            assignments = _load_sv_assignment(assignment_path)

        self.assertEqual(assignments["assigned_code"].tolist(), ["001", "011"])
        self.assertEqual(assignments["majority_code"].tolist(), ["001", "011"])

    def test_viz_supports_generic_and_mei_variants(self):
        self.assertTrue(_viz_supported_for_row({"variant_class": "VAR"}, {}))
        self.assertTrue(_viz_supported_for_row({"variant_class": "MEI"}, {}))

    def test_non_vcf_variant_table_is_not_parsed_as_vcf(self):
        with tempfile.TemporaryDirectory() as td:
            variant_table = Path(td) / "variants.tsv"
            variant_table.write_text("id\nvar1\n", encoding="utf-8")
            variants = pd.DataFrame(
                {
                    "id": pd.Series(["var1"], dtype="string"),
                    "sv_type": pd.Series([""], dtype="string"),
                    "vaf": [pd.NA],
                    "sv_len": [pd.NA],
                }
            )
            with patch("sniffcell.report.report.read_vcf_to_df") as read_vcf:
                result = _backfill_sv_fields_from_manifest_vcf(
                    variants,
                    {"inputs": {"vcf": str(variant_table)}},
                    logging.getLogger("test"),
                )

        read_vcf.assert_not_called()
        self.assertEqual(result.loc[0, "id"], "var1")

    def test_select_high_confidence_defaults(self):
        sv_df = pd.DataFrame(
            {
                "id": ["sv1", "sv2", "sv3", "sv4"],
                "assigned_code": ["10", "", "01", "11"],
                "linked_celltypes": ["A", "A", "B", ""],
                "has_hard_conflict": [False, False, True, False],
                "overlap_pct": [0.8, 0.9, 1.0, 1.0],
                "majority_pct": [0.7, 0.9, 0.95, 0.9],
                "n_overlapped": [6, 8, 10, 10],
            }
        )
        sv_df["has_hard_conflict"] = sv_df["has_hard_conflict"].astype("boolean")

        selected = _select_high_confidence_svs(
            sv_df,
            min_overlap_pct=0.0,
            overlap_filter_mode="gradient",
            overlap_gradient_exponent=0.5,
            min_majority_pct=0.0,
            include_unassigned=False,
            allow_hard_conflict=False,
            max_sv=0,
        )

        self.assertEqual(selected["id"].tolist(), ["sv1"])

    def test_select_high_confidence_with_thresholds_and_limit(self):
        sv_df = pd.DataFrame(
            {
                "id": ["svA", "svB", "svC"],
                "assigned_code": ["10", "11", "01"],
                "linked_celltypes": ["A", "B", "C"],
                "has_hard_conflict": [False, False, False],
                "overlap_pct": [0.4, 0.9, 0.8],
                "majority_pct": [0.99, 0.85, 0.95],
                "n_overlapped": [4, 9, 8],
            }
        )
        sv_df["has_hard_conflict"] = sv_df["has_hard_conflict"].astype("boolean")

        selected = _select_high_confidence_svs(
            sv_df,
            min_overlap_pct=0.5,
            overlap_filter_mode="hard_clip",
            overlap_gradient_exponent=0.5,
            min_majority_pct=0.9,
            include_unassigned=False,
            allow_hard_conflict=False,
            max_sv=1,
        )

        self.assertEqual(selected["id"].tolist(), ["svC"])

    def test_select_high_confidence_gradient_allows_partial_support_with_more_reads(self):
        sv_df = pd.DataFrame(
            {
                "id": ["sv_partial", "sv_drop"],
                "assigned_code": ["10", "10"],
                "linked_celltypes": ["A", "A"],
                "has_hard_conflict": [False, False],
                "overlap_pct": [0.5, 0.5],
                "majority_pct": [1.0, 1.0],
                "n_supporting": [4, 2],
                "n_overlapped": [2, 1],
            }
        )
        sv_df["has_hard_conflict"] = sv_df["has_hard_conflict"].astype("boolean")

        selected = _select_high_confidence_svs(
            sv_df,
            min_overlap_pct=0.8,
            overlap_filter_mode="gradient",
            overlap_gradient_exponent=0.5,
            min_majority_pct=1.0,
            include_unassigned=False,
            allow_hard_conflict=False,
            max_sv=0,
        )

        self.assertEqual(selected["id"].tolist(), ["sv_partial"])
        self.assertEqual(int(selected.iloc[0]["n_overlapped"]), 2)
        self.assertEqual(int(selected.iloc[0]["n_supporting"]), 4)
        self.assertEqual(float(selected.iloc[0]["overlap_required_reads"]), 2.0)
        self.assertAlmostEqual(float(selected.iloc[0]["overlap_required_pct"]), 0.5, places=6)

    def test_select_high_confidence_hard_clip_preserves_fixed_threshold(self):
        sv_df = pd.DataFrame(
            {
                "id": ["sv_partial"],
                "assigned_code": ["10"],
                "linked_celltypes": ["A"],
                "has_hard_conflict": [False],
                "overlap_pct": [0.5],
                "majority_pct": [1.0],
                "n_supporting": [4],
                "n_overlapped": [2],
            }
        )
        sv_df["has_hard_conflict"] = sv_df["has_hard_conflict"].astype("boolean")

        selected = _select_high_confidence_svs(
            sv_df,
            min_overlap_pct=0.8,
            overlap_filter_mode="hard_clip",
            overlap_gradient_exponent=0.5,
            min_majority_pct=1.0,
            include_unassigned=False,
            allow_hard_conflict=False,
            max_sv=0,
        )

        self.assertEqual(len(selected), 0)

    def test_select_high_confidence_keeps_linked_tr_despite_hard_conflict(self):
        sv_df = pd.DataFrame(
            {
                "id": ["tr_keep", "tr_drop"],
                "variant_class": ["TR", "TR"],
                "assigned_code": ["", ""],
                "linked_celltypes": ["A|B", "A"],
                "has_hard_conflict": [True, True],
                "overlap_pct": [0.95, 0.0],
                "majority_pct": [0.8, 0.0],
                "n_supporting": [10, 5],
                "n_overlapped": [8, 0],
            }
        )
        sv_df["has_hard_conflict"] = sv_df["has_hard_conflict"].astype("boolean")

        selected = _select_high_confidence_svs(
            sv_df,
            min_overlap_pct=0.8,
            overlap_filter_mode="gradient",
            overlap_gradient_exponent=0.5,
            min_majority_pct=0.8,
            include_unassigned=False,
            allow_hard_conflict=False,
            max_sv=0,
        )

        self.assertEqual(selected["id"].tolist(), ["tr_keep"])


class TestReportMain(unittest.TestCase):
    def _base_args(self, anno_output: Path, output: str | None = None, reuse: bool = False):
        return SimpleNamespace(
            anno_output=str(anno_output),
            min_overlap_pct=0.0,
            overlap_filter_mode="gradient",
            overlap_gradient_exponent=0.5,
            min_majority_pct=0.0,
            include_unassigned=False,
            allow_hard_conflict=False,
            max_sv=0,
            with_figures=False,
            window=5000,
            max_reads=250,
            format="png",
            reuse_existing_viz=reuse,
            figure_threads=1,
            with_igvviz=False,
            igv_bams=None,
            igv_cmd="igv.sh",
            igv_snapshot_format="png",
            igv_snapshot_width=3600,
            igv_snapshot_height=1600,
            reuse_existing_igvviz=False,
            with_igvreport=False,
            output=output,
        )

    def test_report_main_writes_html_and_selected_table(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            (anno_dir / "anno_run_manifest.json").write_text(json.dumps({"inputs": {}, "outputs": {}}), encoding="utf-8")

            sv_df = pd.DataFrame(
                [
                    {
                        "id": "sv_pass",
                        "assigned_code": "10",
                        "linked_celltypes": "A",
                        "primary_celltype": "A",
                        "has_hard_conflict": False,
                        "overlap_pct": 0.8,
                        "majority_pct": 0.9,
                        "n_supporting": 8,
                        "n_overlapped": 6,
                    },
                    {
                        "id": "sv_drop",
                        "assigned_code": "",
                        "linked_celltypes": "B",
                        "primary_celltype": "B",
                        "has_hard_conflict": False,
                        "overlap_pct": 0.7,
                        "majority_pct": 0.8,
                        "n_supporting": 7,
                        "n_overlapped": 5,
                    },
                ]
            )
            sv_df.to_csv(anno_dir / "sv_assignment.tsv", sep="\t", index=False)

            args = self._base_args(anno_dir)
            args.with_figures = True

            def _fake_viz(viz_args):
                out = Path(viz_args.output)
                out.parent.mkdir(parents=True, exist_ok=True)
                out.write_text("fake figure", encoding="utf-8")

            with patch("sniffcell.report.report.viz_module.viz_main", side_effect=_fake_viz) as mock_viz:
                report_main(args)
                self.assertEqual(mock_viz.call_count, 1)
                self.assertEqual(mock_viz.call_args[0][0].sv_id, "sv_pass")

            report_dir = anno_dir / "report"
            html_path = report_dir / "index.html"
            selected_path = report_dir / "high_confidence_sv.tsv"
            self.assertTrue(html_path.exists())
            self.assertTrue(selected_path.exists())
            text = html_path.read_text(encoding="utf-8")
            self.assertIn("sv_pass", text)
            self.assertNotIn("sv_drop", text)

            selected = pd.read_csv(selected_path, sep="\t")
            self.assertEqual(selected["id"].tolist(), ["sv_pass"])
            self.assertEqual(selected["viz_status"].tolist(), ["rendered"])
            self.assertEqual(selected["overlap_filter_mode"].tolist(), ["gradient"])

    def test_report_main_reuses_existing_viz(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            (anno_dir / "anno_run_manifest.json").write_text(json.dumps({"inputs": {}, "outputs": {}}), encoding="utf-8")

            sv_df = pd.DataFrame(
                [
                    {
                        "id": "sv1",
                        "assigned_code": "10",
                        "linked_celltypes": "A",
                        "primary_celltype": "A",
                        "has_hard_conflict": False,
                        "overlap_pct": 0.9,
                        "majority_pct": 1.0,
                        "n_supporting": 10,
                        "n_overlapped": 9,
                    }
                ]
            )
            sv_df.to_csv(anno_dir / "sv_assignment.tsv", sep="\t", index=False)

            report_dir = anno_dir / "custom_report"
            fig_path = report_dir / "figures" / "sv1.viz.png"
            fig_path.parent.mkdir(parents=True, exist_ok=True)
            fig_path.write_text("existing figure", encoding="utf-8")

            args = self._base_args(anno_dir, output=str(report_dir), reuse=True)
            args.with_figures = True

            with patch("sniffcell.report.report.viz_module.viz_main") as mock_viz:
                report_main(args)
                mock_viz.assert_not_called()

            selected = pd.read_csv(report_dir / "high_confidence_sv.tsv", sep="\t")
            self.assertEqual(selected["viz_status"].tolist(), ["reused"])

    def test_report_main_handles_empty_selection(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            (anno_dir / "anno_run_manifest.json").write_text(json.dumps({"inputs": {}, "outputs": {}}), encoding="utf-8")

            sv_df = pd.DataFrame(
                [
                    {
                        "id": "sv_none",
                        "assigned_code": "",
                        "linked_celltypes": "",
                        "primary_celltype": "",
                        "has_hard_conflict": False,
                        "overlap_pct": 0.1,
                        "majority_pct": 0.2,
                        "n_supporting": 3,
                        "n_overlapped": 1,
                    }
                ]
            )
            sv_df.to_csv(anno_dir / "sv_assignment.tsv", sep="\t", index=False)

            args = self._base_args(anno_dir)

            with patch("sniffcell.report.report.viz_module.viz_main") as mock_viz:
                report_main(args)
                mock_viz.assert_not_called()

            report_dir = anno_dir / "report"
            html_path = report_dir / "index.html"
            selected_path = report_dir / "high_confidence_sv.tsv"
            failed_path = report_dir / "failed_viz.tsv"
            self.assertTrue(html_path.exists())
            self.assertTrue(selected_path.exists())
            self.assertTrue(failed_path.exists())

            text = html_path.read_text(encoding="utf-8")
            self.assertIn("No variants passed the report filters", text)
            selected = pd.read_csv(selected_path, sep="\t")
            self.assertEqual(len(selected), 0)

    def test_report_main_writes_gz_archive_when_output_has_gz_suffix(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            (anno_dir / "anno_run_manifest.json").write_text(json.dumps({"inputs": {}, "outputs": {}}), encoding="utf-8")

            sv_df = pd.DataFrame(
                [
                    {
                        "id": "sv_archive",
                        "assigned_code": "10",
                        "linked_celltypes": "A",
                        "primary_celltype": "A",
                        "has_hard_conflict": False,
                        "overlap_pct": 0.9,
                        "majority_pct": 1.0,
                        "n_supporting": 9,
                        "n_overlapped": 9,
                    }
                ]
            )
            sv_df.to_csv(anno_dir / "sv_assignment.tsv", sep="\t", index=False)

            archive_output = anno_dir / "report_bundle.gz"
            args = self._base_args(anno_dir, output=str(archive_output))
            report_main(args)

            report_dir = anno_dir / "report_bundle"
            self.assertTrue(report_dir.exists())
            self.assertTrue((report_dir / "index.html").exists())
            self.assertTrue(archive_output.exists())
            with tarfile.open(archive_output, mode="r:gz") as tf:
                names = tf.getnames()
            self.assertIn("report_bundle/index.html", names)

    def test_report_main_defaults_to_figureless_and_writes_copy_commands(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            (anno_dir / "anno_run_manifest.json").write_text(json.dumps({"inputs": {}, "outputs": {}}), encoding="utf-8")

            sv_df = pd.DataFrame(
                [
                    {
                        "id": "sv_nofig",
                        "assigned_code": "10",
                        "linked_celltypes": "A",
                        "primary_celltype": "A",
                        "has_hard_conflict": False,
                        "overlap_pct": 0.8,
                        "majority_pct": 0.9,
                        "n_supporting": 8,
                        "n_overlapped": 6,
                    }
                ]
            )
            sv_df.to_csv(anno_dir / "sv_assignment.tsv", sep="\t", index=False)

            args = self._base_args(anno_dir)

            with patch("sniffcell.report.report.viz_module.viz_main") as mock_viz:
                report_main(args)
                mock_viz.assert_not_called()

            report_dir = anno_dir / "report"
            html_path = report_dir / "index.html"
            selected = pd.read_csv(report_dir / "high_confidence_sv.tsv", sep="\t")
            self.assertEqual(selected["viz_status"].tolist(), ["not_rendered"])
            self.assertIn("sniffcell viz --anno_output", str(selected.iloc[0]["viz_command"]))
            text = html_path.read_text(encoding="utf-8")
            self.assertIn("Copy viz command", text)
            self.assertIn("copyVizCommand", text)
            self.assertIn("Interactive Summaries", text)
            self.assertIn("chart-genome-location", text)
            self.assertIn("cdn.plot.ly", text)
            self.assertIn("Variant Review Controls", text)
            self.assertIn("id=\"review-filter\"", text)
            self.assertIn("id=\"celltype-filter\"", text)
            self.assertIn("id=\"svtype-filter\"", text)
            self.assertIn("id=\"hard-conflict-filter\"", text)
            self.assertIn("id=\"numeric-filter-grid\"", text)
            self.assertIn("setSvReview(this)", text)
            self.assertIn("function selectedCelltypeFilter()", text)
            self.assertIn("function selectedSvtypeFilter()", text)
            self.assertIn("function selectedHardConflictFilter()", text)
            self.assertIn("function buildNumericFilterControls()", text)
            self.assertIn("function resetDashboardFilters()", text)
            self.assertNotIn("Reset dashboard filters", text)
            self.assertIn("REVIEW_STORAGE_KEY", text)
            self.assertIn("function persistReviewState()", text)
            self.assertIn("localStorage", text)
            self.assertIn("review-persist-status", text)
            self.assertIn("data-review-status", text)

    def test_report_main_does_not_load_review_status_from_disk_tsv(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            (anno_dir / "anno_run_manifest.json").write_text(json.dumps({"inputs": {}, "outputs": {}}), encoding="utf-8")

            sv_df = pd.DataFrame(
                [
                    {
                        "id": "sv_disk",
                        "assigned_code": "10",
                        "linked_celltypes": "A",
                        "primary_celltype": "A",
                        "has_hard_conflict": False,
                        "overlap_pct": 0.8,
                        "majority_pct": 0.9,
                        "n_supporting": 8,
                        "n_overlapped": 6,
                    }
                ]
            )
            sv_df.to_csv(anno_dir / "sv_assignment.tsv", sep="\t", index=False)
            pd.DataFrame(
                [
                    {
                        "id": "sv_disk",
                        "review_status": "real",
                    }
                ]
            ).to_csv(anno_dir / "report_review.tsv", sep="\t", index=False)

            args = self._base_args(anno_dir)
            report_main(args)

            selected = pd.read_csv(anno_dir / "report" / "high_confidence_sv.tsv", sep="\t")
            self.assertEqual(selected["review_status"].tolist(), ["undecided"])
            html_text = (anno_dir / "report" / "index.html").read_text(encoding="utf-8")
            self.assertIn('data-review-status="undecided"', html_text)
            self.assertIn("review-badge\">undecided</span>", html_text)

    def test_report_main_runs_igvviz_and_records_status(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            (anno_dir / "anno_run_manifest.json").write_text(json.dumps({"inputs": {}, "outputs": {}}), encoding="utf-8")

            sv_df = pd.DataFrame(
                [
                    {
                        "id": "sv_igv",
                        "assigned_code": "10",
                        "linked_celltypes": "A",
                        "primary_celltype": "A",
                        "has_hard_conflict": False,
                        "overlap_pct": 0.8,
                        "majority_pct": 0.9,
                        "n_supporting": 8,
                        "n_overlapped": 6,
                    }
                ]
            )
            sv_df.to_csv(anno_dir / "sv_assignment.tsv", sep="\t", index=False)

            args = self._base_args(anno_dir)
            args.with_igvviz = True
            args.igv_bams = ["a.bam", "b.bam"]
            args.igv_cmd = "igv.sh"
            args.igv_snapshot_format = "png"
            args.reuse_existing_igvviz = False

            def _fake_igvviz(igv_args):
                out_dir = Path(igv_args.output)
                out_dir.mkdir(parents=True, exist_ok=True)
                snap = out_dir / "sv_igv.01.a.igv.png"
                snap.write_text("fake image", encoding="utf-8")
                manifest = out_dir / "sv_igv.igvviz.manifest.json"
                manifest.write_text(
                    json.dumps({"jobs": [{"snapshot": str(snap)}]}),
                    encoding="utf-8",
                )

            with patch("sniffcell.report.report.igvviz_module.igvviz_main", side_effect=_fake_igvviz) as mock_igv:
                report_main(args)
                self.assertEqual(mock_igv.call_count, 1)
                self.assertEqual(mock_igv.call_args[0][0].sv_id, "sv_igv")
                self.assertEqual(mock_igv.call_args[0][0].visibility_window, 5000)
                self.assertEqual(mock_igv.call_args[0][0].phase_tag, "HP")
                self.assertEqual(mock_igv.call_args[0][0].support_tag, "SC")
                self.assertTrue(mock_igv.call_args[0][0].keep_intermediates)
                self.assertFalse(mock_igv.call_args[0][0].batch_only)

            selected = pd.read_csv(anno_dir / "report" / "high_confidence_sv.tsv", sep="\t")
            self.assertEqual(selected["igvviz_status"].tolist(), ["rendered"])
            self.assertIn("sniffcell igvviz --anno_output", str(selected.iloc[0]["igvviz_command"]))

    def test_report_main_embeds_igvviz_snapshots_in_html(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            (anno_dir / "anno_run_manifest.json").write_text(json.dumps({"inputs": {}, "outputs": {}}), encoding="utf-8")

            sv_df = pd.DataFrame(
                [
                    {
                        "id": "sv_img",
                        "assigned_code": "10",
                        "linked_celltypes": "A",
                        "primary_celltype": "A",
                        "has_hard_conflict": False,
                        "overlap_pct": 0.8,
                        "majority_pct": 0.9,
                        "n_supporting": 8,
                        "n_overlapped": 6,
                    }
                ]
            )
            sv_df.to_csv(anno_dir / "sv_assignment.tsv", sep="\t", index=False)

            args = self._base_args(anno_dir)
            args.with_igvviz = True
            args.igv_bams = ["a.bam"]
            args.igv_cmd = "igv.sh"
            args.igv_snapshot_format = "png"
            args.igv_snapshot_width = 3600
            args.igv_snapshot_height = 1600
            args.reuse_existing_igvviz = False

            def _fake_igvviz(igv_args):
                out_dir = Path(igv_args.output)
                out_dir.mkdir(parents=True, exist_ok=True)
                snap = out_dir / "sv_img.01.a.igv.png"
                snap.write_text("fake image", encoding="utf-8")
                manifest = out_dir / "sv_img.igvviz.manifest.json"
                manifest.write_text(
                    json.dumps({"jobs": [{"bam": "/tmp/a.bam", "snapshot": str(snap)}]}),
                    encoding="utf-8",
                )

            with patch("sniffcell.report.report.igvviz_module.igvviz_main", side_effect=_fake_igvviz):
                report_main(args)

            report_dir = anno_dir / "report"
            html_text = (report_dir / "index.html").read_text(encoding="utf-8")
            self.assertIn("sv_img.01.a.igv.png", html_text)
            self.assertIn("<img src=", html_text)
            self.assertIn(".igv-grid{display:flex;flex-direction:column;", html_text)
            self.assertIn("Review controls (IGV):", html_text)
            self.assertIn(".igv-review-layout{display:flex;flex-direction:column;", html_text)

            selected = pd.read_csv(report_dir / "high_confidence_sv.tsv", sep="\t")
            self.assertEqual(selected["igvviz_status"].tolist(), ["rendered"])
            self.assertIn("sv_img.01.a.igv.png", str(selected.iloc[0]["igvviz_snapshots_rel"]))

    def test_report_main_ignores_stale_igvviz_when_disabled(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            (anno_dir / "anno_run_manifest.json").write_text(json.dumps({"inputs": {}, "outputs": {}}), encoding="utf-8")

            sv_df = pd.DataFrame(
                [
                    {
                        "id": "sv_stale_igv",
                        "assigned_code": "10",
                        "linked_celltypes": "A",
                        "primary_celltype": "A",
                        "has_hard_conflict": False,
                        "overlap_pct": 0.8,
                        "majority_pct": 0.9,
                        "n_supporting": 8,
                        "n_overlapped": 6,
                    }
                ]
            )
            sv_df.to_csv(anno_dir / "sv_assignment.tsv", sep="\t", index=False)

            report_dir = anno_dir / "report"
            igvviz_dir = report_dir / "igvviz" / "sv_stale_igv"
            igvviz_dir.mkdir(parents=True, exist_ok=True)
            snap = igvviz_dir / "sv_stale_igv.01.a.igv.png"
            snap.write_text("stale image", encoding="utf-8")
            manifest = igvviz_dir / "sv_stale_igv.igvviz.manifest.json"
            manifest.write_text(
                json.dumps({"jobs": [{"bam": "/tmp/a.bam", "snapshot": str(snap)}]}),
                encoding="utf-8",
            )

            args = self._base_args(anno_dir)
            report_main(args)

            html_text = (report_dir / "index.html").read_text(encoding="utf-8")
            selected = pd.read_csv(report_dir / "high_confidence_sv.tsv", sep="\t")
            manifest_payload = json.loads((report_dir / "report_manifest.json").read_text(encoding="utf-8"))

            self.assertEqual(selected["igvviz_status"].tolist(), ["not_rendered"])
            self.assertNotIn("sv_stale_igv.igvviz.manifest.json", html_text)
            self.assertNotIn("sv_stale_igv.01.a.igv.png", html_text)
            self.assertEqual(int(manifest_payload["counts"]["igvviz_rendered_or_reused"]), 0)

    def test_report_main_runs_igvreport_and_links_alternate_html(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            (anno_dir / "anno_run_manifest.json").write_text(json.dumps({"inputs": {}, "outputs": {}}), encoding="utf-8")

            sv_df = pd.DataFrame(
                [
                    {
                        "id": "sv_alt",
                        "assigned_code": "10",
                        "linked_celltypes": "A",
                        "primary_celltype": "A",
                        "has_hard_conflict": False,
                        "overlap_pct": 0.8,
                        "majority_pct": 0.9,
                        "n_supporting": 8,
                        "n_overlapped": 6,
                    }
                ]
            )
            sv_df.to_csv(anno_dir / "sv_assignment.tsv", sep="\t", index=False)

            args = self._base_args(anno_dir)
            args.with_igvreport = True

            def _fake_igvreport(**kwargs):
                out_dir = Path(kwargs["output_dir"])
                out_dir.mkdir(parents=True, exist_ok=True)
                html_path = out_dir / "index.html"
                html_path.write_text("<html>fake igvreport</html>", encoding="utf-8")
                manifest_path = out_dir / "igvreport_manifest.json"
                manifest_path.write_text(
                    json.dumps(
                        {
                            "status": "rendered",
                            "igvreport_command": "create_report selected_sv_sites.tsv --output index.html",
                        }
                    ),
                    encoding="utf-8",
                )
                return {
                    "status": "rendered",
                    "error": "",
                    "command": "create_report selected_sv_sites.tsv --output index.html",
                    "html_path": str(html_path),
                    "manifest_path": str(manifest_path),
                    "sites_path": str(out_dir / "selected_sv_sites.tsv"),
                    "header_path": str(out_dir / "sniffcell_igvreport_header.html"),
                }

            with patch("sniffcell.report.report.igvreport_module.render_igvreport_bundle", side_effect=_fake_igvreport) as mock_igvreport:
                report_main(args)
                self.assertEqual(mock_igvreport.call_count, 1)
                self.assertEqual(mock_igvreport.call_args.kwargs["window"], 5000)
                self.assertEqual(mock_igvreport.call_args.kwargs["selected_df"]["id"].tolist(), ["sv_alt"])

            report_dir = anno_dir / "report"
            html_text = (report_dir / "index.html").read_text(encoding="utf-8")
            self.assertIn("Alternate IGV Report", html_text)
            self.assertIn("igvreport/index.html", html_text)
            self.assertIn("Copy igvreport command", html_text)

            manifest_payload = json.loads((report_dir / "report_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest_payload["igvreport"]["status"], "rendered")
            self.assertEqual(manifest_payload["counts"]["igvreport_rendered_or_reused"], 1)
            self.assertTrue((report_dir / "igvreport" / "index.html").exists())
            self.assertTrue((report_dir / "igvreport" / "igvreport_manifest.json").exists())

    def test_report_main_backfills_svtype_and_vaf_from_manifest_vcf(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)

            vcf_path = anno_dir / "input.vcf"
            vcf_path.write_text(
                "\n".join(
                    [
                        "##fileformat=VCFv4.2",
                        "##contig=<ID=1>",
                        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">",
                        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">",
                        "##INFO=<ID=VAF,Number=1,Type=Float,Description=\"Variant allele fraction\">",
                        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                        "1\t100\tsv_vcf\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-25;END=125;VAF=0.125",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            (anno_dir / "anno_run_manifest.json").write_text(
                json.dumps({"inputs": {"vcf": str(vcf_path)}, "outputs": {}}),
                encoding="utf-8",
            )

            sv_df = pd.DataFrame(
                [
                    {
                        "id": "sv_vcf",
                        "assigned_code": "10",
                        "linked_celltypes": "A",
                        "primary_celltype": "A",
                        "has_hard_conflict": False,
                        "overlap_pct": 0.8,
                        "majority_pct": 0.9,
                        "n_supporting": 5,
                        "n_overlapped": 4,
                    }
                ]
            )
            sv_df.to_csv(anno_dir / "sv_assignment.tsv", sep="\t", index=False)

            args = self._base_args(anno_dir)
            report_main(args)

            selected = pd.read_csv(anno_dir / "report" / "high_confidence_sv.tsv", sep="\t")
            self.assertEqual(selected["sv_type"].tolist(), ["DEL"])
            self.assertAlmostEqual(float(selected.iloc[0]["vaf"]), 0.125, places=6)
            self.assertEqual(int(selected.iloc[0]["sv_len"]), -25)

            html_text = (anno_dir / "report" / "index.html").read_text(encoding="utf-8")
            self.assertIn('data-sv-type="DEL"', html_text)
            self.assertIn('data-vaf="0.125"', html_text)


if __name__ == "__main__":
    unittest.main()
