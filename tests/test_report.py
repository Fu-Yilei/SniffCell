import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pandas as pd

from sniffcell.report.report import _select_high_confidence_svs, report_main


class TestReportSelection(unittest.TestCase):
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
            min_majority_pct=0.9,
            include_unassigned=False,
            allow_hard_conflict=False,
            max_sv=1,
        )

        self.assertEqual(selected["id"].tolist(), ["svC"])

    def test_select_high_confidence_allows_partial_support_with_strict_majority(self):
        sv_df = pd.DataFrame(
            {
                "id": ["sv_partial", "sv_drop"],
                "assigned_code": ["10", "10"],
                "linked_celltypes": ["A", "A"],
                "has_hard_conflict": [False, False],
                "overlap_pct": [0.8, 0.79],
                "majority_pct": [1.0, 1.0],
                "n_supporting": [5, 5],
                "n_overlapped": [4, 4],
            }
        )
        sv_df["has_hard_conflict"] = sv_df["has_hard_conflict"].astype("boolean")

        selected = _select_high_confidence_svs(
            sv_df,
            min_overlap_pct=0.8,
            min_majority_pct=1.0,
            include_unassigned=False,
            allow_hard_conflict=False,
            max_sv=0,
        )

        self.assertEqual(selected["id"].tolist(), ["sv_partial"])
        self.assertEqual(int(selected.iloc[0]["n_overlapped"]), 4)
        self.assertEqual(int(selected.iloc[0]["n_supporting"]), 5)


class TestReportMain(unittest.TestCase):
    def _base_args(self, anno_output: Path, output: str | None = None, reuse: bool = False):
        return SimpleNamespace(
            anno_output=str(anno_output),
            min_overlap_pct=0.0,
            min_majority_pct=0.0,
            include_unassigned=False,
            allow_hard_conflict=False,
            max_sv=0,
            with_figures=False,
            window=5000,
            max_reads=250,
            format="png",
            export_tables=False,
            reuse_existing_viz=reuse,
            figure_threads=1,
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
            self.assertIn("No SVs passed the report filters", text)
            selected = pd.read_csv(selected_path, sep="\t")
            self.assertEqual(len(selected), 0)

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
            self.assertIn("SV Review Controls", text)
            self.assertIn("id=\"review-filter\"", text)
            self.assertIn("Export review table", text)
            self.assertIn("setSvReview(this)", text)
            self.assertIn("function exportReviewTable()", text)
            self.assertIn("review_status", text)


if __name__ == "__main__":
    unittest.main()
