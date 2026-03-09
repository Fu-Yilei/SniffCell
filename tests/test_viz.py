import unittest
import tempfile
import json
import logging
from pathlib import Path

import pandas as pd

from sniffcell.viz.viz import (
    _assign_ctdmr_label_lanes,
    _collect_large_indels_from_cigar,
    _build_methylation_heatmap_matrix,
    _expand_interval_for_visibility,
    _extend_region_to_first_informative_ctdmr,
    _reference_celltype_mean_columns,
    _resolve_viz_runtime_inputs,
    _summarize_supporting_read_assignments,
)


class TestVizSupportingReadAssignment(unittest.TestCase):
    def test_supporting_reads_assigned_and_unassigned(self):
        assignment_df = pd.DataFrame(
            {
                "chr": ["1", "1", "1"],
                "start": [100, 120, 5000],
                "end": [140, 160, 5060],
                "code_order": ["A|B", "A|B", "A|B"],
                "code": ["10", "xyz", "01"],
                "best_group": ["A", "A", "B"],
                "best_group_leaves": ["A", "A", "B"],
                "is_best_group": [True, False, True],
            },
            index=pd.Index(["r1", "r2", "rX"], name="readname"),
        )
        assignment_df["read_name"] = assignment_df.index.astype(str)

        summary_df, detail_df = _summarize_supporting_read_assignments(
            assignment_df,
            supporting_reads={"r1", "r2", "r3"},
            sv_chrom="1",
            sv_start=1000,
            sv_end=1010,
            window=1000,
            region_start=0,
            region_end=1000,
            assignment_available=True,
        )

        summary = summary_df.set_index("read_name")
        self.assertTrue(bool(summary.loc["r1", "is_assigned"]))
        self.assertEqual(summary.loc["r1", "assignment_status"], "assigned")
        self.assertEqual(summary.loc["r1", "assigned_celltypes"], "A")
        self.assertFalse(bool(summary.loc["r2", "is_assigned"]))
        self.assertEqual(summary.loc["r2", "assignment_status"], "unassigned_unresolved_code")
        self.assertFalse(bool(summary.loc["r3", "is_assigned"]))
        self.assertEqual(summary.loc["r3", "assignment_status"], "unassigned_no_overlap_rows")

        self.assertEqual(len(detail_df), 1)
        self.assertEqual(detail_df.iloc[0]["read_name"], "r1")
        self.assertEqual(detail_df.iloc[0]["assigned_celltypes"], "A")

    def test_supporting_reads_can_remain_assigned_outside_local_plot_window(self):
        assignment_df = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [1200],
                "end": [1245],
                "code_order": ["Neuron|Oligodendrocyte"],
                "code": ["10"],
                "best_group": ["Neuron"],
                "best_group_leaves": ["Neuron"],
                "is_best_group": [True],
            },
            index=pd.Index(["r1"], name="readname"),
        )
        assignment_df["read_name"] = assignment_df.index.astype(str)

        summary_df, detail_df = _summarize_supporting_read_assignments(
            assignment_df,
            supporting_reads={"r1"},
            sv_chrom="1",
            sv_start=2000,
            sv_end=2010,
            window=1000,
            region_start=1900,
            region_end=2100,
            assignment_available=True,
            clip_to_region=False,
        )

        self.assertEqual(len(summary_df), 1)
        self.assertTrue(bool(summary_df.iloc[0]["is_assigned"]))
        self.assertEqual(summary_df.iloc[0]["assignment_status"], "assigned")
        self.assertEqual(summary_df.iloc[0]["assigned_celltypes"], "Neuron")
        self.assertEqual(len(detail_df), 1)


class TestVizHeatmapMatrix(unittest.TestCase):
    def test_heatmap_matrix_limits_and_orders(self):
        methyl_df = pd.DataFrame(
            {
                "read_name": ["r1", "r2", "r3", "r1"],
                "is_assigned": [True, False, True, True],
                "label": ["DMR_A", "DMR_A", "DMR_A", "DMR_B"],
                "chr": ["1", "1", "1", "1"],
                "start": [100, 100, 100, 300],
                "end": [200, 200, 200, 360],
                "mean_methylation": [0.8, 0.2, 0.7, 0.9],
            }
        )

        heat = _build_methylation_heatmap_matrix(methyl_df, max_reads=2, max_dmrs=1)
        self.assertEqual(heat.shape, (2, 1))
        self.assertEqual(list(heat.index), ["r1", "r3"])

    def test_reference_celltype_mean_columns_uses_all_mean_columns(self):
        dmrs = pd.DataFrame(
            {
                "label": ["A", "B"],
                "start": [100, 200],
                "end": [150, 260],
                "mean_T-cell": [0.8, 0.3],
                "mean_B-cell": [0.2, 0.7],
                "mean_best_value": [0.8, 0.7],  # should be excluded
            }
        )
        cols = _reference_celltype_mean_columns(dmrs)
        self.assertEqual(cols, ["mean_T-cell", "mean_B-cell"])


class TestVizManifestResolution(unittest.TestCase):
    def test_resolve_inputs_from_anno_manifest(self):
        with tempfile.TemporaryDirectory() as td:
            anno_dir = Path(td) / "anno_out"
            anno_dir.mkdir(parents=True, exist_ok=True)
            manifest = {
                "inputs": {
                    "bam": "/tmp/in.bam",
                    "vcf": "/tmp/in.vcf.gz",
                    "reference": "/tmp/ref.fa",
                    "bed": "/tmp/dmrs.tsv",
                },
                "runtime": {"window": 12000},
                "outputs": {"reads_classification": str(anno_dir / "reads_classification.tsv")},
            }
            (anno_dir / "anno_run_manifest.json").write_text(
                json.dumps(manifest),
                encoding="utf-8",
            )

            class Args:
                anno_output = str(anno_dir)
                input = None
                vcf = None
                reference = None
                bed = None
                read_assignment = None
                output = None
                format = "png"
                window = 5000
                sv_id = "sv1"
                kanpig_read_names = None

            resolved = _resolve_viz_runtime_inputs(Args(), logging.getLogger("test"))
            self.assertEqual(resolved["bam_path"], "/tmp/in.bam")
            self.assertEqual(resolved["vcf_path"], "/tmp/in.vcf.gz")
            self.assertEqual(resolved["reference_path"], "/tmp/ref.fa")
            self.assertEqual(resolved["bed_path"], "/tmp/dmrs.tsv")
            self.assertEqual(resolved["window"], 12000)
            self.assertTrue(str(resolved["output_path"]).endswith("sv1.viz.png"))


class TestVizLargeIndels(unittest.TestCase):
    def test_collect_large_indels_from_cigar_filters_and_clips(self):
        events = _collect_large_indels_from_cigar(
            cigartuples=[
                (0, 100),  # M
                (1, 55),   # INS
                (0, 20),
                (2, 80),   # DEL
                (0, 30),
                (1, 15),   # too small
            ],
            reference_start=1000,
            region_start=900,
            region_end=1200,
            min_indel_bp=40,
        )

        self.assertEqual(len(events), 2)
        self.assertEqual(events[0]["event_type"], "INS")
        self.assertEqual(events[0]["pos"], 1100)
        self.assertEqual(events[0]["length"], 55)
        self.assertEqual(events[1]["event_type"], "DEL")
        self.assertEqual(events[1]["start"], 1120)
        self.assertEqual(events[1]["end"], 1200)
        self.assertEqual(events[1]["length"], 80)


class TestVizLinkedCtdmrMode(unittest.TestCase):
    def test_expand_interval_for_visibility_enforces_minimum_span(self):
        new_start, new_end = _expand_interval_for_visibility(
            1000,
            1010,
            region_start=0,
            region_end=5000,
            min_span_bp=120,
        )
        self.assertEqual(new_end - new_start, 120)
        self.assertLessEqual(new_start, 1000)
        self.assertGreaterEqual(new_end, 1010)

    def test_extend_region_uses_nearest_informative_ctdmr(self):
        callouts = pd.DataFrame(
            {
                "chr": ["chr1", "chr1"],
                "start": [1200, 3100],
                "end": [1250, 3150],
                "callout_side": ["left", "right"],
                "callout_support_count": [1, 2],
                "callout_distance_bp": [650, 1100],
            }
        )

        new_start, new_end, row = _extend_region_to_first_informative_ctdmr(
            region_start=1900,
            region_end=2900,
            linked_ctdmr_callouts=callouts,
        )

        self.assertEqual(new_start, 700)
        self.assertEqual(new_end, 2900)
        self.assertIsNotNone(row)
        self.assertEqual(int(row["start"]), 1200)

    def test_assign_ctdmr_label_lanes_separates_overlapping_labels(self):
        entries = [
            {"x_center": 100.0, "text": "0.95"},
            {"x_center": 110.0, "text": "0.82"},
            {"x_center": 300.0, "text": "0.14"},
        ]

        lanes, lane_count = _assign_ctdmr_label_lanes(
            entries,
            region_span_bp=1000,
            font_size=12.0,
        )

        self.assertEqual(len(lanes), 3)
        self.assertGreaterEqual(lane_count, 2)
        self.assertNotEqual(lanes[0], lanes[1])
        self.assertEqual(lanes[2], 0)


if __name__ == "__main__":
    unittest.main()
