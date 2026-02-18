import unittest
import tempfile
import json
import logging
from pathlib import Path

import pandas as pd

from sniffcell.viz.viz import (
    _build_methylation_heatmap_matrix,
    _build_reference_celltype_matrix,
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

    def test_reference_celltype_matrix_uses_all_mean_columns(self):
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
        mat = _build_reference_celltype_matrix(dmrs)
        self.assertEqual(list(mat.columns), ["T-cell", "B-cell"])
        self.assertEqual(mat.shape, (2, 2))
        self.assertAlmostEqual(float(mat.iloc[0, 0]), 0.8)


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


if __name__ == "__main__":
    unittest.main()
