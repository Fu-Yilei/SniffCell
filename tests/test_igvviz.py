import json
import logging
import tempfile
import unittest
from pathlib import Path

from sniffcell.viz.igvviz import (
    _build_igv_batch_lines,
    _infer_gene_track,
    _resolve_igvviz_runtime_inputs,
    _split_bam_args,
)


class TestIgvVizHelpers(unittest.TestCase):
    def test_split_bam_args_supports_lists_and_commas(self):
        out = _split_bam_args(["a.bam,b.bam", "c.bam", "a.bam"])
        self.assertEqual(out, ["a.bam", "b.bam", "c.bam"])

    def test_resolve_inputs_from_anno_manifest_when_bam_not_provided(self):
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
                "runtime": {"window": 12345},
            }
            (anno_dir / "anno_run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

            class Args:
                anno_output = str(anno_dir)
                input = None
                vcf = None
                reference = None
                bed = None
                window = 5000
                output = None
                kanpig_read_names = None

            resolved = _resolve_igvviz_runtime_inputs(Args(), logging.getLogger("test"))
            self.assertEqual(resolved["bam_paths"], ["/tmp/in.bam"])
            self.assertEqual(resolved["vcf_path"], "/tmp/in.vcf.gz")
            self.assertEqual(resolved["reference_path"], "/tmp/ref.fa")
            self.assertEqual(resolved["bed_path"], "/tmp/dmrs.tsv")
            self.assertEqual(resolved["window"], 12345)
            self.assertTrue(str(resolved["output_dir"]).endswith("anno_out/igvviz"))

    def test_build_batch_lines_include_phase_sort_and_snapshot(self):
        lines = _build_igv_batch_lines(
            jobs=[
                {
                    "tagged_bam": "/tmp/a.tagged.bam",
                    "ctdmr_track": "/tmp/a.ctdmr.bed",
                    "locus": "chr1:100-200",
                    "snapshot_name": "sv1.a.igv.png",
                }
            ],
            snapshot_dir=Path("/tmp/out"),
            reference_path="/tmp/ref.fa",
            gene_track_path="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqSelect.txt.gz",
            visibility_window=250000,
            phase_tag="HP",
            support_phase_group_tag="SG",
            snapshot_width=2800,
            snapshot_height=1500,
        )
        text = "\n".join(lines)
        self.assertIn("group TAG SG", text)
        self.assertIn("sort TAG HP", text)
        self.assertIn("colorBy BASE_MODIFICATION_2COLOR", text)
        self.assertIn("preference BASEMOD.THRESHOLD 0.7", text)
        self.assertIn("preference BASEMOD.M_COLOR 220,38,38", text)
        self.assertIn("preference BASEMOD.NONE_C_COLOR 65,105,225", text)
        self.assertIn("load https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqSelect.txt.gz", text)
        self.assertIn("snapshot sv1.a.igv.png", text)
        self.assertIn("preference SAM.MAX_VISIBLE_RANGE 250000", text)
        self.assertNotIn("setWindowBounds 50 50 2800 1500", text)

    def test_infer_gene_track_for_hg38_reference(self):
        self.assertEqual(
            _infer_gene_track("/tmp/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"),
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqSelect.txt.gz",
        )
        self.assertIsNone(_infer_gene_track("/tmp/mm39.fa"))


if __name__ == "__main__":
    unittest.main()
