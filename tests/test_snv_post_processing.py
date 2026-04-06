import gzip
import tempfile
import unittest
from pathlib import Path

from sniffcell.discover.snv_post_processing import (
    _resolve_args,
    compare_group_specific_snvs,
    snv_post_processing_main,
)


def _write_gvcf(path: Path, records: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n")
        handle.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">\n")
        handle.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n")
        handle.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n")
        handle.write("##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"AF\">\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for record in records:
            handle.write(record + "\n")


class TestSnvPostProcessing(unittest.TestCase):
    def test_compare_group_specific_snvs(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            group_a = root / "A.pileup.vcf.gz"
            group_b = root / "B.pileup.vcf.gz"
            _write_gvcf(
                group_a,
                [
                    "chr1\t10\t.\tA\tG\t60\tPASS\tP\tGT:GQ:DP:AD:AF\t0/1:30:8:4,4:0.5",
                    "chr1\t20\t.\tC\t.\t20\tRefCall\tP\tGT:GQ:DP:AD:AF\t0/0:20:8:8:1.0",
                    "chr1\t30\t.\tG\tA\t60\tPASS\tP\tGT:GQ:DP:AD:AF\t0/1:30:4:2,2:0.5",
                    "chr1\t40\t.\tT\tC\t60\tPASS\tP\tGT:GQ:DP:AD:AF\t0/1:19:8:4,4:0.5",
                ],
            )
            _write_gvcf(
                group_b,
                [
                    "chr1\t10\t.\tA\t.\t20\tRefCall\tP\tGT:GQ:DP:AD:AF\t0/0:20:9:9:1.0",
                    "chr1\t20\t.\tC\tT\t60\tPASS\tP\tGT:GQ:DP:AD:AF\t0/1:30:9:4,5:0.5",
                    "chr1\t30\t.\tG\tA\t10\tLowQual\tP\tGT:GQ:DP:AD:AF\t0/1:10:7:6,1:0.1429",
                    "chr1\t40\t.\tT\t.\t20\tRefCall\tP\tGT:GQ:DP:AD:AF\t0/0:20:8:8:1.0",
                ],
            )

            summary = compare_group_specific_snvs(
                group_a_label="sample.A",
                group_b_label="sample.B",
                group_a_gvcf=group_a,
                group_b_gvcf=group_b,
                output_dir=root / "out",
                min_dp=5,
                max_dp=50,
                min_dp_absence=5,
                min_gq=21,
                min_other_af=0.9,
            )

            self.assertEqual(summary["group_a_only_count"], 1)
            self.assertEqual(summary["group_b_only_count"], 1)
            self.assertEqual(summary["merged_count"], 2)
            merged = (root / "out" / "snv_changes.tsv").read_text(encoding="utf-8")
            self.assertIn("group_a_only\tchr1\t10\tA\tG\tsample.A\tsample.B", merged)
            self.assertIn("group_b_only\tchr1\t20\tC\tT\tsample.B\tsample.A", merged)
            self.assertNotIn("chr1\t30\tG\tA", merged)
            self.assertNotIn("chr1\t40\tT\tC", merged)
            self.assertIn("\ttarget_alt_ad\t", merged)
            self.assertIn("\tother_alt_ad\t", merged)

    def test_main_runs_end_to_end(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            split_dir = root / "sample1" / "deconv" / "deconv_requested_group_splits"
            split_dir.mkdir(parents=True)
            group_a = root / "A.pileup.vcf.gz"
            group_b = root / "B.pileup.vcf.gz"
            _write_gvcf(
                group_a,
                [
                    "chr1\t10\t.\tA\tG\t60\tPASS\tP\tGT:GQ:DP:AD:AF\t0/1:30:8:4,4:0.5",
                    "chr1\t20\t.\tC\t.\t20\tRefCall\tP\tGT:GQ:DP:AD:AF\t0/0:20:8:8:1.0",
                ],
            )
            _write_gvcf(
                group_b,
                [
                    "chr1\t10\t.\tA\t.\t20\tRefCall\tP\tGT:GQ:DP:AD:AF\t0/0:20:9:9:1.0",
                    "chr1\t20\t.\tC\tT\t60\tPASS\tP\tGT:GQ:DP:AD:AF\t0/1:30:9:4,5:0.5",
                ],
            )

            summary = snv_post_processing_main(
                [
                    "--split-dir", str(split_dir),
                    "--groups", "Neuron,Oligodendrocyte",
                    "--group-a-gvcf", str(group_a),
                    "--group-b-gvcf", str(group_b),
                    "--output-dir", str(root / "run"),
                    "--sample-id", "sample1",
                    "--max-dp", "50",
                    "--min-dp-absence", "5",
                    "--min-other-af", "0.9",
                ]
            )
            self.assertEqual(summary["group_a_only_count"], 1)
            self.assertEqual(summary["group_b_only_count"], 1)
            self.assertTrue(Path(summary["merged_tsv"]).exists())
            self.assertTrue((root / "run" / "summary.json").exists())
            self.assertEqual(summary["params"]["min_dp"], 5)
            self.assertEqual(summary["params"]["max_dp"], 50)
            self.assertEqual(summary["params"]["min_dp_absence"], 5)
            self.assertEqual(summary["params"]["min_gq"], 21)
            self.assertEqual(summary["params"]["min_other_af"], 0.9)

    def test_resolve_args_requires_two_groups(self):
        parser = type(
            "Args",
            (),
            {
                "split_dir": "/tmp/x",
                "groups": "A",
                "group_a_gvcf": "/tmp/a.vcf.gz",
                "group_b_gvcf": "/tmp/b.vcf.gz",
                "output_dir": None,
                "sample_id": None,
                "min_dp": 5,
                "max_dp": 50,
                "min_dp_absence": 15,
                "min_gq": 21,
                "min_other_af": 0.9,
            },
        )
        with self.assertRaises(ValueError):
            _resolve_args(parser)
