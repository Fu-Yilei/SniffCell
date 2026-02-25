import unittest
from pathlib import Path

from sniffcell.sv_discovery import (
    _build_additional_filter_expr,
    _build_sniffles_command,
    _build_svtype_expr,
    _parse_svtypes,
    parse_args,
)


class TestSvDiscoveryHelpers(unittest.TestCase):
    def test_parse_svtypes_normalizes_and_drops_empty(self):
        self.assertEqual(_parse_svtypes("INS,del,, DUP "), ("INS", "DEL", "DUP"))
        self.assertEqual(_parse_svtypes(""), ())

    def test_build_svtype_expr(self):
        self.assertIsNone(_build_svtype_expr(()))
        expr = _build_svtype_expr(("INS", "DEL"))
        self.assertEqual(expr, 'INFO/SVTYPE="INS" || INFO/SVTYPE="DEL"')

    def test_build_additional_filter_expr(self):
        self.assertIsNone(
            _build_additional_filter_expr(
                max_vaf=None,
                max_del_len=None,
                min_support=None,
                max_support=None,
            )
        )

        vaf_only = _build_additional_filter_expr(
            max_vaf=0.2,
            max_del_len=None,
            min_support=None,
            max_support=None,
        )
        self.assertEqual(vaf_only, '((INFO/VAF<=0.2) || (INFO/VAF="." && INFO/AF<=0.2))')

        del_only = _build_additional_filter_expr(
            max_vaf=None,
            max_del_len=500000,
            min_support=None,
            max_support=None,
        )
        self.assertEqual(del_only, '(INFO/SVTYPE!="DEL" || ABS(INFO/SVLEN)<=500000)')

        both = _build_additional_filter_expr(
            max_vaf=0.15,
            max_del_len=1000000,
            min_support=2,
            max_support=20,
        )
        self.assertIn("INFO/VAF<=0.15", both)
        self.assertIn('INFO/SVTYPE!="DEL" || ABS(INFO/SVLEN)<=1000000', both)
        self.assertIn("INFO/SUPPORT>=2", both)
        self.assertIn("INFO/SUPPORT<=20", both)
        self.assertIn("&&", both)

    def test_build_sniffles_command_contains_expected_args(self):
        cmd = _build_sniffles_command(
            sniffles_bin="sniffles",
            input_bam=Path("sample.bam"),
            reference_fa=Path("ref.fa"),
            output_vcf=Path("out.vcf.gz"),
            output_snf=Path("out.snf"),
            threads=12,
            mosaic_af_min=0.01,
            include_germline=True,
            no_qc=True,
            regions_bed=Path("regions.bed"),
            cluster_merge_len=0.02,
            extra_args=("--symbolic",),
        )
        self.assertIn("--mosaic", cmd)
        self.assertIn("--mosaic-include-germline", cmd)
        self.assertIn("--no-qc", cmd)
        self.assertIn("--regions", cmd)
        self.assertIn("regions.bed", cmd)
        self.assertIn("--cluster-merge-len", cmd)
        self.assertIn("0.02", cmd)
        self.assertIn("--symbolic", cmd)
        self.assertIn("--output-rnames", cmd)

    def test_parse_args_defaults_and_toggles(self):
        args = parse_args(["-i", "in.bam", "-r", "ref.fa", "-o", "out"])
        self.assertFalse(args.disable_q100_filter)
        self.assertFalse(args.keep_nonpass)
        self.assertIsNone(args.max_vaf)
        self.assertIsNone(args.max_del_len)
        self.assertIsNone(args.min_support)
        self.assertIsNone(args.max_support)
        self.assertFalse(args.tr_with_qc)
        self.assertFalse(args.tr_exclude_germline)
        self.assertFalse(args.nontr_exclude_germline)
        self.assertEqual(args.nontr_region_mode, "auto")

        args2 = parse_args(
            [
                "-i",
                "in.bam",
                "-r",
                "ref.fa",
                "-o",
                "out",
                "--disable-q100-filter",
                "--keep-nonpass",
                "--max-vaf",
                "0.2",
                "--max-del-len",
                "1000000",
                "--min-support",
                "2",
                "--max-support",
                "20",
                "--tr-with-qc",
                "--nontr-region-mode",
                "postfilter",
                "--tr-exclude-germline",
                "--nontr-exclude-germline",
            ]
        )
        self.assertTrue(args2.disable_q100_filter)
        self.assertTrue(args2.keep_nonpass)
        self.assertEqual(args2.max_vaf, 0.2)
        self.assertEqual(args2.max_del_len, 1000000)
        self.assertEqual(args2.min_support, 2)
        self.assertEqual(args2.max_support, 20)
        self.assertTrue(args2.tr_with_qc)
        self.assertTrue(args2.tr_exclude_germline)
        self.assertTrue(args2.nontr_exclude_germline)
        self.assertEqual(args2.nontr_region_mode, "postfilter")


if __name__ == "__main__":
    unittest.main()
