import csv
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from urllib.parse import unquote

import pysam

from sniffcell.discover.export_variants import export_harmonized_vcfs
from sniffcell.discover.discover import _write_harmonized_variant_summary


HEADER = '''##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##INFO=<ID=TRID,Number=1,Type=String,Description="Repeat ID">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcaller_sample
'''


class ExportVariantsTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.tsv = self.root / "harmonized_variants.tsv"

    def row(self, variant_class="SV", **updates):
        row = dict(chrom="chr1", start="99", end="105", variant_class=variant_class,
                   variant_id="call1", variant_subtype="DEL", category="group_a_only",
                   change_size_bp="5", group_a_alt_reads="3", group_b_alt_reads="0")
        row.update(updates)
        return row

    def write_rows(self, rows):
        with self.tsv.open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=self.row().keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    def source(self, lines, name="caller.vcf"):
        path = self.root / name
        path.write_text(HEADER + "".join(line + "\n" for line in lines))
        return path

    def export(self, source, variant_class="SV", **kwargs):
        return export_harmonized_vcfs(self.tsv, sources={"test": (variant_class, source)},
                                      group_a="T-cell,NK-cell; B-cell", group_b="Neuron", **kwargs)["test"]

    def records(self, result):
        with pysam.VariantFile(result["output"]) as reader:
            return list(reader.fetch("chr1"))

    def test_sv_preserves_caller_fields_and_encodes_group(self):
        self.write_rows([self.row()])
        source = self.source(["chr1\t100\tcall1\tACCCCC\tA\t42\tPASS\tSVTYPE=DEL;SVLEN=-5\tGT:AD\t0|1:7,3"])
        result = self.export(source)
        record = self.records(result)[0]
        self.assertEqual((record.pos, record.ref, record.alts, record.qual), (100, "ACCCCC", ("A",), 42))
        self.assertEqual(record.samples["caller_sample"]["GT"], (0, 1))
        self.assertTrue(record.samples["caller_sample"].phased)
        self.assertEqual(record.samples["caller_sample"]["AD"], (7, 3))
        self.assertEqual(record.info["SVLEN"], -5)
        self.assertEqual(record.info["SC_CHANGE_BP"], 5)
        self.assertEqual(unquote(record.info["SC_TARGET_GROUP"]), "T-cell,NK-cell; B-cell")
        self.assertEqual(result["unmatched_rows"], [])

    def test_tr_matches_catalog_not_vcf_id_and_preserves_locus_evidence(self):
        self.write_rows([self.row("TR", start="100", end="106", variant_id="chr1_100_106",
                                  variant_subtype="expansion_hap1", change_size_bp="80")])
        catalog = self.root / "repeats.bed"
        catalog.write_text("chr1\t100\t106\tID=catalog-id;MOTIFS=CAG;STRUC=<TR>\n")
        source = self.source(["chr1\t100\t.\tACAGCAG\tACAGCAGCAG\t.\t.\tTRID=catalog-id;END=106\tGT:AD\t1/1:0,9"])
        record = self.records(self.export(source, "TR", tr_bed=catalog))[0]
        self.assertEqual(record.id, None)
        self.assertEqual(record.alts, ("ACAGCAGCAG",))
        self.assertEqual(record.info["SC_CHANGE_BP"], 80)
        self.assertNotIn("SVLEN", record.info)

    def test_snv_matches_allele_not_id_and_preserves_multiallelic_gt(self):
        self.write_rows([self.row("SNV", end="100", variant_id="chr1:100:A>G", variant_subtype="A>G",
                                  category="group_b_only", change_size_bp="1")])
        source = self.source(["chr1\t100\t.\tA\tC,G\t30\tPASS\t.\tGT:AD\t0/2:8,0,4"])
        record = self.records(self.export(source, "SNV"))[0]
        self.assertEqual(record.alts, ("C", "G"))
        self.assertEqual(record.samples["caller_sample"]["GT"], (0, 2))
        self.assertEqual(record.info["SC_TARGET_GROUP"], "Neuron")

    def test_shared_duplicate_ids_sorted_and_unmatched_reported(self):
        self.write_rows([self.row(start="199", end="205", category="shared"), self.row(),
                         self.row(variant_id="absent")])
        source = self.source([
            "chr1\t200\tcall1\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=205\tGT\t0/1",
            "chr1\t100\tcall1\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=105\tGT\t0/1",
        ])
        result = self.export(source)
        records = self.records(result)
        self.assertEqual([record.pos for record in records], [100, 200])
        self.assertNotIn("SC_TARGET_GROUP", records[1].info)
        self.assertEqual(result["unmatched_rows"], [3])
        self.assertEqual(result["status"], "partial")

    def test_missing_source_removes_stale_output(self):
        self.write_rows([self.row()])
        stale = self.root / "harmonized_variants.test.vcf.gz"
        stale.write_text("stale")
        result = self.export(None)
        self.assertEqual(result["unmatched_rows"], [1])
        self.assertIsNone(result["output"])
        self.assertFalse(stale.exists())

    def test_empty_input_writes_indexed_header(self):
        self.write_rows([])
        result = self.export(self.source([]))
        self.assertEqual(self.records(result), [])
        self.assertTrue(Path(result["output"] + ".tbi").exists())

    def test_ambiguous_source_records_fail(self):
        self.write_rows([self.row()])
        line = "chr1\t100\tcall1\tA\t<DEL>\t.\t.\tSVTYPE=DEL;END=105\tGT\t0/1"
        with self.assertRaisesRegex(ValueError, "Ambiguous caller"):
            self.export(self.source([line, line]))

    def test_breakend_allele_is_not_reconstructed(self):
        self.write_rows([self.row(variant_subtype="BND", category="unknown", change_size_bp=".")])
        source = self.source(["chr1\t100\tcall1\tA\tA]chr2:200]\t.\t.\tSVTYPE=BND\tGT\t./."])
        record = self.records(self.export(source))[0]
        self.assertEqual(record.alts, ("A]chr2:200]",))
        self.assertNotIn("SC_TARGET_GROUP", record.info)
        self.assertNotIn("SC_CHANGE_BP", record.info)

    def test_multiple_tr_changes_keep_separate_row_annotations(self):
        self.write_rows([self.row("TR", start="100", end="106", variant_subtype="expansion_hap1"),
                         self.row("TR", start="100", end="106", variant_subtype="expansion_hap2")])
        catalog = self.root / "repeats.bed"
        catalog.write_text("chr1\t100\t106\tID=repeat\n")
        source = self.source(["chr1\t100\t.\tACAGCAG\tACAG\t.\t.\tTRID=repeat;END=106\tGT\t0/1"])
        records = self.records(self.export(source, "TR", tr_bed=catalog))
        self.assertEqual([record.info["SC_ROW"] for record in records], [1, 2])
        self.assertEqual([record.info["SC_SUBTYPE"] for record in records], ["expansion_hap1", "expansion_hap2"])

    def test_discover_writes_exports_and_manifest_for_snv_only(self):
        snv = self.root / "snv_changes.tsv"
        snv.write_text("chrom\tpos\tref\talt\tdirection\ttarget_alt_ad\tsnv_pass_for_harmonized\n"
                       "chr1\t100\tA\tG\tgroup_a_only\t3\ttrue\n")
        post = self.root / "snv_post"
        post.mkdir()
        (post / "summary.json").write_text(json.dumps({"merged_tsv": str(snv)}))
        source = self.source(["chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1"])
        ctx = SimpleNamespace(selected_groups=["Neuron", "Oligo"], run_root=self.root,
                              stages=["clair3"], sample_id="donor", dry_run=False,
                              params={"medaka_sample_name_template": "{sample_id}.{group}",
                                      "trgt_sample_name_template": "{sample_id}.{group}"})
        with patch("sniffcell.discover.discover._sv_post_stage_dir", return_value=self.root / "absent"), \
             patch("sniffcell.discover.discover._tr_post_stage_dir", return_value=self.root / "absent"), \
             patch("sniffcell.discover.discover._snv_post_stage_dir", return_value=post), \
             patch("sniffcell.discover.discover._clair3_pileup_output_path", return_value=source):
            self.assertEqual(_write_harmonized_variant_summary(ctx), self.tsv)
        manifest = json.loads((self.root / "harmonized_variants_manifest.json").read_text())
        self.assertEqual(manifest["vcf_exports"]["snv.Neuron"]["exported_records"], 1)
        self.assertEqual(manifest["counts"]["snv_rows"], 1)
        ctx.dry_run = True
        with patch("sniffcell.discover.discover._sv_post_stage_dir", return_value=self.root / "absent"), \
             patch("sniffcell.discover.discover._tr_post_stage_dir", return_value=self.root / "absent"), \
             patch("sniffcell.discover.discover._snv_post_stage_dir", return_value=post), \
             patch("sniffcell.discover.discover.export_harmonized_vcfs") as exporter:
            _write_harmonized_variant_summary(ctx)
            exporter.assert_not_called()


if __name__ == "__main__":
    unittest.main()
