import csv
import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1]/'tools/export_smaht_atlases.py'
spec = importlib.util.spec_from_file_location('export_smaht_atlases', SCRIPT)
exporter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exporter)


class TestSmahtCatalogExport(unittest.TestCase):
    def test_blood_export_preserves_scoring_text_and_corrects_only_metadata(self):
        with tempfile.TemporaryDirectory() as temp:
            source, dest = Path(temp)/'source.gz', Path(temp)/'dest.gz'
            rows = [dict(chr='chr1', start='10', end='20', code_order='B-cell|Granulocyte',
                         modification='5mC', mean_best_value='0.900000000000001',
                         mean_Granulocyte='', unprofiled_celltypes='Granulocyte',
                         atlas_source='/local/private/reference.mmdb',
                         paired_stability_filter='False', paired_donors='0.0')]
            with gzip.open(source, 'wt', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
                writer.writeheader(); writer.writerows(rows)
            report = exporter.export_blood(source, dest, exporter.sha(source))
            with gzip.open(dest, 'rt') as f:
                actual = list(csv.DictReader(f, delimiter='\t'))[0]
            for key in rows[0].keys() - exporter.METADATA_CHANGES:
                self.assertEqual(actual[key], rows[0][key])
            self.assertEqual(actual['atlas_source'], 'reference.mmdb')
            self.assertEqual(actual['paired_donors'], '2')
            self.assertEqual(actual['paired_stability_filter'], 'True')
            self.assertTrue(report['scoring_fields_preserved'])
            with self.assertRaises(FileExistsError):
                exporter.export_blood(source, dest, exporter.sha(source))

    def test_tissue_export_is_lossless_and_rejects_changed_source(self):
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp)/'source.tsv'
            source.write_bytes(b'chr\tstart\tend\tcode_order\nchr1\t1\t10\tA|B\n')
            dest = Path(temp)/'dest.gz'
            exporter.compress_exact(source, dest, exporter.sha(source))
            self.assertEqual(gzip.decompress(dest.read_bytes()), source.read_bytes())
            with self.assertRaisesRegex(ValueError, 'differs from frozen'):
                exporter.compress_exact(source, Path(temp)/'invalid.gz', '0'*64)


if __name__ == '__main__':
    unittest.main()
