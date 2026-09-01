import hashlib
import json
import unittest
from pathlib import Path

import numpy as np
import pysam


FIXTURE = Path(__file__).parent / "data/sh3rf3_dual_sv_tr/inputs"


class TestSH3RF3Fixture(unittest.TestCase):
    def test_fixture_integrity_and_coordinates(self):
        manifest = json.loads((FIXTURE / "fixture_manifest.json").read_text())
        self.assertEqual(manifest["coordinate_shift"], 0)
        self.assertEqual(
            manifest["fixture_coordinates_grch38"], "chr2:109199301-109199876"
        )
        for name, expected in manifest["sha256"].items():
            digest = hashlib.sha256((FIXTURE / name).read_bytes()).hexdigest()
            self.assertEqual(digest, expected, name)

        matrix = np.load(FIXTURE / "atlas.npy")
        self.assertEqual(matrix.shape, (164, 9))
        with pysam.AlignmentFile(FIXTURE / "SH3RF3_example.bam", "rb") as bam:
            self.assertEqual(bam.get_reference_length("chr2"), 242_193_529)
            records = list(bam.fetch("chr2", 109_180_000, 109_225_000))
        self.assertEqual(len(records), 976)
        self.assertEqual(len({record.query_name for record in records}), 914)
        self.assertIn(
            "5d5f11d4-9ec1-4cb9-8d35-ecb128b69613",
            {record.query_name for record in records},
        )


if __name__ == "__main__":
    unittest.main()
