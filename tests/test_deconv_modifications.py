from __future__ import annotations

import array
import tempfile
import unittest
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from sniffcell.anno.anno import _one_dmr
from sniffcell.anno.methyl_matrix import (
    methyl_matrix_from_bam,
    normalize_modification_label,
    wanted_keys_for_modification,
)
from sniffcell.deconv.deconv import _load_ctdmr_bed
from sniffcell.deconv.regions import TargetRegion, select_ctdmrs_for_target
from sniffcell.parse_args import parse_args


class ModificationBamFixture:
    def __init__(self, root: Path, *, partial_cpg_calls: bool = False):
        self.fasta = root / "reference.fa"
        self.fasta.write_text(">chr1\nACGCGTAAAA\n", encoding="utf-8")
        pysam.faidx(str(self.fasta))

        self.bam = root / "dual_modification.bam"
        header = {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 10}],
        }
        with pysam.AlignmentFile(str(self.bam), "wb", header=header) as output:
            if partial_cpg_calls:
                output.write(
                    self._make_read("m_high", [230, 25], mm_tag="C+m,0;C+h,0;")
                )
                output.write(
                    self._make_read("h_high", [25, 230], mm_tag="C+m,0;C+h,0;")
                )
            else:
                output.write(self._make_read("m_high", [230, 230, 25, 25]))
                output.write(self._make_read("h_high", [25, 25, 230, 230]))
        pysam.index(str(self.bam))

    @staticmethod
    def _make_read(
        name: str,
        probabilities: list[int],
        *,
        mm_tag: str = "C+m,0,0;C+h,0,0;",
    ) -> pysam.AlignedSegment:
        read = pysam.AlignedSegment()
        read.query_name = name
        read.query_sequence = "ACGCGT"
        read.flag = 0
        read.reference_id = 0
        read.reference_start = 0
        read.mapping_quality = 60
        read.cigar = ((0, 6),)
        read.query_qualities = pysam.qualitystring_to_array("IIIIII")
        read.set_tag("MM", mm_tag)
        read.set_tag("ML", array.array("B", probabilities))
        return read


def _dmr_row(modification: str) -> dict[str, object]:
    return {
        "chr": "chr1",
        "start": 0,
        "end": 6,
        "best_group": "Neuron",
        "other_group": "Oligodendrocyte",
        "best_dir": "hyper",
        "mean_best_value": 0.90,
        "mean_rest_value": 0.10,
        "code_order": "Neuron|Oligodendrocyte",
        "best_group_leaves": "Neuron",
        "other_group_leaves": "Oligodendrocyte",
        "hyper_group_leaves": "Neuron",
        "hypo_group_leaves": "Oligodendrocyte",
        "modification": modification,
    }


class TestModificationLabels(unittest.TestCase):
    def test_normalizes_supported_aliases(self):
        self.assertEqual(normalize_modification_label("m"), "5mC")
        self.assertEqual(normalize_modification_label("5HMC"), "5hmC")
        self.assertEqual(normalize_modification_label("combined"), "modifiedC")
        self.assertEqual(normalize_modification_label("auto", allow_auto=True), "auto")
        self.assertEqual(
            wanted_keys_for_modification("5hmC"),
            {("C", 0, "h"), ("C", 1, "h")},
        )

    def test_rejects_unknown_catalog_modification(self):
        with self.assertRaisesRegex(ValueError, "Unsupported modification"):
            normalize_modification_label("6mA")


class TestModificationAwareMethylMatrix(unittest.TestCase):
    def test_extracts_5mc_and_5hmc_independently(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            fixture = ModificationBamFixture(Path(temp_dir))
            five_mc = methyl_matrix_from_bam(
                str(fixture.bam), str(fixture.fasta), "chr1", 0, 6, modification="5mC"
            )
            five_hmc = methyl_matrix_from_bam(
                str(fixture.bam), str(fixture.fasta), "chr1", 0, 6, modification="5hmC"
            )
            modified_c = methyl_matrix_from_bam(
                str(fixture.bam), str(fixture.fasta), "chr1", 0, 6, modification="modifiedC"
            )

        self.assertAlmostEqual(float(five_mc.loc[("m_high", -1)].mean()), 230 / 255)
        self.assertAlmostEqual(float(five_mc.loc[("h_high", -1)].mean()), 25 / 255)
        self.assertAlmostEqual(float(five_hmc.loc[("m_high", -1)].mean()), 25 / 255)
        self.assertAlmostEqual(float(five_hmc.loc[("h_high", -1)].mean()), 230 / 255)
        np.testing.assert_allclose(
            modified_c.loc[("m_high", -1)].to_numpy(),
            modified_c.loc[("h_high", -1)].to_numpy(),
        )

    def test_one_dmr_uses_catalog_channel_in_auto_mode_and_records_it(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            fixture = ModificationBamFixture(Path(temp_dir))
            result = _one_dmr(
                (
                    _dmr_row("5hmC"),
                    str(fixture.bam),
                    str(fixture.fasta),
                    "closest_reference_mean",
                    "auto",
                )
            )

        self.assertIsNotNone(result)
        assignments, blocks = result
        self.assertFalse(bool(assignments.loc["m_high", "is_best_group"]))
        self.assertTrue(bool(assignments.loc["h_high", "is_best_group"]))
        self.assertEqual(set(assignments["modification"]), {"5hmC"})
        self.assertEqual(set(assignments["bam_modification"]), {"5hmC"})
        self.assertEqual(blocks.loc[0, "bam_modification"], "5hmC")

    def test_one_dmr_explicit_override_replaces_catalog_channel(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            fixture = ModificationBamFixture(Path(temp_dir))
            result = _one_dmr(
                (
                    _dmr_row("5mC"),
                    str(fixture.bam),
                    str(fixture.fasta),
                    "closest_reference_mean",
                    "5hmC",
                )
            )

        self.assertIsNotNone(result)
        assignments, _ = result
        self.assertTrue(bool(assignments.loc["h_high", "is_best_group"]))
        self.assertEqual(set(assignments["modification"]), {"5mC"})
        self.assertEqual(set(assignments["bam_modification"]), {"5hmC"})

    def test_one_dmr_ignores_cpg_columns_with_no_calls(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            fixture = ModificationBamFixture(Path(temp_dir), partial_cpg_calls=True)
            with warnings.catch_warnings():
                warnings.simplefilter("error", RuntimeWarning)
                result = _one_dmr(
                    (
                        _dmr_row("5mC"),
                        str(fixture.bam),
                        str(fixture.fasta),
                        "closest_reference_mean",
                        "auto",
                    )
                )

        self.assertIsNotNone(result)
        assignments, _ = result
        self.assertTrue(bool(assignments.loc["m_high", "is_best_group"]))
        self.assertFalse(bool(assignments.loc["h_high", "is_best_group"]))


class TestDeconvModificationInterface(unittest.TestCase):
    def test_cli_defaults_to_auto_and_accepts_explicit_5hmc(self):
        base = [
            "deconv",
            "-i",
            "input.bam",
            "-r",
            "ref.fa",
            "-b",
            "ctdmr.tsv",
            "-o",
            "out",
        ]
        self.assertEqual(parse_args(base).bam_modification, "auto")
        self.assertEqual(parse_args(base + ["--bam-modification", "5hmC"]).bam_modification, "5hmC")

    def test_legacy_catalog_defaults_to_modifiedc(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "legacy.tsv"
            pd.DataFrame([_dmr_row("modifiedC")]).drop(columns="modification").to_csv(
                path, sep="\t", index=False
            )
            loaded = _load_ctdmr_bed(str(path))

        self.assertEqual(loaded["modification"].tolist(), ["modifiedC"])

    def test_regional_selection_keeps_all_channels_at_distance_boundary(self):
        catalog = pd.DataFrame(
            [
                {"chr": "chr1", "start": 10, "end": 20, "modification": "5mC"},
                {"chr": "chr1", "start": 10, "end": 20, "modification": "5hmC"},
                {"chr": "chr1", "start": 30, "end": 40, "modification": "5mC"},
                {"chr": "chr1", "start": 30, "end": 40, "modification": "5hmC"},
                {"chr": "chr1", "start": 70, "end": 80, "modification": "5mC"},
                {"chr": "chr1", "start": 70, "end": 80, "modification": "5hmC"},
            ]
        )
        selected = select_ctdmrs_for_target(
            catalog,
            target=TargetRegion("chr1", 50, 60),
            left_ctdmrs=1,
            right_ctdmrs=1,
        )

        self.assertEqual(len(selected), 4)
        self.assertEqual(selected.groupby(["start", "end"]).size().tolist(), [2, 2])


if __name__ == "__main__":
    unittest.main()
