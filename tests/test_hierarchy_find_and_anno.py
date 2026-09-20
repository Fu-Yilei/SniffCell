import unittest

import numpy as np
import pandas as pd

from sniffcell.anno.anno import (
    _assign_target_by_reference_mean,
    _find_dmr_query_interval,
    _resolve_cell_types_and_targets,
)
from sniffcell.find.ctdmr import call_ct_combination_dmrs, call_ct_specific_dmrs
from sniffcell.parse_args import parse_args


class TestFindCombinations(unittest.TestCase):
    def test_default_single_row_call_keeps_effect_and_cpg_requirements(self):
        idx_df = pd.DataFrame({
            "chr": ["chr1"] * 3,
            "start": [100, 10000, 20000],
            "end": [150, 10050, 20050],
            "startCpG": [1, 4, 7],
            "endCpG": [4, 7, 9],
        })
        means = {"A": pd.Series([0.8, 0.58, 0.9]),
                 "B": pd.Series([0.2, 0.2, 0.1])}
        args = parse_args(["find", "-ck", "example", "-o", "unused.tsv"])
        from_cli = call_ct_combination_dmrs(
            idx_df=idx_df, mean_by_group=means,
            diff_threshold=args.diff_threshold, min_rows=args.min_rows,
            min_cpgs=args.min_cpgs, max_gap_bp=args.max_gap_bp, bed_out=None,
        )
        self.assertEqual(from_cli["start"].tolist(), [100])
        self.assertEqual(from_cli["n_rows"].tolist(), [1])
        for caller in (call_ct_combination_dmrs, call_ct_specific_dmrs):
            with self.subTest(caller=caller.__name__):
                actual = caller(idx_df=idx_df, mean_by_group=means, bed_out=None)
                pd.testing.assert_frame_equal(actual, from_cli)
        strict = parse_args(["find", "-ck", "example", "-o", "unused.tsv", "--min_rows", "2"])
        self.assertTrue(call_ct_combination_dmrs(
            idx_df=idx_df, mean_by_group=means, min_rows=strict.min_rows, bed_out=None,
        ).empty)

    def test_call_ct_combination_dmrs_detects_multi_group_patterns(self):
        idx_df = pd.DataFrame(
            {
                "chr": ["1", "1", "1", "1"],
                "start": [100, 200, 300, 400],
                "end": [150, 250, 350, 450],
                "startCpG": [0, 2, 4, 6],
                "endCpG": [2, 4, 6, 8],
            }
        )

        mean_by_group = {
            "A": pd.Series([0.82, 0.80, 0.20, 0.18]),
            "B": pd.Series([0.78, 0.76, 0.22, 0.20]),
            "C": pd.Series([0.20, 0.22, 0.82, 0.80]),
            "D": pd.Series([0.18, 0.20, 0.78, 0.76]),
        }

        dmrs = call_ct_combination_dmrs(
            idx_df=idx_df,
            mean_by_group=mean_by_group,
            diff_threshold=0.40,
            min_rows=2,
            min_cpgs=2,
            min_bp=0,
            direction="both",
            max_gap_bp=100,
            bed_out=None,
        )

        self.assertEqual(len(dmrs), 2)
        first = dmrs.iloc[0]
        second = dmrs.iloc[1]

        self.assertEqual(first["best_group_leaves"], "A|B")
        self.assertEqual(first["other_group_leaves"], "C|D")
        self.assertEqual(first["best_dir"], "hyper")
        self.assertEqual(first["code_order"], "A|B|C|D")

        self.assertEqual(second["best_group_leaves"], "C|D")
        self.assertEqual(second["other_group_leaves"], "A|B")
        self.assertEqual(second["best_dir"], "hyper")
        self.assertEqual(second["code_order"], "A|B|C|D")


class TestAnnoCodeResolution(unittest.TestCase):
    def test_resolve_targets_from_combo_columns(self):
        row = {
            "code_order": "A|B|C|D",
            "best_group_leaves": "A|B",
        }
        cell_types, target_cell_types = _resolve_cell_types_and_targets(
            row=row,
            best_group="A+B",
        )

        self.assertEqual(cell_types, ["A", "B", "C", "D"])
        self.assertEqual(target_cell_types, ["A", "B"])

    def test_resolve_targets_from_legacy_mean_columns(self):
        row = {
            "mean_A": 0.1,
            "mean_B": np.nan,
            "mean_best_value": 0.2,
        }
        cell_types, target_cell_types = _resolve_cell_types_and_targets(
            row=row,
            best_group="A",
        )

        self.assertEqual(cell_types, ["A"])
        self.assertEqual(target_cell_types, ["A"])

    def test_assign_target_by_reference_mean_prefers_closer_mean(self):
        read_means = np.array([0.10, 0.25, 0.70, 0.90], dtype=float)
        mask = _assign_target_by_reference_mean(
            read_means,
            mean_best_value=0.20,
            mean_rest_value=0.75,
        )
        self.assertIsNotNone(mask)
        self.assertEqual(mask.tolist(), [True, True, False, False])

    def test_assign_target_by_reference_mean_returns_none_on_missing_reference(self):
        read_means = np.array([0.10, 0.25, 0.70], dtype=float)
        mask = _assign_target_by_reference_mean(
            read_means,
            mean_best_value=np.nan,
            mean_rest_value=0.75,
        )
        self.assertIsNone(mask)

    def test_find_dmr_query_interval_shifts_left_by_one(self):
        self.assertEqual(_find_dmr_query_interval(10, 20), (9, 20))

    def test_find_dmr_query_interval_preserves_zero_start(self):
        self.assertEqual(_find_dmr_query_interval(0, 20), (0, 20))


if __name__ == "__main__":
    unittest.main()
