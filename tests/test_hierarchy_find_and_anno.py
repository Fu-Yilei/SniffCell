import unittest

import numpy as np
import pandas as pd

from sniffcell.anno.anno import _assign_target_by_reference_mean, _resolve_cell_types_and_targets
from sniffcell.find.ctdmr import call_ct_combination_dmrs


class TestFindCombinations(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
