"""Biological invariants for the partially profiled whole-blood catalog."""
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from sniffcell.anno.anno import _one_dmr
from sniffcell.anno.variant_assignment import (
    _compute_complete_reference_consensus_df,
    _compute_per_read_consensus_df,
    assign_sv_celltypes,
)
from sniffcell.deconv.deconv import _build_read_summary, _prepare_read_assignment_df


CLASSES = ["B-cell", "Granulocyte", "Monocyte", "NK-cell", "T-cell"]
SCHEMA = "|".join(CLASSES)


def evidence(votes):
    rows = []
    for i, (read, modification, candidates) in enumerate(votes):
        code = "".join("1" if cell in candidates else "0" for cell in CLASSES)
        rows.append(dict(readname=read, modification=modification, code=code,
                         code_token=SCHEMA + "::" + code, code_order=SCHEMA,
                         unprofiled_celltypes="" if modification == "modifiedC" else "Granulocyte",
                         best_group="B-cell", best_group_leaves="B-cell",
                         other_group="Granulocyte+Monocyte+NK-cell+T-cell",
                         other_group_leaves="Granulocyte|Monocyte|NK-cell|T-cell",
                         is_best_group=True, chr="chr1", start=10+i*20, end=20+i*20))
    return pd.DataFrame(rows).set_index("readname")


class TestPartialAtlasConsensus(unittest.TestCase):
    def test_missing_granulocyte_assays_preserve_candidate_sets(self):
        cases = [
            ([('5mC', ['B-cell'])], {'B-cell', 'Granulocyte'}),
            ([('5mC', ['B-cell']), ('5hmC', ['T-cell'])], set()),
            ([('modifiedC', ['Granulocyte']), ('5mC', ['B-cell'])], {'Granulocyte'}),
            ([('modifiedC', ['B-cell', 'Granulocyte']), ('5mC', ['T-cell'])], {'B-cell', 'Granulocyte'}),
            ([('modifiedC', ['B-cell', 'T-cell']), ('5mC', ['B-cell'])], {'B-cell'}),
            ([('modifiedC', ['B-cell']), ('modifiedC', ['T-cell']), ('5mC', ['B-cell'])], set()),
            ([('5hmC', ['NK-cell', 'T-cell'])], {'Granulocyte', 'NK-cell', 'T-cell'}),
            ([('modifiedC', ['B-cell']), ('5mC', ['T-cell'])], set()),
            ([('modifiedC', ['Monocyte'])], {'Monocyte'}),
        ]
        # All cases together also exercise reads with only complete profiles
        # inside a batch that has partial profiles for other reads.
        votes = [(str(i), mod, cells) for i, (items, _) in enumerate(cases) for mod, cells in items]
        assignments = evidence(votes)
        for frame in (assignments, _prepare_read_assignment_df(assignments)):
            direct = _compute_per_read_consensus_df(frame).set_index('readname')
            summary = _build_read_summary(frame).set_index('readname')
            for i, (_, expected) in enumerate(cases):
                with self.subTest(case=i):
                    row = summary.loc[str(i)]
                    labels = set(filter(None, row.linked_leaf_celltypes.split('|')))
                    self.assertEqual(labels, expected)
                    self.assertEqual(bool(row.is_mixed), not bool(expected))
                    self.assertEqual(bool(direct.loc[str(i), 'read_is_mixed']), not bool(expected))
                    self.assertTrue(row.assay_consensus_rule)

    def test_variant_annotation_keeps_unknown_granulocyte_possible(self):
        assignments = evidence([('read1', '5mC', ['B-cell']), ('read2', '5hmC', ['B-cell'])])
        variants = pd.DataFrame([dict(chr='chr1', location=20, id='sv1', supporting_reads=['read1', 'read2'])])
        row = assign_sv_celltypes(variants, assignments, use_read_names=True).iloc[0]
        self.assertEqual(set(row.linked_celltypes.split('|')), {'B-cell', 'Granulocyte'})
        self.assertEqual(row.assigned_code, SCHEMA + '::11000')
        self.assertEqual(row.n_overlapped, 2)

    def test_conflicting_mh_reads_have_no_readable_variant_vote(self):
        assignments = evidence([('read1', '5mC', ['B-cell']), ('read1', '5hmC', ['T-cell'])])
        variants = pd.DataFrame([dict(chr='chr1', location=20, id='sv1', supporting_reads=['read1'])])
        row = assign_sv_celltypes(variants, assignments, use_read_names=True).iloc[0]
        self.assertEqual(row.n_mixed_reads, 1)
        self.assertTrue(pd.isna(row.assigned_code))
        self.assertEqual(row.linked_celltypes, '')

    def test_complete_catalog_consensus_is_unchanged(self):
        frame = evidence([('a', 'modifiedC', ['B-cell']), ('b', 'modifiedC', ['Monocyte'])])
        expected = _compute_complete_reference_consensus_df(frame)
        for data in (frame, frame.drop(columns='unprofiled_celltypes')):
            actual = _compute_per_read_consensus_df(data)
            pd.testing.assert_frame_equal(actual[expected.columns], expected)

    def test_one_dmr_never_fills_unknown_reference_with_other_mean(self):
        row = dict(chr='chr1', start=10, end=20, best_group='B-cell',
                   other_group='Monocyte+NK-cell+T-cell', best_dir='hyper',
                   mean_best_value=.9, mean_rest_value=.1,
                   code_order=SCHEMA, best_group_leaves='B-cell',
                   other_group_leaves='Monocyte|NK-cell|T-cell',
                   modification='5mC', unprofiled_celltypes='Granulocyte')
        matrix = pd.DataFrame([[.9, np.nan], [.1, np.nan]],
                              index=pd.Index(['b', 'other'], name='read_name'), columns=[10, 15])
        with patch('sniffcell.anno.anno._get_process_bam'), \
             patch('sniffcell.anno.anno._get_process_fasta'), \
             patch('sniffcell.anno.anno.methyl_matrix_from_bam', return_value=(matrix, [10, 15])):
            assignments, state = _one_dmr((row, 'input.bam', 'reference.fa', 'closest_reference_mean', 'auto'))
        self.assertEqual(assignments.code.tolist(), ['10000', '00111'])
        self.assertTrue(state['Granulocyte_methylation'].isna().all())
        self.assertEqual(set(assignments.unprofiled_celltypes), {'Granulocyte'})


if __name__ == '__main__':
    unittest.main()
