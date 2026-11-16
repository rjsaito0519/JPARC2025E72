"""Independent synthetic checks for fixed-pair PHC comparisons and core fits."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import awkward as ak
import uproot

try:
    from common.lib import phc_peak_study as study
except ModuleNotFoundError:
    import phc_peak_study as study


class CoreFitTests(unittest.TestCase):
    def test_gaussian_with_flat_background_recovers_width_and_displacement(self):
        random = np.random.default_rng(72)
        signal = random.normal(0.08, 0.16, 100000)
        values = np.r_[signal, random.uniform(-2, 2, 20000)]
        result = study.gaussian_core(values)
        self.assertTrue(result['ok'])
        self.assertAlmostEqual(result['sigma'], 0.16, delta=0.003)
        self.assertAlmostEqual(result['mean'], 0.08, delta=0.003)

    def test_low_statistics_is_unresolved(self):
        result = study.gaussian_core(np.linspace(-0.2, 0.2, 20))
        self.assertFalse(result['ok'])
        self.assertEqual('low_statistics', result['reason'])

    def test_nonfinite_events_remain_visible_in_diagnostics(self):
        random = np.random.default_rng(3)
        result = study.gaussian_core(np.r_[random.normal(0, 0.2, 10000), np.nan, np.inf])
        self.assertTrue(result['ok'])
        self.assertEqual(2, result['n_nonfinite'])

    def test_many_pairs_from_one_event_are_not_treated_as_sufficient_statistics(self):
        values = np.random.default_rng(3).normal(0, 0.2, 10000)
        result = study.gaussian_core(values, event_ids=np.ones(len(values)))
        self.assertFalse(result['ok'])
        self.assertEqual(1, result['n_events_window'])


class PairIdentityTests(unittest.TestCase):
    def pairs(self):
        return {'event': np.array([4, 4, 9]), 'bht_seg': np.array([1, 2, 1]),
                'bht_entry': np.array([0, 0, 0]), 'bh2_seg': np.array([0, 0, 0]),
                'bh2_entry': np.array([0, 0, 0]), 'offset': np.zeros(3),
                'raw_btof': np.array([1., 2., 3.]), 'bht_cmt': np.zeros(3),
                'bh2_cmt': np.array([0.2, 0.3, 0.4])}

    def test_alignment_is_by_identity_not_tree_row_order(self):
        first = self.pairs()
        permuted = {name: values[[2, 0, 1]] for name, values in first.items()}
        aligned = study.align_pairs([first, permuted])
        np.testing.assert_equal(study.btof_values(aligned[0]), study.btof_values(aligned[1]))

    def test_missing_pair_is_rejected_instead_of_shrinking_comparison_sample(self):
        first = self.pairs()
        smaller = {name: values[:2] for name, values in first.items()}
        with self.assertRaisesRegex(ValueError, 'membership changed'):
            study.align_pairs([first, smaller])

    def test_changed_uncorrected_times_are_rejected(self):
        first = self.pairs()
        changed = self.pairs()
        changed['raw_btof'][0] += 0.1
        with self.assertRaisesRegex(ValueError, 'Uncorrected times changed'):
            study.align_pairs([first, changed])

    def test_duplicate_pair_identity_is_rejected(self):
        first = self.pairs()
        duplicate = {name: values[[0, 0]] for name, values in first.items()}
        with self.assertRaisesRegex(ValueError, 'Duplicate hit-pair'):
            study.pair_keys(duplicate)

    def test_event_split_does_not_put_pairs_of_one_event_on_opposite_sides(self):
        values = np.repeat(np.arange(1000), 3)
        partition = study.event_partition(values).reshape(-1, 3)
        self.assertTrue(np.all(partition[:, :1] == partition))
        self.assertGreater(np.mean(partition), 0.4)
        self.assertLess(np.mean(partition), 0.6)

    def test_unresolved_channel_keeps_reference_values(self):
        first = self.pairs()
        alternative = self.pairs()
        alternative['bh2_cmt'] += 10
        values, choices = study.channel_mixture([first, alternative], np.ones(3, dtype=bool))
        np.testing.assert_equal(study.btof_values(first), values)
        self.assertTrue(all(row['reason'] == 'keep_reference_unresolved' for row in choices.values()))

    def test_time_offset_is_read_from_time_rows_only(self):
        with tempfile.TemporaryDirectory() as directory:
            parameter = Path(directory)/'HDPRM'
            parameter.write_text('3 0 4 0 2 100 4\n3 0 4 1 2 3 2\n1 0 4 1 2 50 2\n')
            self.assertEqual({4: -6.}, study.bh2_offsets(parameter))

    def test_reconstructed_pairs_reproduce_native_histogram_and_raw_particle_selection(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            parameter = folder/'HDPRM'
            parameter.write_text('3 0 4 1 2 0 1\n')
            root = folder/'hodo.root'
            with uproot.recreate(root) as output:
                output['hodo'] = ak.Array({
                    'run_number': [12, 12, 12], 'event_number': [1, 2, 3], 'beam_flag': [1, 1, 2],
                    'bht_hit_seg': [[2], [2, 3], [2]], 'bh2_hit_seg': [[4], [4], [4]],
                    'bht_mt': [[[1.0]], [[1.2], [1.3]], [[1.5]]],
                    'bht_cmt': [[[1.0]], [[1.2], [1.3]], [[1.5]]],
                    'bh2_mt': [[[1.1]], [[1.1]], [[1.1]]],
                    'bh2_cmt': [[[1.1]], [[1.1]], [[1.1]]]})
                values = np.array([1.1-1.0, 1.1-1.2, 1.1-1.3])
                output['BHT_CBTOF_vs_DeltaE_Pi'] = np.histogram2d(
                    np.ones(3), values, bins=[np.linspace(0, 2, 21), np.linspace(-2, 2, 81)])
            pairs, metadata = study.skim_pairs(root, parameter, 'Pi')
            self.assertTrue(metadata['native_match'])
            self.assertEqual(3, len(pairs['event']))
            np.testing.assert_allclose(values, study.btof_values(pairs))
            self.assertEqual({1, 2}, set(pairs['event']))


if __name__ == '__main__':
    unittest.main()
