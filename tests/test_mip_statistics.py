"""Exercise the ROOT statistics path using generated ADC spectra, without parameter writes."""
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
import uproot

REPO = Path(__file__).resolve().parents[1]
BINARY = REPO / 'bin/BH2_MIP_PREVIEW'


@unittest.skipUnless(BINARY.is_file(), 'build BH2_MIP_PREVIEW first')
class MipStatisticsTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.base = Path(self.directory.name)
        self.edges = np.arange(4097, dtype=float)
        centers = self.edges[:-1]+.5
        pedestal = np.round(1000*np.exp(-.5*((centers-100)/2)**2))
        signal = np.round(25*np.exp(-.5*((centers-250)/15)**2))
        self.counts = pedestal+signal
        with uproot.recreate(self.base/'input.root') as output:
            output['BH2_ADC_seg3U_Pi'] = self.counts, self.edges
            output['BH2_ADC_seg6U_Pi'] = pedestal+30*signal, self.edges

    def preview(self, mode='fixed', lower=230., upper=270., automatic=False):
        command = [str(BINARY), str(self.base/'input.root'), 'Pi', '3', 'U', str(self.base/'result.json')]
        if not automatic:
            settings = self.base/'settings.txt'
            settings.write_text(f'bh2-3-u:Pi mean {mode} {lower} {upper} 1000 250 15 0\n')
            command.append(str(settings))
        result = subprocess.run(command, capture_output=True, text=True, timeout=60)
        self.assertEqual(0, result.returncode, result.stderr[-1000:])
        return json.loads((self.base/'result.json').read_text())

    def expected(self, window, qlo=0., qhi=1.):
        counts = self.counts.reshape(-1, 2).sum(axis=1)
        centers = self.edges[::2][:-1]+1.
        bins = np.arange(len(counts))
        selection = (bins >= int(window[0]//2)) & (bins <= int(window[1]//2))
        weights = counts[selection]
        locations = centers[selection]
        total = weights.sum()
        cumulative = np.cumsum(weights)/total
        before = np.r_[0., cumulative[:-1]]
        selected = np.maximum(0., np.minimum(cumulative, qhi)-np.maximum(before, qlo))*total
        mean = float(np.sum(selected*locations)/selected.sum())
        target = .5*(qlo+qhi)*total
        index = np.searchsorted(np.cumsum(weights), target)
        prior = weights[:index].sum()
        median = float(locations[index]-1.+2.*(target-prior)/weights[index])
        return mean, median, float(total), float(selected.sum())

    def test_fixed_mean_median_counts_and_no_curve(self):
        result = self.preview()
        summary = result['statistics']
        mean, median, total, selected = self.expected([230., 270.])
        self.assertTrue(result['success'])
        self.assertEqual([], result['curve'])
        self.assertAlmostEqual(mean, summary['mean'], places=10)
        self.assertAlmostEqual(median, summary['median'], places=10)
        self.assertEqual(total, summary['entries'])
        self.assertAlmostEqual(selected, summary['selected_entries'])

    def test_automatic_range_keeps_existing_quantile_mean(self):
        result = self.preview('auto')
        summary = result['statistics']
        mean, median, total, selected = self.expected(summary['window'], .2, .5)
        self.assertTrue(result['success'])
        self.assertAlmostEqual(mean, summary['mean'], places=10)
        self.assertAlmostEqual(median, summary['median'], places=10)
        self.assertAlmostEqual(selected, summary['selected_entries'])
        self.assertLess(summary['range'][0], summary['range'][1])

    def test_empty_range_is_not_successful(self):
        result = self.preview(lower=1000., upper=1100.)
        self.assertFalse(result['success'])
        self.assertFalse(result['statistics']['valid'])

    def test_default_low_statistics_exposes_actual_summary(self):
        result = self.preview(automatic=True)
        self.assertEqual('mean', result['function'])
        self.assertEqual(-1, result['status'])
        summary = result['statistics']
        mean, median, _, selected = self.expected(summary['window'], .2, .5)
        self.assertAlmostEqual(mean, summary['mean'], places=10)
        self.assertAlmostEqual(median, summary['median'], places=10)
        self.assertAlmostEqual(selected, summary['selected_entries'])


if __name__ == '__main__':
    unittest.main()
