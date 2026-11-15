"""Reject finite but out-of-histogram residual fits before geometry updates."""
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'calibration/scripts'))
import update_residual


class ResidualHistogramBoundsTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.fit = self.root/'run00012_BLC2_residual_Pi.root'
        self.raw = self.root/'run00012_BcOut.root'
        self.means = np.zeros((8, 2))
        self.widths = np.full((8, 2), .15)
        with uproot.recreate(self.raw) as output:
            for side in 'ab':
                for plane in range(8):
                    output[f'BLC2{side}_Track_Residual_plane{plane}_Pi'] = (np.ones(400), np.linspace(-2.,2.,401))

    def collect(self):
        with uproot.recreate(self.fit) as output:
            output.mktree('tree', {'residual_p0_val': '2 * float64','residual_p1_val': '2 * float64'})
            output['tree'].extend({'residual_p0_val':self.means,'residual_p1_val':self.widths})
        return update_residual.make_dictdata(str(self.fit), input_path=self.raw)

    def test_finite_means_outside_histogram_are_rejected(self):
        self.means[4,1] = 10.97
        with self.assertRaisesRegex(update_residual.ResidualFitInvalid, '1013'):
            self.collect()

    def test_outside_reference_mean_is_rejected(self):
        self.means[0,0] = -10.
        with self.assertRaises(update_residual.ResidualFitInvalid):
            self.collect()

    def test_invalid_width_is_rejected(self):
        for width in (0., -.15, 4.5):
            with self.subTest(width=width):
                self.widths[3,0] = width
                with self.assertRaises(update_residual.ResidualFitInvalid):
                    self.collect()

    def test_missing_histograms_do_not_disable_validation(self):
        with uproot.recreate(self.raw):
            pass
        with self.assertRaisesRegex(update_residual.ResidualFitInvalid, 'cannot verify'):
            self.collect()

    def test_valid_means_keep_existing_correction_convention(self):
        self.means[0,0] = .01
        self.means[3,1] = .03
        result = self.collect()
        self.assertEqual(16,len(result))
        self.assertAlmostEqual(-.02, result['1012'][0])
        self.assertEqual(.15, result['1012'][1])


if __name__ == '__main__':
    unittest.main()
