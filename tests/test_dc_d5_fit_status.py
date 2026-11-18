"""Finite failed native fits must not update DC/D5 calibration."""
import sys
import unittest
from pathlib import Path
from unittest import mock
import numpy as np
sys.path[:0] = [str(Path(__file__).resolve().parents[1] / "calibration/scripts"),
                str(Path(__file__).resolve().parents[1])]
from common.lib import fit_validity
import update_residual
from test_fit_validity_writes import D5PosScoreTests

class NativeFitStatusTests(unittest.TestCase):
    def test_explicit_status_is_integral_and_finite(self):
        for value in (0, 4000, "4000"):
            self.assertTrue(fit_validity.root_fit_status_ok(value))
        for value in (None, 4, .5, float("nan"), float("inf"), "bad"):
            self.assertFalse(fit_validity.root_fit_status_ok(value))

    def residual(self, statuses):
        tree = {"residual_p0_val": np.zeros((8,2)),
                "residual_p1_val": np.full((8,2), .15)}
        if statuses is not None:
            tree["residual_fit_status"] = statuses
        file = mock.MagicMock()
        file.__contains__.return_value = True
        file.__getitem__.return_value.arrays.return_value = tree
        with mock.patch.object(update_residual.uproot, "open", return_value=file):
            return update_residual.make_dictdata("run00001_BLC2_residual_Pi.root")

    def test_residual_success_and_legacy(self):
        for statuses in (None, np.zeros((8,2)), np.full((8,2),4000)):
            self.assertEqual(16,len(self.residual(statuses)))

    def test_residual_finite_failed_fit_rejected(self):
        for value in (4, -1, float("nan"), float("inf")):
            statuses=np.zeros((8,2)); statuses[3,1]=value
            with self.subTest(value=value), self.assertRaises(update_residual.ResidualFitInvalid):
                self.residual(statuses)

    def test_residual_partial_status_rejected(self):
        for statuses in (np.zeros((7,2)), np.zeros((8,1)), np.zeros(8)):
            with self.assertRaises(update_residual.ResidualFitInvalid):
                self.residual(statuses)

class D5NativeStatusTests(D5PosScoreTests):
    BODY = ("dy_a_mm 0\ndy_b_mm 0\ndy_a_err_mm .01\ndy_b_err_mm .02\n"
            "blc1_ab_opposite_shift_s_mm 0\nddxu_mean_percent 0\ndv_mean_mrad 0\n")
    def metadata(self, status=0):
        return "position_fit_ok 1\n" + "".join(
            f"{core}_fit_ok 1\n{core}_fit_status {status}\n" for core in ("dy_a","dy_b","ddxu"))

    def test_finite_failed_fit_rejected(self):
        for status in (4,-1,"nan","inf"):
            self.assertIsNone(self._measure_with_txt(self.BODY+self.metadata(status)))

    def test_partial_or_invalid_aggregate_rejected(self):
        for metadata in ("position_fit_ok 1\n", self.metadata().replace("position_fit_ok 1","position_fit_ok 0")):
            self.assertIsNone(self._measure_with_txt(self.BODY+metadata))

    def test_success_and_info_failure_accepted(self):
        for status in (0,4000):
            self.assertIsNotNone(self._measure_with_txt(
                self.BODY+self.metadata(status)+"dx_a_info_fit_status 4\ndx_a_info_fit_ok 0\n"))

if __name__ == "__main__":
    unittest.main()
