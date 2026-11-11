"""Reject non-finite / failed fits before param write or converge (#51/#52/#53/#54)."""
from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

REPO = Path(__file__).resolve().parents[1]


class FitValidityHelperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import fit_validity

        cls.fv = fit_validity

    def test_all_finite(self):
        self.assertTrue(self.fv.all_finite([1.0, 0.0, -2]))
        self.assertFalse(self.fv.all_finite([1.0, float("nan")]))
        self.assertFalse(self.fv.phc_fit_status_ok(4))
        self.assertTrue(self.fv.phc_fit_status_ok(0))
        self.assertTrue(self.fv.phc_fit_status_ok(None))


class PhcUpdateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        scripts = str(REPO / "calibration" / "scripts")
        if scripts not in sys.path:
            sys.path.insert(0, scripts)
        import update_phc

        cls.phc = update_phc

    def _tree(self):
        # seg0: enough stats, bad status + NaN p0; seg1: no stats
        return {
            "ch": np.array([0, 1]),
            "n_fit": np.array([[1000.0, 1000.0], [0.0, 0.0]]),
            "chi2_ndf": np.array([[1.0, 1.0], [-1.0, -1.0]]),
            "fit_status": np.array([[4, 4], [0, 0]]),
            "p0_val": np.array([[float("nan"), 1.0], [0.0, 0.0]]),
            "p1_val": np.array([[1.0, 1.0], [0.0, 0.0]]),
            "p2_val": np.array([[1.0, 1.0], [0.0, 0.0]]),
        }

    def test_auto_bad_rejects_status_and_nan(self):
        bad = self.phc.auto_bad_segments(self._tree())
        self.assertEqual(bad, [0, 1])

    def test_make_dictdata_no_nan_write_or_mean(self):
        tree = self._tree()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "fake_PHC.root"
            # Bypass ROOT I/O: patch open to return our arrays.
            fake = mock.MagicMock()
            fake.__getitem__.return_value.arrays.return_value = tree
            with mock.patch.object(self.phc.uproot, "open", return_value=fake), \
                 mock.patch.object(self.phc.config, "DETECTOR_I_LIST", {"BHT": 1}):
                data, unused = self.phc.make_dictdata(
                    str(path), "BHT", good_ch_range=None,
                )
                data_h, unused_h = self.phc.make_dictdata(
                    str(path), "BHT", good_ch_range=[0, 1],
                )
        self.assertEqual(data, {})
        self.assertEqual(sorted(unused), [0, 1])
        self.assertEqual(data_h, {})
        self.assertEqual(sorted(unused_h), [0, 1])
        for rows in (data, data_h):
            for vals in rows.values():
                self.assertTrue(all(math.isfinite(float(x)) for x in vals[2:]))


class HdprmUpdateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        scripts = str(REPO / "calibration" / "scripts")
        if scripts not in sys.path:
            sys.path.insert(0, scripts)
        import update_hdprm

        cls.hd = update_hdprm

    def test_ud_counter_skips_nan_tdc(self):
        tree = {
            "ch": np.array([0.0]),
            "tdc_p0_val": np.array([[float("nan")]]),
            "adc_p0_val": np.array([[10.0]]),
            "adc_p1_val": np.array([[2.0]]),
        }
        data = self.hd.ud_counter_rows(tree, cid=2, n_ud=1)
        self.assertIn("2-0-0-0-0", data)
        self.assertNotIn("2-0-0-1-0", data)

    def test_t0_skips_nan(self):
        tree = {
            "ch": np.array([1.0]),
            "offset_p0_val": np.array([[float("nan")]]),
        }
        with mock.patch.object(self.hd, "read_tree", return_value=tree), \
             mock.patch.object(self.hd.config, "DETECTOR_I_LIST", {"BH2": 5}):
            data = self.hd.make_t0_offset_dictdata("/dev/null")
        self.assertEqual(data, {})

    def test_kvc_skips_nan_adc(self):
        tree = {
            "ch": [0],
            "tdc_ok": [1],
            "tdc_p0_val": [[123.0]],
            "adc_p0_val": [[float("nan"), 10.0, 10.0, 10.0, 10.0]],
        }
        data = self.hd.kvc_rows(tree, 6, 5)
        self.assertNotIn("6-0-0-0-0", data)
        self.assertIn("6-0-0-0-1", data)
        self.assertIn("6-0-0-1-4", data)
        for vals in data.values():
            self.assertTrue(all(math.isfinite(float(x)) for x in vals))

    def test_bac_skips_nan_adc(self):
        tree = {
            "ch": np.array([0.0, 1.0]),
            "adc_ok": np.array([1, 1]),
            "tdc_ok": np.array([1]),
            "tdc_p0_val": np.array([[100.0]]),
            "adc_p0_val": np.array([[float("nan")], [5.0]]),
        }
        data = self.hd.bac_rows(tree, cid=3)
        self.assertNotIn("3-0-0-0-0", data)
        self.assertEqual(data.get("3-0-1-0-0"), [5.0])
        self.assertEqual(data.get("3-0-4-1-0")[0], 100.0)

    def test_htof_tdc_skips_inf(self):
        tree = {
            "ch": np.array([0.0]),
            "tdc_p0_val": np.array([[1.0, float("inf"), 3.0]]),
        }
        data = self.hd.htof_tdc_rows(tree, cid=7)
        self.assertIn("7-0-0-1-0", data)
        self.assertNotIn("7-0-0-1-1", data)
        self.assertIn("7-0-0-1-2", data)

    def test_htofprm_skips_nan(self):
        tree = {
            "ch": np.array([0.0]),
            "adc_p0_val": np.array([[1.0, float("nan"), 3.0]]),
            "adc_p1_val": np.array([[10.0, 20.0, 30.0]]),
        }
        fake = mock.MagicMock()
        fake.__getitem__.return_value.arrays.return_value = tree
        with mock.patch.object(self.hd.uproot, "open", return_value=fake), \
             mock.patch.object(self.hd.config, "DETECTOR_I_LIST", {"HTOF": 7}):
            data = self.hd.make_htofprm_dictdata("/dev/null")
        self.assertIn("7-0-0-0-0", data)
        self.assertNotIn("7-0-0-0-1", data)
        self.assertIn("7-0-0-0-2", data)


class ResidualUpdateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        scripts = str(REPO / "calibration" / "scripts")
        if scripts not in sys.path:
            sys.path.insert(0, scripts)
        import update_residual

        cls.res = update_residual

    def test_nan_plane_raises(self):
        tree = {
            "residual_p0_val": np.array([[0.0, float("nan")], [0.0, 0.0]]),
            "residual_p1_val": np.array([[0.1, 0.1], [0.1, 0.1]]),
        }
        fake = mock.MagicMock()
        fake.__contains__.side_effect = lambda k: k == "tree"
        fake.__getitem__.return_value.arrays.return_value = tree
        with mock.patch.object(self.res.uproot, "open", return_value=fake):
            with self.assertRaises(self.res.ResidualFitInvalid):
                self.res.make_dictdata("/tmp/BLC2_residual.root")

    def test_nan_first_also_raises(self):
        tree = {
            "residual_p0_val": np.array([[float("nan"), 0.0], [0.0, 0.0]]),
            "residual_p1_val": np.array([[0.1, 0.1], [0.1, 0.1]]),
        }
        fake = mock.MagicMock()
        fake.__contains__.side_effect = lambda k: k == "tree"
        fake.__getitem__.return_value.arrays.return_value = tree
        with mock.patch.object(self.res.uproot, "open", return_value=fake):
            with self.assertRaises(self.res.ResidualFitInvalid):
                self.res.make_dictdata("/tmp/BLC1_residual.root")


class D5PosScoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        scripts = str(REPO / "calibration" / "scripts")
        if scripts not in sys.path:
            sys.path.insert(0, scripts)
        import run_d5

        cls.d5 = run_d5

    def test_nan_not_converged_score(self):
        r = {"dy_a": 0.0, "dy_b": float("nan"), "dy_err": 0.0, "s": 0.0}
        self.assertGreaterEqual(self.d5.pos_score(r, 0.01), 1.0)
        r2 = {"dy_a": 0.0, "dy_b": 0.0, "dy_err": 0.0, "s": float("nan")}
        self.assertGreaterEqual(self.d5.pos_score(r2, 0.01), 1.0)

    def test_judge_pos_settled_on_nan(self):
        r = {"dy_a": 0.0, "dy_b": float("nan"), "dy_err": 0.0, "s": 0.0}
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            with mock.patch.object(self.d5.config, "OUTPUT_DIR", out), \
                 mock.patch.object(self.d5.calib_anomaly, "remember_check_reason"), \
                 mock.patch.object(self.d5.conv, "record") as rec:
                state = self.d5.judge_pos(1, r, "Pi", 0.01)
        self.assertEqual(state, self.d5.conv.SETTLED)
        rec.assert_called()

    def _write_pos_txt(self, directory, body):
        path = Path(directory) / "d5_pos_review_run00001_Pi.txt"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(body, encoding="utf-8")
        return path

    def _measure_with_txt(self, body):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            self._write_pos_txt(out, body)
            with mock.patch.object(self.d5.runner, "sh", return_value=(0, "")), \
                 mock.patch.object(self.d5, "txt_dir", return_value=out), \
                 mock.patch.object(self.d5.calib_anomaly, "remember_check_reason"):
                return self.d5.measure(1, "Pi")

    def test_measure_rejects_nan_err_b(self):
        body = (
            "dy_a_mm 0\n"
            "dy_b_mm 0\n"
            "dy_a_err_mm 0.01\n"
            "dy_b_err_mm nan\n"
            "blc1_ab_opposite_shift_s_mm 0\n"
            "ddxu_mean_percent 0\n"
            "dv_mean_mrad 0\n"
        )
        self.assertIsNone(self._measure_with_txt(body))

    def test_measure_rejects_nan_err_a(self):
        body = (
            "dy_a_mm 0\n"
            "dy_b_mm 0\n"
            "dy_a_err_mm nan\n"
            "dy_b_err_mm 0.01\n"
            "blc1_ab_opposite_shift_s_mm 0\n"
            "ddxu_mean_percent 0\n"
            "dv_mean_mrad 0\n"
        )
        self.assertIsNone(self._measure_with_txt(body))

    def test_measure_rejects_inf_and_nonpositive_err(self):
        body_inf = (
            "dy_a_mm 0\n"
            "dy_b_mm 0\n"
            "dy_a_err_mm 0.01\n"
            "dy_b_err_mm inf\n"
            "blc1_ab_opposite_shift_s_mm 0\n"
            "ddxu_mean_percent 0\n"
            "dv_mean_mrad 0\n"
        )
        body_zero = (
            "dy_a_mm 0\n"
            "dy_b_mm 0\n"
            "dy_a_err_mm 0.0\n"
            "dy_b_err_mm 0.01\n"
            "blc1_ab_opposite_shift_s_mm 0\n"
            "ddxu_mean_percent 0\n"
            "dv_mean_mrad 0\n"
        )
        self.assertIsNone(self._measure_with_txt(body_inf))
        self.assertIsNone(self._measure_with_txt(body_zero))

    def test_measure_accepts_finite_errs(self):
        body = (
            "dy_a_mm 0\n"
            "dy_b_mm 0\n"
            "dy_a_err_mm 0.01\n"
            "dy_b_err_mm 0.02\n"
            "blc1_ab_opposite_shift_s_mm 0\n"
            "ddxu_mean_percent 0\n"
            "dv_mean_mrad 0\n"
        )
        r = self._measure_with_txt(body)
        self.assertIsNotNone(r)
        self.assertEqual(r["dy_err"], 0.02)
        self.assertLess(self.d5.pos_score(r, 0.01), 1.0)


if __name__ == "__main__":
    unittest.main()
