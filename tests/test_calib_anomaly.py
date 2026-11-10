"""Provisional anomaly checks (#24 Phase B)."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]


class CalibAnomalyTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "out"
        self.out.mkdir()
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import calib_anomaly, config

        self.anomaly = calib_anomaly
        self.cfg = config
        self._old = config.OUTPUT_DIR
        config.OUTPUT_DIR = self.out

    def tearDown(self):
        self.cfg.OUTPUT_DIR = self._old

    def test_ratio_too_high(self):
        self.assertTrue(self.anomaly.ratio_too_high(3, 10, 0.25))
        self.assertFalse(self.anomaly.ratio_too_high(2, 10, 0.25))

    def test_decide_ratio_message(self):
        reason = self.anomaly.decide_ratio(
            "hdprm/BAC", 3, 10, 0.25, ["a", "b", "c"],
        )
        self.assertIsNotNone(reason)
        self.assertIn("3/10", reason)

    def test_remember_and_read(self):
        self.anomaly.remember_check_reason(7, "test reason")
        self.assertEqual(self.anomaly.read_check_reason(7), "test reason")

    def test_bac_skip_counts(self):
        tree = {
            "ch": [0, 1, 2, 3],
            "adc_ok": [1, 0, 0, 1],
            "tdc_ok": [0],
            "tdc_p0_val": [[0.0]],
        }
        bad, total, labels = self.anomaly.bac_skip_counts(tree)
        self.assertEqual(bad, 3)
        self.assertEqual(total, 5)
        self.assertEqual(len(labels), 3)

    def test_drift_missing_file(self):
        reason = self.anomaly.drift_graph_anomaly(Path(self.tmp.name) / "no.root")
        self.assertIsNotNone(reason)
        self.assertIn("0 TGraph", reason)


if __name__ == "__main__":
    unittest.main()
