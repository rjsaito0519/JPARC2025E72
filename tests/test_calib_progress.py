"""progress.jsonl substages and awaiting summary."""
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


class CalibProgressTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "out"
        self.out.mkdir()
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import calib_anomaly, calib_approval, calib_progress, config

        self.progress = calib_progress
        self.approval = calib_approval
        self.anomaly = calib_anomaly
        self.cfg = config
        self._old = config.OUTPUT_DIR
        config.OUTPUT_DIR = self.out

    def tearDown(self):
        self.cfg.OUTPUT_DIR = self._old

    def test_note_and_latest_substage(self):
        self.progress.note_substage(12, "hodo", "hdprm", rc=0, status="done")
        self.progress.note_substage(
            12, "hodo", "hdphc[1]", rc=3, status="stopped",
        )
        latest = self.progress.latest_substage(12)
        self.assertIsNotNone(latest)
        self.assertEqual(latest["stage"], "hdphc[1]")
        self.assertEqual(latest["rc"], 3)
        self.assertEqual(latest["step"], "hodo")

    def test_awaiting_stage_summary(self):
        self.progress.note_substage(
            9, "hodo", "hdphc[1]", rc=3, status="stopped",
        )
        self.anomaly.remember_check_reason(
            9, "phc/BHT anomaly: 28/63 (>= 25%) skipped/flagged [PHC:seg0]",
        )
        st = {
            "step": "hodo",
            "rc": 3,
            "plan": {
                "steps": [
                    {"name": "hodo"},
                    {"name": "dc_chamber"},
                ],
            },
        }
        text = self.approval.awaiting_stage_summary(9, st)
        self.assertIn("plan=[hodo, dc_chamber]", text)
        self.assertIn("last=hodo", text)
        self.assertIn("rc=3", text)
        self.assertIn("stage=hdphc[1]", text)
        self.assertIn("phc/BHT anomaly", text)


if __name__ == "__main__":
    unittest.main()
