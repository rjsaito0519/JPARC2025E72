"""pre-write hold gate (#24 second slice)."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]


class WriteBlockedTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "out"
        self.out.mkdir()
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import calib_approval, config, runner

        self.approval = calib_approval
        self.runner = runner
        self.cfg = config
        self._old = config.OUTPUT_DIR
        config.OUTPUT_DIR = self.out

    def tearDown(self):
        self.cfg.OUTPUT_DIR = self._old

    def _write_state(self, run, payload):
        path = self.approval.dispatch_state_path(run)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload) + "\n", encoding="utf-8")

    def test_no_state_allows_write(self):
        blocked, reason = self.approval.write_blocked(99)
        self.assertFalse(blocked)
        self.assertEqual(reason, "")

    def test_held_blocks(self):
        self._write_state(1, {
            "status": "held",
            "hold": {"reason": "EXIT_CHECK"},
            "hold_requested": False,
        })
        blocked, reason = self.approval.write_blocked(1)
        self.assertTrue(blocked)
        self.assertIn("held", reason)

    def test_hold_requested_blocks(self):
        self._write_state(2, {
            "status": "running",
            "hold_requested": True,
            "hold_reason": "look later",
        })
        blocked, reason = self.approval.write_blocked(2)
        self.assertTrue(blocked)
        self.assertIn("hold_requested", reason)

    def test_awaiting_does_not_block(self):
        self._write_state(3, {
            "status": "awaiting_final_approval",
            "finished": True,
            "hold_requested": False,
        })
        blocked, reason = self.approval.write_blocked(3)
        self.assertFalse(blocked)

    def test_update_param_returns_check_when_held(self):
        self._write_state(4, {
            "status": "held",
            "hold": {"reason": "manual"},
            "hold_requested": False,
        })
        with mock.patch.object(self.runner, "sh") as sh:
            rc = self.runner.update_param(4, "Pi", "hdphc")
        self.assertEqual(rc, self.runner.EXIT_CHECK)
        sh.assert_not_called()

    def test_update_param_runs_when_clear(self):
        with mock.patch.object(
            self.runner, "sh", return_value=(0, ""),
        ) as sh:
            rc = self.runner.update_param(5, "Pi", "hdphc")
        self.assertEqual(rc, 0)
        sh.assert_called_once()

    def test_channel_flag_blocks_update_param(self):
        from common.lib import channel_flags

        channel_flags.append_flag(
            8,
            family="hodo",
            mode="hdprm",
            det="BH2",
            channel="seg1",
            source="human",
            reason="looks bad",
            status="needs_review",
        )
        with mock.patch.object(self.runner, "sh") as sh:
            rc = self.runner.update_param(8, "Pi", "hdprm")
        self.assertEqual(rc, self.runner.EXIT_CHECK)
        sh.assert_not_called()

    def test_do_drift_blocked_when_held(self):
        from calibration.scripts import run_dcprm

        self._write_state(6, {
            "status": "held",
            "hold": {"reason": "manual"},
            "hold_requested": False,
        })
        args = mock.Mock(
            debug=False, force=False, kaon=False, bcin=False, bcout=True,
        )
        with mock.patch.object(
            self.runner, "require_bins", return_value=("BLC_DRFT",),
        ), mock.patch.object(
            self.runner, "input_root", return_value=Path("/tmp/x.root"),
        ), mock.patch.object(
            run_dcprm, "_dc_session", return_value=None,
        ), mock.patch.object(
            run_dcprm, "after_qa_slack",
            side_effect=lambda *a, **k: a[3],
        ), mock.patch.object(self.runner, "sh") as sh:
            rc = run_dcprm.do_drift(6, "Pi", args)
        self.assertEqual(rc, self.runner.EXIT_CHECK)
        sh.assert_not_called()


if __name__ == "__main__":
    unittest.main()
