"""choose_action attach/resume while other lanes run."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]


class ChooseActionParallelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        scripts = str(REPO / "calibration" / "scripts")
        if scripts not in sys.path:
            sys.path.insert(0, scripts)

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "out"
        self.out.mkdir()
        from common.lib import config

        self.cfg = config
        self._old = config.OUTPUT_DIR
        config.OUTPUT_DIR = self.out
        import calib as calib_mod  # noqa: WPS433

        self.calib = calib_mod

    def tearDown(self):
        self.cfg.OUTPUT_DIR = self._old

    def test_bc_running_want_hodo_free_attaches(self):
        st = {
            "status": "running",
            "mode": "lane_coordinator",
            "lanes": {"bcin": "running", "bcout": "running"},
        }
        with mock.patch.object(
            self.calib.cd, "collect_lane_statuses",
            return_value={"bcin": "running", "bcout": "running"},
        ):
            choice = self.calib.choose_action(
                1, st, force=False, interactive=False,
                want_lanes=["hodo"],
            )
        self.assertEqual(choice, self.calib.CHOICE_ATTACH)

    def test_bc_running_want_hodo_held_continues(self):
        st = {
            "status": "running",
            "mode": "lane_coordinator",
            "lanes": {
                "hodo": "held",
                "bcin": "running",
                "bcout": "running",
            },
        }
        with mock.patch.object(
            self.calib.cd, "collect_lane_statuses",
            return_value={
                "hodo": "held",
                "bcin": "running",
                "bcout": "running",
            },
        ):
            choice = self.calib.choose_action(
                1, st, force=False, interactive=False,
                want_lanes=["hodo"],
            )
        self.assertEqual(choice, self.calib.CHOICE_CONTINUE)

    def test_stopped_requested_lane_prompts_before_resuming(self):
        st = {"status": "stopped"}
        for selection in (self.calib.CHOICE_CONTINUE, self.calib.CHOICE_RESTART):
            with self.subTest(selection=selection), mock.patch.object(
                self.calib.cd, "collect_lane_statuses", return_value={"hodo": "stopped"},
            ), mock.patch.object(sys.stdin, "isatty", return_value=True), mock.patch.object(
                self.calib, "prompt_continue_or_restart", return_value=selection,
            ) as prompt:
                choice = self.calib.choose_action(
                    2847, st, force=False, interactive=True, want_lanes=["hodo"],
                )
                self.assertEqual(choice, selection)
                prompt.assert_called_once()

    def test_restart_selection_cannot_disrupt_other_busy_lanes(self):
        with mock.patch.object(self.calib.cd, "collect_lane_statuses",
                               return_value={"hodo": "stopped", "bcin": "running"}), \
                mock.patch.object(sys.stdin, "isatty", return_value=True), \
                mock.patch.object(self.calib, "prompt_continue_or_restart",
                                  return_value=self.calib.CHOICE_RESTART):
            with self.assertRaises(SystemExit):
                self.calib.choose_action(2847, {"status": "running"}, force=False,
                                         interactive=True, want_lanes=["hodo"])

    def test_held_lane_without_waiter_can_restart(self):
        state = {"status": "held", "stop_requested": True}
        with mock.patch.object(self.calib.cd, "collect_lane_statuses", return_value={"hodo": "held"}), \
                mock.patch.object(self.calib.cd, "read_state", return_value=state), \
                mock.patch.object(sys.stdin, "isatty", return_value=True), \
                mock.patch.object(self.calib, "prompt_continue_or_restart", return_value=self.calib.CHOICE_RESTART):
            choice = self.calib.choose_action(3772, state, force=False,
                                             interactive=True, want_lanes=["hodo"])
            self.assertEqual(choice, self.calib.CHOICE_RESTART)
            forced = self.calib.choose_action(3772, state, force=True,
                                              interactive=False, want_lanes=["hodo"])
            self.assertEqual(forced, self.calib.CHOICE_RESTART)

    def test_live_gate_waiter_blocks_restart_even_after_stop_requested(self):
        state = {"status": "held", "stop_requested": True, "gate_waiter": True}
        with mock.patch.object(self.calib.cd, "collect_lane_statuses", return_value={"hodo": "held"}), \
                mock.patch.object(self.calib.cd, "read_state", return_value=state), \
                mock.patch.object(sys.stdin, "isatty", return_value=True), \
                mock.patch.object(self.calib, "prompt_continue_or_restart", return_value=self.calib.CHOICE_RESTART):
            for force in (False, True):
                with self.subTest(force=force), self.assertRaises(SystemExit):
                    self.calib.choose_action(3772, state, force=force,
                                             interactive=True, want_lanes=["hodo"])

    def test_eof_does_not_implicitly_resume(self):
        with mock.patch.object(self.calib, "format_state_summary", return_value="prior"), \
                mock.patch("builtins.input", side_effect=EOFError):
            with self.assertRaises(SystemExit):
                self.calib.prompt_continue_or_restart(2847, {"status": "stopped"})

    def test_want_running_lane_fails(self):
        st = {
            "status": "running",
            "mode": "lane_coordinator",
            "lanes": {"bcin": "running"},
        }
        with mock.patch.object(
            self.calib.cd, "collect_lane_statuses",
            return_value={"bcin": "running"},
        ):
            with self.assertRaises(SystemExit):
                self.calib.choose_action(
                    1, st, force=False, interactive=False,
                    want_lanes=["bcin"],
                )

    def test_continue_coordinator_running_attaches_free(self):
        st = {
            "status": "running",
            "mode": "lane_coordinator",
            "lanes": {"bcin": "running", "bcout": "running"},
        }
        cfg = {"run": 7, "steps": ["hodo_hdprm"]}
        with mock.patch.object(
            self.calib.cd, "collect_lane_statuses",
            return_value={"bcin": "running", "bcout": "running"},
        ), mock.patch.object(self.calib.cd, "do_attach") as attach:
            self.calib.continue_from_state(
                cfg, st, dry_run=True, want_lanes=["hodo"],
            )
        attach.assert_called_once()


if __name__ == "__main__":
    unittest.main()
