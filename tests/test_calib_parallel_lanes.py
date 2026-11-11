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
