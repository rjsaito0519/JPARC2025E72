"""Lane attach scope, Stop≠done, awaiting→deps, single-lane approve."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]


class LaneSemanticsTests(unittest.TestCase):
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

        self.cfg_mod = config
        self._old = config.OUTPUT_DIR
        config.OUTPUT_DIR = self.out
        import calib_dispatch as cd  # noqa: WPS433

        self.cd = cd

    def tearDown(self):
        self.cfg_mod.OUTPUT_DIR = self._old

    def _write_lane(self, run, lane, payload):
        path = self.cd.state_path(run, lane=lane)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2) + "\n")

    def test_do_attach_skips_stopped_starts_only_want(self):
        run = 11
        self._write_lane(run, "bcin", {
            "status": "stopped", "finished": False, "run": run,
        })
        self._write_lane(run, "bcout", {
            "status": "held", "finished": False, "run": run,
            "hold": {"decision": None, "id": "abc"},
        })
        cfg = {
            "run": run,
            "steps": [
                "hodo_hdprm",
                "dc_t0",
                {"dc_t0": {"bcin": True}},
            ],
            "session": {},
            "kaon": True,
            "local": True,
            "auto": False,
            "force": False,
            "max_iter": 10,
        }
        with mock.patch.object(
            self.cd, "run_lane_coordinator",
        ) as coord, mock.patch.object(
            self.cd, "_ensure_slack_bot_for_gates",
        ):
            self.cd.do_attach(cfg, dry_run=True, want_lanes=["hodo"])
        coord.assert_called_once()
        _args, kwargs = coord.call_args
        to_run = _args[2]
        self.assertEqual(list(to_run), ["hodo"])

    def test_wait_external_accepts_awaiting(self):
        run = 12
        self.cd.write_state(run, {
            "status": "running",
            "mode": "lane_coordinator",
            "lanes": {"bcin": "done", "bcout": "done"},
        }, lane="")
        self._write_lane(run, "bcin", {
            "status": "awaiting_final_approval",
            "finished": True,
            "run": run,
        })
        self._write_lane(run, "bcout", {
            "status": "awaiting_final_approval",
            "finished": True,
            "run": run,
        })
        # Would hang forever before the fix; should return immediately.
        self.cd.wait_external_dc_done(
            run,
            {"bcin": ["dc_t0"], "bcout": ["dc_t0"]},
            dry_run=False,
            poll_s=0.01,
        )

    def test_coordinator_outcome_stop_held_not_done(self):
        run = 13
        self._write_lane(run, "bcout", {
            "status": "held",
            "finished": False,
            "hold": {"decision": "stop", "id": "x"},
            "run": run,
        })
        outcome = self.cd._coordinator_outcome_for_lane(run, "bcout", 0)
        self.assertEqual(outcome, "held")

    def test_coordinator_does_not_start_d5_after_stop(self):
        run = 14
        cfg = {
            "run": run,
            "steps": ["dc_t0", "d5_align"],
            "session": {},
            "kaon": True,
            "local": True,
            "auto": False,
            "force": False,
            "max_iter": 10,
        }
        lane_steps = {
            "bcout": ["dc_t0"],
            "d5": [{"d5_align": {"pass": 1}}],
        }
        self._write_lane(run, "bcout", {
            "status": "held",
            "finished": False,
            "hold": {"decision": "stop", "id": "s"},
            "run": run,
        })

        class FakeProc:
            def __init__(self, rc):
                self._rc = rc

            def poll(self):
                return self._rc

            def wait(self, timeout=None):
                return self._rc

            def terminate(self):
                return None

        spawned = []

        def fake_spawn(cfg_arg, lane, steps, dry_run):
            spawned.append(lane)
            return FakeProc(0)

        with mock.patch.object(
            self.cd, "spawn_lane_worker", side_effect=fake_spawn,
        ):
            with self.assertRaises(SystemExit) as caught:
                self.cd.run_lane_coordinator(cfg, dry_run=True, lane_steps=lane_steps)
        self.assertNotIn("d5", spawned)
        top = self.cd.read_state(run, lane="")
        self.assertEqual(top.get("status"), "stopped")
        self.assertEqual((top.get("lanes") or {}).get("bcout"), "held")

    def test_single_lane_approve_reads_lane_state(self):
        run = 15
        self._write_lane(run, "hodo", {
            "status": "awaiting_final_approval",
            "finished": True,
            "run": run,
            "unverified": [],
            "plan": {"steps": [{"name": "hodo_hdprm", "opts": {}}]},
            "session": {"calib_particle": "kaon"},
        })
        # No top-level state.
        cfg = {"run": run}
        with mock.patch.object(
            self.cd.calib_approval, "build_subject",
            return_value={"ok": True, "missing": []},
        ), mock.patch.object(
            self.cd.calib_approval, "subject_missing", return_value=[],
        ), mock.patch.object(
            self.cd.calib_approval, "subject_sha", return_value="abcd" * 4,
        ), mock.patch(
            "common.lib.channel_flags.open_needs_review", return_value=[],
        ):
            with self.assertRaises(SystemExit) as caught:
                self.cd.do_approve(cfg, None, "cli", None, None)
        self.assertEqual(caught.exception.code, 0)
        st = self.cd.read_state(run, lane="hodo")
        self.assertEqual(st.get("status"), "approved")
        self.assertIsNone(self.cd.read_state(run, lane=""))


if __name__ == "__main__":
    unittest.main()
