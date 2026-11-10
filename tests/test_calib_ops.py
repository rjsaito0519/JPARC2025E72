"""Run with python3 -m unittest discover -s tests -v -p 'test_calib_ops.py'."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]


class CalibOpsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "out"
        self.out.mkdir()
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import calib_ops, config

        self.ops = calib_ops
        self.cfg_mod = config
        self._old_out = config.OUTPUT_DIR
        config.OUTPUT_DIR = self.out

    def tearDown(self):
        self.cfg_mod.OUTPUT_DIR = self._old_out

    def test_parse_calib_text(self):
        op, run, err, extras = self.ops.parse_calib_text("status 3774")
        self.assertEqual((op, run, err, extras), ("status", 3774, "", {}))
        op, run, err, extras = self.ops.parse_calib_text("approve 12")
        self.assertEqual((op, run, err), ("approve", 12, ""))
        self.assertEqual(extras, {})
        op, run, err, extras = self.ops.parse_calib_text(
            "approve 12:abcdefghijklmnop",
        )
        self.assertEqual((op, run, err), ("approve", 12, ""))
        self.assertEqual(extras.get("subject_sha"), "abcdefghijklmnop")
        op, run, err, extras = self.ops.parse_calib_text("revoke 12")
        self.assertEqual((op, run, err), ("revoke", 12, ""))
        op, run, err, extras = self.ops.parse_calib_text("hold 12 human flag")
        self.assertEqual(op, "hold")
        self.assertEqual(run, 12)
        self.assertEqual(err, "")
        self.assertEqual(extras.get("reason"), "human flag")
        op, run, err, extras = self.ops.parse_calib_text("decide 9 continue")
        self.assertEqual((op, run, err), ("decide", 9, ""))
        self.assertEqual(extras.get("action"), "continue")
        op, run, err, extras = self.ops.parse_calib_text("proceed 7")
        self.assertEqual((op, run, err), ("proceed", 7, ""))
        self.assertEqual(extras, {})
        op, run, err, extras = self.ops.parse_calib_text("proceed 7:deadbeef")
        self.assertEqual((op, run, err), ("proceed", 7, ""))
        self.assertEqual(extras.get("hold_id"), "deadbeef")
        op, run, err, extras = self.ops.parse_calib_text("retry 7")
        self.assertEqual((op, run, err), ("retry", 7, ""))
        op, run, err, extras = self.ops.parse_calib_text("retry 7:cafebabe")
        self.assertEqual(extras.get("hold_id"), "cafebabe")
        op, run, err, extras = self.ops.parse_calib_text("gate-stop 7")
        self.assertEqual((op, run, err), ("gate-stop", 7, ""))
        op, run, err, extras = self.ops.parse_calib_text("gate-stop 7:aa11bb22")
        self.assertEqual(extras.get("hold_id"), "aa11bb22")
        op, run, err, extras = self.ops.parse_calib_text("hdphc-again 3")
        self.assertEqual((op, run, err, extras), ("hdphc-again", 3, "", {}))
        op, run, err, extras = self.ops.parse_calib_text("phc-accept 2977")
        self.assertEqual((op, run, err, extras), ("phc-accept", 2977, "", {}))
        op, run, err, extras = self.ops.parse_calib_text(
            "phc-range 2977 BHT 15 49",
        )
        self.assertEqual((op, run, err), ("phc-range", 2977, ""))
        self.assertEqual(extras.get("det"), "BHT")
        self.assertEqual(extras.get("first"), 15)
        self.assertEqual(extras.get("last"), 49)
        op, run, err, extras = self.ops.parse_calib_text(
            "flag 12 hdprm BH2 seg1 looks bad",
        )
        self.assertEqual(op, "flag")
        self.assertEqual(run, 12)
        self.assertEqual(extras.get("mode"), "hdprm")
        self.assertEqual(extras.get("channel"), "seg1")
        op, run, err, extras = self.ops.parse_calib_text("clear 12 hdprm BH2 seg1")
        self.assertEqual((op, run), ("clear", 12))
        op, run, err, extras = self.ops.parse_calib_text("hold 12")
        self.assertIsNone(op)
        self.assertIn("hold", err)
        op, run, err, extras = self.ops.parse_calib_text("resume")
        self.assertIsNone(op)
        self.assertIn("usage", err)

    def test_run_dispatch_hold_passes_reason(self):
        yml = self.ops.session_yml(11)
        yml.parent.mkdir(parents=True)
        yml.write_text("run: 11\nsteps:\n  - hodo\n", encoding="utf-8")

        def fake_run(cmd, capture_output=True, text=True, timeout=30.0):
            self.assertEqual(cmd[2], "hold")
            self.assertIn("--reason", cmd)
            self.assertEqual(cmd[cmd.index("--reason") + 1], "need look")
            return mock.Mock(returncode=0, stdout="held\n", stderr="")

        with mock.patch.object(self.ops.subprocess, "run", side_effect=fake_run):
            ok, msg = self.ops.run_dispatch(
                "hold", 11, reason="need look", background_resume=False,
            )
        self.assertTrue(ok)
        self.assertIn("held", msg)

    def test_run_dispatch_hdphc_again_background(self):
        yml = self.ops.session_yml(8)
        yml.parent.mkdir(parents=True)
        yml.write_text(
            "run: 8\nkaon: true\nsteps:\n  - hodo\n", encoding="utf-8",
        )

        class FakeProc:
            pid = 5252

        with mock.patch.object(
            self.ops.subprocess, "Popen", return_value=FakeProc(),
        ) as popen:
            ok, msg = self.ops.run_dispatch("hdphc-again", 8)
        self.assertTrue(ok)
        self.assertIn("hdphc-again started", msg)
        self.assertIn("5252", msg)
        args = popen.call_args[0][0]
        self.assertTrue(str(args[1]).endswith("extra_hdphc_iter.py"))
        self.assertIn("--kaon", args)

    def test_user_channel_allowed(self):
        cfg = {
            "allowed_user_ids": ["U1"],
            "allowed_channel_ids": ["C1"],
        }
        self.assertTrue(self.ops.user_channel_allowed(cfg, "U1", "C1"))
        self.assertFalse(self.ops.user_channel_allowed(cfg, "U2", "C1"))
        self.assertFalse(self.ops.user_channel_allowed(cfg, "U1", "C2"))

    def test_run_dispatch_missing_session(self):
        ok, msg = self.ops.run_dispatch("status", 9)
        self.assertFalse(ok)
        self.assertIn("no session", msg)

    def test_run_dispatch_status_subprocess(self):
        yml = self.ops.session_yml(7)
        yml.parent.mkdir(parents=True)
        yml.write_text("run: 7\nsteps:\n  - hodo\n", encoding="utf-8")

        def fake_run(cmd, capture_output=True, text=True, timeout=30.0):
            self.assertEqual(cmd[2], "status")
            self.assertEqual(Path(cmd[3]), yml)
            return mock.Mock(returncode=0, stdout="status ok\n", stderr="")

        with mock.patch.object(self.ops.subprocess, "run", side_effect=fake_run):
            ok, msg = self.ops.run_dispatch("status", 7, background_resume=False)
        self.assertTrue(ok)
        self.assertIn("status ok", msg)

    def test_run_dispatch_resume_background(self):
        yml = self.ops.session_yml(8)
        yml.parent.mkdir(parents=True)
        yml.write_text("run: 8\nsteps:\n  - hodo\n", encoding="utf-8")

        class FakeProc:
            pid = 4242

        with mock.patch.object(
            self.ops.subprocess, "Popen", return_value=FakeProc(),
        ) as popen:
            ok, msg = self.ops.run_dispatch("resume", 8)
        self.assertTrue(ok)
        self.assertIn("resume started", msg)
        self.assertIn("4242", msg)
        self.assertTrue(popen.called)

    def test_dedup_window(self):
        self.ops.record_op(3, "stop", "U1")
        dup, msg = self.ops.check_dedup(3, "stop", "U1", window_s=60)
        self.assertTrue(dup)
        self.assertIn("dedup", msg)
        dup2, _ = self.ops.check_dedup(3, "status", "U1", window_s=60)
        self.assertFalse(dup2)

    def _write_held_state(self, run, hold_id, source="step_gate"):
        yml = self.ops.session_yml(run)
        yml.parent.mkdir(parents=True, exist_ok=True)
        yml.write_text("run: %d\nsteps:\n  - hodo\n" % run, encoding="utf-8")
        state = yml.parent / "calib_dispatch_state.json"
        state.write_text(json.dumps({
            "run": run,
            "status": "held",
            "hold": {
                "source": source,
                "step": "hodo_hdprm",
                "index": 0,
                "id": hold_id,
                "decision": None,
            },
        }) + "\n", encoding="utf-8")

    def test_proceed_hold_id_mismatch_rejected(self):
        self._write_held_state(21, "currentid01")
        ok, msg = self.ops.run_dispatch(
            "proceed", 21, hold_id="staleid9999", background_resume=False,
        )
        self.assertFalse(ok)
        self.assertIn("stale step gate", msg)

    def test_proceed_hold_id_match_calls_decide(self):
        self._write_held_state(22, "matchid0001")
        calls = []

        def fake_run(cmd, capture_output=True, text=True, timeout=30.0):
            calls.append(cmd[2])
            return mock.Mock(returncode=0, stdout=cmd[2] + " ok\n", stderr="")

        with mock.patch.object(self.ops.subprocess, "run", side_effect=fake_run), \
             mock.patch.object(
                 self.ops.subprocess, "Popen",
                 return_value=mock.Mock(pid=99),
             ) as popen:
            ok, msg = self.ops.run_dispatch(
                "proceed", 22, hold_id="matchid0001", background_resume=True,
            )
        self.assertTrue(ok, msg)
        self.assertEqual(calls, ["decide"])
        self.assertTrue(popen.called)

    def test_proceed_with_gate_waiter_skips_resume(self):
        self._write_held_state(23, "waiterid001")
        state_path = self.out / "calib/run00023/calib_dispatch_state.json"
        st = json.loads(state_path.read_text(encoding="utf-8"))
        st["gate_waiter"] = True
        st["gate_waiter_pid"] = 4242
        state_path.write_text(json.dumps(st) + "\n", encoding="utf-8")
        calls = []

        def fake_run(cmd, capture_output=True, text=True, timeout=30.0):
            calls.append(cmd[2])
            return mock.Mock(returncode=0, stdout=cmd[2] + " ok\n", stderr="")

        with mock.patch.object(self.ops.subprocess, "run", side_effect=fake_run), \
             mock.patch.object(
                 self.ops.subprocess, "Popen",
                 return_value=mock.Mock(pid=99),
             ) as popen:
            ok, msg = self.ops.run_dispatch(
                "proceed", 23, hold_id="waiterid001", background_resume=True,
            )
        self.assertTrue(ok, msg)
        self.assertEqual(calls, ["decide"])
        self.assertFalse(popen.called)
        self.assertIn("no resume", msg)

    def test_proceed_exit_check_waiter_skips_resume(self):
        self._write_held_state(24, "exithold001", source="exit_check")
        state_path = self.out / "calib/run00024/calib_dispatch_state.json"
        st = json.loads(state_path.read_text(encoding="utf-8"))
        st["gate_waiter"] = True
        state_path.write_text(json.dumps(st) + "\n", encoding="utf-8")
        calls = []

        def fake_run(cmd, capture_output=True, text=True, timeout=30.0):
            calls.append(cmd[2])
            return mock.Mock(returncode=0, stdout=cmd[2] + " ok\n", stderr="")

        with mock.patch.object(self.ops.subprocess, "run", side_effect=fake_run), \
             mock.patch.object(
                 self.ops.subprocess, "Popen",
                 return_value=mock.Mock(pid=99),
             ) as popen:
            ok, msg = self.ops.run_dispatch(
                "proceed", 24, hold_id="exithold001", background_resume=True,
            )
        self.assertTrue(ok, msg)
        self.assertEqual(calls, ["decide"])
        self.assertFalse(popen.called)
        self.assertIn("no resume", msg)


class SlackBotConfigTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import local_notify_config as lnc

        self.lnc = lnc
        self.path = Path(self.tmp.name) / "slack_local.yml"

    def test_bot_disabled_and_ok(self):
        self.path.write_text(
            "enabled: true\ntoken: t\nchannel_id: C1\nbot_enabled: false\n",
            encoding="utf-8",
        )
        cfg, st = self.lnc.load_slack_bot_config(self.path)
        self.assertIsNone(cfg)
        self.assertEqual(st, "bot_disabled")

        self.path.write_text(
            "enabled: true\n"
            "token: xoxb-t\n"
            "channel_id: C1\n"
            "bot_enabled: true\n"
            "app_token: xapp-t\n"
            "allowed_user_ids: [U1]\n",
            encoding="utf-8",
        )
        cfg, st = self.lnc.load_slack_bot_config(self.path)
        self.assertEqual(st, "ok")
        self.assertEqual(cfg["allowed_channel_ids"], ["C1"])
        self.assertEqual(cfg["allowed_user_ids"], ["U1"])


class SlackBotHandlerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "out"
        self.out.mkdir()
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        scripts = str(REPO / "calibration" / "scripts")
        if scripts not in sys.path:
            sys.path.insert(0, scripts)
        from common.lib import config
        import calib_slack_bot as bot

        self.bot = bot
        self.cfg_mod = config
        self._old = config.OUTPUT_DIR
        config.OUTPUT_DIR = self.out
        self.bot_cfg = {
            "token": "x",
            "app_token": "y",
            "channel_id": "C1",
            "allowed_user_ids": ["U1"],
            "allowed_channel_ids": ["C1"],
        }

    def tearDown(self):
        self.cfg_mod.OUTPUT_DIR = self._old

    def test_deny_and_dedup_and_dispatch(self):
        ok, msg, eph = self.bot.handle_calib_request(
            self.bot_cfg, "U2", "C1", "status 1",
        )
        self.assertFalse(ok)
        self.assertTrue(eph)
        self.assertIn("not allowed", msg)

        with mock.patch(
            "common.lib.calib_ops.run_dispatch", return_value=(True, "ok"),
        ):
            ok, msg, eph = self.bot.handle_calib_request(
                self.bot_cfg, "U1", "C1", "status 1",
            )
            self.assertTrue(ok)
            self.assertFalse(eph)
            self.assertIn("status", msg)
            ok2, msg2, eph2 = self.bot.handle_calib_request(
                self.bot_cfg, "U1", "C1", "status 1",
            )
        self.assertTrue(ok2)
        self.assertTrue(eph2)
        self.assertIn("dedup", msg2)


if __name__ == "__main__":
    unittest.main()
