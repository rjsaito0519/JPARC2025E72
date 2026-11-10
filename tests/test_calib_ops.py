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
        op, run, err = self.ops.parse_calib_text("status 3774")
        self.assertEqual((op, run, err), ("status", 3774, ""))
        op, run, err = self.ops.parse_calib_text("resume")
        self.assertIsNone(op)
        self.assertIn("usage", err)

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
