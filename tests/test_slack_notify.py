"""Run with python3 -m unittest discover -s tests -v."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

REPO = Path(__file__).resolve().parents[1]


class LocalNotifyConfigTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        lib = self.root / "common/lib"
        lib.mkdir(parents=True)
        (lib / "__init__.py").touch()
        for name in ("local_notify_config.py",):
            shutil.copy(REPO / "common/lib" / name, lib)
        self.cfg_dir = self.root / "calibration/config"
        self.cfg_dir.mkdir(parents=True)
        # Point module CONFIG_DIR at our temp tree via import patch after load
        self.env = {**os.environ, "PYTHONPATH": str(self.root)}

    def py(self, code):
        return subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True, text=True, env=self.env, timeout=30,
            cwd=str(self.root))

    def test_missing_and_disabled(self):
        result = self.py(
            "from pathlib import Path\n"
            "from common.lib import local_notify_config as m\n"
            f"m.SLACK_LOCAL = Path({str(self.cfg_dir / 'slack_local.yml')!r})\n"
            "cfg, st = m.load_slack_config()\n"
            "print(st)\n"
            f"p = Path({str(self.cfg_dir / 'slack_local.yml')!r})\n"
            "p.write_text('enabled: false\\ntoken: x\\nchannel_id: C1\\n')\n"
            "cfg, st = m.load_slack_config(p)\n"
            "print(st)\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip().splitlines(), ["missing", "disabled"])

    def test_slack_ok_and_discord_ok(self):
        slack = self.cfg_dir / "slack_local.yml"
        discord = self.cfg_dir / "discord_local.yml"
        slack.write_text(
            "enabled: true\ntoken: xoxb-test\nchannel_id: C123\n", encoding="utf-8"
        )
        discord.write_text(
            "enabled: true\nwebhook_url: https://example.invalid/hook\n",
            encoding="utf-8",
        )
        result = self.py(
            "from pathlib import Path\n"
            "from common.lib import local_notify_config as m\n"
            f"sc, ss = m.load_slack_config(Path({str(slack)!r}))\n"
            f"dc, ds = m.load_discord_config(Path({str(discord)!r}))\n"
            "print(ss, sc['token'], sc['channel_id'])\n"
            "print(ds, dc['webhook_url'])\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.strip().splitlines()
        self.assertEqual(lines[0], "ok xoxb-test C123")
        self.assertEqual(lines[1], "ok https://example.invalid/hook")


class SlackNotifyUnitTests(unittest.TestCase):
    """In-process tests with stubbed HTTP (no real Slack)."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "out"
        self.out.mkdir()
        self.cfg = Path(self.tmp.name) / "slack_local.yml"
        self.cfg.write_text(
            "enabled: true\ntoken: xoxb-test\nchannel_id: CCHAN\n", encoding="utf-8"
        )
        # Import from the real repo
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import config as cfg_mod
        from common.lib import slack_notify as sn

        self.cfg_mod = cfg_mod
        self.sn = sn
        self._old_out = cfg_mod.OUTPUT_DIR
        cfg_mod.OUTPUT_DIR = self.out

    def tearDown(self):
        self.cfg_mod.OUTPUT_DIR = self._old_out

    def test_dry_run_post(self):
        ok, msg = self.sn.post_hodo_qa(1, "hdprm", dry_run=True, config_path=self.cfg)
        self.assertTrue(ok)
        self.assertIn("dry-run", msg)

    def test_thread_state_reuse(self):
        calls = []

        def fake_api(method, token, payload, as_json=True):
            calls.append(method)
            if method == "chat.postMessage":
                return True, {"ok": True, "ts": "111.222"}, "ok"
            if method == "files.getUploadURLExternal":
                return True, {
                    "ok": True,
                    "upload_url": "http://127.0.0.1/upload",
                    "file_id": "F1",
                }, "ok"
            if method == "files.completeUploadExternal":
                return True, {"ok": True}, "ok"
            return False, {}, f"unexpected {method}"

        pdf = self.out / "img/run00001/run00001_hodo_hdprm_Pi.pdf"
        pdf.parent.mkdir(parents=True, exist_ok=True)
        pdf.write_bytes(b"%PDF-1.4 test")

        with mock.patch.object(self.sn, "_api", side_effect=fake_api), \
             mock.patch("urllib.request.urlopen") as urlopen:
            urlopen.return_value.__enter__.return_value.status = 200
            urlopen.return_value.__enter__.return_value.read.return_value = b""
            ok1, msg1 = self.sn.post_pdf(
                1, "hdprm", pdf, "cap", config_path=self.cfg,
            )
            ok2, msg2 = self.sn.post_pdf(
                1, "t0", pdf, "cap2", particle_suffix="Pi", config_path=self.cfg,
            )
        self.assertTrue(ok1, msg1)
        self.assertTrue(ok2, msg2)
        self.assertEqual(calls.count("chat.postMessage"), 1)
        state = json.loads(self.sn.state_path(1).read_text())
        self.assertEqual(state["thread_ts"], "111.222")
        self.assertIn("hdprm|Pi|run00001_hodo_hdprm_Pi.pdf", state["posts"])

    def test_post_failure_does_not_raise(self):
        def boom(*_a, **_k):
            return False, {}, "Slack chat.postMessage: channel_not_found"

        with mock.patch.object(self.sn, "_api", side_effect=boom):
            ok, msg = self.sn.post_hodo_qa(
                2, "hdprm", config_path=self.cfg,
                pdf_path=self.out / "missing.pdf",
            )
        self.assertFalse(ok)
        self.assertIn("channel_not_found", msg)


class RunHodoSlackRcTests(unittest.TestCase):
    """Slack failure must not change calibration rc from after_qa_merge."""

    def test_after_qa_merge_keeps_rc_on_slack_fail(self):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from calibration.scripts import run_hodo as rh
        from common.lib import runner

        class Args:
            kaon = False
            ftof = False

        with mock.patch.object(
            rh.hodo_qa_pdf, "merge",
            return_value=(True, "merged", Path("/tmp/x.pdf")),
        ), mock.patch.object(
            rh.slack_notify, "try_post_hodo_after_merge",
            return_value=(False, "network down"),
        ), mock.patch.object(rh.runner, "info"), mock.patch.object(
            rh.runner, "warn"
        ) as warn:
            rc = rh.after_qa_merge(1, "hdprm", Args(), runner.EXIT_OK)
        self.assertEqual(rc, runner.EXIT_OK)
        warn.assert_called()
        self.assertIn("Slack", warn.call_args[0][0])


class DiscordConfigNoEnvTests(unittest.TestCase):
    def test_iterate_does_not_read_discord_env(self):
        """DISCORD_WEBHOOK must not appear as a config source in iterate_tpc_calib."""
        text = (REPO / "calibration/scripts/iterate_tpc_calib.py").read_text()
        self.assertNotIn("DISCORD_WEBHOOK", text)
        self.assertIn("load_discord_config", text)
        env_hits = [
            line for line in text.splitlines()
            if "os.environ" in line and "DISCORD" in line
        ]
        self.assertEqual(env_hits, [])


if __name__ == "__main__":
    unittest.main()
