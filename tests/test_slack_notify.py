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

    def test_ensure_owner_only_clears_group_bits(self):
        path = self.cfg_dir / "secret.yml"
        path.write_text("enabled: true\ntoken: t\nchannel_id: C1\n", encoding="utf-8")
        path.chmod(0o644)
        result = self.py(
            "from pathlib import Path\n"
            "from common.lib import local_notify_config as m\n"
            f"p = Path({str(path)!r})\n"
            "m.ensure_owner_only(p)\n"
            "print(oct(p.stat().st_mode & 0o777))\n"
            "_, st = m.load_slack_config(p)\n"
            "print(st)\n"
            "print(oct(p.stat().st_mode & 0o777))\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.strip().splitlines()
        self.assertEqual(lines[0], "0o600")
        self.assertEqual(lines[1], "ok")
        self.assertEqual(lines[2], "0o600")


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

    def test_notify_failure_dedup(self):
        calls = []

        def fake_api(method, token, payload, as_json=True):
            calls.append((method, payload.get("text")))
            if method == "chat.postMessage":
                return True, {"ok": True, "ts": "333.444"}, "ok"
            return False, {}, f"unexpected {method}"

        with mock.patch.object(self.sn, "_api", side_effect=fake_api):
            ok1, msg1, posted1 = self.sn.notify_failure(
                3, "hodo_merge_hdprm_Pi", "[CHECK] merge failed",
                config_path=self.cfg,
            )
            ok2, msg2, posted2 = self.sn.notify_failure(
                3, "hodo_merge_hdprm_Pi", "[CHECK] merge failed",
                config_path=self.cfg,
            )
        self.assertTrue(ok1, msg1)
        self.assertTrue(posted1)
        self.assertTrue(ok2, msg2)
        self.assertFalse(posted2)
        self.assertIn("dedup", msg2)
        # parent create + failure note
        self.assertEqual(sum(1 for m, _ in calls if m == "chat.postMessage"), 2)

    def test_ensure_run_thread_recreates_on_channel_change(self):
        self.sn.save_state(10, {
            "channel_id": "COLD",
            "thread_ts": "1.2",
            "posts": {"k": {"pdf": "old.pdf"}},
            "failures": {"f": {"text": "old"}},
        })
        self.cfg.write_text(
            "enabled: true\ntoken: xoxb-test\nchannel_id: CNEW\n", encoding="utf-8"
        )
        parent_calls = []

        def fake_api(method, token, payload, as_json=True):
            if method == "chat.postMessage":
                parent_calls.append(payload)
                return True, {"ok": True, "ts": "9.9"}, "ok"
            return False, {}, f"unexpected {method}"

        with mock.patch.object(self.sn, "_api", side_effect=fake_api):
            ok, msg, ts = self.sn.ensure_run_thread(
                10, "calib QA run00010", config_path=self.cfg,
            )
        self.assertTrue(ok, msg)
        self.assertEqual(ts, "9.9")
        self.assertEqual(len(parent_calls), 1)
        self.assertEqual(parent_calls[0]["channel"], "CNEW")
        state = json.loads(self.sn.state_path(10).read_text())
        self.assertEqual(state["channel_id"], "CNEW")
        self.assertEqual(state["thread_ts"], "9.9")
        self.assertEqual(state["posts"], {})
        self.assertEqual(state["failures"], {})
        self.assertEqual(len(state["history"]), 1)
        self.assertEqual(state["history"][0]["channel_id"], "COLD")
        self.assertEqual(state["history"][0]["thread_ts"], "1.2")
        self.assertEqual(state["history"][0]["posts"], {"k": {"pdf": "old.pdf"}})
        self.assertEqual(state["history"][0]["failures"], {"f": {"text": "old"}})

    def test_ensure_run_thread_recreates_when_channel_missing(self):
        self.sn.save_state(11, {"thread_ts": "1.2", "posts": {}, "failures": {}})
        parent_calls = []

        def fake_api(method, token, payload, as_json=True):
            if method == "chat.postMessage":
                parent_calls.append(payload)
                return True, {"ok": True, "ts": "2.3"}, "ok"
            return False, {}, f"unexpected {method}"

        with mock.patch.object(self.sn, "_api", side_effect=fake_api):
            ok, msg, ts = self.sn.ensure_run_thread(
                11, "calib QA run00011", config_path=self.cfg,
            )
        self.assertTrue(ok, msg)
        self.assertEqual(ts, "2.3")
        self.assertEqual(len(parent_calls), 1)
        state = json.loads(self.sn.state_path(11).read_text())
        self.assertEqual(state["thread_ts"], "2.3")
        self.assertEqual(state["channel_id"], "CCHAN")
        self.assertEqual(state["history"][0]["thread_ts"], "1.2")

    def test_notify_failure_reposts_after_channel_change(self):
        self.cfg.write_text(
            "enabled: true\ntoken: xoxb-test\nchannel_id: COLD\n", encoding="utf-8"
        )
        seq = {"n": 0}

        def fake_api(method, token, payload, as_json=True):
            if method == "chat.postMessage":
                seq["n"] += 1
                return True, {"ok": True, "ts": f"{seq['n']}.0"}, "ok"
            return False, {}, f"unexpected {method}"

        with mock.patch.object(self.sn, "_api", side_effect=fake_api):
            ok1, msg1, posted1 = self.sn.notify_failure(
                12, "hodo_merge_hdprm_Pi", "[CHECK] merge failed",
                config_path=self.cfg,
            )
            self.assertTrue(ok1, msg1)
            self.assertTrue(posted1)
            self.cfg.write_text(
                "enabled: true\ntoken: xoxb-test\nchannel_id: CNEW\n",
                encoding="utf-8",
            )
            ok2, msg2, posted2 = self.sn.notify_failure(
                12, "hodo_merge_hdprm_Pi", "[CHECK] merge failed",
                config_path=self.cfg,
            )
        self.assertTrue(ok2, msg2)
        self.assertTrue(posted2, msg2)
        state = json.loads(self.sn.state_path(12).read_text())
        self.assertEqual(state["channel_id"], "CNEW")
        self.assertIn("hodo_merge_hdprm_Pi", state["failures"])
        self.assertEqual(len(state["history"]), 1)
        self.assertEqual(state["history"][0]["channel_id"], "COLD")

    def test_notify_failure_skip_without_config(self):
        missing = Path(self.tmp.name) / "no_such_slack.yml"
        ok, msg, posted = self.sn.notify_failure(
            4, "dc_resi_bcout_Pi", "[CHECK] x", config_path=missing,
        )
        self.assertTrue(ok)
        self.assertFalse(posted)
        self.assertIn("slack skip", msg)

    def test_api_timeout_becomes_false_not_raise(self):
        with mock.patch(
            "urllib.request.urlopen", side_effect=TimeoutError("timed out"),
        ):
            ok, body, msg = self.sn._api(
                "chat.postMessage", "xoxb-test", {"channel": "C", "text": "t"},
            )
        self.assertFalse(ok)
        self.assertEqual(body, {})
        self.assertIn("TimeoutError", msg)
        self.assertNotIn("xoxb-", msg)

    def test_post_pdf_timeout_does_not_raise(self):
        pdf = self.out / "img/run00005/x.pdf"
        pdf.parent.mkdir(parents=True, exist_ok=True)
        pdf.write_bytes(b"%PDF")

        def boom_api(*_a, **_k):
            raise TimeoutError("read timed out")

        with mock.patch.object(self.sn, "_api", side_effect=boom_api):
            ok, msg = self.sn.post_pdf(
                5, "hdprm", pdf, "cap", config_path=self.cfg,
            )
        self.assertFalse(ok)
        self.assertIn("TimeoutError", msg)

    def test_save_state_oserror_does_not_raise(self):
        pdf = self.out / "img/run00006/x.pdf"
        pdf.parent.mkdir(parents=True, exist_ok=True)
        pdf.write_bytes(b"%PDF")

        def fake_api(method, token, payload, as_json=True):
            if method == "chat.postMessage":
                return True, {"ok": True, "ts": "1.2"}, "ok"
            if method == "files.getUploadURLExternal":
                return True, {
                    "ok": True,
                    "upload_url": "http://127.0.0.1/upload",
                    "file_id": "F9",
                }, "ok"
            if method == "files.completeUploadExternal":
                return True, {"ok": True}, "ok"
            return False, {}, method

        with mock.patch.object(self.sn, "_api", side_effect=fake_api), \
             mock.patch("urllib.request.urlopen") as urlopen, \
             mock.patch.object(self.sn, "save_state", return_value=False):
            urlopen.return_value.__enter__.return_value.status = 200
            urlopen.return_value.__enter__.return_value.read.return_value = b""
            ok, msg = self.sn.post_pdf(
                6, "hdprm", pdf, "cap", config_path=self.cfg,
            )
        self.assertFalse(ok)
        self.assertIn("state save failed", msg)


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

    def test_after_qa_merge_keeps_rc_when_try_post_times_out(self):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from calibration.scripts import run_hodo as rh
        from common.lib import runner

        class Args:
            kaon = False
            ftof = False

        with mock.patch.object(
            rh.slack_notify, "post_hodo_qa",
            side_effect=TimeoutError("injected"),
        ), mock.patch.object(
            rh.slack_notify, "load_slack_config",
            return_value=({"token": "x", "channel_id": "C"}, "ok"),
        ):
            ok, msg = rh.slack_notify.try_post_hodo_after_merge(1, "hdprm")
        self.assertFalse(ok)
        self.assertIn("TimeoutError", msg)

        with mock.patch.object(
            rh.hodo_qa_pdf, "merge",
            return_value=(True, "merged", Path("/tmp/x.pdf")),
        ), mock.patch.object(
            rh.slack_notify, "try_post_hodo_after_merge",
            return_value=(False, "try_post_hodo_after_merge: TimeoutError"),
        ), mock.patch.object(rh.runner, "info"), mock.patch.object(
            rh.runner, "warn"
        ):
            rc = rh.after_qa_merge(1, "hdprm", Args(), runner.EXIT_OK)
        self.assertEqual(rc, runner.EXIT_OK)


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


class TopSymlinkPolicyTests(unittest.TestCase):
    def test_no_post_merge_calib_set_at_repo_root(self):
        banned = [
            "post_hodo_qa_slack.py",
            "post_dc_qa_slack.py",
            "post_d5_qa_slack.py",
            "merge_hodo_qa_pdf.py",
            "calib_dispatch.py",
            "calib_session.py",
            "set_tdc_window.py",
            "iterate_hodo_chain.py",
        ]
        for name in banned:
            self.assertFalse((REPO / name).exists(), name)
        for name in ("calib.py", "ops.py", "run_hodo.py", "create_runlist.py"):
            self.assertTrue((REPO / name).is_symlink(), name)


if __name__ == "__main__":
    unittest.main()
