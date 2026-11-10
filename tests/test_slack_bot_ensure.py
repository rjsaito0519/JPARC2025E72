"""Unit tests for cluster-wide Slack bot ensure (heartbeat + lock)."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from datetime import datetime, timedelta, timezone
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]


class SlackBotEnsureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))

    def setUp(self):
        self._tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmpdir.cleanup)
        self.root = Path(self._tmpdir.name)
        self.out = self.root / "results"
        self.out.mkdir()
        patcher = mock.patch(
            "common.lib.slack_bot_ensure.config.OUTPUT_DIR", self.out,
        )
        patcher.start()
        self.addCleanup(patcher.stop)
        from common.lib import slack_bot_ensure as sbe

        self.sbe = sbe

    def _write_hb(self, **fields):
        now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        payload = {
            "host": "otherhost",
            "pid": 12345,
            "tmux_session": self.sbe.TMUX_SESSION,
            "started_at": now,
            "updated_at": now,
        }
        payload.update(fields)
        path = self.sbe.heartbeat_path()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload) + "\n", encoding="utf-8")
        return payload

    def test_fresh_remote_heartbeat_skips_spawn(self):
        self._write_hb(host="cw01", pid=99999)
        with mock.patch.object(
            self.sbe.local_notify_config,
            "load_slack_bot_config",
            return_value=({"token": "x"}, "ok"),
        ), mock.patch.object(
            self.sbe, "_spawn_bot_tmux",
        ) as spawn_tmux, mock.patch.object(
            self.sbe, "_spawn_bot_popen",
        ) as spawn_popen:
            ok, msg = self.sbe.ensure_slack_bot()
        self.assertTrue(ok)
        self.assertIn("already running on cw01", msg)
        spawn_tmux.assert_not_called()
        spawn_popen.assert_not_called()

    def test_stale_heartbeat_spawns_once(self):
        old = (
            datetime.now(timezone.utc) - timedelta(seconds=300)
        ).strftime("%Y-%m-%dT%H:%M:%SZ")
        self._write_hb(host="cw01", pid=1, updated_at=old)
        fake_hb = {
            "host": self.sbe._hostname(),
            "pid": 4242,
            "tmux_session": self.sbe.TMUX_SESSION,
            "started_at": self.sbe._utc_now(),
            "updated_at": self.sbe._utc_now(),
        }

        def fake_wait():
            self.sbe.heartbeat_path().write_text(
                json.dumps(fake_hb) + "\n", encoding="utf-8",
            )
            return True

        with mock.patch.object(
            self.sbe.local_notify_config,
            "load_slack_bot_config",
            return_value=({"token": "x"}, "ok"),
        ), mock.patch.object(
            self.sbe, "_tmux_available", return_value=True,
        ), mock.patch.object(
            self.sbe, "_spawn_bot_tmux",
            return_value=(True, "tmux ok"),
        ) as spawn_tmux, mock.patch.object(
            self.sbe, "_wait_for_local_heartbeat", side_effect=fake_wait,
        ), mock.patch.object(
            self.sbe, "BOT_SCRIPT", self.root / "bot.py",
        ):
            (self.root / "bot.py").write_text("# fake\n", encoding="utf-8")
            ok, msg = self.sbe.ensure_slack_bot()
        self.assertTrue(ok, msg)
        self.assertIn("ensured", msg)
        spawn_tmux.assert_called_once()

    def test_bot_disabled_skips(self):
        with mock.patch.object(
            self.sbe.local_notify_config,
            "load_slack_bot_config",
            return_value=(None, "bot_disabled"),
        ), mock.patch.object(
            self.sbe, "_spawn_bot_tmux",
        ) as spawn_tmux:
            ok, msg = self.sbe.ensure_slack_bot()
        self.assertTrue(ok)
        self.assertIn("bot_disabled", msg)
        spawn_tmux.assert_not_called()

    def test_local_dead_pid_is_stale(self):
        now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        data = {
            "host": self.sbe._hostname(),
            "pid": 1,
            "updated_at": now,
        }
        with mock.patch.object(self.sbe, "_pid_alive", return_value=False):
            self.assertFalse(self.sbe.heartbeat_is_fresh(data))


if __name__ == "__main__":
    unittest.main()
