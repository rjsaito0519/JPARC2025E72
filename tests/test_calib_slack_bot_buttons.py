"""FakeApp regression for Slack button handlers (#42)."""
import sys
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "calibration" / "scripts"

BUTTON_ACTION_IDS = (
    "calib_approve_run",
    "calib_revoke_run",
    "calib_step_proceed",
    "calib_step_refit",
    "calib_step_stop",
    "calib_phc_go",
    "calib_phc_finish",
)


class FakeApp:
    def __init__(self, token=None):
        self.token = token
        self.actions = {}
        self.commands = {}
        self.views = {}

    def action(self, action_id):
        def decorator(fn):
            self.actions[action_id] = fn
            return fn
        return decorator

    def command(self, name):
        def decorator(fn):
            self.commands[name] = fn
            return fn
        return decorator

    def view(self, callback_id):
        def decorator(fn):
            self.views[callback_id] = fn
            return fn
        return decorator


class SlackBotButtonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        if str(SCRIPTS) not in sys.path:
            sys.path.insert(0, str(SCRIPTS))

    def test_button_handlers_accept_dict_user(self):
        import calib_slack_bot as bot

        bot_cfg = {
            "token": "xoxb-test",
            "allowed_user_ids": ["U123"],
            "allowed_channel_ids": ["C123"],
        }
        with mock.patch.dict(sys.modules, {"slack_bolt": mock.Mock(App=FakeApp)}):
            # Re-import path: build_app imports App from slack_bolt inside.
            app = bot.build_app(bot_cfg)

        self.assertTrue(set(BUTTON_ACTION_IDS).issubset(set(app.actions)))
        self.assertIn("calib_phc_custom", app.actions)
        self.assertIn("calib_phc_custom_modal", app.views)

        responses = []

        def respond(**kwargs):
            responses.append(kwargs)

        body = {
            "user": {"id": "U123"},
            "channel": {"id": "C123"},
            "actions": [{"value": "7:holdid01"}],
        }

        with mock.patch.object(
            bot, "handle_calib_request",
            return_value=(True, "ok", False),
        ) as handle:
            for action_id in BUTTON_ACTION_IDS:
                responses.clear()
                handler = app.actions[action_id]
                handler(ack=lambda: None, body=body, respond=respond)
                self.assertEqual(len(responses), 1, action_id)
                self.assertTrue(handle.called)
                args = handle.call_args[0]
                self.assertEqual(args[1], "U123", action_id)
                self.assertEqual(args[2], "C123", action_id)
                self.assertFalse(
                    args[1].startswith("{"),
                    f"{action_id} user_id looks like str(dict)",
                )

    def test_success_replaces_buttons_but_failure_preserves_them(self):
        import calib_slack_bot as bot
        with mock.patch.dict(sys.modules, {"slack_bolt": mock.Mock(App=FakeApp)}):
            app = bot.build_app({"token": "test"})
        body = {"user": {"id": "U1"}, "channel": {"id": "C1"},
                "actions": [{"value": "7:hold"}]}
        for success in (False, True):
            with self.subTest(success=success), mock.patch.object(
                bot, "handle_calib_request", return_value=(success, "result", False),
            ), mock.patch.object(bot.runner, "warn"):
                respond = mock.Mock()
                app.actions["calib_step_proceed"](ack=lambda: None, body=body, respond=respond)
                response = respond.call_args.kwargs
                if success:
                    self.assertTrue(response["replace_original"])
                    self.assertIn("選択済み", response["text"])
                    self.assertFalse(any(b["type"] == "actions" for a in response["attachments"] for b in a["blocks"]))
                else:
                    self.assertNotIn("replace_original", response)
                    self.assertEqual(response["response_type"], "ephemeral")

    def test_distinct_lane_clicks_are_both_applied_and_duplicates_do_not_replace_cards(self):
        import tempfile
        import calib_slack_bot as bot
        from common.lib import config
        with tempfile.TemporaryDirectory() as scratch, mock.patch.object(config, "OUTPUT_DIR", Path(scratch)), \
                mock.patch.dict(sys.modules, {"slack_bolt": mock.Mock(App=FakeApp)}):
            cfg = {"token": "test", "allowed_user_ids": ["U1"], "allowed_channel_ids": ["C1"]}
            app = bot.build_app(cfg)
            def body(hold):
                return {"user": {"id": "U1"}, "channel": {"id": "C1"},
                        "actions": [{"value": f"2977:{hold}"}]}
            handler = app.actions["calib_step_proceed"]
            with mock.patch.object(bot.calib_ops, "run_dispatch", return_value=(True, "applied")) as dispatch:
                for hold in ("bcin-gate", "bcout-gate"):
                    reply = mock.Mock()
                    handler(ack=lambda: None, body=body(hold), respond=reply)
                    self.assertTrue(reply.call_args.kwargs["replace_original"])
                self.assertEqual(dispatch.call_count, 2)
                reply = mock.Mock()
                handler(ack=lambda: None, body=body("bcout-gate"), respond=reply)
                self.assertEqual(dispatch.call_count, 2)
                self.assertNotIn("replace_original", reply.call_args.kwargs)
                self.assertEqual(reply.call_args.kwargs["response_type"], "ephemeral")

    def test_failed_request_does_not_suppress_a_valid_retry(self):
        import tempfile
        import calib_slack_bot as bot
        from common.lib import config
        cfg = {"allowed_user_ids": ["U1"], "allowed_channel_ids": ["C1"]}
        with tempfile.TemporaryDirectory() as scratch, mock.patch.object(config, "OUTPUT_DIR", Path(scratch)), \
                mock.patch.object(bot.calib_ops, "run_dispatch", side_effect=[(False, "failed"), (True, "applied")]) as dispatch:
            self.assertFalse(bot.handle_calib_request(cfg, "U1", "C1", "proceed 2977:hold")[0])
            self.assertTrue(bot.handle_calib_request(cfg, "U1", "C1", "proceed 2977:hold")[0])
            self.assertEqual(dispatch.call_count, 2)


if __name__ == "__main__":
    unittest.main()
