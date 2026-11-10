"""FakeApp regression for Slack button handlers (#42)."""
import sys
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "calibration" / "scripts"

ACTION_IDS = (
    "calib_approve_run",
    "calib_revoke_run",
    "calib_step_proceed",
    "calib_step_refit",
    "calib_step_stop",
    "calib_phc_accept",
)


class FakeApp:
    def __init__(self, token=None):
        self.token = token
        self.actions = {}
        self.commands = {}

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

        self.assertEqual(set(app.actions), set(ACTION_IDS))

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
            for action_id in ACTION_IDS:
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


if __name__ == "__main__":
    unittest.main()
