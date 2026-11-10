"""Awaiting Slack Block Kit (#25 Phase E)."""
import sys
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]


class AwaitingBlocksTests(unittest.TestCase):
    def setUp(self):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import calib_approval

        self.approval = calib_approval

    def test_blocks_have_approve_revoke(self):
        blocks = self.approval.awaiting_action_blocks(42, "abcdefghijklmnop")
        self.assertEqual(len(blocks), 2)
        actions = blocks[1]["elements"]
        ids = {e["action_id"] for e in actions}
        self.assertEqual(ids, {"calib_approve_run", "calib_revoke_run"})
        self.assertEqual(actions[0]["value"], "42:abcdefghijklmnop")
        self.assertEqual(actions[1]["value"], "42")

    def test_step_gate_blocks_have_three_actions(self):
        blocks = self.approval.step_gate_action_blocks(
            7, "hodo", 0, hold_id="abc12def",
        )
        self.assertEqual(len(blocks), 2)
        actions = blocks[1]["elements"]
        ids = [e["action_id"] for e in actions]
        self.assertEqual(
            ids,
            ["calib_step_proceed", "calib_step_refit", "calib_step_stop"],
        )
        self.assertEqual({e["value"] for e in actions}, {"7:abc12def"})

    def test_notify_step_gate_key_includes_hold_id(self):
        from common.lib import slack_notify

        keys = []

        def fake_try_notify(run, key, text, *, blocks=None):
            keys.append(key)
            return True, "ok"

        with mock.patch.object(
            slack_notify, "try_notify_failure", side_effect=fake_try_notify,
        ):
            ok, _msg = self.approval.notify_step_gate_slack(
                7, "hodo_hdprm", 0, hold_id="deadbeef",
            )
        self.assertTrue(ok)
        self.assertEqual(keys, ["step_gate:0:hodo_hdprm:deadbeef"])


if __name__ == "__main__":
    unittest.main()
