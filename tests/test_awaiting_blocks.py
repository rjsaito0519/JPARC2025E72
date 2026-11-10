"""Awaiting Slack Block Kit (#25 Phase E)."""
import sys
import unittest
from pathlib import Path

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


if __name__ == "__main__":
    unittest.main()
