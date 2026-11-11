"""Review notices explain why a PDF is absent and identify the held lane."""
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.lib import calib_approval, slack_notify


class ExitCheckNoticeTests(unittest.TestCase):
    def test_stale_input_notice_says_fit_and_pdf_are_missing(self):
        hold = {
            "id": "current-hold", "step": "dc_t0", "index": 0, "lane": "bcin",
            "reason": "step dc_t0 EXIT_CHECK: stale_input: decode again; HDPRM changed",
        }
        with mock.patch.object(slack_notify, "try_notify_failure", return_value=(True, "posted")) as notify:
            calib_approval.notify_exit_check_slack(2977, hold)
        text = notify.call_args.args[2]
        self.assertIn("[bcin]", text)
        self.assertIn("fit 未実行", text)
        self.assertIn("未生成", text)
        self.assertIn("Refit", text)
        blocks = notify.call_args.kwargs["blocks"]
        self.assertIn("stale_input", blocks[0]["text"]["text"])
        buttons = blocks[1]["elements"]
        self.assertNotIn("calib_step_proceed", [button["action_id"] for button in buttons])
        self.assertEqual(buttons[0]["value"], "2977:current-hold")

    def test_delivery_failure_is_distinguished_from_missing_fit(self):
        detail = calib_approval.exit_check_detail("qa_delivery_failed: network down")
        self.assertIn("network down", detail)
        self.assertIn("更新済み", detail)
        self.assertNotIn("fit 未実行", detail)


if __name__ == "__main__":
    unittest.main()
