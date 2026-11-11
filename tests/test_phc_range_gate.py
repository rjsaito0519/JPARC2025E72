"""PHC range proposal helpers."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]


class PhcRangeGateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "out"
        self.out.mkdir()
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import config, phc_range_gate

        self.gate = phc_range_gate
        self.cfg = config
        self._old = config.OUTPUT_DIR
        config.OUTPUT_DIR = self.out

    def tearDown(self):
        self.cfg.OUTPUT_DIR = self._old

    def test_propose_from_good(self):
        self.assertEqual(
            self.gate.propose_from_good([1, 2, 3, 4, 5], [1, 5]),
            (2, 4),
        )
        self.assertIsNone(self.gate.propose_from_good([1, 2], [1, 2]))

    def test_upsert_phc_conf(self):
        conf = Path(self.tmp.name) / "phc_conf.py"
        conf.write_text(
            'limits_dict = {\n    "00001_Pi_bht" : [10, 20],\n}\n',
            encoding="utf-8",
        )
        key = self.gate.upsert_phc_conf(
            1, "Pi", "BHT", 12, 18, path=conf,
        )
        self.assertEqual(key, "00001_Pi_bht")
        text = conf.read_text(encoding="utf-8")
        self.assertIn('"00001_Pi_bht" : [12, 18]', text)
        key2 = self.gate.upsert_phc_conf(
            2, "Pi", "BH2", 2, 9, path=conf,
        )
        self.assertEqual(key2, "00002_Pi_bh2")
        text = conf.read_text(encoding="utf-8")
        self.assertIn('"00002_Pi_bh2" : [ 2,  9]', text)
        keys = self.gate.conf_keys_on_disk(path=conf)
        self.assertIn("00001_Pi_bht", keys)
        self.assertIn("00002_Pi_bh2", keys)

    def test_write_and_read_proposal(self):
        payload = {
            "run": 7,
            "suffix": "Pi",
            "proposals": [{"det": "BHT", "first": 15, "last": 49}],
        }
        path = self.gate.write_proposal(7, payload)
        self.assertTrue(path.is_file())
        got = self.gate.read_proposal(7)
        self.assertEqual(got["suffix"], "Pi")
        self.assertEqual(got["proposals"][0]["first"], 15)

    def test_notify_held_skips_phc_and_step_gate(self):
        from common.lib import calib_approval

        with mock.patch(
            "common.lib.slack_notify.try_notify_failure",
            return_value=(True, "notified"),
        ) as notify:
            ok, msg = calib_approval.notify_held_slack(
                3, {"source": "step_gate", "reason": "ok", "step": "hodo"},
            )
            self.assertTrue(ok)
            self.assertIn("step_gate", msg)
            notify.assert_not_called()

            ok, msg = calib_approval.notify_held_slack(
                3,
                {
                    "source": "exit_check",
                    "reason": "step hodo EXIT_CHECK (rc=3): other anomaly",
                    "step": "hodo",
                    "index": 0,
                },
            )
            self.assertTrue(ok)
            self.assertIn("exit_check", msg)
            notify.assert_not_called()

            ok, msg = calib_approval.notify_exit_check_slack(
                3,
                {
                    "source": "exit_check",
                    "reason": "step hodo EXIT_CHECK (rc=3): phc_range_proposal: BHT",
                    "step": "hodo",
                    "index": 0,
                    "id": "abc",
                },
            )
            self.assertTrue(ok)
            self.assertIn("phc_range_proposal", msg)
            notify.assert_not_called()

            ok, msg = calib_approval.notify_exit_check_slack(
                3,
                {
                    "source": "exit_check",
                    "reason": "step hodo EXIT_CHECK (rc=3): other anomaly",
                    "step": "hodo",
                    "index": 0,
                    "id": "def",
                },
            )
            self.assertTrue(ok)
            self.assertEqual(msg, "notified")
            notify.assert_called_once()

    def test_phc_blocks_include_custom(self):
        rows = [
            {
                "det": "BHT",
                "first": 1,
                "last": 10,
                "bad": 2,
                "total": 5,
                "current_first": 15,
                "current_last": 49,
                "anomaly": True,
            },
            {
                "det": "T1",
                "first": None,
                "last": None,
                "bad": None,
                "total": None,
                "current_first": None,
                "current_last": None,
                "anomaly": False,
            },
        ]
        blocks = self.gate.phc_range_action_blocks(9, rows)
        text = blocks[0]["text"]["text"]
        self.assertIn("current [15, 49]", text)
        self.assertIn("suggest [1, 10]", text)
        self.assertIn("suggest (none)", text)
        actions = blocks[1]["elements"]
        ids = [e["action_id"] for e in actions]
        self.assertEqual(ids, ["calib_phc_accept", "calib_phc_custom"])
        view = self.gate.phc_custom_modal_view(9, rows)
        self.assertEqual(view["callback_id"], "calib_phc_custom_modal")
        self.assertEqual(view["private_metadata"], "9")
        det_el = view["blocks"][0]["element"]
        self.assertEqual(det_el["type"], "static_select")
        self.assertEqual(
            [o["value"] for o in det_el["options"]], ["BHT", "T1"],
        )
        self.assertEqual(det_el["initial_option"]["value"], "BHT")
        self.assertEqual(
            view["blocks"][1]["element"]["placeholder"]["text"], "1",
        )

    def test_pending_and_missing_ask_list(self):
        self.gate.write_proposal(
            11,
            {
                "run": 11,
                "suffix": "Pi",
                "proposals": [
                    {"det": "BHT", "first": 15, "last": 49},
                    {"det": "BH2", "first": None, "last": None},
                    {"det": "T1", "first": 0, "last": 0},
                ],
            },
        )
        self.assertTrue(self.gate.pending_phc_ask(11))
        conf = Path(self.tmp.name) / "phc_conf.py"
        conf.write_text("limits_dict = {\n}\n", encoding="utf-8")
        rows = self.gate.read_proposal(11)["proposals"]
        missing = self.gate.missing_proposal_dets(11, "Pi", rows, path=conf)
        self.assertEqual(missing, ["BHT", "BH2", "T1"])
        self.gate.upsert_phc_conf(11, "Pi", "BHT", 15, 49, path=conf)
        self.gate.upsert_phc_conf(11, "Pi", "T1", 0, 0, path=conf)
        missing2 = self.gate.missing_proposal_dets(11, "Pi", rows, path=conf)
        self.assertEqual(missing2, ["BH2"])


if __name__ == "__main__":
    unittest.main()
