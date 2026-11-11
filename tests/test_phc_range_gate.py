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

    def test_needs_range_ask_skips_t1(self):
        self.assertTrue(self.gate.needs_range_ask("BHT"))
        self.assertFalse(self.gate.needs_range_ask("T1"))
        self.assertEqual(
            [o["value"] for o in self.gate._det_select_options(None)],
            ["BHT", "BH2", "CVC"],
        )

    def test_clear_proposal(self):
        self.gate.write_proposal(
            21,
            {
                "run": 21,
                "suffix": "Pi",
                "proposals": [{"det": "BHT", "first": 1, "last": 10}],
            },
        )
        self.assertTrue(self.gate.pending_phc_ask(21))
        self.assertTrue(self.gate.clear_proposal(21))
        self.assertFalse(self.gate.pending_phc_ask(21))
        self.assertFalse(self.gate.clear_proposal(21))

    def test_range_confirmed_set_and_clear(self):
        self.assertFalse(self.gate.is_range_confirmed(31))
        path = self.gate.set_range_confirmed(31)
        self.assertTrue(path.is_file())
        self.assertTrue(self.gate.is_range_confirmed(31))
        self.assertTrue(self.gate.clear_range_confirmed(31))
        self.assertFalse(self.gate.is_range_confirmed(31))
        self.assertFalse(self.gate.clear_range_confirmed(31))

    def test_clear_phc_ask_state_drops_proposal_and_confirmed(self):
        self.gate.write_proposal(
            32,
            {
                "run": 32,
                "suffix": "Pi",
                "proposals": [{"det": "BHT", "first": 1, "last": 10}],
            },
        )
        self.gate.set_range_confirmed(32)
        self.gate.clear_phc_ask_state(32)
        self.assertFalse(self.gate.pending_phc_ask(32))
        self.assertFalse(self.gate.is_range_confirmed(32))

    def test_apply_phc_go_sets_confirmed(self):
        from common.lib import calib_ops

        conf = Path(self.tmp.name) / "phc_conf.py"
        conf.write_text("limits_dict = {\n}\n", encoding="utf-8")
        self.gate.write_proposal(
            33,
            {
                "run": 33,
                "suffix": "Pi",
                "proposals": [{"det": "BHT", "first": 1, "last": 10}],
            },
        )
        self.gate.upsert_phc_conf(33, "Pi", "BHT", 1, 10, path=conf)
        with mock.patch.object(
            self.gate, "missing_proposal_dets", return_value=[],
        ), mock.patch.object(
            calib_ops, "_phc_refit_after_conf", return_value=(True, "refit ok"),
        ):
            ok, msg = calib_ops.apply_phc_go(33, background_resume=False)
        self.assertTrue(ok)
        self.assertIn("phc-go", msg)
        self.assertFalse(self.gate.pending_phc_ask(33))
        self.assertTrue(self.gate.is_range_confirmed(33))

    def test_clear_rounds_keeps_confirmed_by_default(self):
        from common.lib import hodo_hdphc_rounds

        self.gate.set_range_confirmed(34)
        with mock.patch.object(hodo_hdphc_rounds.config, "OUTPUT_DIR", self.out):
            hodo_hdphc_rounds.clear_rounds(34)
        self.assertTrue(self.gate.is_range_confirmed(34))
        self.assertTrue(
            (self.out / "calib" / "run00034" / "hdphc_rounds").is_dir()
        )

    def test_clear_rounds_optional_reset_confirmed(self):
        from common.lib import hodo_hdphc_rounds

        self.gate.set_range_confirmed(34)
        with mock.patch.object(hodo_hdphc_rounds.config, "OUTPUT_DIR", self.out):
            hodo_hdphc_rounds.clear_rounds(34, reset_range_confirmed=True)
        self.assertFalse(self.gate.is_range_confirmed(34))

    def test_run_hdphc_round1_clear_journal_keeps_confirmed(self):
        """Go → refit round1 must not drop session range_confirmed (#46)."""
        scripts = str(REPO / "calibration" / "scripts")
        if scripts not in sys.path:
            sys.path.insert(0, scripts)
        import iterate_hodo_chain as chain
        from types import SimpleNamespace

        self.gate.set_range_confirmed(36)
        args = SimpleNamespace(
            run=36, kaon=True, local=True, ftof=False, force=False,
        )
        with mock.patch.object(
            chain, "stage_prepare_and_decode", return_value=Path("/dev/null"),
        ), mock.patch.object(
            chain, "run_hodo_ok",
        ), mock.patch.object(
            chain.decode_run, "run_decode",
        ), mock.patch.object(
            chain.hodo_btof_monitor, "hodo_root", return_value=Path("/dev/null"),
        ), mock.patch.object(
            chain.hodo_btof_monitor, "score_btof",
            return_value={"width": 1.0, "n_entries": 1, "ok": True},
        ), mock.patch.object(
            chain.hodo_hdphc_rounds, "snapshot_round",
            return_value={"width": 1.0, "n_entries": 1, "ok": True},
        ), mock.patch.object(
            chain.hodo_qa_pdf, "source_paths", return_value=[],
        ):
            chain.run_hdphc_one_round(args, 1, clear_journal=True)
        self.assertTrue(self.gate.is_range_confirmed(36))

    def test_clear_dispatch_state_clears_confirmed(self):
        scripts = str(REPO / "calibration" / "scripts")
        if scripts not in sys.path:
            sys.path.insert(0, scripts)
        import calib_dispatch as cd

        self.gate.write_proposal(
            35,
            {
                "run": 35,
                "suffix": "Pi",
                "proposals": [{"det": "BHT", "first": 1, "last": 5}],
            },
        )
        self.gate.set_range_confirmed(35)
        cd.clear_dispatch_state(35)
        self.assertFalse(self.gate.pending_phc_ask(35))
        self.assertFalse(self.gate.is_range_confirmed(35))

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

        self.gate.write_proposal(
            3,
            {
                "run": 3,
                "suffix": "Pi",
                "proposals": [{"det": "BHT", "first": 1, "last": 10}],
            },
        )
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

    def test_exit_check_posts_buttons_when_phc_reason_stale(self):
        from common.lib import calib_approval

        # Stale reason string, no proposal on disk → Proceed/Refit/Stop.
        self.assertFalse(self.gate.pending_phc_ask(41))
        with mock.patch(
            "common.lib.slack_notify.try_notify_failure",
            return_value=(True, "notified"),
        ) as notify:
            ok, msg = calib_approval.notify_exit_check_slack(
                41,
                {
                    "source": "exit_check",
                    "reason": (
                        "step hodo_hdphc EXIT_CHECK (rc=3): "
                        "phc_range_proposal: trial BHT,BH2,T1"
                    ),
                    "step": "hodo_hdphc",
                    "index": 3,
                    "id": "stalehold",
                },
            )
        self.assertTrue(ok)
        self.assertEqual(msg, "notified")
        notify.assert_called_once()
        _args, kwargs = notify.call_args
        blocks = kwargs.get("blocks") or []
        action_ids = [
            el.get("action_id")
            for b in blocks if b.get("type") == "actions"
            for el in b.get("elements") or []
        ]
        self.assertEqual(
            action_ids,
            ["calib_step_proceed", "calib_step_refit", "calib_step_stop"],
        )

    def test_phc_blocks_include_custom_and_go(self):
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
                "det": "BH2",
                "first": 2,
                "last": 8,
                "bad": 1,
                "total": 4,
                "current_first": None,
                "current_last": None,
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
        self.assertIn("*range* `15`–`49`", text)
        self.assertIn("_hint_ `2`–`8`", text)
        self.assertIn("_range unset_", text)
        self.assertIn(":large_green_circle:", text)
        self.assertIn(":large_orange_circle:", text)
        self.assertIn(":red_circle:", text)
        self.assertNotIn("unused", text)
        self.assertNotIn("Go で hdphc", text)
        actions = blocks[1]["elements"]
        ids = [e["action_id"] for e in actions]
        self.assertEqual(ids, ["calib_phc_custom", "calib_phc_go"])
        view = self.gate.phc_custom_modal_view(9, rows)
        self.assertEqual(view["callback_id"], "calib_phc_custom_modal")
        self.assertEqual(view["private_metadata"], "9")
        det_el = view["blocks"][0]["element"]
        self.assertEqual(det_el["type"], "static_select")
        self.assertEqual(
            [o["value"] for o in det_el["options"]], ["BHT", "BH2", "T1"],
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
