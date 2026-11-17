"""PHC keeps updating through worse rounds, then restores the global best."""
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "calibration/scripts"))
from common.lib import config, fit_validity, hodo_hdphc_rounds as rounds, runner
import calib_dispatch as cd
import calib_session
import iterate_hodo_chain as iterator
import run_hodo


def comparable_score(width, n_entries=100, **fields):
    return {
        "width": width, "ok": True, "n_entries": n_entries,
        "scoring_context": {"version": 1, "particle": "Pi", "ftof": False,
                            "required_detectors": ["BHT"]},
        "per_det": {"BHT": {"segments": [
            {"histogram": "BHT_seg0U_CBTOF_vs_DeltaE_Pi", "width": width,
             "n_entries": n_entries, "n_total": n_entries, "axes_sha": "same"}
        ], "width": width, "n_seg": 1, "n_entries": n_entries}},
        **fields,
    }




class AdaptivePhcTests(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.root = Path(self.scratch.name)
        for name, value in (("OUTPUT_DIR", self.root / "out"), ("PARAM_DIR", self.root / "param")):
            override = mock.patch.object(config, name, value)
            override.start()
            self.addCleanup(override.stop)

    def history(self, widths):
        return [comparable_score(width, round=i) for i, width in enumerate(widths, 1)]

    def test_waits_for_three_rounds_without_new_best(self):
        widths = [3.0, 2.0, 2.4, 2.5]
        self.assertFalse(rounds.adaptive_progress(self.history(widths))["finished"])
        summary = rounds.adaptive_progress(self.history(widths + [2.6]))
        self.assertTrue(summary["finished"])
        self.assertEqual(summary["best_round"], 2)
        self.assertEqual(summary["reason"], "no_improvement")

    def test_new_best_resets_patience(self):
        summary = rounds.adaptive_progress(self.history([3.0, 2.0, 2.2, 1.9, 2.0]))
        self.assertFalse(summary["finished"])
        self.assertEqual(summary["without_improvement"], 1)
        self.assertEqual(summary["best_round"], 4)

    def test_improving_series_stops_at_explicit_cap(self):
        summary = rounds.adaptive_progress(self.history([3.0, 2.0, 1.0]), max_rounds=3)
        self.assertTrue(summary["finished"])
        self.assertEqual(summary["reason"], "round_limit")
        self.assertEqual(summary["best_round"], 3)

    def test_worse_round_params_feed_next_round_and_best_is_restored_last(self):
        parameter = rounds.hdphc_param_path(1, "Pi")
        parameter.parent.mkdir(parents=True)
        parameter.write_text("initial")
        root_file = self.root / "run00001_Hodo.root"
        root_file.write_bytes(b"input")
        seen_inputs = []
        widths = [3.0, 2.0, 2.4, 2.5, 2.6]

        def fit(*args):
            seen_inputs.append(parameter.read_text())
            parameter.write_text(f"round{len(seen_inputs)}")

        args = SimpleNamespace(run=1, kaon=False, ftof=False, force=False,
                               local=False, hdphc_iter=10, hdphc_patience=3)
        scores = [comparable_score(width, n_entries=1000) for width in widths]
        with mock.patch.object(iterator, "stage_prepare_and_decode", return_value=self.root / "runlist.yml"), \
                mock.patch.object(iterator, "run_hodo_ok", side_effect=fit), \
                mock.patch.object(iterator.decode_run, "run_decode"), \
                mock.patch.object(iterator.hodo_btof_monitor, "hodo_root", return_value=root_file), \
                mock.patch.object(iterator.hodo_btof_monitor, "score_btof", side_effect=scores), \
                mock.patch.object(iterator.hodo_qa_pdf, "source_paths", return_value=[]):
            for index in range(1, 6):
                self.assertEqual(iterator.run_hdphc_one_round(args, index, clear_journal=index == 1), runner.EXIT_OK)
            self.assertEqual(iterator.run_hdphc_one_round(args, 6, clear_journal=False), runner.EXIT_SKIP)
        self.assertEqual(seen_inputs, ["initial", "round1", "round2", "round3", "round4"])
        self.assertEqual(parameter.read_text(), "round5")
        chosen = rounds.pick_round(rounds.load_journal(1), "auto")
        rounds.restore_round(1, chosen, "Pi")
        self.assertEqual(chosen, 2)
        self.assertEqual(parameter.read_text(), "round2")

    def test_later_round_requires_completed_predecessor_before_decode(self):
        args = SimpleNamespace(run=1, kaon=False, hdphc_iter=10, hdphc_patience=3)
        with mock.patch.object(iterator, "stage_prepare_and_decode") as decode:
            with self.assertRaises(SystemExit) as stopped:
                iterator.run_hdphc_one_round(args, 2, clear_journal=False)
            self.assertEqual(stopped.exception.code, runner.EXIT_CHECK)
            decode.assert_not_called()

    def completed_series(self):
        parameter = rounds.hdphc_param_path(1, "Pi")
        parameter.parent.mkdir(parents=True, exist_ok=True)
        for index, width in enumerate((2.0, 1.0, 1.2), 1):
            parameter.write_text(f"round{index}")
            rounds.snapshot_round(1, index, "Pi", comparable_score(width))
        return parameter

    def test_manual_finish_skips_remaining_updates_then_restores_best(self):
        parameter = self.completed_series()
        record = rounds.request_finish(1, hold_id="gate", user_id="U1")
        self.assertEqual(record["best_round"], 2)
        args = SimpleNamespace(run=1, kaon=False, hdphc_iter=10)
        with mock.patch.object(iterator, "stage_prepare_and_decode") as decode:
            self.assertEqual(iterator.run_hdphc_one_round(args, 4, clear_journal=False), runner.EXIT_SKIP)
            self.assertEqual(iterator.run_hdphc_one_round(args, 5, clear_journal=False), runner.EXIT_SKIP)
        decode.assert_not_called()
        self.assertEqual(parameter.read_text(), "round3")
        adopt_args = SimpleNamespace(run=1, kaon=False, ftof=False, force=False, local=True, adopt="1")
        with mock.patch.object(iterator.config, "ANALYZER_DIR", self.root / "analyzer"), \
                mock.patch.object(iterator, "stage_prepare_and_decode", return_value=self.root / "runlist.yml"), \
                mock.patch.object(iterator.decode_run, "run_decode"), \
                mock.patch.object(iterator, "run_hodo_ok") as qa, \
                mock.patch.object(iterator.hodo_btof_monitor, "score_btof",
                                  return_value=comparable_score(1.0)):
            iterator.run_hdphc_adopt_stage(adopt_args)
        self.assertEqual(parameter.read_text(), "round2")
        self.assertIn("--qa-only", qa.call_args.args[1])
        self.assertIn("--adopted-round 2", qa.call_args.args[1])

    def test_refit_invalidates_manual_finish_request(self):
        parameter = self.completed_series()
        rounds.request_finish(1, hold_id="gate")
        parameter.write_text("refit3")
        rounds.snapshot_round(1, 3, "Pi", comparable_score(.5))
        self.assertFalse(rounds.finish_requested(1, rounds.load_journal(1)))

    def test_finish_operation_requires_successful_current_round_and_preserves_final_approval(self):
        from common.lib import calib_ops
        self.completed_series()
        session = Path(config.OUTPUT_DIR)/"calib/run00001/calib_session.yml"
        session.write_text("run: 1\nsteps: [hodo_hdphc]\n")
        hold = cd.make_hold("step_gate", "completed", "hodo_hdphc", 2, 0)
        state = {"run":1,"status":"held","step":"hodo_hdphc","index":2,"rc":0,
                 "opts":{"hdphc_round":3},"hold":hold,"gate_waiter":False}
        with mock.patch.dict(os.environ, {"E72_CALIB_LANE":"hodo"}):
            cd.write_state(1, state)
            with mock.patch.object(calib_ops, "run_dispatch", return_value=(True,"resume started")) as resume:
                ok, message = calib_ops.apply_phc_finish(1, hold_id=hold["id"], user_id="U1")
            self.assertTrue(ok, message)
            self.assertEqual(resume.call_args.kwargs["lane"], "hodo")
            self.assertNotEqual(cd.read_state(1)["status"], "approved")
            state["rc"]=3
            cd.write_state(1,state)
            self.assertFalse(calib_ops.apply_phc_finish(1,hold_id=hold["id"])[0])

    def test_refit_invalidates_later_round_scores(self):
        path = rounds.journal_path(1)
        path.parent.mkdir(parents=True)
        rows = self.history([3.0, 2.0, 1.0]) + [comparable_score(2.5, round=2)]
        path.write_text("".join(json.dumps(row) + "\n" for row in rows))
        journal = rounds.load_journal(1)
        self.assertEqual([row["round"] for row in journal], [1, 2])
        self.assertEqual(rounds.pick_round(journal, "auto"), 2)

    def test_improve_only_failure_is_usable_but_primary_failures_are_rejected(self):
        self.assertTrue(fit_validity.phc_fit_status_ok(4000))
        for status in (4, 4070, 4140, 4210, -1):
            self.assertFalse(fit_validity.phc_fit_status_ok(status))
        self.assertFalse(fit_validity.all_finite([float("nan"), 1.0, 2.0]))

    def test_default_plan_exposes_cap_and_allows_completed_rounds_to_skip(self):
        steps = calib_session.expand_coarse_steps(["hodo"])
        phc = [step["hodo_hdphc"] for step in steps if isinstance(step, dict) and "hodo_hdphc" in step]
        self.assertEqual(len(phc), 10)
        self.assertTrue(all(step["allow_skip"] for step in phc))
        cfg = {"run": 1, "kaon": False, "force": False, "max_iter": 10, "steps": steps}
        command = cd.step_cmd("hodo_hdphc", phc[0], cfg, False)
        self.assertIn("--hdphc-iter 10", command)

    def test_only_explicitly_allowed_completion_skips_count_as_success(self):
        accepted = [{"index": 0, "rc": 4, "opts": {"allow_skip": True}},
                    {"index": 1, "rc": 0, "opts": {}}]
        self.assertEqual(cd.overall_rc(accepted), runner.EXIT_OK)
        self.assertEqual(cd.overall_rc([{"index": 0, "rc": 4, "opts": {}}]), runner.EXIT_SKIP)

    def test_no_usable_fit_is_review_not_an_adaptive_skip(self):
        args = SimpleNamespace(kaon=False, ftof=False, force=False, debug=False,
                               qa_only=False, skip_qa=True)
        from common.lib import fit_overrides
        with mock.patch.dict(os.environ, {"E72_CALIB_LANE": ""}), \
                mock.patch.object(runner, "require_bins", return_value=[]), \
                mock.patch.object(runner, "input_root", return_value=self.root / "input.root"), \
                mock.patch.object(run_hodo.qa_session, "begin", return_value={"meta": {}}), \
                mock.patch.object(run_hodo.qa_session, "attach_params"), \
                mock.patch.object(fit_overrides, "write_flat_mip_left", return_value=self.root / "override"), \
                mock.patch.object(runner, "update_param", return_value=runner.EXIT_SKIP), \
                mock.patch.object(run_hodo, "after_qa_merge", side_effect=lambda run, mode, args, rc, **kwargs: rc):
            self.assertEqual(run_hodo.fit_and_update(1, "Pi", args, "hdphc", "PHC"), runner.EXIT_CHECK)

    def test_completed_external_lane_finalizes_without_restarting_any_work(self):
        cfg = {"run": 1, "steps": ["hodo_hdprm", "dc_t0"], "kaon": False,
               "force": False, "local": False, "max_iter": 10, "auto": False}
        with mock.patch.dict(os.environ, {"E72_CALIB_LANE": ""}):
            cd.write_state(1, {"status": "stopped", "mode": "lane_coordinator",
                               "lanes": {"hodo": "stopped", "bcout": "done"}}, lane="")
            for lane in ("hodo", "bcout"):
                cd.write_state(1, {"status": "awaiting_final_approval", "finished": True}, lane=lane)
            with mock.patch.object(cd, "run_lane_coordinator") as spawn, \
                    mock.patch.object(cd, "mark_awaiting_final") as finalize:
                with self.assertRaises(SystemExit) as completed:
                    cd.do_resume(cfg, False)
        self.assertEqual(completed.exception.code, runner.EXIT_OK)
        spawn.assert_not_called()
        self.assertEqual(finalize.call_args.args[1]["status"], "awaiting_final_approval")
        self.assertTrue(finalize.call_args.args[1]["finished"])


if __name__ == "__main__":
    unittest.main()
