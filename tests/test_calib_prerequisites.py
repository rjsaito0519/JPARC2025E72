"""Unexecuted fits and unsent QA cannot be approved past a step boundary."""
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "calibration/scripts"))
from common.lib import calib_anomaly, calib_approval, calib_ops, config, qa_session, runner
from common.lib.calib_state_io import write_json
import calib_dispatch as cd


class PrerequisiteTests(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.output = Path(self.scratch.name)
        output_patch = mock.patch.object(config, "OUTPUT_DIR", self.output)
        output_patch.start()
        self.addCleanup(output_patch.stop)
        env_patch = mock.patch.dict(os.environ, {"E72_CALIB_LANE": "bcin"})
        env_patch.start()
        self.addCleanup(env_patch.stop)

    def held(self, reason):
        hold = cd.make_hold("exit_check", reason, "dc_t0", 0, 3)
        cd.write_state(2977, {"run": 2977, "status": "held", "hold": hold})
        return hold

    def test_stale_fit_cannot_be_skipped_by_cli(self):
        hold = self.held("stale_input: HDPRM changed")
        with self.assertRaises(SystemExit) as rejected:
            cd.apply_decide(2977, "continue", hold_id=hold["id"])
        self.assertEqual(rejected.exception.code, runner.EXIT_FAIL)
        self.assertEqual(cd.read_state(2977)["status"], "held")
        cd.apply_decide(2977, "refit", hold_id=hold["id"])
        self.assertEqual(cd.read_state(2977)["hold"]["decision"], "refit")

    def test_stale_fit_cannot_be_skipped_by_slack(self):
        hold = self.held("stale_input: HDPRM changed")
        session = self.output / "calib/run02977/calib_session.yml"
        session.write_text("run: 2977\nsteps: [dc_t0]\n")
        success, message = calib_ops.run_dispatch("proceed", 2977, hold_id=hold["id"])
        self.assertFalse(success)
        self.assertIn("fit 未実行", message)

    def sealed_qa(self):
        exec_id = "20261111T000000Z_test"
        directory = qa_session.manifest_dir(2977)
        pdf = directory / exec_id / "qa.pdf"
        pdf.parent.mkdir(parents=True)
        pdf.write_bytes(b"%PDF current QA")
        digest = qa_session.file_sha256(pdf)
        write_json(directory / f"{exec_id}.json", {
            "run": 2977, "exec_id": exec_id,
            "outputs": [{"path": str(pdf), "name": "qa.pdf", "sha": digest}],
        })
        return exec_id, pdf, digest

    def test_resend_of_exact_manifest_unblocks_delivery_failure(self):
        exec_id, pdf, digest = self.sealed_qa()
        hold = self.held(f"qa_delivery_failed: timeout; qa_exec_id={exec_id}")
        self.assertTrue(calib_approval.continue_block_reason(2977, hold))
        write_json(self.output / "calib/run02977/slack_thread.json", {"posts": {
            "qa": {"exec_id": "older", "pdf": pdf.name, "pdf_sha": digest},
        }})
        self.assertTrue(calib_approval.continue_block_reason(2977, hold))
        write_json(self.output / "calib/run02977/slack_thread.json", {"posts": {
            "qa": {"exec_id": exec_id, "pdf": pdf.name, "pdf_sha": digest},
        }})
        self.assertEqual(calib_approval.continue_block_reason(2977, hold), "")
        cd.apply_decide(2977, "continue", hold_id=hold["id"])
        self.assertEqual(cd.read_state(2977)["status"], "stopped")

    def test_missing_or_modified_pdf_never_counts_as_delivered(self):
        exec_id, pdf, digest = self.sealed_qa()
        hold = self.held(f"qa_delivery_failed: timeout; qa_exec_id={exec_id}")
        write_json(self.output / "calib/run02977/slack_thread.json", {"posts": {
            "qa": {"exec_id": exec_id, "pdf": pdf.name, "pdf_sha": digest},
        }})
        pdf.write_bytes(b"different QA")
        self.assertTrue(calib_approval.continue_block_reason(2977, hold))
        pdf.unlink()
        self.assertTrue(calib_approval.continue_block_reason(2977, hold))

    def test_unsealed_qa_failure_requires_refit(self):
        hold = self.held("qa_delivery_failed: missing PDF")
        with self.assertRaises(SystemExit):
            cd.apply_decide(2977, "continue", hold_id=hold["id"])

    def simulate_step(self, failures):
        cfg = {"run": 2977, "steps": ["dc_t0", "dc_drift"], "kaon": False,
               "force": False, "local": False, "max_iter": 10, "auto": False}
        attempts = []

        def child(command, **kwargs):
            attempts.append(command)
            if len(attempts) <= failures:
                calib_anomaly.remember_check_reason(2977, "stale_input: changed during decode")
                return runner.EXIT_CHECK, ""
            self.assertEqual(calib_anomaly.read_check_reason(2977), "")
            return runner.EXIT_OK, ""

        with mock.patch.object(cd, "step_cmd", side_effect=lambda name, *args: name), \
                mock.patch.object(runner, "sh", side_effect=child), \
                mock.patch.object(cd, "_ensure_slack_bot_for_gates"), \
                mock.patch.object(calib_approval, "notify_step_gate_slack", return_value=(True, "posted")), \
                mock.patch.object(calib_approval, "notify_exit_check_slack", return_value=(True, "posted")), \
                mock.patch.object(cd, "wait_held_decision", return_value="stop"):
            with self.assertRaises(SystemExit):
                cd.run_steps(cfg, False, 0, [], "test")
        return attempts, cd.read_state(2977)

    def test_stale_input_retries_same_stage_then_waits_for_approval(self):
        attempts, state = self.simulate_step(failures=1)
        self.assertEqual(attempts, ["dc_t0", "dc_t0"])
        self.assertEqual(state["status"], "held")
        self.assertEqual(state["hold"]["source"], "step_gate")
        self.assertEqual(state["step"], "dc_t0")

    def test_stale_retries_are_bounded_and_do_not_start_next_stage(self):
        attempts, state = self.simulate_step(failures=10)
        self.assertEqual(attempts, ["dc_t0"] * 3)
        self.assertEqual(state["hold"]["source"], "exit_check")
        self.assertIn("stale_input", state["hold"]["reason"])


if __name__ == "__main__":
    unittest.main()
