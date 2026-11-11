"""Approval, delivery and state races in parallel calibration lanes (no Slack I/O)."""
import json
import os
import sys
import tempfile
import threading
import time
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "calibration/scripts"))

from common.lib import calib_anomaly, calib_ops, config, runner, slack_notify
from common.lib.calib_state_io import write_json
import calib
import calib_dispatch as cd


class ParallelSlackTests(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.output = Path(self.scratch.name)
        output_patch = mock.patch.object(config, "OUTPUT_DIR", self.output)
        output_patch.start()
        self.addCleanup(output_patch.stop)
        env_patch = mock.patch.dict(os.environ, {"E72_CALIB_LANE": ""})
        env_patch.start()
        self.addCleanup(env_patch.stop)

    def lane_hold(self, waiter):
        root = self.output / "calib/run02977"
        root.mkdir(parents=True)
        (root / "calib_session.yml").write_text("run: 2977\nsteps: [hodo_hdphc]\n")
        write_json(root / "calib_dispatch_state.json", {
            "status": "running", "mode": "lane_coordinator",
        })
        lane = root / "lanes/hodo/calib_dispatch_state.json"
        write_json(lane, {
            "status": "held", "gate_waiter": waiter,
            "hold": {"id": "phc-hold", "source": "exit_check",
                     "step": "hodo_hdphc", "lane": "hodo", "index": 2, "rc": 3},
        })
        return lane

    def test_phc_go_targets_live_lane_without_spawning_resume(self):
        lane = self.lane_hold(True)
        with mock.patch.object(calib_ops.subprocess, "Popen") as spawn:
            success, message = calib_ops._phc_refit_after_conf(
                2977, background_resume=True, timeout_s=30, user_id="reviewer",
            )
        self.assertTrue(success, message)
        spawn.assert_not_called()
        state = json.loads(lane.read_text())
        self.assertEqual(state["hold"]["decision"], "refit")
        self.assertEqual(state["hold"]["user_id"], "reviewer")
        journal = (lane.parents[2] / "hold_journal.jsonl").read_text()
        self.assertIn('"hold_id": "phc-hold"', journal)
        self.assertIn('"user_id": "reviewer"', journal)

    def test_retry_without_waiter_resumes_only_selected_lane(self):
        self.lane_hold(False)
        with mock.patch.object(
            calib_ops.subprocess, "Popen", return_value=mock.Mock(pid=42),
        ) as spawn:
            success, message = calib_ops.run_dispatch(
                "retry", 2977, hold_id="phc-hold", user_id="reviewer",
            )
        self.assertTrue(success, message)
        command = spawn.call_args.args[0]
        self.assertEqual(command[-2:], ["--lane", "hodo"])

    def test_cli_continue_leaves_work_to_live_waiter(self):
        self.lane_hold(True)
        with mock.patch.object(cd, "do_resume") as resume:
            calib._resume_held_lane({"run": 2977}, "hodo", False)
        resume.assert_not_called()

    def test_successful_step_refit_resumes_same_index(self):
        state = {"index": 1, "rc": 0, "hold": {"decision": "refit"}}
        self.assertEqual(cd.resume_start_index(state), 1)
        state["hold"]["decision"] = "continue"
        self.assertEqual(cd.resume_start_index(state), 2)

    def test_gate_is_published_with_waiter_before_buttons(self):
        with mock.patch.dict(os.environ, {"E72_CALIB_LANE": "hodo"}):
            hold = cd.make_hold("step_gate", "review", "hodo_hdprm", 0, 0)
            state = cd.enter_held(2977, {}, hold, "hold")
        self.assertTrue(state["gate_waiter"])
        self.assertEqual(state["gate_waiter_pid"], os.getpid())

    def test_parallel_first_posts_share_one_parent_thread(self):
        start = threading.Barrier(2)
        replies = []
        calls = []

        def api(*args, **kwargs):
            calls.append(args)
            time.sleep(0.03)
            return True, {"ts": str(len(calls))}, "ok"

        def ensure():
            start.wait(timeout=5)
            replies.append(slack_notify.ensure_run_thread(2977, "QA"))

        with mock.patch.object(slack_notify, "load_slack_config", return_value=(
            {"token": "test", "channel_id": "channel"}, "ok",
        )), mock.patch.object(slack_notify, "_api", side_effect=api):
            workers = [threading.Thread(target=ensure) for _ in range(2)]
            for worker in workers:
                worker.start()
            for worker in workers:
                worker.join(timeout=5)
                self.assertFalse(worker.is_alive())
        self.assertEqual(len(calls), 1)
        self.assertEqual([reply[2] for reply in replies], ["1", "1"])

    def test_parallel_pdf_and_review_preserve_both_receipts(self):
        start = threading.Barrier(2)
        replies = []
        pdf = self.output / "qa.pdf"
        pdf.write_bytes(b"%PDF test")

        def post():
            start.wait(timeout=5)
            replies.append(slack_notify.post_pdf(2977, "t0", pdf, "QA"))

        def notify():
            start.wait(timeout=5)
            replies.append(slack_notify.notify_failure(2977, "gate", "Review"))

        with mock.patch.object(slack_notify, "load_slack_config", return_value=(
            {"token": "test", "channel_id": "channel"}, "ok",
        )), mock.patch.object(slack_notify, "_api", return_value=(True, {"ts": "1"}, "ok")), \
                mock.patch.object(slack_notify, "_upload_pdf", return_value=(True, "uploaded")):
            workers = [threading.Thread(target=post), threading.Thread(target=notify)]
            for worker in workers:
                worker.start()
            for worker in workers:
                worker.join(timeout=5)
                self.assertFalse(worker.is_alive())
        self.assertTrue(all(reply[0] for reply in replies))
        state = slack_notify.load_state(2977)
        self.assertIn("t0|Pi|qa.pdf", state["posts"])
        self.assertIn("gate", state["failures"])

    def test_atomic_publication_never_exposes_empty_json(self):
        path = self.output / "state.json"
        write_json(path, {"index": 0})
        finished = threading.Event()
        errors = []

        def read():
            while not finished.is_set():
                try:
                    json.loads(path.read_text())
                except Exception as error:
                    errors.append(error)

        reader = threading.Thread(target=read)
        reader.start()
        try:
            for index in range(100):
                write_json(path, {"index": index})
        finally:
            finished.set()
            reader.join(timeout=5)
        self.assertEqual(errors, [])

    def test_stale_automated_input_holds_without_terminal_prompt(self):
        root = self.output / "run02977_BcOut.root"
        root.touch()
        with mock.patch.dict(os.environ, {"E72_CALIB_DISPATCH": "1"}), \
                mock.patch.object(runner, "stale_params", return_value=["HDPRM: changed"]), \
                mock.patch.object(runner, "check_input") as interactive_check:
            with self.assertRaises(SystemExit) as stopped:
                runner.checked_root(root, "Pi", False)
        self.assertEqual(stopped.exception.code, runner.EXIT_CHECK)
        interactive_check.assert_not_called()
        self.assertIn("stale_input", calib_anomaly.read_check_reason(2977))

    def test_lane_check_reasons_remain_separate(self):
        for lane in ("hodo", "bcin"):
            with mock.patch.dict(os.environ, {"E72_CALIB_LANE": lane}):
                calib_anomaly.remember_check_reason(2977, lane)
        for lane in ("hodo", "bcin"):
            with mock.patch.dict(os.environ, {"E72_CALIB_LANE": lane}):
                self.assertEqual(calib_anomaly.read_check_reason(2977), lane)

    def test_drift_graph_count_reads_uproot_classname_property(self):
        root = self.output / "drift.root"
        root.touch()
        graph_file = mock.MagicMock()
        graph_file.keys.return_value = ["graph", "hist"]
        graph_file.__getitem__.side_effect = lambda key: SimpleNamespace(
            classname="TGraph" if key == "graph" else "TH1D",
        )
        with mock.patch("uproot.open") as open_root:
            open_root.return_value.__enter__.return_value = graph_file
            self.assertEqual(calib_anomaly.count_drift_graphs(root), 1)


if __name__ == "__main__":
    unittest.main()
