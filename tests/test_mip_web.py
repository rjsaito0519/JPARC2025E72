"""MIP Web persistence, review and replay boundaries, without ROOT or Slack writes."""
import copy
import importlib.util
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
SPEC = importlib.util.spec_from_file_location("mip_web", REPO / "calibration/mip_web/mip_web.py")
web = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(web)


class MipWebTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        root = Path(self.directory.name)
        for name in ("OUTPUT_DIR", "WORK_DIR", "PARAM_DIR", "DATA_DIR"):
            patcher = patch.object(web.config, name, root / name)
            patcher.start()
            self.addCleanup(patcher.stop)
        self.app = web.App(12)
        self.source = "input-v1"
        self.app.fingerprint = lambda: self.source
        self.app.gate = lambda: (_ for _ in ()).throw(ValueError("no gate"))
        self.settings = {"parameters": [12000., 250., 30., 4.], "range": [180., 340.],
                         "function": "auto", "range_mode": "auto"}
        self.result = {"success": True, "status": 0, "covariance": 3, "ndf": 50,
                       "parameters": [13000., 255., 28., 3.], "range": [200., 320.],
                       "function": "gaussian", "curve": [[200., 2.], [255., 100.]], "origin": "test"}
        self.app.execute_preview = lambda *args, **kwargs: copy.deepcopy(self.result)
        self.app.automatic = lambda *args: {**copy.deepcopy(self.result), **copy.deepcopy(self.settings)}
        self.request = {"fingerprint": self.source, "key": "bh2-3-u", "settings": self.settings}

    def staged(self):
        fitted = self.app.fit(self.request)
        return self.app.stage({**self.request, "receipt": fitted["receipt"]})

    def test_fit_alone_does_not_stage_or_write(self):
        self.app.fit(self.request)
        self.assertFalse(self.app.state()["staged"])
        self.assertFalse(web.fit_overrides.overrides_path(12).exists())

    def test_success_receipt_required(self):
        with self.assertRaisesRegex(ValueError, "成功"):
            self.app.stage({**self.request, "receipt": "forged"})
        self.result["success"] = False
        fitted = self.app.fit(self.request)
        self.assertNotIn("receipt", fitted)
        self.assertFalse(self.app.state()["staged"])

    def test_restart_restores_stage(self):
        original = self.staged()
        restored = web.App(12)
        restored.fingerprint = lambda: self.source
        restored.gate = self.app.gate
        self.assertEqual(original["staged"], restored.state()["staged"])

    def test_source_change_discards_stage_and_receipts(self):
        self.staged()
        self.source = "input-v2"
        current = self.app.state()
        self.assertFalse(current["staged"])
        self.assertIn("破棄", current["notice"])
        with self.assertRaisesRegex(ValueError, "変わり"):
            self.app.fit(self.request)

    def test_failed_refit_preserves_previous_stage(self):
        original = self.staged()["staged"]
        self.result["success"] = False
        self.app.fit(self.request)
        self.assertEqual(original, self.app.state()["staged"])

    def test_remove_affects_only_one_staged_channel(self):
        self.staged()
        second = {**self.request, "key": "bh2-4-d"}
        result = self.app.fit(second)
        self.app.stage({**second, "receipt": result["receipt"]})
        current = self.app.stage({**self.request, "remove": True})
        self.assertEqual(["bh2-4-d"], list(current["staged"]))
        self.assertFalse(web.fit_overrides.overrides_path(12).exists())

    def test_review_contains_only_changed_fields(self):
        self.staged()
        self.result["parameters"] = list(self.settings["parameters"])
        rows = self.app.review(self.request)["rows"]
        labels = [change[0] for change in rows[0]["changes"]]
        self.assertNotIn("function", labels)
        self.assertNotIn("range_mode", labels)
        self.assertIn("p0", labels)
        self.assertIn("下限", labels)

    def install_gate(self):
        state_path = self.app.base / "lanes/hodo/calib_dispatch_state.json"
        hold = {"id": "hold-1", "step": "hodo_hdprm", "rc": 0}
        state = {"status": "held", "hold": hold}
        web.write_json(state_path, state)
        self.app.gate = lambda: (state_path, state, hold)
        return state_path, hold

    def test_write_changes_only_staged_channel_and_particle(self):
        self.staged()
        self.install_gate()
        old = {"fit": {"T1": {"hdprm": {"t1-0-u": {"mip_left": 123}}},
                       "BH2": {"hdprm": {"bh2-8-d": {"mip_left": 150},
                                          "bh2-3-u": {"mip_fit": {"K": self.settings}}}}},
               "calib_locks": {"BH2": {"hdprm": ["existing-lock"]}}}
        web.fit_overrides.save(12, old)
        current = self.app.state()
        with patch("common.lib.calib_ops.run_dispatch", return_value=(True, "accepted")) as replay:
            result = self.app.write({**self.request, "revision": current["revision"], "hold_id": "hold-1"})
        replay.assert_called_once_with("retry", 12, hold_id="hold-1", background_resume=True)
        self.assertEqual("running", result["job"]["status"])
        updated = web.fit_overrides.load(12)
        self.assertEqual(old["fit"]["T1"], updated["fit"]["T1"])
        self.assertEqual(old["calib_locks"], updated["calib_locks"])
        self.assertEqual(old["fit"]["BH2"]["hdprm"]["bh2-8-d"], updated["fit"]["BH2"]["hdprm"]["bh2-8-d"])
        self.assertEqual(self.settings, updated["fit"]["BH2"]["hdprm"]["bh2-3-u"]["mip_fit"]["K"])
        self.assertIn("Pi", updated["fit"]["BH2"]["hdprm"]["bh2-3-u"]["mip_fit"])

    def test_stale_review_and_gate_rejected_before_write(self):
        current = self.staged()
        self.install_gate()
        for revision, hold in ((current["revision"]-1, "hold-1"), (current["revision"], "wrong")):
            with self.assertRaises(ValueError):
                self.app.write({**self.request, "revision": revision, "hold_id": hold})
        self.assertFalse(web.fit_overrides.overrides_path(12).exists())

    def test_replay_start_failure_keeps_written_settings(self):
        current = self.staged()
        self.install_gate()
        with patch("common.lib.calib_ops.run_dispatch", return_value=(False, "failed")):
            result = self.app.write({**self.request, "revision": current["revision"], "hold_id": "hold-1"})
        self.assertEqual("failed", result["job"]["status"])
        self.assertTrue(web.fit_overrides.overrides_path(12).is_file())
        self.assertIn("Pi", web.fit_overrides.load(12)["fit"]["BH2"]["hdprm"]["bh2-3-u"]["mip_fit"])

    def test_running_job_blocks_edit_and_duplicate_write(self):
        current = self.staged()
        path, _ = self.install_gate()
        doc = self.app.document()
        doc["job"] = {"status": "running", "state_path": str(path), "hold_id": "hold-1"}
        web.write_json(self.app.staging_path, doc)
        with self.assertRaisesRegex(ValueError, "再実行中"):
            self.app.fit(self.request)
        with self.assertRaisesRegex(ValueError, "再実行中"):
            self.app.write({**self.request, "revision": current["revision"]})

    def test_fresh_hold_finishes_job_and_clears_staging(self):
        self.staged()
        path, _ = self.install_gate()
        doc = self.app.document()
        doc["job"] = {"status": "running", "state_path": str(path), "hold_id": "hold-1"}
        web.write_json(self.app.staging_path, doc)
        web.write_json(path, {"status": "held", "hold": {"id": "hold-2", "step": "hodo_hdprm", "rc": 0}})
        result = self.app.state()
        self.assertEqual("done", result["job"]["status"])
        self.assertFalse(result["staged"])

    def test_invalid_values_never_reach_root(self):
        for changes in ({"parameters": [1., 2., float("nan"), 0.]}, {"range": [300., 100.]},
                        {"function": "shell"}, {"parameters": [1., 2., 0., 0.]}):
            with self.assertRaises(ValueError):
                self.app.fit({**self.request, "settings": {**self.settings, **changes}})

    def saved_fit(self, *, status=4000, success=True):
        self.app.root.parent.mkdir(parents=True, exist_ok=True)
        self.app.root.write_bytes(b"input fixture")
        output = Path(web.config.OUTPUT_DIR) / "root/run00012/run00012_BH2_HDPRM_Pi.root"
        output.parent.mkdir(parents=True)
        with web.uproot.recreate(output) as root_file:
            tree = root_file.mktree("mip_fit", {"seg": "int32", "side": "int32", "parameters": "4 * float64",
                "range": "2 * float64", "function": "string", "success": "bool", "status": "int32",
                "covariance": "int32", "ndf": "int32"})
            tree.extend({"seg": [3], "side": [0], "parameters": [[12300., 251., 32., 2.]],
                "range": [[190., 330.]], "function": ["landau"], "success": [success], "status": [status],
                "covariance": [3], "ndf": [45]})
        return output

    def test_saved_root_function_metadata_is_used(self):
        self.saved_fit()
        self.app.execute_preview = lambda *args, **kwargs: {"curve": self.result["curve"]}
        restored = web.App.automatic(self.app, "bh2-3-u", self.source)
        self.assertEqual("landau", restored["function"])
        self.assertEqual([12300., 251., 32., 2.], restored["parameters"])
        self.assertEqual("前回実行のFit結果", restored["origin"])

    def test_later_decode_does_not_replace_saved_fit(self):
        output = self.saved_fit()
        os.utime(output, (1000, 1000))
        self.app.execute_preview = lambda *args, **kwargs: {"curve": self.result["curve"]} if kwargs.get("curve_only") else self.fail("implicit refit")
        restored = web.App.automatic(self.app, "bh2-3-u", self.source)
        self.assertEqual([12300., 251., 32., 2.], restored["parameters"])

    def test_unfitted_estimate_is_not_drawn_as_fit(self):
        self.saved_fit(status=-1, success=False)
        self.app.execute_preview = lambda *args, **kwargs: self.fail("estimate must not become a curve")
        restored = web.App.automatic(self.app, "bh2-3-u", self.source)
        self.assertEqual([], restored["curve"])
        self.assertIn("Fit未実施", restored["origin"])

    def test_missing_saved_result_requires_hdprm_without_refit(self):
        self.app.execute_preview = lambda *args, **kwargs: self.fail("implicit refit")
        with self.assertRaisesRegex(ValueError, "HDPRM"):
            web.App.automatic(self.app, "bh2-3-u", self.source)

    def test_overview_and_refit_use_histogram_saved_with_result(self):
        output = self.saved_fit()
        counts = web.np.arange(8, dtype=float)
        edges = web.np.arange(9, dtype=float)
        with web.uproot.update(output) as root_file:
            for segment in range(web.config.NUM_OF_CH["BH2"]):
                for side in ("U", "D"):
                    root_file[f"BH2_ADC_seg{segment}{side}_Pi"] = counts, edges
        self.assertEqual(output, self.app.fit_input())
        result = self.app.overview({"fingerprint": self.source, "page": 0})
        self.assertEqual([1., 5., 9., 13.], result["channels"][0]["counts"])

    def test_legacy_result_without_histograms_requires_refit(self):
        self.saved_fit()
        with self.assertRaisesRegex(ValueError, "元のADC"):
            self.app.fit_input()

    def test_statistics_mode_stages_and_replays_mean_method(self):
        self.request["settings"] = {**self.settings, "function": "mean", "range_mode": "fixed"}
        self.result.update(function="mean", covariance=-1, ndf=0,
                           statistics={"valid": True, "mean": 255., "median": 253.})
        staged = self.staged()
        self.assertEqual("mean", staged["staged"]["bh2-3-u"]["settings"]["function"])
        self.install_gate()
        with patch("common.lib.calib_ops.run_dispatch", return_value=(True, "accepted")):
            self.app.write({**self.request, "revision": staged["revision"], "hold_id": "hold-1"})
        line = web.fit_overrides.flat_mip_fit_path(12).read_text()
        self.assertIn("bh2-3-u:Pi mean fixed", line)

    def test_empty_statistics_cannot_be_staged(self):
        self.result.update(function="mean", statistics={"valid": False})
        response = self.app.fit({**self.request, "settings": {**self.settings, "function": "mean"}})
        self.assertFalse(response["success"])
        self.assertNotIn("receipt", response)

    def test_mean_initial_method_is_preserved(self):
        result = {**self.result, "function": "mean"}
        self.assertEqual("mean", self.app.initial_settings(result)["function"])

    def test_unacknowledged_replay_does_not_wait_forever(self):
        self.staged()
        path, _ = self.install_gate()
        doc = self.app.document()
        doc["job"] = {"status": "running", "state_path": str(path), "hold_id": "hold-1", "started_at": 1.}
        web.write_json(self.app.staging_path, doc)
        result = self.app.state()
        self.assertEqual("failed", result["job"]["status"])
        self.assertIn("開始を確認できません", result["job"]["message"])


if __name__ == "__main__":
    unittest.main()
