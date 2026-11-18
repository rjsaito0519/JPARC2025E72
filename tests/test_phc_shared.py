"""Shared PHC must recover shape without mixing channel offsets or readout sides."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import uproot

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO), str(REPO / "calibration/scripts")]
from common.lib import phc_shared as shared, config
import update_phc


class SharedShapeTests(unittest.TestCase):
    def samples(self, amplitude=1.8, shift=-.2):
        x = np.linspace(.625, 1.375, 16)
        return {ch: (x, amplitude * shared.basis(x, shift) + c, np.full(len(x), 100.))
                for ch, c in [(0, -.4), (1, .7), (2, .1)]}

    def test_shape_is_not_biased_by_different_channel_offsets(self):
        fit = shared.fit_shape(self.samples())
        self.assertTrue(fit["ok"], fit)
        self.assertAlmostEqual(fit["amplitude"], 1.8, places=3)
        self.assertAlmostEqual(fit["shift"], -.2, places=3)
        self.assertAlmostEqual(fit["constants"]["0"], -.4, places=4)
        self.assertAlmostEqual(fit["constants"]["1"], .7, places=4)

    def test_one_or_no_donors_is_explicitly_unresolved(self):
        for samples in [{}, {0: self.samples()[0]}]:
            fit = shared.fit_shape(samples)
            self.assertFalse(fit["ok"])
            self.assertEqual(fit["reason"], "fewer_than_two_source_segments")

    def test_constant_energy_cannot_identify_shape(self):
        fit = shared.fit_shape({ch: (np.ones(10), np.ones(10), np.ones(10)) for ch in [0, 1]})
        self.assertFalse(fit["ok"])

    def test_robust_fit_limits_sparse_timing_outliers(self):
        samples = self.samples()
        x, y, w = samples[0]
        samples[0] = (np.r_[x, .8], np.r_[y, 15.], np.r_[w, 5.])
        fit = shared.fit_shape(samples)
        self.assertTrue(fit["ok"], fit)
        self.assertLess(abs(fit["constants"]["0"] + .4), .02)

    def test_native_type1_sign_and_pivot_conversion(self):
        fit = shared.fit_shape(self.samples())
        a, b = fit["amplitude"], fit["shift"]
        c = fit["constants"]["1"]
        x = np.linspace(.7, 1.3, 20)
        p2 = a / np.sqrt(abs(1. - b)) + c
        np.testing.assert_allclose(-a / np.sqrt(abs(x - b)) + p2, a * shared.basis(x, b) + c)


class NativeHistogramTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.fit = self.directory / "fit.root"
        self.native = self.directory / "native.root"
        tree = {"ch": np.arange(4, dtype=np.int32),
                "n_fit": np.array([[1000., 1000.], [1000., 1000.], [20., 20.], [0., 0.]]),
                "chi2_ndf": np.ones((4, 2)), "fit_status": np.zeros((4, 2), dtype=np.int32)}
        for j, value in enumerate([2., -.2, 2.]):
            tree[f"p{j}_val"] = np.full((4, 2), value)
        # One healthy side remains individual even when its other side is bad.
        tree["n_fit"][2, 1] = 1000.
        with uproot.recreate(self.fit) as root:
            root.mktree("tree", tree)
        edges_x = np.linspace(-.5, 4.5, 101)
        edges_y = np.linspace(-10., 10., 401)
        rng = np.random.default_rng(20261118)
        with uproot.recreate(self.native) as root:
            for det in ["BHT", "BH2", "CVC"]:
                for side in [0, 1]:
                    for ch in [0, 1, 2]:
                        n = 30000 if ch < 2 else 40
                        x = rng.uniform(.61, 1.39, n)
                        y = (1.6 + .5 * side) * shared.basis(x, -.2) + .3 * ch + rng.normal(0, .12, n)
                        counts, _, _ = np.histogram2d(x, y, bins=(edges_x, edges_y))
                        root[shared.histogram_name(det, ch, side, "Pi")] = (counts, edges_x, edges_y)

    def update(self, det="BHT", hand=None, source=True):
        report = {}
        rows, unused = update_phc.make_dictdata(
            self.fit, det, good_ch_range=hand,
            source_root=self.native if source else None, suffix="Pi", report=report)
        return rows, unused, report

    def test_preserves_good_individual_side_and_fills_only_bad_side(self):
        rows, unused, report = self.update()
        self.assertEqual(rows["1-0-2-1"], [1, 3, 2., -.2, 2.])
        self.assertNotEqual(rows["1-0-2-0"], [1, 3, 2., -.2, 2.])
        self.assertEqual(report["counts"]["individual"], 5)
        self.assertEqual(unused, [2, 3])
        self.assertEqual(report["counts"]["shared_prediction"], 2)

    def test_detector_and_side_adapters_do_not_mix(self):
        for det in ["BHT", "BH2", "CVC"]:
            _, _, report = self.update(det)
            self.assertTrue(report["models"]["U"]["ok"])
            self.assertTrue(report["models"]["D"]["ok"])
            self.assertGreater(report["models"]["D"]["amplitude"], report["models"]["U"]["amplitude"] + .2)

    def test_manual_range_restricts_sources_and_is_not_expanded(self):
        rows, _, report = self.update(hand=[0, 0])
        self.assertEqual(set(rows), {"1-0-0-0", "1-0-0-1"})
        self.assertEqual(report["range_origin"], "manual")
        self.assertEqual(report["models"]["U"]["sources"], [0])
        self.assertEqual(report["counts"]["unresolved"], 6)

    def test_missing_native_data_never_falls_back_to_parameter_mean(self):
        rows, _, report = self.update(source=False)
        self.assertNotIn("1-0-2-0", rows)
        self.assertEqual(report["counts"]["unresolved"], 3)

    def test_summary_and_plot_contain_actual_fallback_decisions(self):
        _, _, det = self.update()
        report = {"run": 1, "suffix": "Pi", "decision": "rows_prepared", "detectors": [det]}
        pdf = shared.render_pdf(report, self.native, self.directory / "qa.pdf", "test")
        self.assertTrue(pdf.read_bytes().startswith(b"%PDF"))
        self.assertIn("未解決", shared.summary_line(det))
        self.assertEqual(det["channels"][-1]["mode"], "shared_prediction")

    def test_slack_caption_uses_the_sealed_manifest_report(self):
        from common.lib import slack_notify
        _, _, detector = self.update()
        manifest = {"meta": {"phc_shared": {"detectors": [detector]}}}
        with mock.patch.object(slack_notify, "post_pdf", return_value=(True, "mocked")) as post:
            success, _ = slack_notify.post_hodo_qa(1, "hdphc", pdf_path=self.directory / "book.pdf", manifest=manifest)
        self.assertTrue(success)
        caption = post.call_args.args[3]
        self.assertIn("BHT:", caption)
        self.assertIn("共有", caption)
        self.assertIn("未解決", caption)

    def test_report_roundtrip_is_strict_json(self):
        _, _, det = self.update()
        with mock.patch.object(config, "OUTPUT_DIR", self.directory):
            path = shared.write_report(1, "Pi", {"detectors": [det]})
        self.assertEqual(json.loads(path.read_text())["detectors"][0]["counts"], det["counts"])


class ReferenceClockTests(unittest.TestCase):
    def test_reference_clock_is_fixed_without_changing_relative_shapes(self):
        previous = {"3-0-5-0": [1, 3, 2., -.2, 2.], "3-0-5-1": [1, 3, 2., -.2, 2.]}
        rows = {"1-0-0-0": [1, 3, 1., -.2, 1.],
                "3-0-5-0": [1, 3, 2., -.2, 1.6], "3-0-5-1": [1, 3, 2., -.2, 1.6],
                "3-0-6-0": [1, 3, 3., -.4, 2.5]}
        report = {"detector": "BH2", "channels": [],
                  "models": {"U": {"sources": [5]}, "D": {"sources": [5]}}}
        record = shared.preserve_bh2_reference(rows, previous, [report])
        self.assertEqual(record["segment"], 5)
        self.assertAlmostEqual(record["common_shift_ns"], .4)
        self.assertAlmostEqual(rows["3-0-5-0"][4], 2.)
        self.assertAlmostEqual(rows["3-0-6-0"][4], 2.9)
        self.assertEqual(rows["3-0-6-0"][2:4], [3., -.4])
        self.assertEqual(rows["1-0-0-0"], [1, 3, 1., -.2, 1.])

    def test_uncalibrated_type_zero_is_a_valid_zero_correction_reference(self):
        previous = {"3-0-5-0": [0, 0], "3-0-5-1": [0, 3, 0., 0., 0.]}
        rows = {"1-0-0-0": [1, 3, 1., -.2, 1.],
                "3-0-5-0": [1, 3, 2., -.2, 1.6], "3-0-5-1": [1, 3, 2., -.2, 1.6]}
        report = {"detector": "BH2", "channels": [],
                  "models": {"U": {"sources": [5]}, "D": {"sources": [5]}}}
        record = shared.preserve_bh2_reference(rows, previous, [report])
        self.assertEqual(record["previous_mean_correction_ns"], 0.)
        value = rows["3-0-5-0"]
        self.assertAlmostEqual(-value[2] / np.sqrt(abs(1. - value[3])) + value[4], 0.)

    def test_sparse_reference_retains_bh2_and_allows_bht_update(self):
        rows = {"1-0-0-0": [1, 3, 1., -.2, 1.], "3-0-5-0": [1, 3, 2., -.2, 1.6]}
        report = {"detector": "BH2",
                  "channels": [{"segment": 5, "side": "U", "mode": "individual", "counts": 400.}],
                  "models": {"U": {"sources": [5]}, "D": {"sources": [6]}}}
        record = shared.preserve_bh2_reference(rows, {}, [report])
        self.assertEqual(record["status"], "retained")
        self.assertEqual(set(rows), {"1-0-0-0"})
        self.assertEqual(report["counts"]["unresolved"], 1)
        self.assertEqual(report["channels"][0]["candidate_mode"], "individual")
        self.assertEqual(report["channels"][0]["reason"], "no_usable_clock_reference")
        self.assertIn("前回値保持", shared.summary_line(report))

    def test_single_detector_update_has_no_clock_adjustment(self):
        rows = {"3-0-5-0": [1, 3, 2., -.2, 1.6]}
        self.assertIsNone(shared.preserve_bh2_reference(rows, {}, []))

    def test_round_report_is_sealed_and_checked_before_restore(self):
        from common.lib import hodo_hdphc_rounds as rounds
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with mock.patch.object(config, "OUTPUT_DIR", root / "results"), mock.patch.object(config, "PARAM_DIR", root / "param"):
                parameter = rounds.hdphc_param_path(1, "Pi")
                parameter.parent.mkdir(parents=True)
                parameter.write_text("first")
                shared.write_report(1, "Pi", {"decision": "rows_prepared", "detectors": []})
                shared.mark_written(1, "Pi", parameter)
                entry = rounds.snapshot_round(1, 1, "Pi", {"ok": True, "width": .2})
                parameter.write_text("second")
                rounds.restore_round(1, 1, "Pi")
                self.assertEqual(parameter.read_text(), "first")
                sealed = rounds.round_dir(1, 1) / entry["shared_phc"]["name"]
                sealed.write_text("corrupted")
                parameter.write_text("keep")
                with self.assertRaisesRegex(ValueError, "report snapshot"):
                    rounds.restore_round(1, 1, "Pi")
                self.assertEqual(parameter.read_text(), "keep")


class SessionChoiceTests(unittest.TestCase):
    def setUp(self):
        from common.lib import phc_range_gate, calib_ops
        self.gate, self.ops = phc_range_gate, calib_ops
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.directory = Path(self.tmp.name)
        self.conf = self.directory / "phc_conf.py"
        self.conf.write_text("limits_dict = {}\n")
        for patch in [mock.patch.object(config, "OUTPUT_DIR", self.directory),
                      mock.patch.object(self.gate, "PHC_CONF_PATH", self.conf),
                      mock.patch.object(self.ops, "_phc_refit_after_conf", return_value=(True, "refit"))]:
            patch.start(); self.addCleanup(patch.stop)

    def proposal(self):
        self.gate.write_proposal(1, {"run": 1, "suffix": "Pi",
                                    "proposals": [{"det": "BHT", "first": 2, "last": 10}]})

    def test_go_uses_proposal_without_writing_manual_file(self):
        self.proposal()
        before = self.conf.read_bytes()
        success, message = self.ops.apply_phc_go(1, background_resume=False)
        self.assertTrue(success, message)
        self.assertEqual(self.gate.session_range(1, "Pi", "BHT"), [2, 10])
        self.assertIsNone(self.gate.session_range(1, "K", "BHT"))
        self.assertFalse(self.gate.is_range_confirmed(1, suffix="K"))
        self.assertEqual(self.conf.read_bytes(), before)
        self.gate.clear_range_confirmed(1)
        self.assertIsNone(self.gate.session_range(1, "Pi", "BHT"))

    def test_manual_edit_always_overrides_accepted_automatic_range(self):
        import update_param
        from types import SimpleNamespace
        self.proposal()
        self.ops.apply_phc_go(1, background_resume=False)
        self.conf.write_text('limits_dict = {"00001_Pi_bht": [4, 8]}\n')
        self.assertEqual(update_param.phc_range(SimpleNamespace(run=1, suffix="Pi"), "BHT"), (4, 8))

    def test_failed_refit_restores_prior_session_selection(self):
        self.gate.set_range_confirmed(1, suffix="Pi", ranges={"BHT": [3, 9]})
        before = self.gate.confirmed_path(1).read_bytes()
        self.proposal()
        with mock.patch.object(self.ops, "_phc_refit_after_conf", return_value=(False, "failed")):
            success, _ = self.ops.apply_phc_go(1, background_resume=False)
        self.assertFalse(success)
        self.assertEqual(self.gate.confirmed_path(1).read_bytes(), before)
        self.assertTrue(self.gate.pending_phc_ask(1))

    def test_out_of_detector_proposal_cannot_be_adopted(self):
        self.gate.write_proposal(1, {"run": 1, "suffix": "Pi",
                                    "proposals": [{"det": "BHT", "first": 0, "last": 999}]})
        with mock.patch.object(self.gate, "update_phc_ask_slack"):
            success, _ = self.ops.apply_phc_go(1, background_resume=False)
        self.assertFalse(success)
        self.assertFalse(self.gate.is_range_confirmed(1))


if __name__ == "__main__":
    unittest.main()
