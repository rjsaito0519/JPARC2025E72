"""PHC selection must compare the same populated histograms."""
import copy
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common.lib import hodo_btof_monitor as monitor
from common.lib import hodo_hdphc_rounds as rounds
from common.lib.phc_comparison import compare_candidates


def candidate(round_i, widths=(1., 2., 3.)):
    segments = [
        {"histogram": f"BHT_seg{i}U_CBTOF_vs_DeltaE_Pi", "width": width,
         "n_entries": 100., "n_total": 110., "axes_sha": "same"}
        for i, width in enumerate(widths)]
    return {"round": round_i, "width": float(np.median(widths)), "ok": True,
            "n_entries": len(segments)*100.,
            "scoring_context": {"version": 1, "particle": "Pi", "ftof": False,
                                "required_detectors": ["BHT"]},
            "per_det": {"BHT": {"segments": segments}}}


class ComparisonTests(unittest.TestCase):
    def test_missing_wide_segment_cannot_win(self):
        initial = candidate(1)
        missing = candidate(2, (.1,))
        self.assertEqual(rounds.pick_round([initial, missing], "auto"), 1)
        self.assertFalse(compare_candidates([initial, missing])[1]["comparable"])

    def test_detector_loss_and_empty_initial_detector_require_review(self):
        initial = candidate(1)
        initial["scoring_context"]["required_detectors"].append("BH2")
        self.assertFalse(compare_candidates([initial])[0]["comparable"])
        with self.assertRaisesRegex(ValueError, "no comparable"):
            rounds.pick_round([initial], "auto")

    def test_lost_detector_cannot_make_the_remaining_narrow_detector_win(self):
        initial = candidate(1, (.1,))
        initial["scoring_context"]["required_detectors"].append("BH2")
        initial["per_det"]["BH2"] = copy.deepcopy(initial["per_det"]["BHT"])
        initial["per_det"]["BH2"]["segments"][0].update(histogram="BH2_segment", width=3.)
        lost = copy.deepcopy(initial)
        lost.update(round=2, width=.1)
        del lost["per_det"]["BH2"]
        self.assertEqual(rounds.pick_round([initial, lost], "auto"), 1)
        self.assertEqual(compare_candidates([initial, lost])[1]["missing_segments"], 1)

    def test_statistics_loss_cannot_win_even_with_same_flow_fraction(self):
        initial = candidate(1)
        fewer = candidate(2, (.1, .2, .3))
        for segment in fewer["per_det"]["BHT"]["segments"]:
            segment["n_entries"] /= 2
            segment["n_total"] /= 2
        self.assertEqual(rounds.pick_round([initial, fewer], "auto"), 1)

    def test_increased_flow_cannot_win_even_if_normal_statistics_unchanged(self):
        initial = candidate(1)
        leaky = candidate(2, (.1, .2, .3))
        leaky["per_det"]["BHT"]["segments"][0]["n_total"] = 200.
        self.assertEqual(rounds.pick_round([initial, leaky], "auto"), 1)

    def test_time_flow_increase_cannot_be_hidden_by_less_energy_flow(self):
        initial = candidate(1)
        changed = candidate(2, (.1, .2, .3))
        initial["per_det"]["BHT"]["segments"][0].update(
            n_tof_flow=0., n_de_flow=10., n_both_flow=0.)
        changed["per_det"]["BHT"]["segments"][0].update(
            n_tof_flow=10., n_de_flow=0., n_both_flow=0.)
        self.assertEqual(rounds.pick_round([initial, changed], "auto"), 1)

    def test_added_segments_do_not_change_population_score(self):
        initial = candidate(1)
        extra = candidate(2, (1., 2., 3., .01, .01, .01, .01))
        comparisons = compare_candidates([initial, extra])
        self.assertTrue(comparisons[1]["comparable"])
        self.assertEqual(comparisons[1]["width"], 2.)
        self.assertEqual(rounds.pick_round([initial, extra], "auto"), 1)

    def test_valid_improvement_drives_both_adoption_and_patience(self):
        journal = [candidate(1), candidate(2, (.5, 1., 1.5)),
                   candidate(3, (.7, 1.2, 1.7))]
        self.assertEqual(rounds.pick_round(journal, "auto"), 2)
        self.assertEqual(rounds.adaptive_progress(journal)["best_round"], 2)
        self.assertEqual(rounds.adaptive_progress(journal)["without_improvement"], 1)

    def test_rejected_candidates_do_not_claim_convergence(self):
        initial = candidate(1)
        missing = [candidate(i, (.1,)) for i in range(2, 6)]
        progress = rounds.adaptive_progress([initial] + missing)
        self.assertFalse(progress["finished"])
        self.assertTrue(progress["review_required"])
        self.assertEqual(progress["without_improvement"], 0)

    def test_legacy_journal_requires_review_and_still_allows_explicit_manual_choice(self):
        journal = [{"round": 1, "width": .1, "ok": True}]
        with self.assertRaisesRegex(ValueError, "no comparable"):
            rounds.pick_round(journal, "auto")
        self.assertEqual(rounds.pick_round(journal, "1"), 1)
        self.assertIsNone(rounds.adaptive_progress(journal)["best_width"])

    def test_particle_and_axes_change_are_not_comparable(self):
        for changed_field in ("particle", "axes"):
            with self.subTest(changed_field=changed_field):
                initial = candidate(1)
                changed = candidate(2, (.1, .2, .3))
                if changed_field == "particle":
                    changed["scoring_context"]["particle"] = "K"
                else:
                    changed["per_det"]["BHT"]["segments"][0]["axes_sha"] = "changed"
                self.assertFalse(compare_candidates([initial, changed])[1]["comparable"])

    def test_invalid_flow_statistics_require_review(self):
        for count in (None, float("nan"), 0., -1.):
            with self.subTest(count=count):
                initial = candidate(1)
                initial["per_det"]["BHT"]["segments"][0]["n_total"] = count
                self.assertFalse(compare_candidates([initial])[0]["comparable"])

    def test_approval_context_explains_comparison_failure(self):
        from common.lib import calib_approval, config
        journal = [candidate(1), candidate(2, (.1,))]
        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.object(rounds, "load_journal", return_value=journal), \
                    mock.patch.object(config, "OUTPUT_DIR", Path(directory)):
                message = calib_approval.gate_context(1, "hodo_hdphc", 1, lane="hodo")
        self.assertIn("要確認（自動採用の対象外）", message)
        self.assertIn("欠損 2", message)
        self.assertIn("第1ラウンド", message)

    def test_summary_contains_reason_and_same_segment_statistics(self):
        with tempfile.TemporaryDirectory() as directory:
            journal = [candidate(1), candidate(2, (.1,))]
            with mock.patch.object(rounds, "load_journal", return_value=journal), \
                    mock.patch.object(rounds, "rounds_dir", return_value=Path(directory)), \
                    mock.patch.object(rounds, "summary_path", return_value=Path(directory)/"summary.txt"):
                report = rounds.write_summary(1)
            self.assertIn("REVIEW", report.read_text())
            self.assertIn("reference 100.0", report.read_text())
            self.assertTrue((Path(directory)/"comparison.json").is_file())


class MonitorTests(unittest.TestCase):
    def test_root_cycles_not_double_counted_and_context_is_recorded(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"test.root"
            name = "BHT_seg0U_CBTOF_vs_DeltaE_Pi"
            with uproot.recreate(path) as root:
                for count in (1., 2.):
                    root[name] = (np.ones((2, 3))*count,
                                  np.array([0., 1., 2.]), np.array([-1., 0., 1., 2.]))
            with mock.patch.object(monitor.config, "hodo_detectors", return_value=["BHT"]):
                score = monitor.score_btof(1, root_path=path)
            self.assertEqual(score["per_det"]["BHT"]["n_seg"], 1)
            segment = score["per_det"]["BHT"]["segments"][0]
            self.assertEqual(segment["n_entries"], 12.)
            self.assertEqual(segment["n_tof_flow"], 0.)
            self.assertEqual(segment["n_de_flow"], 0.)
            self.assertTrue(segment["axes_sha"])
            self.assertTrue(compare_candidates([dict(score, round=1)])[0]["comparable"])

    def test_axis_flow_directions_are_reported_separately(self):
        values = np.zeros((4, 4))
        values[1:-1, 1:-1] = 1
        values[1, 0] = 3
        values[0, 1] = 5
        values[0, 0] = 7
        histogram = mock.Mock()
        histogram.values.side_effect = lambda flow=False: values if flow else values[1:-1, 1:-1]
        histogram.axes = [mock.Mock(), mock.Mock()]
        for axis in histogram.axes:
            axis.edges.return_value = np.array([0., 1., 2.])
        root = mock.Mock()
        root.keys.return_value = ["BHT_seg0U_CBTOF_vs_DeltaE_Pi;1"]
        root.__getitem__ = mock.Mock(return_value=histogram)
        with mock.patch.object(monitor.uproot, "open", return_value=root), \
                mock.patch.object(monitor.Path, "is_file", return_value=True), \
                mock.patch.object(monitor.config, "hodo_detectors", return_value=["BHT"]):
            score = monitor.score_btof(1)
        segment = score["per_det"]["BHT"]["segments"][0]
        self.assertEqual((segment["n_entries"], segment["n_total"]), (4, 19))
        self.assertEqual((segment["n_tof_flow"], segment["n_de_flow"], segment["n_both_flow"]), (3,5,7))
        root.close.assert_called_once()
