"""Unattended PHC selection retains fit validation and never fabricates approval."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "calibration/scripts"))
import update_param as update
from common.lib import calib_auto, channel_flags, config, phc_range_gate


class AutomaticCalibrationTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.output = Path(self.directory.name)
        self.patch(mock.patch.object(config, "OUTPUT_DIR", self.output))
        self.patch(mock.patch.dict("os.environ", {"E72_CALIB_DISPATCH": "1", "E72_CALIB_AUTO": "1"}))
        self.patch(mock.patch.object(update.phc_conf, "limits_dict", {}))
        self.patch(mock.patch.object(update, "fresh_outputs", return_value={"BHT": self.output / "fit.root"}))
        self.patch(mock.patch.object(update.calib_anomaly, "decide_ratio", return_value="unused segment ratio high"))
        self.flags = self.patch(mock.patch.object(update.calib_anomaly, "dump_candidates"))
        self.values = {"1-0-2-0": [0., 1., 2.]}
        self.make = self.patch(mock.patch.object(update.update_phc, "make_dictdata", return_value=(self.values, [0, 1, 3])))
        fake = mock.MagicMock()
        fake.__getitem__.return_value.arrays.return_value = {"ch": np.array([0, 1, 2, 3])}
        self.patch(mock.patch("uproot.open", return_value=fake))
        self.args = SimpleNamespace(run=1, suffix="Pi", ftof=False)

    def patch(self, patcher):
        self.addCleanup(patcher.stop)
        return patcher.start()

    def collect(self):
        return update.collect_phc(self.args, ["BHT"], [])

    def test_auto_selects_valid_span_and_leaves_no_confirmation(self):
        result, states = self.collect()
        self.assertEqual(result, self.values)
        self.assertEqual(states, [])
        self.assertEqual(self.make.call_args.kwargs["good_ch_range"], [2, 2])
        self.assertFalse(phc_range_gate.pending_phc_ask(1))
        self.assertFalse(phc_range_gate.is_range_confirmed(1))
        journal = self.output / "calib/run00001/auto_decisions.jsonl"
        record = json.loads(journal.read_text().splitlines()[-1])
        self.assertEqual(record["source"], "auto")
        self.assertEqual(record["selection"], "valid_fit_span")

    def test_auto_prioritizes_configured_hand_range(self):
        update.phc_conf.limits_dict["00001_Pi_bht"] = [1, 3]
        result, states = self.collect()
        self.assertEqual(result, self.values)
        self.assertFalse(states)
        self.make.assert_called_once()
        self.assertFalse(phc_range_gate.pending_phc_ask(1))
        self.assertFalse(phc_range_gate.is_range_confirmed(1))

    def test_no_valid_fit_stops_instead_of_adopting(self):
        self.make.return_value = ({}, [0, 1, 2, 3])
        with self.assertRaises(SystemExit) as stopped:
            self.collect()
        self.assertEqual(stopped.exception.code, 1)
        self.assertFalse(phc_range_gate.pending_phc_ask(1))

    def test_nonfinite_selected_fit_is_rejected(self):
        self.make.side_effect = [(self.values, [0, 1, 3]), ({"1-0-2-0": [float("nan")]}, [])]
        with self.assertRaises(SystemExit):
            self.collect()

    def test_manual_mode_keeps_the_existing_range_confirmation(self):
        update.phc_conf.limits_dict["00001_Pi_bht"] = [1, 3]
        with mock.patch.dict("os.environ", {"E72_CALIB_AUTO": "0"}):
            result, states = self.collect()
        self.assertEqual(result, {})
        self.assertTrue(states)
        self.assertTrue(phc_range_gate.pending_phc_ask(1))

    def test_auto_skip_warning_does_not_create_a_write_block(self):
        with mock.patch.object(update.calib_anomaly, "kvc_skip_counts", return_value=(3, 8, ["seg5", "seg6", "seg7"])):
            self.assertIsNone(update._hdprm_skip_anomaly(self.args, "KVC", self.output / "fit.root"))
        self.flags.assert_not_called()
        self.assertFalse(channel_flags.open_needs_review(1))

    def test_restart_retires_only_automatic_findings(self):
        for channel, source in (("old-auto", "auto"), ("human-note", "human")):
            channel_flags.append_flag(1, family="hodo", mode="hdprm", det="KVC", channel=channel, source=source)
        self.assertEqual(calib_auto.retire_previous_auto_flags(1), 1)
        self.assertEqual([row["channel"] for row in channel_flags.open_needs_review(1)], ["human-note"])
