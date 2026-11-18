import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import awkward as ak
import numpy as np
import uproot
from common.lib import beam_selection,runner

class BeamSelectionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        requirement=self.root/"requirements.yml"
        requirement.write_text("2514:\n  Pi:\n    bac: veto\n")
        change=patch.object(beam_selection,"REQUIREMENTS",requirement);change.start();self.addCleanup(change.stop)

    def source(self, tdcs, flags=None, user="BAC_TDC 0 10", embedded=True):
        path=self.root/"input.root"
        with uproot.recreate(path) as file:
            if embedded:file["param/USER"]=user
            file["hodo"]={"beam_flag":np.asarray(flags if flags is not None else [1]*len(tdcs),dtype=np.int32),
                          "bac_tdc_u":ak.Array(tdcs)}
        return path

    def test_unlisted_run_does_not_require_new_metadata(self):
        self.assertIsNone(beam_selection.check(1,"Pi",self.root/"missing.root"))

    def test_veto_rejects_current_bac_positive_pi_tags(self):
        root=self.source([[[5]],[[6]]])
        result=beam_selection.check(2514,"Pi",root)
        self.assertFalse(result["ok"])
        self.assertEqual(2,result["violations"])

    def test_native_strict_window_boundaries_and_other_particles_are_respected(self):
        root=self.source([[[0,10]],[[5]],[[12]]],flags=[1,2,1])
        result=beam_selection.check(2514,"Pi",root)
        self.assertTrue(result["ok"],result)
        self.assertEqual(2,result["selected_events"])

    def test_no_selected_particles_cannot_pass_validation(self):
        result=beam_selection.check(2514,"Pi",self.source([[[5]]],flags=[2]))
        self.assertFalse(result["ok"])
        self.assertEqual("no_selected_particle",result["reason"])

    def test_missing_or_ambiguous_window_is_not_guessed(self):
        for text in ["other 0 10","BAC_TDC 0 10\nBAC_TDC 0 20"]:
            result=beam_selection.check(2514,"Pi",self.source([[[12]]],user=text))
            self.assertFalse(result["ok"])
            self.assertEqual("beam_selection_unavailable",result["reason"])

    def test_external_user_handles_root_without_parameter_macros(self):
        root=self.source([[[12]]],embedded=False)
        user=self.root/"USER";user.write_text("BAC_TDC 0 10\n")
        self.assertTrue(beam_selection.check(2514,"Pi",root,parameter_files={"USER":user})["ok"])

    def test_requirement_failure_stops_check_before_parameter_write(self):
        root=self.source([[[5]]])
        with patch("common.lib.calib_anomaly.remember_check_reason") as record:
            with self.assertRaises(SystemExit) as result:
                beam_selection.require(2514,"Pi",root)
        self.assertEqual(runner.EXIT_CHECK,result.exception.code)
        self.assertTrue(record.call_args.args[1].startswith("beam_selection_mismatch:"))

    def test_explicit_require_mode_uses_the_same_native_window(self):
        beam_selection.REQUIREMENTS.write_text("2514:\n  Pi:\n    bac: require\n")
        self.assertTrue(beam_selection.check(2514,"Pi",self.source([[[5]]]))["ok"])
        self.assertFalse(beam_selection.check(2514,"Pi",self.source([[[12]]]))["ok"])

if __name__=="__main__":unittest.main()
