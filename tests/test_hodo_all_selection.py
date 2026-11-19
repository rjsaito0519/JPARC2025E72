"""All histogram input remains separate from parameter particle and full-population QA."""
import copy
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import awkward as ak
import numpy as np
import uproot
sys.path[:0]=[str(Path(__file__).resolve().parents[1]),str(Path(__file__).resolve().parent)]
from common.lib import hodo_selection as selection,hodo_btof_monitor as monitor,phc_quality,phc_comparison,slack_notify
import test_phc_core_adoption as core
import test_phc_quality as fixtures

class AllSelectionTests(unittest.TestCase):
    def setUp(self):
        directory=tempfile.TemporaryDirectory();self.addCleanup(directory.cleanup)
        self.directory=Path(directory.name);self.config=self.directory/"selection.yml"
        self.config.write_text("1:\n  Pi:\n    histograms: All\n    peak_window_ns: [-1.5, 1.5]\n")
        change=patch.object(selection,"SELECTIONS",self.config);change.start();self.addCleanup(change.stop)

    def test_only_declared_run_and_particle_changes(self):
        self.assertEqual("BH2_ADC_seg0U",selection.histogram_name("BH2_ADC_seg0U",1,"Pi"))
        self.assertEqual("BH2_ADC_seg0U_Pi",selection.histogram_name("BH2_ADC_seg0U",2,"Pi"))
        self.assertEqual("BH2_ADC_seg0U_K",selection.histogram_name("BH2_ADC_seg0U",1,"K"))

    def test_environment_cannot_leak_into_normal_runs(self):
        base={"E72_CALIB_HIST_SELECTION":"All","E72_CALIB_PHC_TIME_WINDOW":"-1 1","E72_CALIB_T0_TIME_WINDOW":"0 3","KEEP":"yes"}
        self.assertEqual({"KEEP":"yes"},selection.environment(2,"Pi",base))
        self.assertEqual("All",selection.environment(1,"Pi",base)["E72_CALIB_HIST_SELECTION"])

    def test_invalid_selection_is_not_silently_ignored(self):
        for text in ("1: []","1:\n  Pi:\n    histograms: Pi\n    peak_window_ns: [-1,1]",
                     "1:\n  Pi:\n    histograms: All\n    peak_window_ns: [1,2]",
                     "1:\n  Pi:\n    histograms: All\n    peak_window_ns: [true,2]"):
            self.config.write_text(text)
            with self.subTest(text=text),self.assertRaises(ValueError):selection.selection(1,"Pi")

    def root(self):
        path=self.directory/"data.root";row=fixtures.segment()
        with uproot.recreate(path) as root:
            root["param/HDPRM"]="same";root["param/USER"]="same"
            root["BHT_seg0U_CBTOF_vs_DeltaE"]=(np.asarray(row["tof_counts"][1:-1])[None,:],np.array([0.,3.]),np.array(fixtures.EDGES))
            root[row["histogram"]]=(np.ones((1,100)),np.array([0.,3.]),np.array(fixtures.EDGES))
            root["hodo"]={"event_number":np.array([1,2,3],dtype=np.uint32),"beam_flag":np.array([1,2,3],dtype=np.int32),
                          "bht_hit_seg":ak.Array([[0],[0],[0]]),"bht_mt":ak.Array([[[0.]],[[0.]],[[0.]]]),
                          "bht_cmt":ak.Array([[[0.]],[[0.]],[[0.]]]),"bht_de_u":ak.Array([[1.],[1.],[1.]])}
        return path

    def test_all_does_not_fall_back_to_pi_and_includes_unknown_events(self):
        with patch.object(monitor.config,"hodo_detectors",return_value=["BHT"]):
            result=monitor.initial_candidate(1,root_path=self.root())
        self.assertTrue(result["ok"])
        self.assertEqual("BHT_seg0U_CBTOF_vs_DeltaE",result["per_det"]["BHT"]["segments"][0]["histogram"])
        self.assertEqual(3,result["per_det"]["BHT"]["time_validity"]["n_events"])
        self.assertEqual("All",result["scoring_context"]["input_selection"]["histograms"])
        self.assertEqual("Pi",result["scoring_context"]["particle"])

    def test_missing_all_histograms_are_not_replaced_by_pi(self):
        self.assertEqual([],monitor._pick_keys(["BHT_seg0U_CBTOF_vs_DeltaE_Pi"],"BHT","CBTOF","All"))

    def test_slack_caption_declares_input_and_output_separately(self):
        with patch.object(slack_notify,"post_pdf",return_value=(True,"posted")) as post:
            slack_notify.post_hodo_qa(1,"hdprm",pdf_path=self.directory/"qa.pdf",
                                     manifest={"meta":selection.metadata(1,"Pi")})
        caption=post.call_args.args[3]
        self.assertIn("入力: All",caption);self.assertIn("出力param: Pi",caption)

    def test_raw_fit_window_is_detector_specific_and_recorded(self):
        self.config.write_text("1:\n  Pi:\n    histograms: All\n    peak_window_ns: [-1.5, 1.5]\n"
                               "    fit_time_windows_ns:\n      BHT: [-1.5, 1.5]\n      BH2: [-4, -1]\n")
        self.assertEqual("-4.0 -1.0",selection.environment(1,"Pi",detector="BH2")["E72_CALIB_PHC_TIME_WINDOW"])
        self.assertEqual([-4.,-1.],selection.metadata(1,"Pi")["input_fit_time_windows_ns"]["BH2"])

    def test_t0_window_is_declared_and_validated(self):
        original = self.config.read_text()
        self.config.write_text(original+"    t0_fit_window_ns: [0, 3]\n")
        env=selection.environment(1,"Pi")
        self.assertEqual("0.0 3.0",env["E72_CALIB_T0_TIME_WINDOW"])
        self.assertEqual([0.,3.],selection.metadata(1,"Pi")["input_t0_fit_window_ns"])
        for bounds in ["[1, 1]", "[0, .nan]", "[true, 3]"]:
            self.config.write_text(original+"    t0_fit_window_ns: "+bounds+"\n")
            with self.subTest(bounds=bounds),self.assertRaises(ValueError):selection.selection(1)

    def test_all_phc_requires_declared_raw_fit_windows(self):
        with self.assertRaisesRegex(ValueError,"declared raw-time fit windows"):
            selection.require_sources(1,"Pi","hdphc",self.root())

    def test_stale_pi_fit_output_is_rejected_even_when_its_name_matches(self):
        path=self.directory/"run00001_T0_Offset_Pi.root"
        with uproot.recreate(path) as output:
            output["other"]="old Pi result"
        with self.assertRaisesRegex(ValueError,"input selection differs"):
            selection.require_fit(path,1,"Pi")
        with uproot.update(path) as output:
            output["calib_input_selection"]="All"
            output["calib_parameter_particle"]="Pi"
        selection.require_fit(path,1,"Pi")

    def test_missing_all_source_stops_before_fitting(self):
        path=self.root()
        with self.assertRaisesRegex(ValueError,"All histograms missing"):
            selection.require_sources(1,"Pi","t0",path)
        selection.require_sources(2,"Pi","t0",path)

    def test_web_uses_all_source_but_pi_override_namespace(self):
        import test_mip_web
        web=test_mip_web.web
        with patch.object(web.config,"WORK_DIR",self.directory),patch.object(web.config,"OUTPUT_DIR",self.directory):
            app=web.App(1,"Pi","hdprm")
            self.assertEqual("BH2_ADC_seg0U",app.histogram_name(0,"U"))
            self.assertEqual("Pi",app.particle)
            t0=web.App(1,"Pi","t0")
            self.assertEqual("T0_seg0_TimeOffset",t0.histogram_name(0,"U"))

class AllPeakQualityTests(unittest.TestCase):
    def sample(self,sigma=.3):
        value=core.score(sigma=sigma)
        row=value["per_det"]["BHT"]["segments"][0]
        background=fixtures.segment(sigma=.2,center=-3.,count=3000)
        row["tof_counts"]=[a+b for a,b in zip(row["tof_counts"],background["tof_counts"])]
        row["n_total"]+=3000;row["n_entries"]+=3000
        value["per_det"]["BHT"]["peak_window"]=[-1.5,1.5]
        value["scoring_context"]["input_selection"]={"histograms":"All","peak_window_ns":[-1.5,1.5]}
        value["input_integrity"]=copy.deepcopy(phc_quality.integrity_from_score(value))
        return value

    def test_main_peak_width_keeps_full_counts_and_background(self):
        value=self.sample()
        metrics=phc_quality.aggregate(value["per_det"]["BHT"])
        self.assertEqual(11000,metrics["n_total"]);self.assertGreater(metrics["peak_entries"],7900)
        self.assertLess(metrics["half68"],.4);self.assertGreater(metrics["tail_fraction"],.25)

    def test_narrowing_the_peak_is_a_valid_update(self):
        initial=self.sample();new=self.sample(.24);new["round"]=1
        result=phc_comparison.compare_candidates([initial,new])[-1]
        self.assertTrue(result["comparable"],result["reasons"])

    def test_narrower_appearance_cannot_hide_peak_population_loss(self):
        initial=self.sample();new=copy.deepcopy(initial);new["round"]=1
        row=new["per_det"]["BHT"]["segments"][0]
        index=int(np.argmax(row["tof_counts"]));row["tof_counts"][index]-=300;row["tof_counts"][30]+=300
        result=phc_comparison.compare_candidates([initial,new])[-1]
        self.assertFalse(result["comparable"])
        self.assertTrue(any("main-peak population loss" in text for text in result["reasons"]))

    def test_changing_peak_window_invalidates_comparison(self):
        initial=self.sample();new=copy.deepcopy(initial);new["round"]=1
        new["per_det"]["BHT"]["peak_window"]=[-1.,1.]
        result=phc_comparison.compare_candidates([initial,new])[-1]
        self.assertFalse(result["comparable"])

if __name__=="__main__":unittest.main()
