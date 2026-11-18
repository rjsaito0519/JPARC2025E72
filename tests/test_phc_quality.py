"""Quality ranking must keep the same hits while tolerating bounded bin migration."""
import copy
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common.lib import phc_quality as quality, phc_comparison, hodo_hdphc_rounds

EDGES=np.linspace(-10.,10.,101).tolist()

def segment(name="BHT_seg0U_CBTOF_vs_DeltaE_Pi", sigma=.3, count=8000, center=0.):
    x=np.random.default_rng(17).normal(center,sigma,count)
    normal,_=np.histogram(x,EDGES)
    counts=[0.,*normal.astype(float).tolist(),0.]
    metrics=quality.profile_metrics(counts,EDGES)
    return {"histogram":name,"width":metrics["half68"],"n_entries":float(count),
            "n_total":float(count),"n_tof_flow":0.,"n_de_flow":0.,"n_both_flow":0.,
            "axes_sha":"same","tof_counts":counts}

def score(rows=None):
    rows=rows or [segment()]
    total=sum(row["n_total"] for row in rows)
    value={"round":1,"ok":True,"width":.3,"n_entries":total,
           "scoring_context":{"version":2,"particle":"Pi","ftof":False,
                              "required_detectors":["BHT"],"primary_detectors":["BHT"],
                              "quality_policy":quality.policy()},
           "per_det":{"BHT":{"segments":rows,"tof_edges":EDGES,
                            "time_validity":{"n_cmt":int(total),"n_events":int(total),
                                             "raw_hit_signature":"same-raw-hits","nonfinite_keys":[]}}}}
    value["input_integrity"]=copy.deepcopy(quality.integrity_from_score(value))
    return value

def candidate(reference):
    value=copy.deepcopy(reference);value["round"]=2
    return value

class QualityTests(unittest.TestCase):
    def compare(self, first, second):
        return phc_comparison.compare_candidates([first,second])[-1]

    def move_flow(self, value, count):
        row=value["per_det"]["BHT"]["segments"][0]
        index=int(np.argmax(row["tof_counts"]))
        row["tof_counts"][index]-=count;row["tof_counts"][-1]+=count
        row["n_entries"]-=count;row["n_tof_flow"]+=count

    def test_first_update_cannot_make_a_bad_center_its_own_reference(self):
        first=score()
        first["per_det"]["BHT"]["segments"]=[segment(center=.3)]
        result=phc_comparison.compare_candidates([first])[0]
        self.assertFalse(result["comparable"])
        self.assertTrue(any("before first update" in reason and "center" in reason
                            for reason in result["reasons"]))

    def test_first_update_cannot_make_new_tails_its_own_reference(self):
        first=score()
        row=first["per_det"]["BHT"]["segments"][0]
        row["tof_counts"][50]-=20
        row["tof_counts"][65]+=20
        result=phc_comparison.compare_candidates([first])[0]
        self.assertFalse(result["comparable"])
        self.assertTrue(any("before first update" in reason and "tail" in reason
                            for reason in result["reasons"]))

    def test_one_of_eight_thousand_axis_migrations_is_allowed(self):
        first=score();second=candidate(first);self.move_flow(second,1)
        result=self.compare(first,second)
        self.assertTrue(result["comparable"],result["reasons"])
        self.assertEqual(result["allowed_flow_changes"][0]["delta"],1)

    def test_large_flow_increase_is_rejected_even_if_width_is_small(self):
        first=score();second=candidate(first);self.move_flow(second,20)
        result=self.compare(first,second)
        self.assertFalse(result["comparable"])
        self.assertTrue(any("allowance" in text for text in result["reasons"]))

    def test_fill_loss_cannot_be_hidden_as_a_small_axis_change(self):
        first=score();second=candidate(first)
        row=second["per_det"]["BHT"]["segments"][0]
        row["n_total"]-=1;row["n_entries"]-=1;row["tof_counts"][50]-=1
        result=self.compare(first,second)
        self.assertFalse(result["comparable"])
        self.assertTrue(any("fill total" in text for text in result["reasons"]))

    def test_adc_population_change_is_rejected(self):
        first=score();second=candidate(first)
        second["per_det"]["BHT"]["segments"][0]["n_de_flow"]=1
        result=self.compare(first,second)
        self.assertFalse(result["comparable"])
        self.assertTrue(any("ADC-domain" in text for text in result["reasons"]))

    def test_new_nonfinite_hit_cannot_be_cancelled_by_recovering_another(self):
        first=score()
        first["per_det"]["BHT"]["time_validity"]["nonfinite_keys"]=[[1,0]]
        first["input_integrity"]=quality.integrity_from_score(first)
        second=candidate(first);second["per_det"]["BHT"]["time_validity"]["nonfinite_keys"]=[[2,0]]
        result=self.compare(first,second)
        self.assertFalse(result["comparable"])
        self.assertTrue(any("non-finite" in text for text in result["reasons"]))

    def test_changed_raw_times_reject_an_apparently_better_profile(self):
        first=score();second=candidate(first)
        second["per_det"]["BHT"]["time_validity"]["raw_hit_signature"]="changed"
        result=self.compare(first,second)
        self.assertFalse(result["comparable"])
        self.assertTrue(any("uncorrected timing" in text for text in result["reasons"]))

    def test_missing_fixed_segment_is_rejected(self):
        first=score();second=candidate(first);second["per_det"]["BHT"]["segments"]=[]
        result=self.compare(first,second)
        self.assertFalse(result["comparable"])
        self.assertTrue(any("missing/empty" in text for text in result["reasons"]))

    def test_narrow_peak_farther_from_zero_is_rejected(self):
        first=score();second=candidate(first)
        second["per_det"]["BHT"]["segments"]=[segment(sigma=.15,center=.5)]
        result=self.compare(first,second)
        self.assertFalse(result["comparable"])
        self.assertTrue(any("center moved" in text for text in result["reasons"]))

    def test_narrow_core_with_more_tail_is_rejected(self):
        first=score();second=candidate(first)
        row=segment(sigma=.2);middle=int(np.argmax(row["tof_counts"]))
        row["tof_counts"][middle]-=40;row["tof_counts"][62]+=40
        second["per_det"]["BHT"]["segments"]=[row]
        result=self.compare(first,second)
        self.assertFalse(result["comparable"])
        self.assertTrue(any("tail fraction" in text for text in result["reasons"]))

    def test_high_statistics_btof_controls_rank_instead_of_sparse_channel_width(self):
        main=segment()
        sparse=segment(name="BHT_seg1U_CBTOF_vs_DeltaE_Pi",count=3,sigma=.001)
        first=score([main,sparse]);second=candidate(first)
        new_main=segment(sigma=.2)
        sparse=copy.deepcopy(sparse);sparse["tof_counts"]=[0.]*102
        for location in [40,50,60]:sparse["tof_counts"][location]=1.
        sparse["width"]=2.
        second["per_det"]["BHT"]["segments"]=[new_main,sparse]
        result=self.compare(first,second)
        self.assertTrue(result["comparable"],result["reasons"])
        self.assertEqual(hodo_hdphc_rounds.pick_round([first,second],"auto"),2)

    def test_missing_initial_integrity_is_not_auto_approved(self):
        first=score();first.pop("input_integrity")
        result=phc_comparison.compare_candidates([first])[0]
        self.assertFalse(result["comparable"])
        self.assertTrue(any("initial decode" in text for text in result["reasons"]))

    def test_policy_cannot_change_mid_series(self):
        first=score();second=candidate(first)
        second["scoring_context"]["quality_policy"]["flow_fraction"]=.01
        self.assertFalse(self.compare(first,second)["comparable"])

    def test_secondary_detector_without_data_is_unmeasured_and_keeps_integrity_checks(self):
        first=score()
        first["scoring_context"]["required_detectors"].append("T1")
        empty=segment(name="T1_seg0U_CTOF_vs_DeltaE_Pi",count=0)
        empty["width"]=None
        first["per_det"]["T1"]={"segments":[empty],"tof_edges":EDGES,
            "time_validity":{"n_cmt":0,"n_events":0,"raw_hit_signature":"empty","nonfinite_keys":[]}}
        first["input_integrity"]=quality.integrity_from_score(first)
        second=candidate(first)
        result=self.compare(first,second)
        self.assertTrue(result["comparable"],result["reasons"])
        self.assertIn("no_evaluable",result["quality"]["T1"]["status"])
        second["per_det"]["T1"]["time_validity"]["raw_hit_signature"]="changed"
        self.assertFalse(self.compare(first,second)["comparable"])

    def test_empty_profile_does_not_break_global_aggregation(self):
        first=score()
        empty=segment(name="BHT_seg2U_CBTOF_vs_DeltaE_Pi",count=0)
        first["per_det"]["BHT"]["segments"].append(empty)
        self.assertGreater(quality.aggregate(first["per_det"]["BHT"])["n_total"],0)

class MonitorInputTests(unittest.TestCase):
    def test_real_root_records_initial_raw_hits_and_validates_adoption_profiles(self):
        import tempfile
        import awkward as ak
        import uproot
        from common.lib import hodo_btof_monitor as monitor
        values=np.random.default_rng(72).normal(0,.3,2000)
        edges=np.linspace(-10,10,101)
        bins,_,_=np.histogram2d(np.ones(len(values)),values,bins=(np.linspace(-.5,4.5,101),edges))
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"run00001_Hodo.root"
            with uproot.recreate(path) as root:
                for prefix,tag in [("BHT","CBTOF"),("T0","CBTOF"),("T1","CTOF")]:
                    root[f"{prefix}_seg0U_{tag}_vs_DeltaE_Pi"]=(bins,np.linspace(-.5,4.5,101),edges)
                fields={"event_number":np.arange(len(values),dtype=np.uint32),"beam_flag":np.ones(len(values),dtype=np.int32)}
                for detector in ["bht","bh2","t1"]:
                    fields[detector+"_hit_seg"]=ak.Array([[0.] for _ in values])
                    times=np.zeros(len(values)) if detector=="bh2" else -values
                    fields[detector+"_mt"]=ak.Array([[[float(v)]] for v in times])
                    fields[detector+"_cmt"]=ak.Array([[[float(v)]] for v in times])
                root["hodo"]=fields
            measured=monitor.score_btof(1,root_path=path)
            self.assertEqual(measured["scoring_context"]["version"],2)
            self.assertEqual(measured["per_det"]["BHT"]["time_validity"]["n_cmt"],2000)
            self.assertEqual(len(measured["per_det"]["BHT"]["time_validity"]["raw_hit_signature"]),64)
            measured["round"]=1
            measured["input_integrity"]=monitor.initial_integrity(1,root_path=path)
            result=phc_comparison.compare_candidates([measured])[0]
            self.assertTrue(result["comparable"],result["reasons"])
            fresh=monitor.score_btof(1,root_path=path)
            self.assertEqual(hodo_hdphc_rounds.score_differences(measured,fresh),[])
            fresh["per_det"]["BHT"]["segments"][0]["tof_counts"][50]+=1
            self.assertTrue(any("TOF profile" in text for text in hodo_hdphc_rounds.score_differences(measured,fresh)))


if __name__=="__main__":
    unittest.main()
