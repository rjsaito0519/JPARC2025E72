"""Core peak adoption keeps a valid initial correction without hiding data loss."""
import copy
import json
import unittest
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import awkward as ak
import uproot

import test_phc_quality as quality_fixtures
import test_phc_adaptive as adaptive
from common.lib import phc_quality, phc_comparison, hodo_hdphc_rounds as rounds, hodo_btof_monitor as monitor, slack_notify

def score(round_index=0, sigma=.3, center=0.):
    value=quality_fixtures.score([quality_fixtures.segment(sigma=sigma,center=center)])
    value["round"]=round_index
    value["scoring_context"].update(version=3,quality_policy=phc_quality.policy(version=3),
                                    fixed_parameters={"HDPRM":"same","USER":"same"})
    value["input_integrity"]=copy.deepcopy(phc_quality.integrity_from_score(value))
    return value

def update(initial, *, sigma=.3, center=0., round_index=1):
    value=copy.deepcopy(initial);value["round"]=round_index
    value["per_det"]["BHT"]["segments"]=[quality_fixtures.segment(sigma=sigma,center=center)]
    return value

class CoreAdoptionTests(unittest.TestCase):
    setUp=adaptive.AdaptivePhcTests.setUp

    def test_small_tail_change_is_not_data_loss_and_is_reported(self):
        initial=score();new=update(initial)
        row=new["per_det"]["BHT"]["segments"][0]
        row["tof_counts"][50]-=20;row["tof_counts"][65]+=20
        result=phc_comparison.compare_candidates([initial,new])[-1]
        self.assertTrue(result["comparable"],result["reasons"])
        self.assertTrue(any("small tail" in x for x in result["warnings"]))

    def test_a_rare_segment_can_move_to_flow_without_being_silently_dropped(self):
        initial=score()
        rare=quality_fixtures.segment(name="BHT_seg1U_CBTOF_vs_DeltaE_Pi",count=1)
        initial["per_det"]["BHT"]["segments"].append(rare)
        initial["input_integrity"]=copy.deepcopy(phc_quality.integrity_from_score(initial))
        new=copy.deepcopy(initial);new["round"]=1
        row=new["per_det"]["BHT"]["segments"][1]
        row["tof_counts"]=[0.]*len(row["tof_counts"]);row["tof_counts"][-1]=1.
        row.update(n_entries=0.,n_tof_flow=1.,width=None)
        result=phc_comparison.compare_candidates([initial,new])[-1]
        self.assertTrue(result["comparable"],result["reasons"])
        self.assertEqual(2,result["required_segments"])
        self.assertEqual(0,result["missing_segments"])
        self.assertEqual(8001,result["quality"]["BHT"]["current"]["n_total"])

    def test_external_parameter_fallback_and_runlist_resolution(self):
        root=self.root/"without_params.root"
        self.write_root(root,include_params=False)
        hdprm=self.root/"HDPRM";hdprm.write_text("gain and time0\n")
        user=self.root/"USER";user.write_text("window\n")
        conf=self.root/"hodo.conf";conf.write_text(f"HDPRM: {hdprm}\nUSER: {user}\n")
        yml=self.root/"hodo.yml";yml.write_text(f"DEFAULT:\n  conf: {conf}\nRUN:\n  1: {{}}\n")
        files=monitor.parameters_from_runlist(yml)
        self.assertEqual(hdprm.resolve(),files["HDPRM"])
        with patch.object(monitor.config,"hodo_detectors",return_value=["BHT"]):
            initial=monitor.score_btof(1,root_path=root,scoring_version=3,parameter_files=files)
            initial["round"]=0;initial["input_integrity"]=phc_quality.integrity_from_score(initial)
            self.assertTrue(phc_comparison.compare_candidates([initial])[0]["comparable"])
            user.write_text("changed window\n")
            new=monitor.score_btof(1,root_path=root,scoring_version=3,parameter_files=files)
            new["round"]=1
        self.assertFalse(phc_comparison.compare_candidates([initial,new])[-1]["comparable"])

    def test_old_v2_series_resumes_with_its_existing_policy(self):
        iterator=adaptive.iterator
        initial=quality_fixtures.score()
        parameter=rounds.hdphc_param_path(1,"Pi");parameter.parent.mkdir(parents=True)
        parameter.write_text("existing")
        rounds.snapshot_round(1,1,"Pi",initial)
        args=SimpleNamespace(run=1,kaon=False,ftof=False,force=False,local=True,hdphc_iter=4,hdphc_patience=3)
        with patch.object(iterator,"stage_prepare_and_decode"), \
             patch.object(monitor,"initial_candidate") as new_series, \
             patch.object(monitor,"score_btof",return_value=initial) as scorer, \
             patch.object(iterator.decode_run,"run_decode"), \
             patch.object(iterator.hodo_qa_pdf,"source_paths",return_value=[]), \
             patch.object(iterator,"run_hodo_ok"):
            iterator.run_hdphc_one_round(args,2,clear_journal=False)
        new_series.assert_not_called()
        self.assertEqual(2,scorer.call_args.kwargs["scoring_version"])
        self.assertEqual(phc_quality.DEFAULT_POLICY,scorer.call_args.kwargs["quality_policy"])

    def test_large_tail_change_remains_a_blocker(self):
        initial=score();new=update(initial)
        row=new["per_det"]["BHT"]["segments"][0]
        row["tof_counts"][50]-=80;row["tof_counts"][65]+=80
        result=phc_comparison.compare_candidates([initial,new])[-1]
        self.assertFalse(result["comparable"])

    def test_small_flow_is_visible_but_large_flow_still_rejected(self):
        helper=quality_fixtures.QualityTests()
        initial=score();new=update(initial)
        helper.move_flow(new,20)
        result=phc_comparison.compare_candidates([initial,new])[-1]
        self.assertTrue(result["comparable"],result["reasons"])
        self.assertTrue(result["warnings"])
        helper.move_flow(new,60)
        self.assertFalse(phc_comparison.compare_candidates([initial,new])[-1]["comparable"])

    def test_tiny_width_gain_keeps_initial_correction(self):
        initial=score()
        self.assertEqual(0,rounds.pick_round([initial,update(initial,sigma=.299)],"auto"))

    def test_material_peak_gain_adopts_update(self):
        initial=score()
        self.assertEqual(1,rounds.pick_round([initial,update(initial,sigma=.26)],"auto"))

    def test_center_gain_can_justify_comparable_width(self):
        initial=score(center=.2)
        self.assertEqual(1,rounds.pick_round([initial,update(initial,center=0.)],"auto"))

    def test_center_gain_does_not_justify_a_much_wider_peak(self):
        initial=score(center=.2)
        self.assertEqual(0,rounds.pick_round([initial,update(initial,sigma=.4,center=0.)],"auto"))

    def test_nonfinite_hit_is_never_a_soft_warning(self):
        initial=score();new=update(initial,sigma=.2)
        new["per_det"]["BHT"]["time_validity"]["nonfinite_keys"]=[[1,0]]
        result=phc_comparison.compare_candidates([initial,new])[-1]
        self.assertFalse(result["comparable"])
        self.assertTrue(any("non-finite" in x for x in result["reasons"]))

    def test_changed_adc_or_user_parameters_are_not_comparable(self):
        for name in ["HDPRM","USER"]:
            initial=score();new=update(initial,sigma=.2)
            new["scoring_context"]["fixed_parameters"][name]="changed"
            result=phc_comparison.compare_candidates([initial,new])[-1]
            self.assertFalse(result["comparable"])

    def test_rejected_rounds_count_toward_patience_after_baseline(self):
        initial=score();history=[initial]
        for index in range(1,4):
            new=update(initial,sigma=.4,round_index=index)
            history.append(new)
        progress=rounds.adaptive_progress(history)
        self.assertTrue(progress["finished"])
        self.assertEqual(0,progress["best_round"])
        self.assertEqual(3,progress["without_improvement"])
        self.assertEqual(3,progress["completed"])

    def test_initial_snapshot_restores_params_and_invalidates_unmatched_report(self):
        parameter=rounds.hdphc_param_path(1,"Pi");parameter.parent.mkdir(parents=True)
        parameter.write_text("initial")
        initial=score()
        rounds.snapshot_round(1,0,"Pi",initial)
        parameter.write_text("updated")
        from common.lib import phc_shared
        phc_shared.write_report(1,"Pi",{"decision":"written","param_sha":"other","detectors":[]})
        rounds.restore_round(1,0,"Pi")
        self.assertEqual("initial",parameter.read_text())
        self.assertEqual("retained_initial",json.loads(phc_shared.report_path(1,"Pi").read_text())["decision"])

    def test_initial_snapshot_does_not_attach_invalid_old_report(self):
        from common.lib import phc_shared
        parameter=rounds.hdphc_param_path(1,"Pi");parameter.parent.mkdir(parents=True)
        parameter.write_text("initial")
        report=phc_shared.report_path(1,"Pi");report.parent.mkdir(parents=True,exist_ok=True)
        report.write_text("[]")
        entry=rounds.snapshot_round(1,0,"Pi",score())
        self.assertNotIn("shared_phc",entry)
        with self.assertRaises(ValueError):
            rounds.snapshot_round(1,1,"Pi",score(1))

    def test_failed_report_restore_cannot_leave_old_validation_marked_successful(self):
        from common.lib import phc_shared,qa_session
        parameter=rounds.hdphc_param_path(1,"Pi");parameter.parent.mkdir(parents=True)
        parameter.write_text("initial")
        phc_shared.write_report(1,"Pi",{"decision":"written","param_sha":qa_session.file_sha256(parameter),"detectors":[]})
        rounds.snapshot_round(1,0,"Pi",score())
        parameter.write_text("updated")
        validation=rounds.rounds_dir(1)/"adoption_validation.json"
        validation.write_text(json.dumps({"round":1,"matches":True}))
        with patch.object(phc_shared,"write_report",side_effect=OSError("publication failed")):
            with self.assertRaises(OSError):
                rounds.restore_round(1,0,"Pi")
        result=json.loads(validation.read_text())
        self.assertIsNone(result["matches"])
        self.assertEqual("pending_decode",result["status"])

    def test_first_update_snapshots_initial_values_before_fitting(self):
        iterator=adaptive.iterator
        parameter=rounds.hdphc_param_path(1,"Pi");parameter.parent.mkdir(parents=True)
        parameter.write_text("initial")
        args=SimpleNamespace(run=1,kaon=False,ftof=False,force=False,local=True,hdphc_iter=4,hdphc_patience=3)
        initial=score()
        with patch.object(iterator,"stage_prepare_and_decode"), \
             patch.object(monitor,"initial_candidate",return_value=initial), \
             patch.object(monitor,"score_btof",return_value=update(initial,sigma=.2)), \
             patch.object(iterator.decode_run,"run_decode"), \
             patch.object(iterator.hodo_qa_pdf,"source_paths",return_value=[]), \
             patch.object(iterator,"run_hodo_ok",side_effect=lambda *args:parameter.write_text("updated")):
            iterator.run_hdphc_one_round(args,1,clear_journal=True)
        history=rounds.load_journal(1)
        self.assertEqual([0,1],[r["round"] for r in history])
        self.assertEqual("initial",(rounds.round_dir(1,0)/parameter.name).read_text())
        self.assertEqual("updated",parameter.read_text())

    def test_first_round_retry_keeps_the_original_baseline(self):
        iterator=adaptive.iterator
        parameter=rounds.hdphc_param_path(1,"Pi");parameter.parent.mkdir(parents=True)
        initial=score();parameter.write_text("original");rounds.snapshot_round(1,0,"Pi",initial)
        parameter.write_text("partly-updated")
        seen=[]
        args=SimpleNamespace(run=1,kaon=False,ftof=False,force=False,local=True,hdphc_iter=4,hdphc_patience=3,only="hdphc")
        def fit(*args):
            seen.append(parameter.read_text());parameter.write_text("successful-retry")
            from common.lib import phc_shared,qa_session
            phc_shared.write_report(1,"Pi",{"decision":"written","param_sha":qa_session.file_sha256(parameter),"detectors":[]})
        with patch.dict("os.environ",{"E72_CALIB_STEP_RETRY":"1"}), \
             patch.object(iterator,"stage_prepare_and_decode"), \
             patch.object(monitor,"initial_candidate",return_value=initial) as initial_read, \
             patch.object(monitor,"score_btof",return_value=update(initial,sigma=.2)), \
             patch.object(iterator.decode_run,"run_decode"), \
             patch.object(iterator.hodo_qa_pdf,"source_paths",return_value=[]), \
             patch.object(iterator,"run_hodo_ok",side_effect=fit):
            iterator.run_hdphc_one_round(args,1,clear_journal=True)
        self.assertEqual(["original"],seen)
        initial_read.assert_called_once()
        self.assertEqual("original",(rounds.round_dir(1,0)/parameter.name).read_text())

    def test_hold_history_is_saved_before_the_held_state_is_published(self):
        order=[]
        with patch.object(adaptive.cd,"append_hold_journal",side_effect=lambda *args:order.append("journal")), \
             patch.object(adaptive.cd,"write_state",side_effect=lambda *args:order.append("state")):
            adaptive.cd.enter_held(1,{"run":1},{"source":"exit_check","id":"hold"},"hold")
        self.assertEqual(["journal","state"],order)

    def test_source_changes_between_rounds_stop_before_another_fit(self):
        iterator=adaptive.iterator
        parameter=rounds.hdphc_param_path(1,"Pi");parameter.parent.mkdir(parents=True)
        initial=score();parameter.write_text("initial");rounds.snapshot_round(1,0,"Pi",initial)
        parameter.write_text("first-update");rounds.snapshot_round(1,1,"Pi",update(initial,sigma=.28))
        changed=copy.deepcopy(initial)
        changed["scoring_context"]["fixed_parameters"]["HDPRM"]="changed"
        args=SimpleNamespace(run=1,kaon=False,ftof=False,force=False,local=True,hdphc_iter=4,hdphc_patience=3)
        with patch.object(iterator,"stage_prepare_and_decode"), \
             patch.object(monitor,"initial_candidate",return_value=changed), \
             patch.object(iterator,"run_hodo_ok") as fit:
            with self.assertRaises(SystemExit) as stopped:
                iterator.run_hdphc_one_round(args,2,clear_journal=False)
        self.assertEqual(iterator.runner.EXIT_CHECK,stopped.exception.code)
        fit.assert_not_called()
        self.assertEqual("first-update",parameter.read_text())

    def test_preround_audit_allows_corrected_times_and_profiles_to_evolve(self):
        initial=score();new=update(initial,sigma=.2)
        new["per_det"]["BHT"]["time_validity"]["nonfinite_keys"]=[[1,0]]
        self.assertEqual(phc_quality.stable_input(initial),phc_quality.stable_input(new))
        new["per_det"]["BHT"]["time_validity"]["raw_hit_signature"]="changed"
        self.assertNotEqual(phc_quality.stable_input(initial),phc_quality.stable_input(new))

    def test_initial_only_rollback_generates_qa_without_a_zero_round_count(self):
        iterator=adaptive.iterator
        parameter=rounds.hdphc_param_path(1,"Pi");parameter.parent.mkdir(parents=True)
        initial=score()
        initial["per_det"]["BHT"].update(width=.3,n_entries=8000.,n_seg=1)
        parameter.write_text("initial");rounds.snapshot_round(1,0,"Pi",initial)
        parameter.write_text("partial-update")
        args=SimpleNamespace(run=1,kaon=False,ftof=False,force=False,local=True,adopt="0")
        with patch.object(iterator,"stage_prepare_and_decode"), \
             patch.object(monitor,"score_btof",return_value=initial), \
             patch.object(iterator,"run_hodo_ok") as qa:
            iterator.run_hdphc_adopt_stage(args)
        self.assertEqual("initial",parameter.read_text())
        self.assertIn("--adopted-round 0",qa.call_args.args[1])
        self.assertNotIn("--hdphc-rounds 0",qa.call_args.args[1])
        validation=json.loads((rounds.rounds_dir(1)/"adoption_validation.json").read_text())
        self.assertTrue(validation["matches"])
        self.assertEqual("retained_initial",validation["decision"])

    def test_finish_request_accepts_initial_snapshot_and_real_rounds(self):
        parameter=rounds.hdphc_param_path(1,"Pi");parameter.parent.mkdir(parents=True)
        initial=score();parameter.write_text("initial");rounds.snapshot_round(1,0,"Pi",initial)
        for index in [1,2]:
            parameter.write_text("updated")
            rounds.snapshot_round(1,index,"Pi",update(initial,sigma=.4,round_index=index))
        request=rounds.request_finish(1,hold_id="review")
        self.assertEqual(0,request["best_round"])
        self.assertEqual(2,request["completed_round"])
        self.assertEqual("更新前の値",rounds.round_label(0))

    def test_retained_caption_never_claims_round_zero_was_fitted(self):
        with patch.object(slack_notify,"post_pdf",return_value=(True,"ok")) as post:
            slack_notify.post_hodo_qa(1,"hdphc",manifest={"meta":{"adopted_round":0}})
        caption=post.call_args.args[3]
        self.assertIn("更新前の値を維持",caption)
        self.assertNotIn("第0",caption)

    def test_adopted_caption_shows_the_selected_quality_and_small_change_note(self):
        with patch.object(slack_notify,"post_pdf",return_value=(True,"ok")) as post:
            slack_notify.post_hodo_qa(1,"hdphc",manifest={"meta":{"adopted_round":1,
                "phc_quality_lines":["BHT: 中心 +0.001 ns / 幅 0.210 ns"],
                "phc_quality_warnings":["small tail"]}})
        caption=post.call_args.args[3]
        self.assertIn("BHT: 中心",caption)
        self.assertIn("許容範囲内",caption)

    def write_root(self,path,energy=1.,user="same",include_params=True):
        row=quality_fixtures.segment()
        normal=np.asarray(row["tof_counts"][1:-1],dtype=float)
        with uproot.recreate(path) as root:
            if include_params:
                root["param/HDPRM"]="fixed HDPRM"
                root["param/USER"]=user
            root[row["histogram"]]=(normal[None,:],np.array([0.,3.]),np.asarray(quality_fixtures.EDGES))
            root["hodo"]={"event_number":np.array([1,2,3],dtype=np.uint32),
                          "beam_flag":np.ones(3,dtype=np.int32),
                          "bht_hit_seg":ak.Array([[0],[0],[0]]),
                          "bht_mt":ak.Array([[[0.]],[[0.]],[[0.]]]),
                          "bht_cmt":ak.Array([[[0.]],[[0.]],[[0.]]]),
                          "bht_de_u":ak.Array([[energy],[1.],[1.]])}

    def test_nonfinite_adc_can_be_fingerprinted_without_losing_validity(self):
        root=self.root/"nan_adc.root"
        self.write_root(root,energy=float("nan"))
        with patch.object(monitor.config,"hodo_detectors",return_value=["BHT"]):
            initial=monitor.score_btof(1,root_path=root,scoring_version=3)
        self.assertTrue(initial["per_det"]["BHT"]["time_validity"])

    def test_missing_initial_input_audit_stops_before_fit(self):
        iterator=adaptive.iterator
        initial=score();initial["scoring_context"]["fixed_parameters"]["USER"]=None
        args=SimpleNamespace(run=1,kaon=False,ftof=False,force=False,local=True,hdphc_iter=4,hdphc_patience=3)
        with patch.object(iterator,"stage_prepare_and_decode"), \
             patch.object(monitor,"initial_candidate",return_value=initial), \
             patch.object(iterator,"run_hodo_ok") as fit:
            with self.assertRaises(SystemExit) as stopped:
                iterator.run_hdphc_one_round(args,1,clear_journal=True)
        self.assertEqual(iterator.runner.EXIT_CHECK,stopped.exception.code)
        fit.assert_not_called()

    def test_native_energy_and_parameter_fingerprints_detect_changes(self):
        paths=[self.root/f"input{i}.root" for i in range(3)]
        self.write_root(paths[0]);self.write_root(paths[1],energy=1.1);self.write_root(paths[2],user="changed")
        with patch.object(monitor.config,"hodo_detectors",return_value=["BHT"]):
            values=[monitor.score_btof(1,root_path=path,scoring_version=3) for path in paths]
        first=values[0];first["round"]=0;first["input_integrity"]=phc_quality.integrity_from_score(first)
        for changed in values[1:]:
            changed["round"]=1
            result=phc_comparison.compare_candidates([first,changed])[-1]
            self.assertFalse(result["comparable"],result)
        self.assertTrue(values[0]["per_det"]["BHT"]["time_validity"])
        self.assertNotEqual(values[0]["per_det"]["BHT"]["time_validity"]["raw_hit_signature"],
                            values[1]["per_det"]["BHT"]["time_validity"]["raw_hit_signature"])

if __name__=="__main__":unittest.main()
