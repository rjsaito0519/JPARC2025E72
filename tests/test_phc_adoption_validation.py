"""Adoption must use the recorded parameter snapshot and reproduce its score."""
import copy
import json
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import test_phc_adaptive as fixtures
iterator=fixtures.iterator
rounds=fixtures.rounds


class PhcAdoptionValidationTests(unittest.TestCase):
    setUp=fixtures.AdaptivePhcTests.setUp
    completed_series=fixtures.AdaptivePhcTests.completed_series

    def adopt(self,score):
        self.completed_series()
        stale=self.root/'analyzer/runmanager/runlist/e72/hodo_run00001.yml'
        stale.parent.mkdir(parents=True)
        stale.write_text('old runlist from another stage\n')
        args=SimpleNamespace(run=1,kaon=False,ftof=False,force=False,local=True,adopt='auto')
        with patch.object(iterator.config,'ANALYZER_DIR',self.root/'analyzer'), \
             patch.object(iterator,'stage_prepare_and_decode',return_value=stale) as prepare, \
             patch.object(iterator.hodo_btof_monitor,'score_btof',return_value=score), \
             patch.object(iterator,'run_hodo_ok') as qa:
            try:iterator.run_hdphc_adopt_stage(args)
            finally:
                prepare.assert_called_once_with(1,False,'hdphc',local=True,dry_run=False)
                qa.assert_called_once()

    def test_no_eligible_candidate_sends_qa_and_stops_check_without_restoring(self):
        live=self.completed_series()
        args=SimpleNamespace(run=1,kaon=False,ftof=False,force=False,local=True,adopt='auto')
        with patch.object(rounds,'compare_candidates',return_value=[{'round':1,'comparable':False}]), \
             patch.object(rounds,'write_summary'), \
             patch.object(rounds,'restore_round') as restore, \
             patch.object(iterator,'stage_prepare_and_decode') as decode, \
             patch.object(iterator,'run_hodo_ok') as qa:
            with self.assertRaises(SystemExit) as stopped:
                iterator.run_hdphc_adopt_stage(args)
            self.assertEqual(iterator.runner.EXIT_CHECK,stopped.exception.code)
            restore.assert_not_called()
            decode.assert_not_called()
            qa.assert_called_once()
        self.assertEqual('round3',live.read_text())

    def test_correct_adoption_regenerates_conf_and_records_verified_score(self):
        self.adopt(fixtures.comparable_score(1.))
        validation=json.loads((rounds.rounds_dir(1)/'adoption_validation.json').read_text())
        self.assertTrue(validation['matches'])
        self.assertEqual(2,validation['round'])

    def test_fresh_decode_disagrees_with_recorded_candidate_stops_check(self):
        with self.assertRaises(SystemExit) as stopped:
            self.adopt({'width':1.66,'ok':True,'n_entries':100})
        self.assertEqual(3,stopped.exception.code)
        validation=json.loads((rounds.rounds_dir(1)/'adoption_validation.json').read_text())
        self.assertFalse(validation['matches'])
        self.assertTrue(validation['differences'])

    def test_unchanged_median_does_not_hide_changed_detector_score(self):
        expected={'ok':True,'width':1.,'n_entries':100,
                  'per_det':{'BHT':{'width':.8,'n_entries':40,'n_seg':5}}}
        actual=copy.deepcopy(expected)
        actual['per_det']['BHT']['width']=.9
        self.assertTrue(rounds.score_differences(expected,actual))
        self.assertFalse(rounds.score_differences(expected,copy.deepcopy(expected)))

    def test_changed_channel_cannot_hide_behind_equal_detector_median(self):
        expected={'ok':True,'width':1.,'n_entries':100,'per_det':{'BH2':{
            'width':1.,'n_entries':100,'n_seg':1,
            'segments':[{'histogram':'segment','width':1.,'n_entries':100,'n_total':110}]}}}
        actual=copy.deepcopy(expected)
        actual['per_det']['BH2']['segments'][0]['n_total']=120
        self.assertTrue(rounds.score_differences(expected,actual))

    def test_empty_detector_remains_unavailable_without_false_mismatch(self):
        expected={'ok':True,'width':1.,'n_entries':100,'per_det':{'T1':{
            'width':float('nan'),'n_entries':0,'n_seg':0,
            'segments':[{'histogram':'empty','width':None,'n_entries':0,'n_total':20}]}}}
        self.assertFalse(rounds.score_differences(expected,copy.deepcopy(expected)))
        actual=copy.deepcopy(expected)
        actual['per_det']['T1']['segments'][0]['n_total']=21
        self.assertTrue(rounds.score_differences(expected,actual))

    def test_restoring_another_round_invalidates_previous_validation(self):
        self.completed_series()
        path=rounds.rounds_dir(1)/'adoption_validation.json'
        path.write_text(json.dumps({'round':3,'matches':True}))
        rounds.restore_round(1,2,'Pi')
        result=json.loads(path.read_text())
        self.assertIsNone(result['matches'])
        self.assertEqual('pending_decode',result['status'])
        self.assertEqual(2,result['round'])

    def test_restoring_old_snapshot_does_not_hide_change_with_old_mtime(self):
        live=self.completed_series()
        snapshot=rounds.round_dir(1,2)/live.name
        rounds.os.utime(snapshot,(1,1))
        rounds.restore_round(1,2,'Pi')
        self.assertGreater(live.stat().st_mtime,snapshot.stat().st_mtime)

    def test_corrupt_snapshot_does_not_replace_live_parameters(self):
        live=self.completed_series()
        snapshot=rounds.round_dir(1,2)/live.name
        snapshot.write_text('corrupted snapshot\n')
        with self.assertRaises(ValueError):rounds.restore_round(1,2,'Pi')
        self.assertEqual('round3',live.read_text())

    def test_failed_publication_keeps_live_parameters_intact(self):
        live=self.completed_series()
        with patch.object(rounds.os,'replace',side_effect=OSError('publication failed')):
            with self.assertRaises(OSError):rounds.restore_round(1,2,'Pi')
        self.assertEqual('round3',live.read_text())


if __name__=='__main__':unittest.main()
