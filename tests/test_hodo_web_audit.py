"""Regression cases found when auditing Hodo review boundaries and persistence."""
import copy
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import test_t0_web as fixtures
web=fixtures.web


class WebAuditTests(unittest.TestCase):
    setUp=fixtures.T0WebTests.setUp
    saved_input=fixtures.T0WebTests.saved_input
    gate=fixtures.T0WebTests.gate
    staged=fixtures.T0WebTests.staged
    final_review=fixtures.T0WebTests.final_review

    def parent_final(self, step):
        web.write_json(self.app.base/'calib_dispatch_state.json',{
            'status':'awaiting_final_approval','step':step,'rc':0,'finished':True,'index':0,'calib_particle':'pion'})

    def test_active_lane_cannot_fall_back_to_old_parent_final(self):
        path,_=self.gate()
        web.write_json(path,{'status':'running','step':'hodo_t0'})
        self.parent_final('hodo_t0')
        with self.assertRaises(ValueError):self.app.gate()

    def test_different_lane_stage_cannot_use_parent_adc_review(self):
        self.gate(step='hodo_t0')
        self.parent_final('hodo_hdprm')
        with self.assertRaises(ValueError):web.App(12).gate()

    def test_approved_lane_cannot_fall_back_to_old_parent_review(self):
        self.final_review(status='approved')
        self.parent_final('hodo_t0')
        with self.assertRaises(ValueError):self.app.gate()

    def test_lane_auto_setting_is_authoritative(self):
        path,_=self.gate()
        (path.parent/'session.yml').write_text('auto: true\n')
        with self.assertRaises(ValueError):self.app.gate()

    def test_lane_auto_state_cannot_use_manual_session(self):
        path,state=self.gate()
        state["approval_mode"]="auto"
        web.write_json(path,state)
        with self.assertRaises(ValueError):self.app.gate()

    def replay_job(self, status, step):
        path=self.final_review(status=status,step=step)
        doc=self.app.document()
        doc['job']={'status':'running','state_path':str(path),'hold_id':'old-hold','started_at':1.}
        web.write_json(self.app.staging_path,doc)
        return self.app.state()['job']

    def test_quick_slack_approval_is_still_successful_web_replay(self):
        self.assertEqual('done',self.replay_job('approved','hodo_t0')['status'])

    def test_different_completed_stage_is_not_our_replay_success(self):
        self.assertEqual('failed',self.replay_job('awaiting_final_approval','hodo_hdphc')['status'])

    def test_invalid_new_settings_do_not_replace_saved_document(self):
        original={'fit':{},'calib_locks':{}}
        web.fit_overrides.save(12,original)
        path=web.fit_overrides.overrides_path(12)
        before=path.read_bytes()
        invalid=copy.deepcopy(original)
        invalid['fit']={'BH2':{'t0':{'bh2-3-t0':{'gaussian_fit':{'Pi':{**self.settings,'parameters':[1.,0.,float('nan')]}}}}}}
        with self.assertRaises(ValueError):web.fit_overrides.save(12,invalid)
        self.assertEqual(before,path.read_bytes())

    def test_broken_existing_settings_are_not_silently_replaced(self):
        path=web.fit_overrides.overrides_path(12)
        path.parent.mkdir(parents=True,exist_ok=True)
        for original in ('{broken','[]','{"fit": []}'):
            path.write_text(original)
            with self.assertRaises(ValueError):web.fit_overrides.set_mip_left(12,'BH2',3,'U',120.)
            self.assertEqual(original,path.read_text())

    def test_t0_flat_publication_is_atomic(self):
        path=web.fit_overrides.write_flat_t0_fit(12)
        before=path.read_bytes()
        with patch.object(web.fit_overrides.os,'replace',side_effect=OSError('publish failed')):
            with self.assertRaises(OSError):web.fit_overrides.write_flat_t0_fit(12)
        self.assertEqual(before,path.read_bytes())

    def test_curve_cannot_return_success_after_source_changes(self):
        def preview(*args, **kwargs):
            self.app.fingerprint=lambda:'changed-source'
            return {'curve':[[0.,1.]]}
        self.app.execute_preview=preview
        with self.assertRaisesRegex(ValueError,'変わり'):
            self.app.fit(self.request,curve_only=True)

    def test_bad_fit_row_does_not_hide_all_channels(self):
        import numpy as np
        import uproot
        path=self.app.fit_output()
        columns=uproot.open(path)['t0_fit'].arrays(library='np')
        columns['parameters'][3][0]=np.nan
        with uproot.update(path) as output:
            schema={'seg':'int32','parameters':'3 * float64','range':'2 * float64',
                    'status':'int32','covariance':'int32','ndf':'int32','success':'bool'}
            output.mktree('t0_fit',schema).extend(columns)
        result=self.app.overview({'fingerprint':'source-t0','page':0})
        bad=next(channel['automatic'] for channel in result['channels'] if channel['key']=='bh2-3-t0')
        self.assertEqual(12,len(result['channels']))
        self.assertFalse(bad['success'])
        self.assertTrue(bad['unusable'])
        self.assertFalse(bad['curve'])
        self.assertTrue(np.all(np.isfinite(bad['parameters'])))

    def test_legacy_parent_approval_refuses_reopened_active_lane(self):
        self.check_approval_guard(legacy=True)

    def test_final_approval_refuses_reopened_active_lane(self):
        self.check_approval_guard()

    def check_approval_guard(self, legacy=False):
        sys.path.insert(0,str(web.PROJECT/'calibration/scripts'))
        import calib_dispatch as dispatch
        path,_=self.gate()
        web.write_json(path,{'status':'running','step':'hodo_t0'})
        parent={'status':'awaiting_final_approval','finished':True,'rc':0,
                'plan':{'steps':[{'name':'hodo_t0','opts':{}}]}}
        if legacy:
            parent.pop('plan')
            parent['step']='hodo_t0'
        with patch.object(dispatch,'resolve_approval_state',return_value=(parent,'')), \
             patch.object(dispatch.calib_approval,'build_subject',return_value={}) as subject, \
             patch.object(dispatch.calib_approval,'subject_missing',return_value=[]), \
             patch.object(dispatch.calib_approval,'subject_sha',return_value='sha'), \
             patch.object(dispatch,'write_state'), \
             patch.object(dispatch.calib_approval,'append_approval_journal'):
            with self.assertRaises(SystemExit) as stopped:
                dispatch.do_approve({'run':12},None,'cli',None)
        self.assertEqual(1,stopped.exception.code)
        subject.assert_not_called()


if __name__=='__main__':unittest.main()
