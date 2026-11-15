"""T0 tuning stages, saved inputs, overrides and hold boundaries."""
import copy
import importlib.util
import json
import sys
import tempfile
import threading
import urllib.request
import urllib.error
import unittest
from types import SimpleNamespace
from pathlib import Path
from unittest.mock import patch

import numpy as np
import uproot

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO))
spec = importlib.util.spec_from_file_location('t0_web',REPO/'calibration/mip_web/mip_web.py')
web=importlib.util.module_from_spec(spec)
spec.loader.exec_module(web)


class T0WebTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        root=Path(self.directory.name)
        for name in ('OUTPUT_DIR','WORK_DIR','PARAM_DIR','DATA_DIR'):
            patcher=patch.object(web.config,name,root/name)
            patcher.start()
            self.addCleanup(patcher.stop)
        self.app=web.App(12,mode='t0')
        self.app.fingerprint=lambda:'source-t0'
        self.settings={'parameters':[20.,-.3,.1],'range':[-.7,.1],'function':'gaussian','range_mode':'fixed'}
        self.result={'parameters':[21.,-.31,.09],'range':[-.7,.1],'function':'gaussian',
                     'success':True,'status':0,'covariance':3,'ndf':20,'curve':[[-.3,80.]]}
        self.request={'fingerprint':'source-t0','key':'bh2-3-t0','settings':self.settings}
        self.app.execute_preview=lambda *args,**kwargs:copy.deepcopy(self.result)
        self.saved_input()

    def saved_input(self):
        output=self.app.fit_output()
        output.parent.mkdir(parents=True,exist_ok=True)
        with uproot.recreate(output) as root:
            for segment in range(15):
                root[f'T0_seg{segment}_TimeOffset_Pi']=np.ones(2000),np.linspace(-10.,10.,2001)
            tree=root.mktree('t0_fit',{'seg':'int32','parameters':'3 * float64','range':'2 * float64',
                'status':'int32','covariance':'int32','ndf':'int32','success':'bool'})
            tree.extend({'seg':list(range(15)),'parameters':[[20.,-.3,.1]]*15,'range':[[-.7,.1]]*15,
                'status':[0]*15,'covariance':[3]*15,'ndf':[20]*15,'success':[True]*15})

    def gate(self,step='hodo_t0',particle='pion',auto=False):
        path=self.app.base/'lanes/hodo/calib_dispatch_state.json'
        hold={'id':'t0-hold','step':step,'rc':0}
        state={'status':'held','hold':hold,'calib_particle':particle}
        web.write_json(path,state)
        (self.app.base/'calib_session.yml').write_text(f'auto: {str(auto).lower()}\n')
        return path,hold

    def staged(self):
        result=self.app.fit(self.request)
        return self.app.stage({**self.request,'receipt':result['receipt']})

    def test_channels_and_staging_are_separate_from_adc(self):
        adc=web.App(12)
        self.assertEqual(15,len(self.app.channel_keys()))
        self.assertEqual((3,'T0'),self.app.channel('bh2-3-t0'))
        self.assertNotEqual(adc.staging_path,self.app.staging_path)
        with self.assertRaises(ValueError):self.app.channel('bh2-3-u')
        with self.assertRaises(ValueError):adc.channel('bh2-3-t0')

    def test_saved_gaussian_and_negative_axis_are_used(self):
        self.assertEqual([20.,-.3,.1],self.app.automatic('bh2-3-t0','source-t0')['parameters'])
        overview=self.app.overview({'fingerprint':'source-t0','page':0})
        self.assertEqual(12,len(overview['channels']))
        self.assertEqual(2,overview['pages'])
        self.assertEqual(666,len(overview['channels'][0]['counts']))
        self.assertEqual(-10.,overview['channels'][0]['bounds'][0])
        self.assertEqual('gaussian',self.app.edit(self.request)['initial']['function'])

    def test_invalid_parameters_and_ranges_never_reach_preview(self):
        self.app.execute_preview=lambda *args,**kwargs:self.fail('invalid settings reached ROOT')
        for change in ({'function':'mean'},{'parameters':[20.,0.,0.]},{'parameters':[20.,float('nan'),.1]},
                       {'parameters':[20.,0.,.1,0.]},{'range':[-11.,0.]},{'range':[0.,10.]},{'range':[1.,-1.]}):
            with self.assertRaises(ValueError):
                self.app.fit({**self.request,'settings':{**self.settings,**change}})

    def test_t0_hold_required_and_auto_rejected(self):
        self.gate(step='hodo_hdprm')
        with self.assertRaises(ValueError):self.app.gate()
        self.gate()
        self.assertEqual('t0-hold',self.app.gate()[2]['id'])
        self.gate(auto=True)
        with self.assertRaises(ValueError):self.app.gate()
        self.gate(particle='kaon')
        with self.assertRaises(ValueError):self.app.gate()

    def test_write_replays_t0_and_preserves_adc_and_other_particle(self):
        original={'fit':{'BH2':{'hdprm':{'bh2-3-u':{'mip_left':140.}},
                      't0':{'bh2-3-t0':{'gaussian_fit':{'K':self.settings}}}}},'calib_locks':{'BH2':{'hdprm':['lock']}}}
        web.fit_overrides.save(12,original)
        staged=self.staged()
        self.gate()
        with patch('common.lib.calib_ops.run_dispatch',return_value=(True,'accepted')) as replay:
            result=self.app.write({**self.request,'revision':staged['revision'],'hold_id':'t0-hold'})
        replay.assert_called_once_with('retry',12,hold_id='t0-hold',background_resume=True)
        self.assertEqual('running',result['job']['status'])
        saved=web.fit_overrides.load(12)
        self.assertEqual(original['fit']['BH2']['hdprm'],saved['fit']['BH2']['hdprm'])
        self.assertEqual(original['calib_locks'],saved['calib_locks'])
        self.assertEqual(self.settings,saved['fit']['BH2']['t0']['bh2-3-t0']['gaussian_fit']['K'])
        self.assertIn('Pi',saved['fit']['BH2']['t0']['bh2-3-t0']['gaussian_fit'])
        self.assertIn('bh2-3-t0:Pi gaussian fixed',web.fit_overrides.flat_t0_fit_path(12).read_text())

    def test_failed_fit_cannot_stage(self):
        self.result['success']=False
        self.assertNotIn('receipt',self.app.fit(self.request))
        with self.assertRaises(ValueError):self.app.stage({**self.request,'receipt':'fake'})

    def test_wrong_hold_does_not_write(self):
        staged=self.staged()
        self.gate()
        with self.assertRaises(ValueError):
            self.app.write({**self.request,'revision':staged['revision'],'hold_id':'wrong-stage'})
        self.assertFalse(web.fit_overrides.overrides_path(12).exists())

    def test_t0_runner_passes_saved_conditions_to_cpp(self):
        sys.path.insert(0,str(REPO/'calibration/scripts'))
        import run_hodo
        web.fit_overrides.save(12,{'fit':{'BH2':{'t0':{'bh2-3-t0':{'gaussian_fit':{'Pi':self.settings}}}}},'calib_locks':{}})
        args=SimpleNamespace(force=False,kaon=False,debug=True)
        with patch.object(run_hodo.runner,'require_bins',return_value=('T0_Offset',)), \
             patch.object(run_hodo.runner,'input_root',return_value=Path('/input.root')), \
             patch.object(run_hodo.qa_session,'begin',return_value={}), \
             patch.object(run_hodo.qa_session,'attach_params'), \
             patch.object(run_hodo.runner,'sh') as command, \
             patch.object(run_hodo.runner,'update_param',return_value=0), \
             patch.object(run_hodo,'after_qa_merge',return_value=0):
            run_hodo.do_t0(12,'Pi',args)
        override=Path(command.call_args.kwargs['env']['E72_T0_FIT_OVERRIDE'])
        self.assertIn('bh2-3-t0:Pi gaussian fixed',override.read_text())

    def final_review(self, status='awaiting_final_approval', step='hodo_t0'):
        path,_=self.gate()
        web.write_json(path,{'status':status,'step':step,'index':1,'rc':0,'finished':True,'calib_particle':'pion'})
        return path

    def test_pending_final_stage_can_be_reviewed_without_state_changes(self):
        path=self.final_review()
        original=path.read_bytes()
        _,_,context=self.app.gate()
        self.assertTrue(context['id'].startswith('final_review:'))
        self.assertEqual(original,path.read_bytes())

    def test_final_review_write_reopens_only_on_explicit_write(self):
        staged=self.staged()
        path=self.final_review()
        context=self.app.gate()[2]
        with patch('common.lib.calib_ops.run_dispatch',return_value=(True,'accepted')) as replay:
            result=self.app.write({**self.request,'revision':staged['revision'],'hold_id':context['id']})
        state=json.loads(path.read_text())
        self.assertEqual('held',state['status'])
        self.assertFalse(state['finished'])
        self.assertEqual('web_refit_review',state['hold']['source'])
        self.assertEqual('hodo_t0',state['hold']['step'])
        self.assertEqual(1,state['hold']['index'])
        self.assertNotEqual(context['id'],state['hold']['id'])
        self.assertEqual(state['hold']['id'],result['job']['hold_id'])
        replay.assert_called_once_with('retry',12,hold_id=state['hold']['id'],background_resume=True)

    def test_changed_final_review_context_cannot_be_written(self):
        staged=self.staged()
        path=self.final_review()
        context=self.app.gate()[2]
        state=json.loads(path.read_text())
        state['updated_at']='changed'
        web.write_json(path,state)
        with self.assertRaises(ValueError):
            self.app.write({**self.request,'revision':staged['revision'],'hold_id':context['id']})
        self.assertEqual(state,json.loads(path.read_text()))
        self.assertFalse(web.fit_overrides.overrides_path(12).exists())

    def test_approved_or_different_final_stage_cannot_be_reopened(self):
        for status,step in [('approved','hodo_t0'),('awaiting_final_approval','hodo_hdprm')]:
            self.final_review(status,step)
            with self.assertRaises(ValueError):self.app.gate()

    def test_slack_t0_link_selects_t0_screen(self):
        from common.lib import calib_approval
        text=json.dumps(calib_approval.step_gate_action_blocks(12,'hodo_t0',1,hold_id='t0-hold'))
        self.assertIn('mode=t0',text)
        self.assertIn('BH2 T0 offset',text)

    def test_http_routing_and_tokens_are_separate(self):
        adc=web.App(12)
        adc.fingerprint=lambda:'source-adc'
        server=web.ThreadingHTTPServer(('127.0.0.1',0),web.handler(adc,self.app))
        worker=threading.Thread(target=server.serve_forever,daemon=True)
        worker.start()
        base=f'http://127.0.0.1:{server.server_port}'
        try:
            with urllib.request.urlopen(base+'/state?mode=t0') as response:
                self.assertEqual('t0',json.load(response)['mode'])
            request=urllib.request.Request(base+'/edit?mode=t0',data=b'{}',headers={'X-MIP-Token':adc.token})
            with self.assertRaises(urllib.error.HTTPError) as error:
                urllib.request.urlopen(request)
            self.assertEqual(403,error.exception.code)
            with self.assertRaises(urllib.error.HTTPError) as error:
                urllib.request.urlopen(base+'/state?mode=unknown')
            self.assertEqual(409,error.exception.code)
        finally:
            server.shutdown()
            server.server_close()
            worker.join()

    def test_t0_replay_returns_to_t0_gate(self):
        self.staged()
        path,_=self.gate()
        doc=self.app.document()
        doc['job']={'status':'running','state_path':str(path),'hold_id':'old-hold'}
        web.write_json(self.app.staging_path,doc)
        result=self.app.state()
        self.assertEqual('done',result['job']['status'])
        self.assertFalse(result['staged'])


if __name__=='__main__':unittest.main()
