"""Review cards follow current PDFs and retire only after a successful choice."""
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO));sys.path.insert(0,str(REPO/'calibration/scripts'))
from common.lib import config, calib_approval, phc_range_gate, slack_notify
import run_hodo

class ReviewFlowTests(unittest.TestCase):
    def setUp(self):
        self.scratch=tempfile.TemporaryDirectory();self.addCleanup(self.scratch.cleanup)
        self.directory=Path(self.scratch.name)
        patch=mock.patch.object(config,'OUTPUT_DIR',self.directory);patch.start();self.addCleanup(patch.stop)
        patch=mock.patch.dict(os.environ,{'E72_CALIB_LANE':'hodo'});patch.start();self.addCleanup(patch.stop)

    def test_refit_range_card_is_new_for_each_pdf_and_updates_its_own_card(self):
        phc_range_gate.write_proposal(1,{'run':1,'proposals':[{'det':'BHT','first':1,'last':10}]})
        with mock.patch.object(phc_range_gate,'sync_proposal_currents',side_effect=lambda run:phc_range_gate.read_proposal(run)), \
                mock.patch.object(slack_notify,'try_notify_failure',return_value=(True,'posted')) as post, \
                mock.patch.object(slack_notify,'try_update_failure_message',return_value=(True,'updated')) as update:
            for result in ('date_first','date_second'):
                phc_range_gate.notify_pending_phc_ask(1,qa_exec_id=result)
                self.assertEqual(post.call_args.args[1],'phc_range_ask:'+result)
                self.assertIn(result.split('_')[-1],post.call_args.kwargs['blocks'][0]['text']['text'])
                phc_range_gate.update_phc_ask_slack(1)
                self.assertEqual(update.call_args.args[1],post.call_args.args[1])

    def test_pdf_is_posted_before_range_buttons_and_failure_has_no_range_buttons(self):
        args=SimpleNamespace(skip_qa=True,kaon=False,ftof=False,qa_only=False)
        for merge_ok,post_ok in ((True,True),(True,False),(False,False)):
            events=[]
            with self.subTest(merge_ok=merge_ok,post_ok=post_ok), \
                    mock.patch.object(phc_range_gate,'pending_phc_ask',return_value=True), \
                    mock.patch.object(run_hodo.hodo_qa_pdf,'merge',return_value=(merge_ok,'merge result',Path('qa.pdf'))), \
                    mock.patch.object(slack_notify,'try_post_hodo_after_merge',side_effect=lambda *a,**k:(events.append('pdf') or (post_ok,'result'))), \
                    mock.patch.object(phc_range_gate,'notify_pending_phc_ask',side_effect=lambda *a,**k:(events.append('buttons') or (True,'posted'))), \
                    mock.patch.object(slack_notify,'try_notify_failure',return_value=(True,'mocked')):
                run_hodo.after_qa_merge(1,'hdphc',args,3)
            self.assertEqual(events,['pdf','buttons'] if post_ok else (['pdf'] if merge_ok else []))

    def test_card_identifies_round_next_step_and_current_pdf(self):
        base=self.directory/'calib/run00001';lane=base/'lanes/hodo';lane.mkdir(parents=True)
        state={'step':'hodo_hdphc','index':0,'opts':{'hdphc_round':2},'plan':{'steps':[
            {'name':'hodo_hdphc','opts':{'hdphc_round':2}},
            {'name':'hodo_hdphc','opts':{'hdphc_round':3}}]}}
        (lane/'calib_dispatch_state.json').write_text(json.dumps(state))
        (base/'slack_thread.json').write_text(json.dumps({'posts':{'test':{'mode':'hdphc','pdf':'qa.pdf','exec_id':'date_current'}}}))
        text=calib_approval.step_gate_action_blocks(1,'hodo_hdphc',0)[0]['text']['text']
        for expected in ('第2ラウンド','第3ラウンド','qa.pdf','current','確認待ち'):
            self.assertIn(expected,text)
        state['hold']={'reason':'stale_input: changed'}
        (lane/'calib_dispatch_state.json').write_text(json.dumps(state))
        self.assertNotIn('qa.pdf',calib_approval.gate_context(1,'hodo_hdphc',0))

if __name__=='__main__':unittest.main()
