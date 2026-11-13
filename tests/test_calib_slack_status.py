"""Slack parent status mirrors durable lane state without changing calibration."""
import json
import sys
import tempfile
import threading
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from common.lib import calib_slack_status as status, config, slack_notify
from common.lib.calib_state_io import locked_path, write_json

class ParentStatusTests(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.output = Path(self.scratch.name)
        patch = mock.patch.object(config, 'OUTPUT_DIR', self.output)
        patch.start(); self.addCleanup(patch.stop)
        self.directory = self.output / 'calib/run02977'
        self.directory.mkdir(parents=True)
        self.cfg = {'token': 'test', 'channel_id': 'C1'}
        patch = mock.patch.object(slack_notify, 'load_slack_config', return_value=(self.cfg, 'ok'))
        patch.start(); self.addCleanup(patch.stop)

    def lane(self, name, **fields):
        payload = {'status': 'running', 'step': 'hodo_hdphc',
                   'opts': {'hdphc_round': 3}, 'updated_at': '2026-11-13T06:00:00Z', **fields}
        write_json(self.directory / f'lanes/{name}/calib_dispatch_state.json', payload)
        return payload

    def thread(self, **fields):
        write_json(slack_notify.state_path(2977), {'thread_ts': '100.001', 'channel_id': 'C1', 'posts': {}, **fields})

    def test_user_stopped_successful_round_is_not_waiting_or_complete(self):
        self.lane('hodo', status='held', rc=0, hold={'decision': 'stop'})
        text = status.summary_text(2977)
        self.assertIn('停止中（ユーザー操作）', text)
        self.assertIn('第3ラウンド', text)
        self.assertNotIn('確認待ち', text)
        self.assertIn('11/13 15:00 JST', text)

    def test_actual_lanes_override_old_coordinator_and_d5_skipped_pass_number(self):
        write_json(self.directory/'calib_dispatch_state.json', {'status':'running', 'lanes':{'d5':'running'}})
        self.lane('d5', status='awaiting_final_approval', step='d5_pass', opts={'pass':10})
        (self.directory/'d5_converged.txt').write_text('pass=2\n')
        text = status.summary_text(2977)
        self.assertIn('第2周回で収束', text)
        self.assertIn('最終承認待ち', text)
        self.assertNotIn('第10', text)
        self.assertNotIn('実行中', text)
        self.assertNotIn('✅ 承認済み', text)

    def test_final_approval_covers_completed_lanes_but_not_restarted_work(self):
        write_json(self.directory/'calib_dispatch_state.json', {'status':'approved'})
        self.lane('hodo', status='awaiting_final_approval')
        self.assertIn('✅ 承認済み', status.summary_text(2977))
        self.lane('hodo', status='running')
        self.assertNotIn('✅ 承認済み', status.summary_text(2977))

    def test_other_finished_lanes_do_not_hide_stopped_hodo(self):
        self.lane('hodo', status='stopped')
        self.lane('bcin', status='awaiting_final_approval', step='dc_reso', opts={})
        text = status.summary_text(2977)
        self.assertIn('停止中・処理完了のレーンあり', text.splitlines()[0])
        self.assertIn('BcIn / BLC1', text)

    def test_parent_update_uses_original_ts_and_deduplicates(self):
        self.lane('hodo')
        self.thread()
        path = self.directory/'lanes/hodo/calib_dispatch_state.json'
        before = path.read_bytes()
        with mock.patch.object(slack_notify,'_api',return_value=(True,{},'ok')) as api:
            self.assertEqual(status.refresh_parent(2977),(True,'updated'))
            self.assertEqual(status.refresh_parent(2977),(True,'unchanged'))
        api.assert_called_once()
        self.assertEqual(api.call_args.args[0],'chat.update')
        self.assertEqual(api.call_args.args[2]['ts'],'100.001')
        self.assertEqual(path.read_bytes(),before)

    def test_state_change_refreshes_the_same_parent(self):
        self.lane('hodo'); self.thread()
        with mock.patch.object(slack_notify, '_api', return_value=(True, {}, 'ok')) as api:
            status.refresh_parent(2977)
            self.lane('hodo', status='held', rc=0, hold={'decision':None})
            status.refresh_parent(2977)
        self.assertEqual(api.call_count, 2)
        self.assertTrue(all(call.args[0] == 'chat.update' for call in api.call_args_list))
        self.assertTrue(all(call.args[2]['ts'] == '100.001' for call in api.call_args_list))
        self.assertIn('確認待ち', api.call_args.args[2]['text'])

    def test_failure_retries_and_no_thread_is_created(self):
        self.lane('hodo')
        with mock.patch.object(slack_notify,'_api') as api:
            self.assertTrue(status.refresh_parent(2977)[0])
        api.assert_not_called()
        self.assertFalse(slack_notify.state_path(2977).exists())
        self.thread()
        with mock.patch.object(slack_notify,'_api',return_value=(False,{},'rate limit')):
            self.assertEqual(status.refresh_parent(2977),(False,'rate limit'))
        self.assertNotIn('header_status',slack_notify.load_state(2977))
        with mock.patch.object(slack_notify,'_api',return_value=(True,{},'ok')) as api:
            self.assertTrue(status.refresh_parent(2977)[0])
        api.assert_called_once()

    def test_concurrent_pdf_receipt_is_preserved_and_not_blocked_by_network(self):
        self.lane('hodo'); self.thread()
        finished = threading.Event()
        def upload():
            with locked_path(slack_notify.state_path(2977)):
                state = slack_notify.load_state(2977)
                state['posts']['current'] = {'pdf':'new.pdf'}
                slack_notify.save_state(2977,state)
            finished.set()
        def api(*args):
            worker = threading.Thread(target=upload)
            worker.start()
            self.assertTrue(finished.wait(2), 'parent update blocked PDF state lock during network call')
            worker.join()
            return True,{},'ok'
        with mock.patch.object(slack_notify,'_api',side_effect=api):
            self.assertTrue(status.refresh_parent(2977)[0])
        self.assertEqual(slack_notify.load_state(2977)['posts']['current']['pdf'],'new.pdf')

    def test_old_channel_is_never_updated(self):
        self.lane('hodo'); self.thread(channel_id='OLD')
        with mock.patch.object(slack_notify,'_api') as api:
            self.assertIn('current channel',status.refresh_parent(2977)[1])
        api.assert_not_called()

if __name__ == '__main__': unittest.main()
