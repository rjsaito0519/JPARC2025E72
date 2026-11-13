"""Readable Slack blocks preserve the reviewed result and pinned actions."""
import sys
import unittest
from pathlib import Path

REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO))
from common.lib import calib_slack_status, slack_notify

class MessageLayoutTests(unittest.TestCase):
    def test_parent_uses_header_and_short_lane_columns(self):
        text='\n'.join([
            '*E72 Calibration · run02977* — ⏸ 停止中・処理完了のレーンあり',
            '*Hodo*：⏸ 停止中｜Hodo PHC 第1ラウンド',
            '*BcIn / BLC1*：☑ 処理完了・最終承認待ち｜DC 分解能補正',
            '*BcOut / BLC2*：☑ 処理完了・最終承認待ち｜DC 分解能補正',
            '*D5*：☑ 処理完了・最終承認待ち｜第2周回で収束',
            '状態の最終更新: 11/13 15:00 JST',
        ])
        blocks=calib_slack_status.summary_blocks(2977,text)
        self.assertEqual(blocks[0]['type'],'header')
        self.assertEqual(blocks[0]['text']['type'],'plain_text')
        fields=next(block['fields'] for block in blocks if 'fields' in block)
        self.assertEqual(len(fields),4)
        self.assertIn('停止中',fields[0]['text'])
        self.assertIn('第2周回で収束',fields[1]['text'])
        self.assertTrue(any(block['type']=='context' for block in blocks))

    def test_review_actions_stay_last_with_their_original_hold_id(self):
        actions={'type':'actions','block_id':'original','elements':[
            {'type':'button','action_id':'calib_step_proceed','value':'2977:hold-id'}]}
        original=[{'type':'section','text':{'type':'mrkdwn','text':
            '*run02977* [hodo] *Hodo PHC 第2ラウンド*\n'
            '対象PDF: `qa.pdf`（結果ID: current）\n'
            '今回の幅: 0.6000\n'
            '【確認待ち】更新・評価が完了しました。\n'
            '次へ進むと: 第3ラウンド'}},actions]
        rendered=slack_notify.review_blocks(original)
        self.assertEqual(rendered[0]['type'],'header')
        self.assertNotIn('*',rendered[0]['text']['text'])
        self.assertEqual(rendered[-1],actions)
        body=next(block['text']['text'] for block in rendered if block['type']=='section')
        self.assertTrue(body.startswith('【確認待ち】'))
        reference=next(block for block in rendered if block['type']=='context')
        self.assertIn('qa.pdf',reference['elements'][0]['text'])
        self.assertIn('current',reference['elements'][0]['text'])
        self.assertEqual(original[-1]['elements'][0]['value'],'2977:hold-id')
        self.assertEqual(slack_notify.review_blocks(rendered),rendered)

if __name__=='__main__':unittest.main()
