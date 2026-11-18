
"""Manual PHC choices remain identical between Python, Slack and persisted edits."""
import ast
import os
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import patch

from common.lib import config, phc_range_gate as gate


class ManualPhcSettingsTests(unittest.TestCase):
    def setUp(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.folder = Path(folder.name)
        self.conf = self.folder/'phc_conf.py'

    def write(self, text):
        self.conf.write_text(text, encoding='utf-8')

    def test_manual_python_syntax_and_comments(self):
        self.write("# '00001_Pi_bht': [0, 99]\nlimits_dict = {\n"
                   "    '00001_Pi_bht': [\n        12, # first\n        18,\n    ],\n"
                   "    \"00001_Pi_bh2\": (3, 9),\n}\n")
        self.assertEqual((12, 18), gate.read_hand_range(1, 'Pi', 'BHT', path=self.conf))
        self.assertEqual((3, 9), gate.read_hand_range(1, 'Pi', 'BH2', path=self.conf))
        self.assertEqual({'00001_Pi_bht', '00001_Pi_bh2'},
                         gate.conf_keys_on_disk(path=self.conf))

    def test_commented_range_is_not_manual(self):
        self.write("# \"00001_Pi_bht\": [1, 10]\nlimits_dict = {}\n")
        self.assertIsNone(gate.read_hand_range(1, 'Pi', 'BHT', path=self.conf))
        self.assertEqual(['BHT'], gate.missing_proposal_dets(
            1, 'Pi', [{'det': 'BHT'}], path=self.conf))

    def test_explicit_edit_preserves_comments_and_other_dictionary(self):
        self.write("limits_dict = {'00001_Pi_bht': [1, 10]} # keep this\n"
                   "other = {'untouched': 7}\n")
        gate.upsert_phc_conf(1, 'Pi', 'BHT', 4, 8, path=self.conf)
        gate.upsert_phc_conf(2, 'K', 'CVC', 2, 5, path=self.conf)
        text = self.conf.read_text()
        self.assertIn('# keep this', text)
        self.assertIn("other = {'untouched': 7}", text)
        self.assertEqual((4, 8), gate.read_hand_range(1, 'Pi', 'BHT', path=self.conf))
        self.assertEqual((2, 5), gate.read_hand_range(2, 'K', 'CVC', path=self.conf))
        self.assertEqual(1, text.count("'00001_Pi_bht'"))
        ast.parse(text)

    def test_inline_comment_without_comma_remains_valid(self):
        self.write("limits_dict = {\n    '00001_Pi_bht': [1, 8] # keep\n}\n")
        gate.upsert_phc_conf(2, 'Pi', 'BH2', 3, 9, path=self.conf)
        self.assertIn('# keep', self.conf.read_text())
        self.assertEqual((3, 9), gate.read_hand_range(2, 'Pi', 'BH2', path=self.conf))

    def test_unicode_before_value_preserves_exact_source_positions(self):
        self.write("limits_dict = {'説明': [1, 2], '00001_Pi_bht': [3, 8]}\n")
        gate.upsert_phc_conf(1, 'Pi', 'BHT', 4, 7, path=self.conf)
        self.assertIn("'説明': [1, 2]", self.conf.read_text())
        self.assertEqual((4, 7), gate.read_hand_range(1, 'Pi', 'BHT', path=self.conf))

    def test_invalid_settings_are_not_silently_replaced(self):
        for text in (
            "limits_dict = {\n",
            "limits_dict = {'00001_Pi_bht': [8, 2]}\n",
            "limits_dict = {'00001_Pi_bht': [1, 2], '00001_Pi_bht': [3, 4]}\n",
            "limits_dict = {'00001_Pi_bht': [True, 9]}\n",
        ):
            with self.subTest(text=text):
                self.write(text)
                with self.assertRaises(ValueError):
                    gate.read_hand_range(1, 'Pi', 'BHT', path=self.conf)
                with self.assertRaises(ValueError):
                    gate.upsert_phc_conf(1, 'Pi', 'BHT', 1, 9, path=self.conf)
                self.assertEqual(text, self.conf.read_text())

    def test_dynamic_updates_do_not_silently_use_an_old_declaration(self):
        self.write("limits_dict = {'00001_Pi_bht': [1, 8]}\n"
                   "limits_dict['00001_Pi_bht'] = [3, 7]\n")
        before = self.conf.read_bytes()
        with self.assertRaises(ValueError):
            gate.read_hand_range(1, 'Pi', 'BHT', path=self.conf)
        with self.assertRaises(ValueError):
            gate.upsert_phc_conf(1, 'Pi', 'BHT', 2, 6, path=self.conf)
        self.assertEqual(before, self.conf.read_bytes())

    def test_failed_atomic_save_keeps_original_settings(self):
        self.write("limits_dict = {'00001_Pi_bht': [1, 8]}\n")
        before = self.conf.read_bytes()
        with patch.object(gate.os, 'replace', side_effect=OSError('save failed')):
            with self.assertRaises(OSError):
                gate.upsert_phc_conf(1, 'Pi', 'BHT', 2, 7, path=self.conf)
        self.assertEqual(before, self.conf.read_bytes())
        self.assertEqual([], list(self.folder.glob('*.tmp')))

    def test_noop_preserves_file_timestamp_and_mode(self):
        self.write('limits_dict = {"00001_Pi_bht" : [12, 18]}\n')
        os.chmod(self.conf, 0o640)
        before = self.conf.stat()
        gate.upsert_phc_conf(1, 'Pi', 'BHT', 12, 18, path=self.conf)
        self.assertEqual(before.st_mtime_ns, self.conf.stat().st_mtime_ns)
        gate.upsert_phc_conf(1, 'Pi', 'BHT', 13, 18, path=self.conf)
        self.assertEqual(0o640, self.conf.stat().st_mode & 0o777)

    def test_concurrent_explicit_updates_keep_both_users_changes(self):
        self.write('limits_dict = {}\n')
        def edit(run):
            gate.upsert_phc_conf(run, 'Pi', 'BHT', 10, 20, path=self.conf)
        with ThreadPoolExecutor(max_workers=6) as pool:
            list(pool.map(edit, range(1, 19)))
        for run in range(1, 19):
            self.assertEqual((10, 20), gate.read_hand_range(run, 'Pi', 'BHT', path=self.conf))

    def test_proposal_and_confirmation_do_not_persist_automatic_ranges(self):
        self.write("limits_dict = {'00001_Pi_bht': [12, 18]}\n")
        before = self.conf.read_bytes()
        with patch.object(config, 'OUTPUT_DIR', self.folder/'out'), \
                patch.object(gate, 'PHC_CONF_PATH', self.conf):
            gate.write_proposal(1, {'run': 1, 'suffix': 'Pi',
                                   'proposals': [{'det': 'BHT', 'first': 0, 'last': 50}]})
            gate.set_range_confirmed(1)
            self.assertEqual((12, 18), gate.read_hand_range(1, 'Pi', 'BHT'))
            self.assertEqual([], gate.missing_proposal_dets(1, 'Pi', [{'det': 'BHT'}]))
            gate.clear_range_confirmed(1)
        self.assertEqual(before, self.conf.read_bytes())


if __name__ == '__main__':
    unittest.main()
