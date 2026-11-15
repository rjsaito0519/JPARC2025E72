"""A newly merged ROOT can still contain chunks decoded with different parameters."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from common.lib import decode_run


class DecodeInputChangesTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.directory = Path(temporary.name)
        self.parameter = self.directory/'UserParam'
        self.parameter.write_text('calibration v1\n')
        self.conf = self.directory/'analyzer.conf'
        self.conf.write_text('USER: UserParam\n')
        self.root = self.directory/'run00012_Hodo.root'
        self.runlist = self.directory/'hodo.yml'
        self.runlist.write_text(f'WORKDIR: {self.directory}\nDEFAULT:\n  conf: analyzer.conf\nRUN:\n  run00012:\n    root: {self.root}\n')
        stat_dir = self.directory/'stat'
        stat_dir.mkdir()
        stat_patch = patch.object(decode_run, 'STAT_DIR', stat_dir)
        stat_patch.start()
        self.addCleanup(stat_patch.stop)
        (stat_dir/'hodo.json').write_text(json.dumps({'run00012':{'stat':'DONE'}}))
        self.calls = 0

    def decode(self, change=False):
        self.calls += 1
        if change:
            self.parameter.write_text(f'calibration v{self.calls+1}\n')
        self.root.write_text(self.parameter.read_text())
        return 0, ''

    def test_change_before_merge_is_detected_even_with_done_and_fresh_root(self):
        def myrun(*args, **kwargs):
            return self.decode(change=self.calls==0)
        with patch.object(decode_run.runner,'sh',side_effect=myrun):
            self.assertEqual(0,decode_run.run_decode(self.runlist))
        self.assertEqual(2,self.calls)
        self.assertEqual(self.parameter.read_text(),self.root.read_text())

    def test_continuously_changing_inputs_stop_before_fit(self):
        with patch.object(decode_run.runner,'sh',side_effect=lambda *args,**kwargs:self.decode(change=True)):
            with self.assertRaises(SystemExit):
                decode_run.run_decode(self.runlist)
        self.assertEqual(3,self.calls)

    def test_stable_inputs_only_decode_once(self):
        with patch.object(decode_run.runner,'sh',side_effect=lambda *args,**kwargs:self.decode()) as shell:
            self.assertEqual(0,decode_run.run_decode(self.runlist,local=True))
        self.assertEqual(1,self.calls)
        self.assertIn('--local',shell.call_args.args[0])

    def test_dry_run_does_not_read_inputs_or_start_jobs(self):
        self.runlist.unlink()
        with patch.object(decode_run.runner,'sh') as shell,patch.object(decode_run,'calibration_inputs') as snapshot:
            self.assertEqual(0,decode_run.run_decode(self.runlist,dry_run=True))
        shell.assert_not_called()
        snapshot.assert_not_called()


if __name__ == '__main__':
    unittest.main()
