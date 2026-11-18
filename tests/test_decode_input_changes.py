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

    def d5_inputs(self):
        self.conf.write_text(
            'USER: UserParam\nHDPRM: HodoParam\nHDPHC: HodoPHCParam\n'
            'DCGEO: DCGeomParam\nD5MTX: D5Matrix\n'
        )
        for name in ('HodoParam', 'HodoPHCParam', 'DCGeomParam', 'D5Matrix'):
            (self.directory/name).write_text('v1\n')
        text = self.runlist.read_text()
        if '  bin: ./bin/D5Tracking\n' not in text:
            text = text.replace('DEFAULT:\n', 'DEFAULT:\n  bin: ./bin/D5Tracking\n')
        self.runlist.write_text(text)

    def test_d5_hodo_parameter_changes_do_not_repeat_an_independent_decode(self):
        self.d5_inputs()
        def myrun(*args, **kwargs):
            for name in ('HodoParam', 'HodoPHCParam'):
                (self.directory/name).write_text('v2\n')
            return self.decode()
        with patch.object(decode_run.runner, 'sh', side_effect=myrun):
            self.assertEqual(0, decode_run.run_decode(self.runlist))
        self.assertEqual(1, self.calls)

    def test_d5_tracking_inputs_still_require_a_fresh_decode(self):
        for name in ('UserParam', 'DCGeomParam', 'D5Matrix'):
            with self.subTest(parameter=name):
                self.d5_inputs()
                self.calls = 0
                def myrun(*args, **kwargs):
                    if self.calls == 0:
                        (self.directory/name).write_text('v2\n')
                    return self.decode()
                with patch.object(decode_run.runner, 'sh', side_effect=myrun):
                    self.assertEqual(0, decode_run.run_decode(self.runlist))
                self.assertEqual(2, self.calls)

    def test_mixed_d5_and_hodo_runlist_retains_hodo_dependencies(self):
        self.d5_inputs()
        with self.runlist.open('a') as stream:
            stream.write('  hodo:\n    bin: ./bin/Hodoscope\n'
                         '    conf: analyzer.conf\n    root: hodo.root\n')
        before = decode_run.calibration_inputs(self.runlist)
        (self.directory/'HodoPHCParam').write_text('v2\n')
        after = decode_run.calibration_inputs(self.runlist)
        self.assertNotEqual(before, after)

    def test_dry_run_does_not_read_inputs_or_start_jobs(self):
        self.runlist.unlink()
        with patch.object(decode_run.runner,'sh') as shell,patch.object(decode_run,'calibration_inputs') as snapshot:
            self.assertEqual(0,decode_run.run_decode(self.runlist,dry_run=True))
        shell.assert_not_called()
        snapshot.assert_not_called()


if __name__ == '__main__':
    unittest.main()
