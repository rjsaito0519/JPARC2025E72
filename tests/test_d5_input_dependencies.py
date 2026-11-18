"""D5 freshness guards track DC/beam inputs without serializing Hodo PHC."""
import os
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

from common.lib import calib_anomaly, input_check, runner


class D5InputDependenciesTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.directory = Path(directory.name)
        self.root = self.directory/'run00012_D5_Pi.root'
        self.root.write_text('decoded')
        self.now = time.time()
        os.utime(self.root, (self.now, self.now))
        self.parameters = {}
        keys = (*input_check.D5_TRACKING_PARAM_KEYS, 'HDPRM', 'HDPHC')
        for key in keys:
            path = self.directory/key
            path.write_text('v1')
            os.utime(path, (self.now-20, self.now-20))
            self.parameters[key] = path
        conf = self.directory/'analyzer.conf'
        conf.write_text(''.join(f'{key}: {path}\n'
                                for key, path in self.parameters.items()))
        (self.directory/'d5.yml').write_text(
            f'RUN:\n  d5:\n    conf: {conf}\n    root: {self.root}\n'
        )
        context = patch.object(input_check, 'RUNLIST_DIR', self.directory)
        context.start()
        self.addCleanup(context.stop)

    def touch_parameter(self, key):
        path = self.parameters[key]
        path.write_text('v2')
        os.utime(path, (self.now+20, self.now+20))

    def checked(self):
        with patch.dict(os.environ, {'E72_CALIB_DISPATCH': '1'}), \
                patch.object(calib_anomaly, 'remember_check_reason'):
            return runner.checked_root(
                self.root, 'Pi', False,
                keys=input_check.D5_TRACKING_PARAM_KEYS,
            )

    def test_hodo_update_does_not_make_d5_root_stale(self):
        self.touch_parameter('HDPRM')
        self.touch_parameter('HDPHC')
        self.assertEqual(self.root, self.checked())

    def test_actual_d5_inputs_still_block_fit(self):
        for key in input_check.D5_TRACKING_PARAM_KEYS:
            with self.subTest(parameter=key):
                for path in self.parameters.values():
                    os.utime(path, (self.now-20, self.now-20))
                self.touch_parameter(key)
                with self.assertRaises(SystemExit) as exit_result:
                    self.checked()
                self.assertEqual(runner.EXIT_CHECK, exit_result.exception.code)

    def test_default_hodo_guard_is_preserved(self):
        self.touch_parameter('HDPHC')
        self.assertTrue(input_check.stale_params(self.root, 'Pi'))
        self.assertFalse(input_check.stale_params(
            self.root, 'Pi', keys=input_check.D5_TRACKING_PARAM_KEYS,
        ))


if __name__ == '__main__':
    unittest.main()
