"""Run with python3 -m unittest discover -s tests -v (requires PyYAML).

The dispatcher and runner are exercised unchanged in an isolated tree.
Only site configuration, color, and unused ROOT-facing helpers are stubbed.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]


class DispatchTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.scripts = self.root / "calibration/scripts"
        self.scripts.mkdir(parents=True)
        lib = self.root / "common/lib"
        lib.mkdir(parents=True)
        (lib / "__init__.py").touch()
        (lib / "config.py").write_text(
            f"from pathlib import Path\nOUTPUT_DIR = Path({str(self.root / 'out')!r})\n")
        (lib / "input_check.py").write_text("def check_input(*a): return True\n")
        (lib / "param.py").write_text("def backup_file(*a): pass\n")
        (self.root / "termcolor.py").write_text("def colored(text, *a, **kw): return text\n")
        shutil.copy(REPO / "calibration/scripts/calib_dispatch.py", self.scripts)
        shutil.copy(REPO / "common/lib/runner.py", lib)
        self.state = self.root / "out/calib/run00001/calib_dispatch_state.json"

    def run_chain(self, child, steps="  - hodo\n  - dc_resi\n", dry=False):
        (self.scripts / "iterate_hodo_chain.py").write_text(child)
        (self.scripts / "iterate_dc_resi.py").write_text(
            f"from pathlib import Path\nPath({str(self.root / 'next')!r}).touch()\n")
        cfg = self.root / "config.yml"
        cfg.write_text("run: 1\nsteps:\n" + steps)
        return subprocess.run(
            [sys.executable, str(self.scripts / "calib_dispatch.py"), "start", str(cfg)]
            + (["--dry-run"] if dry else []), capture_output=True, text=True,
            env={**os.environ, "PYTHONPATH": str(self.root)}, timeout=15)

    def test_nonzero_exit_records_stop(self):
        for rc in (1, 2, 3, 4, 127, 137):
            with self.subTest(rc=rc):
                result = self.run_chain(f"raise SystemExit({rc})\n")
                self.assertEqual(result.returncode, rc, result.stderr)
                state = json.loads(self.state.read_text())
                self.assertEqual(state['status'], 'stopped')
                self.assertEqual(state['steps'][0]['rc'], rc)
                self.assertFalse(state['finished'])
                self.assertFalse((self.root / 'next').exists())

    def test_interruption_records_stop(self):
        runner = self.root / 'common/lib/runner.py'
        original = runner.read_text()
        for exception in ('KeyboardInterrupt', 'OSError'):
            with self.subTest(exception=exception):
                runner.write_text(original + f'\ndef sh(*a, **kw): raise {exception}\n')
                result = self.run_chain('pass\n')
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertEqual(json.loads(self.state.read_text())['status'], 'stopped')
                self.assertFalse((self.root / 'next').exists())

    def test_success_and_allowed_skip(self):
        for rc in (0, 4):
            with self.subTest(rc=rc):
                result = self.run_chain(f'raise SystemExit({rc})\n',
                                        '  - hodo: {allow_skip: true}\n  - dc_resi\n')
                self.assertEqual(result.returncode, rc, result.stderr)
                self.assertTrue((self.root / 'next').exists())
                self.assertEqual(json.loads(self.state.read_text())['status'], 'finished')

    def test_dry_run_has_no_state(self):
        result = self.run_chain('raise SystemExit(1)\n', dry=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(self.state.exists())
        self.assertFalse((self.root / 'next').exists())


if __name__ == '__main__':
    unittest.main()
