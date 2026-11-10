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
        self.cfg = self.root / "config.yml"
        self.env = {**os.environ, "PYTHONPATH": str(self.root)}

    def write_cfg(self, steps="  - hodo\n  - dc_resi\n"):
        self.cfg.write_text("run: 1\nsteps:\n" + steps)

    def dispatch(self, command, dry=False):
        return subprocess.run(
            [sys.executable, str(self.scripts / "calib_dispatch.py"), command, str(self.cfg)]
            + (["--dry-run"] if dry else []),
            capture_output=True, text=True, env=self.env, timeout=15)

    def run_chain(self, child, steps="  - hodo\n  - dc_resi\n", dry=False):
        (self.scripts / "iterate_hodo_chain.py").write_text(child)
        (self.scripts / "iterate_dc_resi.py").write_text(
            f"from pathlib import Path\nPath({str(self.root / 'next')!r}).touch()\n")
        self.write_cfg(steps)
        return self.dispatch("start", dry=dry)

    def test_nonzero_exit_records_stop(self):
        for rc in (1, 2, 3, 4, 127, 137):
            with self.subTest(rc=rc):
                result = self.run_chain(f"raise SystemExit({rc})\n")
                self.assertEqual(result.returncode, rc, result.stderr)
                state = json.loads(self.state.read_text())
                self.assertEqual(state['status'], 'stopped')
                self.assertEqual(state['steps'][0]['rc'], rc)
                self.assertFalse(state['finished'])
                self.assertFalse(state.get('stop_requested'))
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

    def test_stop_request_pauses_before_next_step(self):
        state_path = str(self.state)
        child = f"""
import json
from pathlib import Path
p = Path({state_path!r})
st = json.loads(p.read_text())
st['stop_requested'] = True
p.write_text(json.dumps(st) + '\\n')
"""
        result = self.run_chain(child)
        self.assertEqual(result.returncode, 0, result.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'stopped')
        self.assertEqual(state['step'], 'hodo')
        self.assertFalse(state['finished'])
        self.assertFalse(state.get('stop_requested'))
        self.assertFalse((self.root / 'next').exists())

    def test_resume_retries_failed_step(self):
        first = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(first.returncode, 1, first.stderr)
        self.assertEqual(json.loads(self.state.read_text())['status'], 'stopped')
        (self.scripts / "iterate_hodo_chain.py").write_text('pass\n')
        (self.scripts / "iterate_dc_resi.py").write_text(
            f"from pathlib import Path\nPath({str(self.root / 'next')!r}).touch()\n")
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'finished')
        self.assertTrue((self.root / 'next').exists())
        self.assertEqual([e['index'] for e in state['steps']], [0, 0, 1])

    def test_resume_after_cooperative_stop_continues_next(self):
        state_path = str(self.state)
        child = f"""
import json
from pathlib import Path
p = Path({state_path!r})
st = json.loads(p.read_text())
st['stop_requested'] = True
p.write_text(json.dumps(st) + '\\n')
"""
        first = self.run_chain(child)
        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertFalse((self.root / 'next').exists())
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'finished')
        self.assertTrue((self.root / 'next').exists())
        self.assertEqual([e['index'] for e in state['steps']], [0, 1])

    def test_resume_rejects_finished_and_not_started(self):
        self.write_cfg()
        missing = self.dispatch("resume")
        self.assertEqual(missing.returncode, 1, missing.stderr)
        ok = self.run_chain('pass\n')
        self.assertEqual(ok.returncode, 0, ok.stderr)
        self.assertEqual(json.loads(self.state.read_text())['status'], 'finished')
        again = self.dispatch("resume")
        self.assertEqual(again.returncode, 1, again.stderr)

    def test_stop_cli_sets_flag_and_ignores_finished(self):
        running = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(running.returncode, 1)
        flag = self.dispatch("stop")
        self.assertEqual(flag.returncode, 0, flag.stderr)
        self.assertTrue(json.loads(self.state.read_text())['stop_requested'])
        (self.scripts / "iterate_hodo_chain.py").write_text('pass\n')
        (self.scripts / "iterate_dc_resi.py").write_text('pass\n')
        done = self.dispatch("start")
        self.assertEqual(done.returncode, 0, done.stderr)
        ignored = self.dispatch("stop")
        self.assertEqual(ignored.returncode, 0, ignored.stderr)
        self.assertEqual(json.loads(self.state.read_text())['status'], 'finished')

    def test_start_records_plan(self):
        result = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(result.returncode, 1, result.stderr)
        plan = json.loads(self.state.read_text())['plan']
        self.assertEqual([e['name'] for e in plan['steps']], ['hodo', 'dc_resi'])
        self.assertEqual(plan['calib_particle'], 'pion')
        self.assertFalse(plan['kaon'])
        self.assertFalse(plan['force'])
        self.assertEqual(plan['max_iter'], 10)

    def test_resume_rejects_steps_change(self):
        first = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(first.returncode, 1, first.stderr)
        before = self.state.read_text()
        (self.scripts / "iterate_dc_chamber.py").write_text('pass\n')
        (self.scripts / "iterate_dc_resi.py").write_text(
            f"from pathlib import Path\nPath({str(self.root / 'next')!r}).touch()\n")
        self.write_cfg('  - dc_chamber\n  - dc_resi\n')
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 1, resumed.stdout + resumed.stderr)
        self.assertIn('plan mismatch', resumed.stdout)
        self.assertEqual(self.state.read_text(), before)
        self.assertFalse((self.root / 'next').exists())

    def test_resume_rejects_particle_change(self):
        first = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(first.returncode, 1, first.stderr)
        before = self.state.read_text()
        self.cfg.write_text('run: 1\nkaon: true\nsteps:\n  - hodo\n  - dc_resi\n')
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 1, resumed.stdout + resumed.stderr)
        self.assertIn('plan mismatch', resumed.stdout)
        self.assertEqual(self.state.read_text(), before)

    def test_resume_rejects_missing_plan(self):
        first = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(first.returncode, 1, first.stderr)
        state = json.loads(self.state.read_text())
        del state['plan']
        self.state.write_text(json.dumps(state) + '\n')
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 1, resumed.stdout + resumed.stderr)
        self.assertIn('no plan fingerprint', resumed.stdout)
        self.assertEqual(json.loads(self.state.read_text())['status'], 'stopped')
        self.assertNotIn('plan', json.loads(self.state.read_text()))


if __name__ == '__main__':
    unittest.main()
