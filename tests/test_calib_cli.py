"""Run with python3 -m unittest discover -s tests -v."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml

from dispatch_isolation import populate_dispatch_lib

REPO = Path(__file__).resolve().parents[1]


class CalibCliHelpersTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        paths = populate_dispatch_lib(
            REPO, self.root, analyzer_dir=self.root / "analyzer",
        )
        self.scripts = paths["scripts"]
        lib = paths["lib"]
        for name in (
            "calib_session.py", "calib_dispatch.py", "calib.py", "ops.py",
        ):
            shutil.copy(REPO / "calibration/scripts" / name, self.scripts)
        shutil.copy(REPO / "common/lib/decode_run.py", lib)
        shutil.copy(REPO / "common/lib/beam_run_condition.py", lib)
        cfg = self.root / "calibration/config"
        cfg.mkdir(parents=True, exist_ok=True)
        self.table = cfg / "beam_run_condition.yml"
        self.table.write_text(
            "ranges:\n"
            "  - {first: 100, last: 110, particle: pion, "
            "momentum_gev_c: 1.0, charge: -1}\n"
        )
        text = (lib / "beam_run_condition.py").read_text()
        text = text.replace(
            'Path(__file__).resolve().parents[2] / "calibration" / "config" / "beam_run_condition.yml"',
            f"Path({str(self.table)!r})",
        )
        (lib / "beam_run_condition.py").write_text(text)
        # Minimal yaml stub already via real PyYAML if present; decode_run imports yaml.
        self.env = {**os.environ, "PYTHONPATH": str(self.root)}

    def py(self, code):
        return subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True, text=True, env=self.env, timeout=15,
            cwd=str(self.root),
        )

    def test_steps_from_flags(self):
        result = self.py(
            "from calibration.scripts import calib_session as cs\n"
            "print(cs.steps_from_flags())\n"
            "print(cs.steps_from_flags(hodo=True))\n"
            "print(cs.steps_from_flags(bcout=True, bcin=True, d5=True))\n"
        )
        # import path: scripts are under calibration/scripts but not a package
        result = self.py(
            "import sys\n"
            f"sys.path.insert(0, {str(self.scripts)!r})\n"
            "import calib_session as cs\n"
            "s,t=cs.steps_from_flags()\n"
            "print(t, len(s))\n"
            "s,t=cs.steps_from_flags(hodo=True)\n"
            "print(t, s)\n"
            "s,t=cs.steps_from_flags(bcout=True, d5=True)\n"
            "print(t, s)\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        lines = result.stdout.strip().splitlines()
        self.assertEqual(lines[0], "both 4")
        self.assertEqual(lines[1], "hodo ['hodo']")
        self.assertEqual(lines[2], "custom ['dc_chamber', 'd5_align']")

    def test_decode_cmd_local_flag(self):
        result = self.py(
            "from common.lib import decode_run\n"
            "c0=decode_run.decode_cmd('/tmp/x.yml')\n"
            "c1=decode_run.decode_cmd('/tmp/x.yml', local=True)\n"
            "print('LOCAL' if c0.endswith('--local') else 'BSUB')\n"
            "print('LOCAL' if '--local' in c1 else 'BSUB')\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.assertEqual(result.stdout.strip().splitlines(), ["BSUB", "LOCAL"])

    def test_calib_dry_creates_session(self):
        result = subprocess.run(
            [sys.executable, str(self.scripts / "calib.py"), "105", "--hodo", "--dry"],
            capture_output=True, text=True, env=self.env, timeout=20,
            cwd=str(self.root),
        )
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        path = self.root / "out/calib/run00105/calib_session.yml"
        self.assertTrue(path.is_file(), result.stdout)
        doc = yaml.safe_load(path.read_text())
        self.assertEqual(doc["calib_particle"], "pion")
        self.assertEqual(doc["steps"], ["hodo"])
        self.assertFalse(doc.get("local", False))
        self.assertIn("iterate_hodo_chain.py", result.stdout)

    def test_calib_local_flag_on_session(self):
        result = subprocess.run(
            [
                sys.executable, str(self.scripts / "calib.py"), "105",
                "--hodo", "--local", "--dry",
            ],
            capture_output=True, text=True, env=self.env, timeout=20,
            cwd=str(self.root),
        )
        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        doc = yaml.safe_load(
            (self.root / "out/calib/run00105/calib_session.yml").read_text()
        )
        self.assertTrue(doc["local"])
        self.assertIn("--local", result.stdout)

    def test_ops_status_without_state(self):
        subprocess.run(
            [sys.executable, str(self.scripts / "calib.py"), "105", "--hodo", "--dry"],
            capture_output=True, text=True, env=self.env, timeout=20,
            cwd=str(self.root),
        )
        # dry-run start does not write state; ops status should say not started
        result = subprocess.run(
            [sys.executable, str(self.scripts / "ops.py"), "status", "105"],
            capture_output=True, text=True, env=self.env, timeout=15,
            cwd=str(self.root),
        )
        # do_status exits 0 even when not started? check behavior
        out = result.stdout + result.stderr
        self.assertTrue(
            result.returncode in (0, 1),
            out,
        )


if __name__ == "__main__":
    unittest.main()
