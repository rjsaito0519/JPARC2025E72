"""Run with python3 -m unittest discover -s tests -v.

Exercises iterate_dc_chamber / fit-loop helpers with stubbed run_dcprm and decode.
"""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]


class DcChamberTests(unittest.TestCase):
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
            "from pathlib import Path\n"
            f"OUTPUT_DIR = Path({str(self.root / 'out')!r})\n"
            f"DATA_DIR = Path({str(self.root / 'data')!r})\n"
            f"ANALYZER_DIR = Path({str(self.root / 'analyzer')!r})\n"
            f"PARAM_DIR = Path({str(self.root / 'param')!r})\n"
            "SUB_DIR = 'e72'\n"
        )
        (lib / "input_check.py").write_text(
            "def stale_params(*a, **k): return []\n"
            "def check_input(*a, **k): return True\n"
        )
        (lib / "param.py").write_text("def backup_file(*a): pass\n")
        (lib / "decode_run.py").write_text(
            "from pathlib import Path\n"
            "def run_decode(runlist, dry_run=False, local=False):\n"
            f"    Path({str(self.root / 'decoded')!r}).write_text("
            "str(runlist) + (' dry' if dry_run else '') + '\\n', encoding='utf-8')\n"
        )
        (lib / "beam_run_condition.py").write_text(
            "from pathlib import Path\n"
            "def user_param_paths(run):\n"
            "    return [Path('/no/such')]\n"
            "def resolve(run):\n"
            "    return {'momentum_gev_c': 1.0}\n"
        )
        shutil.copy(REPO / "common/lib/calib_runlist.py", lib)
        (self.root / "termcolor.py").write_text(
            "def colored(text, *a, **kw): return text\n")
        shutil.copy(REPO / "calibration/scripts/iterate_dc_resi.py", self.scripts)
        shutil.copy(REPO / "calibration/scripts/iterate_dc_reso.py", self.scripts)
        shutil.copy(REPO / "calibration/scripts/iterate_dc_chamber.py", self.scripts)
        shutil.copy(REPO / "common/lib/runner.py", lib)
        # Stub create_runlist used by calib_runlist
        create = self.scripts / "create_runlist.py"
        create.write_text(
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "raise SystemExit(0)\n"
        )
        create.chmod(0o755)
        myrun = self.root / "analyzer/runmanager/myrun.py"
        myrun.parent.mkdir(parents=True)
        myrun.write_text("#!/usr/bin/env python3\n")
        myrun.chmod(0o755)
        runlist = self.root / "analyzer/runmanager/runlist/e72/bcout_run00001.yml"
        runlist.parent.mkdir(parents=True)
        runlist.write_text("run: 1\n")
        (self.root / "data").mkdir()
        (self.root / "data/run00001_BcOut.root").write_bytes(b"root")
        self.env = {**os.environ, "PYTHONPATH": str(self.root)}
        self.runlist = runlist

    def write_dcprm(self, body):
        path = self.scripts / "run_dcprm.py"
        path.write_text(body)
        path.chmod(0o755)

    def run_chamber(self, *extra):
        return subprocess.run(
            [sys.executable, str(self.scripts / "iterate_dc_chamber.py"), "1", *extra],
            capture_output=True, text=True, env=self.env, timeout=20,
            cwd=str(self.root))

    def test_dry_run_prints_full_chain(self):
        self.write_dcprm("raise SystemExit(0)\n")
        result = self.run_chamber("--dry-run")
        self.assertEqual(result.returncode, 0, result.stderr)
        out = result.stdout + result.stderr
        self.assertIn("create_runlist", out)
        self.assertIn("Pi_tdc", out)
        self.assertIn("Pi_drift", out)
        self.assertIn("Pi_resi", out)
        self.assertIn(" t0 ", out)
        self.assertIn(" drift ", out)
        self.assertIn(" resi ", out)
        self.assertIn(" reso ", out)

    def test_t0_nonzero_stops_before_drift(self):
        log = self.root / "dcprm_log.txt"
        self.write_dcprm(
            "import sys\n"
            f"log = open({str(log)!r}, 'a')\n"
            "log.write(sys.argv[2] + '\\n'); log.close()\n"
            "raise SystemExit(1 if sys.argv[2] == 't0' else 0)\n"
        )
        result = self.run_chamber()
        self.assertEqual(result.returncode, 1, result.stderr)
        modes = log.read_text().split()
        self.assertEqual(modes, ["t0"])
        # create_runlist stage still decoded before t0
        self.assertTrue((self.root / "decoded").exists())

    def test_success_runs_all_modes_once_when_converged(self):
        log = self.root / "dcprm_log.txt"
        self.write_dcprm(
            "import sys\n"
            f"open({str(log)!r}, 'a').write(sys.argv[2] + '\\n')\n"
            "raise SystemExit(0)\n"
        )
        result = self.run_chamber()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(log.read_text().split(), ["t0", "drift", "resi", "reso"])
        self.assertTrue((self.root / "decoded").exists())


if __name__ == "__main__":
    unittest.main()
