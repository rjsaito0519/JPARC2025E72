"""Dry-run / window-gate checks for iterate_hodo_chain (tempfile isolation)."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]


class HodoChainRunlistTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.scripts = self.root / "calibration/scripts"
        self.scripts.mkdir(parents=True)
        (self.root / "common").mkdir(parents=True)
        (self.root / "common/__init__.py").touch()
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
        (self.root / "termcolor.py").write_text(
            "def colored(text, *a, **kw): return text\n"
        )
        (lib / "decode_run.py").write_text(
            "def run_decode(runlist, dry_run=False, local=False):\n"
            "    print(f'decode {runlist} dry={dry_run} local={local}')\n"
            "    return 0\n"
        )
        for name in (
            "hodo_btof_monitor.py", "hodo_hdphc_rounds.py", "hodo_qa_pdf.py",
        ):
            (lib / name).write_text("pass\n")
        (lib / "beam_run_condition.py").write_text(
            "from pathlib import Path\n"
            f"_UP = Path({str(self.root / 'param/USER/e72')!r})\n"
            "def user_param_paths(run):\n"
            "    return [_UP / f'UserParam_run{run:05d}']\n"
            "def resolve(run):\n"
            "    return {'momentum_gev_c': 1.0}\n"
        )
        shutil.copy(REPO / "common/lib/calib_runlist.py", lib)
        shutil.copy(REPO / "common/lib/runner.py", lib)
        shutil.copy(REPO / "calibration/scripts/iterate_hodo_chain.py", self.scripts)
        (self.scripts / "run_hodo.py").write_text("raise SystemExit(0)\n")
        (self.scripts / "create_runlist.py").write_text("raise SystemExit(0)\n")
        myrun = self.root / "analyzer/runmanager/myrun.py"
        myrun.parent.mkdir(parents=True)
        myrun.write_text("#!/usr/bin/env python3\n")
        myrun.chmod(0o755)
        user_dir = self.root / "param/USER/e72"
        user_dir.mkdir(parents=True)
        tpl = self.root / "param/USER/UserParam_e72_20251104"
        tpl.write_text(
            "BHT_TDC 100 200\nBHT_TDC_Pi 100 150\nBHT_TDC_K 150 200\n"
            "BH2_TDC 10 20\nBAC_TDC 10 20\nHTOF_TDC 10 20\n"
            "KVC_TDC 10 20\nT1_TDC 10 20\n"
        )
        self.tpl = tpl
        self.user = user_dir / "UserParam_run00007"
        self.env = {**os.environ, "PYTHONPATH": str(self.root)}

    def run_chain(self, *extra):
        return subprocess.run(
            [sys.executable, str(self.scripts / "iterate_hodo_chain.py"),
             "7", *extra],
            capture_output=True, text=True, env=self.env, timeout=20,
            cwd=str(self.root),
        )

    def test_dry_run_window_unset_stops_check_message(self):
        self.user.write_text(self.tpl.read_text())
        result = self.run_chain("--dry-run")
        self.assertEqual(result.returncode, 0, result.stderr)
        out = result.stdout + result.stderr
        self.assertIn("0", out)  # suffix 0
        self.assertIn("EXIT_CHECK", out)
        self.assertIn("create_runlist", out)

    def test_dry_run_ready_prints_stage_order(self):
        text = self.tpl.read_text().replace("T1_TDC 10 20", "T1_TDC 30 40")
        self.user.write_text(text)
        result = self.run_chain("--dry-run", "--skip-hdphc")
        self.assertEqual(result.returncode, 0, result.stderr)
        out = result.stdout + result.stderr
        self.assertIn("Pi_hdprm", out)
        self.assertIn("Pi_t0", out)
        self.assertIn("hdprm", out)
        # hdphc skipped
        self.assertNotIn("Pi_hdphc", out)

    def test_window_unset_exits_check(self):
        self.user.write_text(self.tpl.read_text())
        result = self.run_chain("--skip-window-bootstrap")
        self.assertEqual(result.returncode, 3, result.stderr + result.stdout)
        out = result.stdout + result.stderr
        self.assertIn("TDC window", out)


if __name__ == "__main__":
    unittest.main()
