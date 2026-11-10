"""Unit tests for common.lib.calib_runlist (tempfile isolation)."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]


class CalibRunlistTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        lib = self.root / "common/lib"
        lib.mkdir(parents=True)
        (lib / "__init__.py").touch()
        (lib / "config.py").write_text(
            "from pathlib import Path\n"
            f"PARAM_DIR = Path({str(self.root / 'param')!r})\n"
            f"ANALYZER_DIR = Path({str(self.root / 'analyzer')!r})\n"
            "SUB_DIR = 'e72'\n"
        )
        (lib / "runner.py").write_text(
            "EXIT_OK=0\nEXIT_FAIL=1\nEXIT_CHECK=3\n"
            "def info(m): print(m)\n"
            "def warn(m): print(m)\n"
            "def step(m): print(m)\n"
            "def fail(m):\n"
            "    raise SystemExit(f'FAIL:{m}')\n"
            "def sh(cmd, ok=(0,)):\n"
            "    return 0, ''\n"
        )
        (lib / "beam_run_condition.py").write_text(
            "from pathlib import Path\n"
            f"_UP = Path({str(self.root / 'param/USER/e72')!r})\n"
            "def user_param_paths(run):\n"
            "    return [_UP / f'UserParam_run{run:05d}']\n"
            "def resolve(run):\n"
            "    return {'momentum_gev_c': 0.735}\n"
        )
        shutil.copy(REPO / "common/lib/calib_runlist.py", lib)
        scripts = self.root / "calibration/scripts"
        scripts.mkdir(parents=True)
        (scripts / "create_runlist.py").write_text("raise SystemExit(0)\n")
        user_dir = self.root / "param/USER/e72"
        user_dir.mkdir(parents=True)
        tpl_dir = self.root / "param/USER"
        self.tpl = tpl_dir / "UserParam_e72_20251104"
        self.tpl.write_text(
            "BHT_TDC 100 200\n"
            "BHT_TDC_Pi 100 150\n"
            "BHT_TDC_K 150 200\n"
            "BH2_TDC 10 20\n"
            "BAC_TDC 10 20\n"
            "HTOF_TDC 10 20\n"
            "KVC_TDC 10 20\n"
            "T1_TDC 10 20\n"
        )
        self.user = user_dir / "UserParam_run00001"
        self.env = {**os.environ, "PYTHONPATH": str(self.root)}

    def py(self, code):
        return subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True, text=True, env=self.env, timeout=15,
            cwd=str(self.root),
        )

    def test_suffixes(self):
        result = self.py(
            "from common.lib import calib_runlist as cr\n"
            "print(cr.hodo_suffix(False, 'hdprm'))\n"
            "print(cr.hodo_suffix(True, 't0'))\n"
            "print(cr.hodo_suffix(False, 0))\n"
            "print(cr.dc_suffix(False, 'tdc'))\n"
            "print(cr.d5_suffix(True))\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            result.stdout.strip().splitlines(),
            ["Pi_hdprm", "K_t0", "0", "Pi_tdc", "K"],
        )

    def test_windows_match_template_not_ready(self):
        self.user.write_text(self.tpl.read_text())
        result = self.py(
            "from common.lib import calib_runlist as cr\n"
            "ready, path, detail = cr.tdc_windows_ready(1)\n"
            "print(ready)\n"
            "print(detail)\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.strip().splitlines()
        self.assertEqual(lines[0], "False")
        self.assertIn("match template", lines[1])

    def test_windows_differ_ready(self):
        text = self.tpl.read_text().replace("BH2_TDC 10 20", "BH2_TDC 11 21")
        self.user.write_text(text)
        result = self.py(
            "from common.lib import calib_runlist as cr\n"
            "ready, path, detail = cr.tdc_windows_ready(1)\n"
            "print(ready)\n"
            "print(detail)\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.strip().splitlines()
        self.assertEqual(lines[0], "True")
        self.assertIn("differ", lines[1])

    def test_missing_userparam_not_ready(self):
        result = self.py(
            "from common.lib import calib_runlist as cr\n"
            "ready, path, detail = cr.tdc_windows_ready(1)\n"
            "print(ready)\n"
            "print(detail)\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.strip().splitlines()
        self.assertEqual(lines[0], "False")
        self.assertIn("missing", lines[1].lower())

    def test_create_runlist_cmd_order(self):
        result = self.py(
            "from common.lib import calib_runlist as cr\n"
            "cmd = cr.create_runlist_cmd('--hodo', 2977, 'Pi_hdprm')\n"
            "print(cmd)\n"
            "cmd_d5 = cr.create_runlist_cmd('--d5', 2580, 'K', mom=0.735)\n"
            "print(cmd_d5)\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.strip().splitlines()
        cmd = lines[0]
        self.assertIn("--hodo", cmd)
        self.assertLess(cmd.index("--hodo"), cmd.index("2977"))
        self.assertLess(cmd.index("2977"), cmd.index("Pi_hdprm"))
        self.assertIn("--mom", lines[1])
        self.assertIn("0.735", lines[1])


if __name__ == "__main__":
    unittest.main()
