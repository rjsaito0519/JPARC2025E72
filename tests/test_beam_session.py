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


class BeamSessionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        cfg = self.root / "calibration/config"
        cfg.mkdir(parents=True)
        self.table = cfg / "beam_run_condition.yml"
        self.table.write_text(
            "ranges:\n"
            "  - {first: 100, last: 110, particle: pion, momentum_gev_c: 1.0, charge: -1}\n"
            "  - {first: 200, last: 200, particle: kaon, momentum_gev_c: 0.9, charge: -1}\n"
        )
        paths = populate_dispatch_lib(
            REPO, self.root, analyzer_dir=self.root / "analyzer",
        )
        self.scripts = paths["scripts"]
        lib = paths["lib"]
        shutil.copy(REPO / "common/lib/beam_run_condition.py", lib)
        shutil.copy(REPO / "calibration/scripts/calib_session.py", self.scripts)
        # Force table path to the fixture file.
        text = (lib / "beam_run_condition.py").read_text()
        text = text.replace(
            'Path(__file__).resolve().parents[2] / "calibration" / "config" / "beam_run_condition.yml"',
            f"Path({str(self.table)!r})",
        )
        (lib / "beam_run_condition.py").write_text(text)
        self.env = {**os.environ, "PYTHONPATH": str(self.root)}

    def py(self, code):
        return subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True, text=True, env=self.env, timeout=15,
            cwd=str(self.root))

    def session_cmd(self, *args):
        return subprocess.run(
            [sys.executable, str(self.scripts / "calib_session.py"), *args],
            capture_output=True, text=True, env=self.env, timeout=15,
            cwd=str(self.root))

    def test_lookup_table_hit_and_miss(self):
        result = self.py(
            "from common.lib import beam_run_condition as brc\n"
            "h=brc.resolve(105); m=brc.resolve(9999)\n"
            "print(h['particle'], h['momentum_gev_c'], h['momentum_source'])\n"
            "print(m['momentum_gev_c'], m['momentum_source'])\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.strip().splitlines()
        self.assertEqual(lines[0], "pion 1.0 run_table")
        self.assertEqual(lines[1], "None none")

    def test_user_param_overrides_momentum(self):
        path = self.root / "param/USER/e72/UserParam_run00105"
        path.parent.mkdir(parents=True)
        path.write_text("BeamMom 1.234\nBeamParticle 1\n")
        result = self.py(
            "from common.lib import beam_run_condition as brc\n"
            "h=brc.resolve(105)\n"
            "print(h['particle'], h['particle_source'], "
            "h['momentum_gev_c'], h['momentum_source'])\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "kaon user_param 1.234 user_param")

    def test_init_fails_without_momentum(self):
        result = self.session_cmd(
            "init", "9999", "--calib-particle", "pion", "--target", "hodo")
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn("momentum unknown", result.stderr + result.stdout)

    def test_init_with_momentum_and_write_table(self):
        result = self.session_cmd(
            "init", "9999", "--calib-particle", "pion", "--target", "hodo",
            "--momentum", "0.8", "--beam-particle", "pion", "--charge", "-1",
            "--write-table")
        self.assertEqual(result.returncode, 0, result.stderr)
        doc = yaml.safe_load(
            (self.root / "out/calib/run09999/calib_session.yml").read_text())
        self.assertEqual(doc["momentum_gev_c"], 0.8)
        self.assertEqual(doc["momentum_source"], "user_input")
        self.assertTrue(doc["table_updated"])
        table = yaml.safe_load(self.table.read_text())
        self.assertTrue(any(r["first"] == 9999 for r in table["ranges"]))

    def test_write_table_conflict_in_range(self):
        result = self.session_cmd(
            "init", "105", "--calib-particle", "pion", "--target", "hodo",
            "--momentum", "2.0", "--beam-particle", "pion", "--charge", "-1",
            "--write-table")
        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn("inside range", result.stderr + result.stdout)

    def test_dispatch_reads_session_kaon(self):
        (self.scripts / "iterate_hodo_chain.py").write_text("raise SystemExit(0)\n")
        session = self.root / "session.yml"
        session.write_text(
            "run: 1\ncalib_particle: kaon\nkaon: true\nmax_iter: 2\n"
            "momentum_gev_c: 0.9\nmomentum_source: run_table\n"
            "steps:\n  - hodo\n"
        )
        result = subprocess.run(
            [sys.executable, str(self.scripts / "calib_dispatch.py"),
             "start", str(session), "--dry-run"],
            capture_output=True, text=True, env=self.env, timeout=15)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--kaon", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
