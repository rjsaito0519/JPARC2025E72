"""Run with python3 -m unittest discover -s tests -v."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

REPO = Path(__file__).resolve().parents[1]


class D5QaPathTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        lib = self.root / "common/lib"
        lib.mkdir(parents=True)
        (lib / "__init__.py").touch()
        self.out = self.root / "out"
        (lib / "config.py").write_text(
            "from pathlib import Path\n"
            f"OUTPUT_DIR = Path({str(self.out)!r})\n"
        )
        for name in ("d5_qa_slack.py", "local_notify_config.py", "slack_notify.py",
                     "hodo_qa_pdf.py"):
            shutil.copy(REPO / "common/lib" / name, lib)
        self.env = {**os.environ, "PYTHONPATH": str(self.root)}

    def py(self, code):
        return subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True, text=True, env=self.env, timeout=30,
            cwd=str(self.root))

    def test_path_names(self):
        result = self.py(
            "from common.lib import d5_qa_slack as m\n"
            "print(m.source_path(2580, 'pos_review', kaon=True).name)\n"
            "print(m.source_path(2489, 'plane_resi_qa').name)\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.strip().splitlines()
        self.assertEqual(lines[0], "d5_pos_review_run02580_K.pdf")
        self.assertEqual(lines[1], "d5_plane_resi_qa_run02489_Pi.pdf")

    def test_dry_run(self):
        result = self.py(
            "from common.lib import d5_qa_slack as m\n"
            "ok, msg = m.post_d5_qa(1, 'pos_review', dry_run=True)\n"
            "print(ok); print(msg)\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("True", result.stdout)
        self.assertIn("dry-run", result.stdout)


class RunD5SlackRcTests(unittest.TestCase):
    def test_after_qa_slack_does_not_raise_on_fail(self):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        scripts = str(REPO / "calibration" / "scripts")
        if scripts not in sys.path:
            sys.path.insert(0, scripts)
        from calibration.scripts import run_d5 as rd

        with mock.patch.object(
            rd.d5_qa_slack, "try_post_after_d5",
            return_value=(False, "network down"),
        ), mock.patch.object(rd.runner, "info"), mock.patch.object(
            rd.runner, "warn"
        ) as warn:
            rd.after_qa_slack(1, "pos_review", "K")
        warn.assert_called()
        self.assertIn("Slack", warn.call_args[0][0])


if __name__ == "__main__":
    unittest.main()
