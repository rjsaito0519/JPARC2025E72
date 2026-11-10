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


class DcQaPathTests(unittest.TestCase):
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
        for name in ("dc_qa_slack.py", "local_notify_config.py", "slack_notify.py",
                     "hodo_qa_pdf.py", "qa_session.py"):
            shutil.copy(REPO / "common/lib" / name, lib)
        self.env = {**os.environ, "PYTHONPATH": str(self.root)}

    def py(self, code):
        return subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True, text=True, env=self.env, timeout=30,
            cwd=str(self.root))

    def test_bcout_and_bcin_names(self):
        result = self.py(
            "from common.lib import dc_qa_slack as m\n"
            "print([p.name for p in m.source_paths(2580, 't0')])\n"
            "print([p.name for p in m.source_paths(2580, 'drift', bcin=True, kaon=True)])\n"
            "print([p.name for p in m.source_paths(2580, 'reso', all_chambers=True)])\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.strip().splitlines()
        self.assertEqual(lines[0], "['run02580_BLC2_TDC_Pi.pdf']")
        self.assertEqual(lines[1], "['run02580_BLC1_DRIFT_K.pdf']")
        self.assertEqual(
            lines[2],
            "['run02580_BLC1_pull_Pi.pdf', 'run02580_BLC2_pull_Pi.pdf']",
        )

    def test_dry_run(self):
        result = self.py(
            "from common.lib import dc_qa_slack as m\n"
            "ok, msg = m.post_dc_qa(1, 'resi', dry_run=True)\n"
            "print(ok); print(msg)\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("True", result.stdout)
        self.assertIn("dry-run", result.stdout)


class RunDcprmSlackRcTests(unittest.TestCase):
    def test_after_qa_slack_keeps_rc(self):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from calibration.scripts import run_dcprm as rd
        from common.lib import runner

        class Args:
            kaon = False
            bcin = False
            all = False

        with mock.patch.object(
            rd.dc_qa_slack, "try_post_after_dcprm",
            return_value=(False, "network down"),
        ), mock.patch.object(rd.runner, "info"), mock.patch.object(
            rd.runner, "warn"
        ) as warn:
            rc = rd.after_qa_slack(1, "t0", Args(), runner.EXIT_OK)
        self.assertEqual(rc, runner.EXIT_OK)
        warn.assert_called()
        self.assertIn("Slack", warn.call_args[0][0])


class DcSealWithoutSlackTests(unittest.TestCase):
    def test_seals_when_slack_disabled(self):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import config, dc_qa_slack, qa_session

        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        out = root / "out"
        param = root / "param"
        img = out / "img/run00001"
        img.mkdir(parents=True)
        pdir = param / "DCTDC" / "e72"
        pdir.mkdir(parents=True)
        pfile = pdir / "DCTdcParam_run00001_Pi"
        pfile.write_bytes(b"tdc")
        pdf = img / "run00001_BLC2_TDC_Pi.pdf"
        pdf.write_bytes(b"%PDF-before")
        old = {
            "OUTPUT_DIR": config.OUTPUT_DIR,
            "PARAM_DIR": config.PARAM_DIR,
            "SUB_DIR": config.SUB_DIR,
        }
        config.OUTPUT_DIR = out
        config.PARAM_DIR = param
        config.SUB_DIR = "e72"
        self.addCleanup(lambda: [setattr(config, k, v) for k, v in old.items()])
        session = qa_session.begin(
            1, "dc_t0", [pdf],
            input_path=None,
            meta={"family": "dc", "mode": "t0", "particle": "pion", "bcin": False},
        )
        # Simulate fit rewriting the PDF.
        pdf.write_bytes(b"%PDF-after-fit")
        with mock.patch.object(
            dc_qa_slack, "load_slack_config", return_value=({}, "disabled"),
        ), mock.patch.object(dc_qa_slack, "post_dc_qa") as post:
            ok, msg = dc_qa_slack.try_post_after_dcprm(1, "t0", session=session)
        self.assertTrue(ok, msg)
        self.assertIn("slack skip", msg)
        post.assert_not_called()
        man = out / "calib/run00001/qa_manifests" / f"{session['exec_id']}.json"
        self.assertTrue(man.is_file())
        sealed = out / "calib/run00001/qa_manifests" / session["exec_id"] / pdf.name
        self.assertTrue(sealed.is_file())


if __name__ == "__main__":
    unittest.main()
