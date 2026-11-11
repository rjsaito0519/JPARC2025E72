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
                     "hodo_qa_pdf.py", "qa_session.py", "calib_state_io.py"):
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
    def setUp(self):
        from common.lib import config
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        output_patch = mock.patch.object(config, "OUTPUT_DIR", Path(self.scratch.name))
        output_patch.start()
        self.addCleanup(output_patch.stop)

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


class D5IterationQaTests(unittest.TestCase):
    def setUp(self):
        sys.path.insert(0, str(REPO / "calibration/scripts"))
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.directory = Path(self.scratch.name)

    def test_bcin_updates_repeat_until_convergence(self):
        from types import SimpleNamespace
        import iterate_d5_align as align
        from common.lib import config, runner
        runlist = self.directory / "bcin.yml"
        runlist.touch()
        args = SimpleNamespace(run=2977, kaon=False, force=False,
                               bcin_max_iter=4, local=True)
        with mock.patch.object(config, "OUTPUT_DIR", self.directory), \
                mock.patch.object(align.iterate_dc_resi, "MYRUN", runlist), \
                mock.patch.object(align.iterate_dc_resi, "RUN_DCPRM", runlist), \
                mock.patch.object(runner, "sh", side_effect=[(2, ""), (2, ""), (0, "")]) as fit, \
                mock.patch.object(align.decode_run, "run_decode") as decode:
            self.assertEqual(align.refine_bcin(args, runlist, 2), runner.EXIT_OK)
        self.assertEqual(fit.call_count, 3)
        self.assertEqual(decode.call_count, 2)
        self.assertTrue(all(call.kwargs["local"] for call in decode.call_args_list))

    def test_bcin_reaching_limit_requires_review(self):
        from types import SimpleNamespace
        import iterate_d5_align as align
        from common.lib import config, runner
        executable = self.directory / "program"
        executable.touch()
        args = SimpleNamespace(run=2977, kaon=False, force=False,
                               bcin_max_iter=2, local=True)
        with mock.patch.object(config, "OUTPUT_DIR", self.directory), \
                mock.patch.object(align.iterate_dc_resi, "MYRUN", executable), \
                mock.patch.object(align.iterate_dc_resi, "RUN_DCPRM", executable), \
                mock.patch.object(runner, "sh", return_value=(2, "")), \
                mock.patch.object(align.decode_run, "run_decode") as decode:
            self.assertEqual(align.refine_bcin(args, executable, 1), runner.EXIT_CHECK)
        self.assertEqual(decode.call_count, 1)

    def test_missing_or_empty_layer_histograms_fail_validation(self):
        import numpy as np
        import uproot
        import run_d5 as d5
        paths = [self.directory / name for name in ("d5.root", "bcin.root", "bcout.root")]
        with uproot.recreate(paths[0]) as output:
            for det in ("BLC1a", "BLC1b", "BLC2a", "BLC2b"):
                for kind in ("D5Fit", "D5Excl", "D5Pull"):
                    output[f"D5WireResi_{det}_{kind}_vs_plane_Pi"] = (
                        np.ones((8, 2)), np.arange(9), np.arange(3),
                    )
        for path, detectors in ((paths[1], ("BLC1a", "BLC1b")),
                                (paths[2], ("BLC2a", "BLC2b"))):
            with uproot.recreate(path) as output:
                for det in detectors:
                    for plane in range(8):
                        output[f"{det}_Track_Residual_plane{plane}_Pi"] = (
                            np.ones(2), np.arange(3),
                        )
        with mock.patch.object(d5, "d5_root", return_value=paths[0]), \
                mock.patch.object(d5, "chamber_root", side_effect=lambda run, det, s:
                                  paths[1] if det == "BcIn" else paths[2]):
            self.assertEqual(d5.layer_qa_input_error(2977, "Pi"), "")
            with uproot.update(paths[0]) as output:
                values = np.ones((8, 2))
                values[3] = 0
                output["D5WireResi_BLC1a_D5Fit_vs_plane_Pi"] = (
                    values, np.arange(9), np.arange(3),
                )
            self.assertIn("empty layer histogram", d5.layer_qa_input_error(2977, "Pi"))

    def test_layer_pdf_failure_prevents_geometry_judgement_and_write(self):
        from types import SimpleNamespace
        import run_d5 as d5
        from common.lib import runner
        args = SimpleNamespace(runs=[2977], force=False, debug=False, tol=.01)
        with mock.patch.object(d5, "check_pos_inputs"), \
                mock.patch.object(d5.qa_session, "begin", return_value={}), \
                mock.patch.object(d5.qa_session, "attach_params"), \
                mock.patch.object(d5.d5_qa_slack, "pass_metadata", return_value={}), \
                mock.patch.object(d5, "measure", return_value={"dy_a":.1,"dy_b":.1,"dy_err":.01,"s":.1}), \
                mock.patch.object(d5, "after_qa_slack", return_value=True), \
                mock.patch.object(d5, "layer_qa", return_value=False) as layers, \
                mock.patch.object(d5, "write_summary"), \
                mock.patch.object(d5, "judge_pos") as judge, \
                mock.patch.object(d5, "write_ofs") as write:
            with self.assertRaises(SystemExit) as stopped:
                d5.do_pos(args, "Pi")
            self.assertEqual(stopped.exception.code, runner.EXIT_CHECK)
        layers.assert_called_once_with(2977, "Pi")
        judge.assert_not_called()
        write.assert_not_called()


if __name__ == "__main__":
    unittest.main()
