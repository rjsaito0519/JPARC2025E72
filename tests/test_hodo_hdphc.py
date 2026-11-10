"""hdphc BTOF monitor, round adopt, and run_hodo --skip-qa."""
import math
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from common.lib import hodo_btof_monitor, hodo_hdphc_rounds  # noqa: E402


class RobustWidthTests(unittest.TestCase):
    def test_half_width_uniform_box(self):
        # Flat weights on centers 0..100 → q16≈16, q84≈84, half-width≈34
        centers = np.arange(0.0, 101.0, 1.0)
        weights = np.ones_like(centers)
        width, n = hodo_btof_monitor.robust_half_width(centers, weights)
        self.assertAlmostEqual(width, 34.0, places=0)
        self.assertEqual(n, float(len(centers)))

    def test_empty_is_nan(self):
        width, n = hodo_btof_monitor.robust_half_width(
            np.array([1.0]), np.array([0.0]),
        )
        self.assertTrue(math.isnan(width))
        self.assertEqual(n, 0.0)


class PickRoundTests(unittest.TestCase):
    def test_auto_picks_min_width(self):
        journal = [
            {"round": 1, "width": 3.0, "ok": True},
            {"round": 2, "width": 1.0, "ok": True},
            {"round": 3, "width": 2.0, "ok": True},
        ]
        self.assertEqual(hodo_hdphc_rounds.pick_round(journal, "auto"), 2)

    def test_manual_adopt(self):
        journal = [
            {"round": 1, "width": 3.0, "ok": True},
            {"round": 2, "width": 1.0, "ok": True},
            {"round": 3, "width": 2.0, "ok": True},
        ]
        self.assertEqual(hodo_hdphc_rounds.pick_round(journal, "3"), 3)

    def test_auto_skips_bad_ok(self):
        journal = [
            {"round": 1, "width": 0.5, "ok": False},
            {"round": 2, "width": 2.0, "ok": True},
        ]
        self.assertEqual(hodo_hdphc_rounds.pick_round(journal, "auto"), 2)


class SkipQaTests(unittest.TestCase):
    def test_skip_qa_bypasses_merge(self):
        sys.path.insert(0, str(REPO / "calibration" / "scripts"))
        import run_hodo  # noqa: WPS433

        args = SimpleNamespace(skip_qa=True, kaon=False, ftof=False)
        with mock.patch.object(run_hodo, "hodo_qa_pdf") as qa:
            rc = run_hodo.after_qa_merge(1, "hdphc", args, 0)
        self.assertEqual(rc, 0)
        qa.merge.assert_not_called()


class DispatchHodoOptsTests(unittest.TestCase):
    def test_hdphc_iter_and_adopt_in_cmd(self):
        sys.path.insert(0, str(REPO / "calibration" / "scripts"))
        # Import after stubbing? calib_dispatch imports yaml and runner at load.
        import calib_dispatch as cd  # noqa: WPS433

        cfg = {"run": 1, "kaon": False, "force": False, "max_iter": 10}
        cmd = cd.step_cmd(
            "hodo", {"hdphc_iter": 3, "adopt": "auto"}, cfg, dry_run=False,
        )
        self.assertIn("--hdphc-iter 3", cmd)
        self.assertIn("--adopt auto", cmd)

    def test_normalize_adopt_int(self):
        import calib_dispatch as cd  # noqa: WPS433

        name, opts = cd.normalize_step({"hodo": {"adopt": 2, "hdphc_iter": 4}}, 0)
        self.assertEqual(name, "hodo")
        self.assertEqual(opts["adopt"], 2)
        self.assertEqual(opts["hdphc_iter"], 4)


class SnapshotRestoreTests(unittest.TestCase):
    def test_snapshot_and_restore(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out = root / "out"
            param = root / "param" / "HDPHC" / "e72"
            param.mkdir(parents=True)
            live = param / "HodoPHCParam_run00001_Pi"
            live.write_text("round1\n", encoding="utf-8")
            with mock.patch.object(hodo_hdphc_rounds.config, "OUTPUT_DIR", out), \
                 mock.patch.object(hodo_hdphc_rounds.config, "PARAM_DIR", root / "param"), \
                 mock.patch.object(hodo_hdphc_rounds.config, "SUB_DIR", "e72"), \
                 mock.patch.object(
                     hodo_hdphc_rounds.qa_session, "file_sha256", return_value="abc",
                 ):
                hodo_hdphc_rounds.clear_rounds(1)
                e1 = hodo_hdphc_rounds.snapshot_round(
                    1, 1, "Pi", {"width": 3.0, "n_entries": 10, "ok": True},
                )
                live.write_text("round2\n", encoding="utf-8")
                e2 = hodo_hdphc_rounds.snapshot_round(
                    1, 2, "Pi", {"width": 1.0, "n_entries": 10, "ok": True},
                )
                live.write_text("round3\n", encoding="utf-8")
                hodo_hdphc_rounds.snapshot_round(
                    1, 3, "Pi", {"width": 2.0, "n_entries": 10, "ok": True},
                )
                journal = hodo_hdphc_rounds.load_journal(1)
                chosen = hodo_hdphc_rounds.pick_round(journal, "auto")
                self.assertEqual(chosen, 2)
                hodo_hdphc_rounds.restore_round(1, chosen, "Pi")
                self.assertEqual(live.read_text(encoding="utf-8"), "round2\n")
                self.assertEqual(e1["width"], 3.0)
                self.assertEqual(e2["width"], 1.0)


if __name__ == "__main__":
    unittest.main()
