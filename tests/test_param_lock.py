"""Shared param RMW: exclusive_param_lock prevents lost chamber updates."""
from __future__ import annotations

import sys
import tempfile
import threading
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


class ParamLockTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import param as param_mod

        cls.param = param_mod

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "DCTdcParam_run00001_K"
        self.path.write_text("1\t0\n2\t0\n")

    def test_unlocked_lost_update_sequence(self):
        """Issue #47 reproduction: A then B write from the same stale load."""
        pf_a = self.param.ParamFile(self.path)
        pf_b = self.param.ParamFile(self.path)
        pf_a.update({"1": ["10"]}, key_len=1, start_col=1)
        self.path.write_text(pf_a.serialized_text())
        pf_b.update({"2": ["20"]}, key_len=1, start_col=1)
        self.path.write_text(pf_b.serialized_text())
        data = self.param.ParamFile(self.path).get_data_dict(key_len=1)
        self.assertEqual(data["2"], ["20"])
        self.assertEqual(data["1"], ["0"])  # A's change lost

    def test_locked_keeps_both_writers(self):
        ready_a, ready_b = threading.Event(), threading.Event()
        go = threading.Event()
        errors = []
        lock = self.param.exclusive_param_lock
        ParamFile = self.param.ParamFile
        path = self.path

        def writer(key, value, ready):
            try:
                ready.set()
                go.wait(timeout=5)
                with lock(path):
                    pf = ParamFile(path)
                    pf.update({str(key): [str(value)]}, key_len=1, start_col=1)
                    path.write_text(pf.serialized_text())
            except Exception as exc:  # noqa: BLE001 — surface in main thread
                errors.append(exc)

        threads = [
            threading.Thread(target=writer, args=(1, 10, ready_a)),
            threading.Thread(target=writer, args=(2, 20, ready_b)),
        ]
        for thread in threads:
            thread.start()
        self.assertTrue(ready_a.wait(5) and ready_b.wait(5))
        go.set()
        for thread in threads:
            thread.join(timeout=5)
        self.assertEqual(errors, [])
        data = ParamFile(path).get_data_dict(key_len=1)
        self.assertEqual(data["1"], ["10"])
        self.assertEqual(data["2"], ["20"])


if __name__ == "__main__":
    unittest.main()
