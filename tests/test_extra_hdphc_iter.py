"""Unit tests for extra_hdphc_iter next-round numbering."""
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]


class NextRoundTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "out"
        self.out.mkdir()
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        scripts = str(REPO / "calibration" / "scripts")
        if scripts not in sys.path:
            sys.path.insert(0, scripts)
        from common.lib import config
        import extra_hdphc_iter as extra

        self.extra = extra
        self.cfg = config
        self._old = config.OUTPUT_DIR
        config.OUTPUT_DIR = self.out

    def tearDown(self):
        self.cfg.OUTPUT_DIR = self._old

    def test_next_round_empty_is_one(self):
        self.assertEqual(self.extra.next_round_index(1), 1)

    def test_next_round_after_journal(self):
        path = (
            self.out / "calib" / "run00001" / "hdphc_rounds" / "journal.jsonl"
        )
        path.parent.mkdir(parents=True)
        path.write_text(
            json.dumps({"round": 1, "width": 3.0}) + "\n"
            + json.dumps({"round": 2, "width": 1.0}) + "\n",
            encoding="utf-8",
        )
        with mock.patch.object(
            self.extra.hodo_hdphc_rounds, "load_journal",
            return_value=[{"round": 1}, {"round": 2}],
        ):
            self.assertEqual(self.extra.next_round_index(1), 3)


if __name__ == "__main__":
    unittest.main()
