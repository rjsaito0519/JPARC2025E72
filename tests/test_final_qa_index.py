"""final_qa_index (#25 Phase F)."""
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


class FinalQaIndexTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "results"
        (self.out / "img" / "run00001").mkdir(parents=True)
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import config, final_qa_index

        self.idx = final_qa_index
        self.cfg = config
        self._old = config.OUTPUT_DIR
        config.OUTPUT_DIR = self.out

    def tearDown(self):
        self.cfg.OUTPUT_DIR = self._old

    def test_lists_known_pdfs(self):
        d = self.out / "img" / "run00001"
        (d / "run00001_hodo_hdprm_Pi.pdf").write_bytes(b"%PDF")
        (d / "d5_pos_review_run00001_Pi.pdf").write_bytes(b"%PDF")
        (d / "noise.txt").write_text("x")
        rows = self.idx.list_qa_pdfs(1)
        names = {r["name"] for r in rows}
        self.assertIn("run00001_hodo_hdprm_Pi.pdf", names)
        self.assertIn("d5_pos_review_run00001_Pi.pdf", names)
        self.assertNotIn("noise.txt", names)


if __name__ == "__main__":
    unittest.main()
