"""BTOF0 QA path helpers and hdphc merge prepend order."""
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parents[1]


class HodoBtofQaTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name)
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import config, hodo_btof_qa, hodo_qa_pdf

        self.qa = hodo_btof_qa
        self.pdf = hodo_qa_pdf
        self.cfg = config
        self._old = config.OUTPUT_DIR
        config.OUTPUT_DIR = self.out

    def tearDown(self):
        self.cfg.OUTPUT_DIR = self._old

    def test_btof0_pdf_path(self):
        path = self.qa.btof0_pdf_path(12, kaon=False)
        self.assertEqual(path.name, "run00012_hodo_btof0_Pi.pdf")
        self.assertTrue(str(path).endswith("img/run00012/run00012_hodo_btof0_Pi.pdf"))

    def test_prepend_puts_btof0_first(self):
        btof = self.qa.btof0_pdf_path(3, kaon=False)
        btof.parent.mkdir(parents=True, exist_ok=True)
        btof.write_bytes(b"%PDF-1.4\n")
        det = self.out / "img" / "run00003" / "run00003_BHT_PHC_Pi.pdf"
        det.write_bytes(b"%PDF-1.4\n")
        with mock.patch.object(
            self.qa, "ensure_btof0_pdf",
            return_value=(True, "reuse", btof),
        ):
            sources, note = self.pdf.prepend_btof0_source(
                3, [det], kaon=False,
            )
        self.assertEqual(sources[0], btof)
        self.assertEqual(sources[1], det)
        self.assertIn("reuse", note)

    def test_draw_missing_root(self):
        ok, msg, path = self.qa.draw_btof0_pdf(99, kaon=False)
        self.assertFalse(ok)
        self.assertIsNone(path)
        self.assertIn("missing", msg.lower())


if __name__ == "__main__":
    unittest.main()
