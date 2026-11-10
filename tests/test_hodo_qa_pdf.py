"""Run with python3 -m unittest discover -s tests -v."""
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]


def _make_pdf(path, pages):
    """Minimal multi-page PDF via ghostscript."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # one blank page repeated
    cmd = [
        "gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite",
        f"-sOutputFile={path}",
    ]
    for _ in range(pages):
        cmd.append("-c")
        cmd.append("showpage")
    subprocess.run(cmd, check=True, timeout=30)


class HodoQaPdfTests(unittest.TestCase):
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
            "HODO_DETECTORS = {\n"
            "  'hdprm': {'beam': ['BHT', 'BH2'], 'ftof': ['BHT', 'CVC']},\n"
            "  'hdphc': {'beam': ['BHT'], 'ftof': ['BHT', 'CVC']},\n"
            "}\n"
            "def hodo_detectors(mode, ftof):\n"
            "  return HODO_DETECTORS[mode]['ftof' if ftof else 'beam']\n"
        )
        shutil.copy(REPO / "common/lib/hodo_qa_pdf.py", lib)
        self.env = {**os.environ, "PYTHONPATH": str(self.root)}

    def py(self, code):
        return subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True, text=True, env=self.env, timeout=60,
            cwd=str(self.root))

    def test_missing_pdf_fails(self):
        result = self.py(
            "from common.lib import hodo_qa_pdf as m\n"
            "ok, msg, _ = m.merge(1, 'hdprm')\n"
            "print(ok); print(msg)\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.strip().splitlines()
        self.assertEqual(lines[0], "False")
        self.assertIn("missing PDF", lines[1])

    def test_merge_order_and_page_count(self):
        img = self.out / "img/run00001"
        _make_pdf(img / "run00001_BHT_HDPRM_Pi.pdf", 2)
        _make_pdf(img / "run00001_BH2_HDPRM_Pi.pdf", 3)
        result = self.py(
            "from common.lib import hodo_qa_pdf as m\n"
            "ok, msg, out = m.merge(1, 'hdprm')\n"
            "print(ok); print(out); print(m.pdf_page_count(out)); print(msg)\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.strip().splitlines()
        self.assertEqual(lines[0], "True")
        self.assertTrue(lines[1].endswith("run00001_hodo_hdprm_Pi.pdf"))
        self.assertEqual(lines[2], "5")
        # originals kept
        self.assertTrue((img / "run00001_BHT_HDPRM_Pi.pdf").is_file())
        self.assertTrue((img / "run00001_BH2_HDPRM_Pi.pdf").is_file())

    def test_ftof_detector_set(self):
        img = self.out / "img/run00002"
        _make_pdf(img / "run00002_BHT_HDPRM_Pi.pdf", 1)
        _make_pdf(img / "run00002_CVC_HDPRM_Pi.pdf", 1)
        result = self.py(
            "from common.lib import hodo_qa_pdf as m\n"
            "srcs = m.source_paths(2, 'hdprm', False, True)\n"
            "print([p.name for p in srcs])\n"
            "ok, msg, out = m.merge(2, 'hdprm', ftof=True)\n"
            "print(ok, out.name if out else None)\n"
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("['run00002_BHT_HDPRM_Pi.pdf', 'run00002_CVC_HDPRM_Pi.pdf']",
                      result.stdout)
        self.assertIn("run00002_hodo_hdprm_ftof_Pi.pdf", result.stdout)


if __name__ == "__main__":
    unittest.main()
