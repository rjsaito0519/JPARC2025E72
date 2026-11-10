"""fit_overrides (#23)."""
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


class FitOverridesTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name) / "results"
        self.out.mkdir()
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import config, fit_overrides

        self.fo = fit_overrides
        self.cfg = config
        self._old = config.OUTPUT_DIR
        config.OUTPUT_DIR = self.out

    def tearDown(self):
        self.cfg.OUTPUT_DIR = self._old

    def test_roundtrip_mip_and_flat(self):
        self.fo.set_mip_left(3, "BH2", 6, "U", 125.0)
        doc = self.fo.load(3)
        self.assertEqual(
            doc["fit"]["BH2"]["hdprm"]["bh2-6-u"]["mip_left"], 125.0,
        )
        flat = self.fo.flat_override_path(3).read_text()
        self.assertIn("bh2-6-u 125", flat)

    def test_calib_lock(self):
        self.fo.set_calib_lock(4, "BH2", "1-0-6-0-0", True)
        self.assertIn("1-0-6-0-0", self.fo.calib_lock_keys(4, "BH2"))
        self.fo.set_calib_lock(4, "BH2", "1-0-6-0-0", False)
        self.assertNotIn("1-0-6-0-0", self.fo.calib_lock_keys(4, "BH2"))


if __name__ == "__main__":
    unittest.main()
