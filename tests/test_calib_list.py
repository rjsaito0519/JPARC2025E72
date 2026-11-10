"""caliblist DEFAULT+RUN loader and particle aliases."""
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


class CalibListTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        sys.path.insert(0, str(REPO / "calibration" / "scripts"))
        from common.lib import calib_list
        import calib_session as cs

        self.calib_list = calib_list
        self.cs = cs

    def test_particle_aliases(self):
        self.assertEqual(self.cs.normalize_calib_particle("pi"), "pion")
        self.assertEqual(self.cs.normalize_calib_particle("K"), "kaon")
        self.assertEqual(self.cs.normalize_beam_particle("p"), "proton")

    def test_d5_forces_bc(self):
        steps, target = self.cs.steps_from_flags(d5=True)
        self.assertEqual(target, "chamber")
        self.assertEqual(len(steps), 3)

    def test_load_caliblist(self):
        path = Path(self.tmp.name) / "batch.yml"
        path.write_text(
            "DEFAULT:\n"
            "  calib_particle: pi\n"
            "  hodo: true\n"
            "  local: false\n"
            "RUN:\n"
            "  10:\n"
            "  11:\n"
            "    calib_particle: k\n"
            "    d5: true\n"
            "    local: true\n",
            encoding="utf-8",
        )
        entries = self.calib_list.load_caliblist(path)
        self.assertEqual(len(entries), 2)
        run0, opts0 = entries[0]
        self.assertEqual(run0, 10)
        self.assertEqual(opts0["calib_particle"], "pi")
        self.assertTrue(opts0["hodo"])
        self.assertFalse(opts0["d5"])
        run1, opts1 = entries[1]
        self.assertEqual(run1, 11)
        self.assertEqual(opts1["calib_particle"], "k")
        self.assertTrue(opts1["d5"])
        self.assertTrue(opts1["bcin"])
        self.assertTrue(opts1["bcout"])
        self.assertTrue(opts1["local"])


if __name__ == "__main__":
    unittest.main()
