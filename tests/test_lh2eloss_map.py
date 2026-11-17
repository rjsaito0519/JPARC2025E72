"""lh2eloss_map.py: map I/O, the enum codes shared with the analyzer's Lh2ElossMan, strict reading."""
import re
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "analysis" / "scripts" / "lh2eloss"))
import lh2eloss_map as M  # noqa: E402

try:
    sys.path.insert(0, str(REPO))
    from common.lib.paths import ANALYZER_DIR
except ImportError:
    ANALYZER_DIR = None

GEO = dict(TargetRadius=40.0, TargetHalfHeight=50.0, TargetCenterX=0.0, TargetCenterY=0.0, TargetCenterZ=-143.0,
           VesselRadius=48.5, VesselHalfHeight=50.0, HolderInnerRadius=53.5, HolderRadius=55.0,
           HolderWindowHalfPhi=58.0, HolderWindowHalfHeight=50.0, HolderBottomY=-119.5)


def toy_grids(particle=M.PROTON):
    rng = np.random.default_rng(1)
    grids = {}
    for rg, ex, ho in M.REQUIRED:
        L = np.array([0.0]) if rg == M.GAP else np.array([0.0, 10.0, 50.0])
        p, c = np.array([100.0, 200.0, 400.0]), np.array([0.0, 0.5, 1.0])
        grids[(particle, rg, ex, ho)] = dict(p=p, L=L, c=c, loss=rng.uniform(0, 20, (3, len(L), 3)).round(6))
    return grids


class Lh2ElossMapTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "map"

    def tearDown(self):
        self.tmp.cleanup()

    def write_text(self, text):
        self.path.write_text(text)
        return str(self.path)

    def test_round_trip(self):
        grids = toy_grids()
        M.write_map(str(self.path), grids, GEO, ["toy"])
        g2, geo2 = M.read_map(str(self.path))
        self.assertEqual(geo2, GEO)
        self.assertEqual(set(g2), set(grids))
        for k in grids:
            np.testing.assert_array_equal(g2[k]["loss"], grids[k]["loss"])

    def test_rejects_broken_maps(self):
        M.write_map(str(self.path), toy_grids(), GEO, ["toy"])
        good = self.path.read_text()
        broken = {
            "string value": good.replace("particle 0", "particle proton", 1),
            "label out of range": good.replace("holder   1", "holder   2", 1),
            "missing key": good.replace("HolderBottomY", "# HolderBottomY", 1),
            "unknown key": "Foo 1\n" + good,
            "missing grid": good[:good.rindex("\nparticle")],
            "negative loss": re.sub(r"(\n100 0 )[0-9.]+", r"\1-1", good, count=1),
        }
        for what, text in broken.items():
            with self.subTest(what):
                with self.assertRaises(M.MapError):
                    M.read_map(self.write_text(text))

    def test_classify_regions(self):
        v = np.array([[0, 0, -143], [0, 0, -143 + 51], [0, 0, -143 + 51], [0, 0, -143 + 60]], float)
        d = np.array([[0, 0, 1], [0, 0, 1], [0, 0, -1], [0, 0, 1]], float)
        t = M.classify(GEO, v, d)
        self.assertEqual(t["region"].tolist(), [M.VESSEL, M.GAP, M.GAP_CROSS, M.OUTSIDE])
        self.assertAlmostEqual(t["L"][0], 40.0)
        self.assertAlmostEqual(t["L"][2], 80.0)
        self.assertEqual(t["holder"][0], M.WINDOW)
        self.assertEqual(M.classify(GEO, v[:1], [[1, 0, 0]])["holder"][0], M.G10)

    @unittest.skipIf(ANALYZER_DIR is None, "common/lib/paths.py not configured")
    def test_codes_match_analyzer_enums(self):
        hh = Path(ANALYZER_DIR) / "include" / "Lh2ElossMan.hh"
        if not hh.exists():
            self.skipTest(f"{hh} not found")
        enums = dict(re.findall(r"\b(k[A-Z]\w*)\s*=\s*(\d+)", hh.read_text()))
        expect = dict(kProton=M.PROTON, kPion=M.PION, kKaon=M.KAON, kVessel=M.VESSEL, kGap=M.GAP,
                      kGapCross=M.GAP_CROSS, kSide=M.SIDE, kCap=M.CAP, kWindow=M.WINDOW, kG10=M.G10)
        for name, code in expect.items():
            self.assertEqual(int(enums[name]), code, name)

    @unittest.skipIf(ANALYZER_DIR is None, "common/lib/paths.py not configured")
    def test_production_map_loads(self):
        path = Path(ANALYZER_DIR) / "param" / "ELOSS" / "Lh2ElossMap_e72"
        if not path.exists():
            self.skipTest(f"{path} not found")
        grids, geo = M.read_map(str(path))
        self.assertIn((M.PROTON, M.VESSEL, M.SIDE, M.WINDOW), grids)


if __name__ == "__main__":
    unittest.main()
