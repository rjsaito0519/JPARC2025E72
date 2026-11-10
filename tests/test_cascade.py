"""Unverified cascade helpers (#24 Phase D)."""
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


class CascadeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        scripts = REPO / "calibration" / "scripts"
        if str(scripts) not in sys.path:
            sys.path.insert(0, str(scripts))
        import calib_dispatch as disp
        cls.disp = disp

    def test_hodo_invalidates_dc_and_d5(self):
        steps = ["hodo", "dc_chamber", "d5_align"]
        self.assertEqual(
            self.disp.dependent_indices(steps, 0), [1, 2],
        )

    def test_dc_invalidates_d5_only(self):
        steps = ["hodo", "dc_chamber", "d5_align"]
        self.assertEqual(self.disp.dependent_indices(steps, 1), [2])

    def test_d5_has_no_dependents(self):
        steps = ["hodo", "dc_chamber", "d5_align"]
        self.assertEqual(self.disp.dependent_indices(steps, 2), [])

    def test_bcin_chamber_invalidates_d5(self):
        steps = ["hodo", "dc_chamber", {"dc_chamber": {"bcin": True}}, "d5_align"]
        self.assertEqual(self.disp.dependent_indices(steps, 2), [3])

    def test_fine_hodo_invalidates_dc_and_d5(self):
        steps = [
            "hodo_hdprm", "hodo_t0", "dc_t0", {"d5_pass": {"pass": 1}},
        ]
        self.assertEqual(self.disp.dependent_indices(steps, 0), [2, 3])


if __name__ == "__main__":
    unittest.main()
