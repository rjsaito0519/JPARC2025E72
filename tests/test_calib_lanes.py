"""Lane split and d5 readiness."""
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


class CalibLanesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if str(REPO) not in sys.path:
            sys.path.insert(0, str(REPO))
        from common.lib import calib_lanes

        cls.lanes = calib_lanes

    def test_split_and_ready(self):
        steps = [
            "hodo_hdprm",
            "hodo_t0",
            "dc_t0",
            {"dc_t0": {"bcin": True}},
            {"d5_pass": {"pass": 1, "allow_skip": True}},
        ]
        split = self.lanes.split_steps_by_lane(steps)
        self.assertEqual(list(split), ["hodo", "bcout", "bcin", "d5"])
        self.assertEqual(len(split["hodo"]), 2)
        self.assertEqual(
            self.lanes.parallel_ready_lanes(split),
            ["hodo", "bcout", "bcin"],
        )
        status = {"hodo": "running", "bcout": "done", "bcin": "pending", "d5": "pending"}
        self.assertFalse(self.lanes.d5_ready(status, split))
        status["bcin"] = "done"
        self.assertTrue(self.lanes.d5_ready(status, split))

    def test_hodo_bcin_only(self):
        split = self.lanes.split_steps_by_lane(
            ["hodo", {"dc_chamber": {"bcin": True}}],
        )
        self.assertEqual(set(split), {"hodo", "bcin"})
        self.assertTrue(
            self.lanes.d5_ready({"hodo": "running", "bcin": "done"}, split),
        )

    def test_busy_and_overlap(self):
        busy = self.lanes.busy_lanes(
            {"hodo": "held", "bcout": "pending", "bcin": "running"},
        )
        self.assertEqual(busy, {"hodo", "bcin"})
        self.assertEqual(
            self.lanes.overlap_lanes(["bcout", "bcin"], busy),
            ["bcin"],
        )
        self.assertEqual(
            self.lanes.overlap_lanes(["bcout"], busy),
            [],
        )
        self.assertEqual(
            self.lanes.requested_lanes(
                ["hodo_hdprm", {"dc_t0": {"bcin": True}}],
            ),
            ["hodo", "bcin"],
        )

    def test_running_held_free_split(self):
        status = {
            "hodo": "held",
            "bcout": "running",
            "bcin": "running",
        }
        self.assertEqual(
            self.lanes.running_lanes(status), {"bcout", "bcin"},
        )
        self.assertEqual(self.lanes.held_lanes(status), {"hodo"})
        self.assertEqual(
            self.lanes.free_requested_lanes(["hodo", "d5"], status),
            ["d5"],
        )
        self.assertEqual(
            self.lanes.free_requested_lanes(["hodo"], status),
            [],
        )
        self.assertEqual(
            self.lanes.free_requested_lanes(
                ["hodo"], {"bcout": "running", "bcin": "running"},
            ),
            ["hodo"],
        )


if __name__ == "__main__":
    unittest.main()
