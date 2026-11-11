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

    def test_pool_local_cap(self):
        alive = set()
        local_of = {}
        peak_local = {"n": 0}
        spawn_order = []
        age = {}

        class Handle:
            def __init__(self, run):
                self.run = run

        def spawn(run, opts):
            local = bool(opts.get("local"))
            local_of[run] = local
            spawn_order.append(run)
            alive.add(run)
            age[run] = 0
            n_local = sum(1 for r in alive if local_of.get(r))
            peak_local["n"] = max(peak_local["n"], n_local)
            return Handle(run)

        def is_alive(handle):
            # Age every alive run once per outer-style sweep: finish bsub
            # immediately; finish locals after they have been alive one reap.
            for r in list(alive):
                age[r] = age.get(r, 0) + 1
                if not local_of.get(r) or age[r] > 1:
                    alive.discard(r)
            return handle.run in alive

        entries = [
            (1, {"local": True}),
            (2, {"local": True}),
            (3, {"local": True}),
            (4, {"local": False}),
            (5, {"local": True}),
        ]
        agg, finished = self.calib_list.run_caliblist_pool(
            entries,
            max_parallel_local=2,
            spawn=spawn,
            poll_s=0.0,
            is_alive=is_alive,
            wait_rc=lambda h: 0,
        )
        self.assertEqual(agg, 0)
        self.assertEqual(len(finished), 5)
        self.assertLessEqual(peak_local["n"], 2)
        self.assertEqual(set(spawn_order), {1, 2, 3, 4, 5})
        self.assertLess(spawn_order.index(4), spawn_order.index(5))


if __name__ == "__main__":
    unittest.main()
