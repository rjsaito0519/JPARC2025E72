"""Run with python3 -m unittest discover -s tests -v (requires PyYAML).

The dispatcher and runner are exercised unchanged in an isolated tree.
Only site configuration, color, and unused ROOT-facing helpers are stubbed.
"""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest

from dispatch_isolation import populate_dispatch_lib

REPO = Path(__file__).resolve().parents[1]


class DispatchTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        paths = populate_dispatch_lib(REPO, self.root)
        self.scripts = paths["scripts"]
        self.data = paths["data"]
        self.decode = paths["decode"]
        self.scratch = paths["scratch"]
        self.state = self.root / "out/calib/run00001/calib_dispatch_state.json"
        self.cfg = self.root / "config.yml"
        self.env = {
            **os.environ,
            "PYTHONPATH": str(self.root),
            "CALIB_GATE_WAIT_POLL_S": "0.05",
        }

    def _sha(self, path):
        import hashlib
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]

    def _token(self, path):
        path = Path(path)
        st = path.stat()
        return f"{path.name}:{st.st_size}:{st.st_mtime_ns}"

    def _write_sealed_manifest(
        self, *, exec_id, family, mode, decode_path, param_paths,
        bcin=False, finished_at="2026-01-01T00:00:00Z",
    ):
        run = 1
        tag = f"run{run:05d}"
        man_dir = self.root / "out/calib" / tag / "qa_manifests"
        man_dir.mkdir(parents=True, exist_ok=True)
        seal_dir = man_dir / exec_id
        seal_dir.mkdir(parents=True, exist_ok=True)
        pdf = seal_dir / f"{exec_id}.pdf"
        pdf.write_bytes(b"%PDF " + exec_id.encode())
        meta = {"family": family, "mode": mode, "particle": "pion"}
        if family == "dc":
            meta["bcin"] = bool(bcin)
            meta["all"] = False
        written = []
        for path in param_paths:
            path = Path(path)
            written.append({
                "name": path.name,
                "path": str(path),
                "sha": self._sha(path),
                "role": "written",
            })
        (man_dir / f"{exec_id}.json").write_text(json.dumps({
            "exec_id": exec_id,
            "run": run,
            "meta": meta,
            "finished_at": finished_at,
            "input_token": self._token(decode_path),
            "outputs": [{
                "name": pdf.name,
                "path": str(pdf),
                "sha": self._sha(pdf),
            }],
            "params": [],
            "params_written": written,
        }) + "\n")

    def seed_subject_files(self, *, hodo=True, dc=True, d5=False):
        """Minimal param/decode/QA files so approve subject has no missing.

        Matches default plan modes: hodo (t0/hdprm/hdphc) + dc_resi (resi/bcout).
        """
        run = 1
        suffix = "Pi"
        tag = f"run{run:05d}"
        self.data.mkdir(parents=True, exist_ok=True)
        if hodo:
            paths = {}
            for cat, prefix, ext in (
                ("HDPRM", "HodoParam", ""),
                ("HDPHC", "HodoPHCParam", ""),
            ):
                pdir = self.root / "param" / cat / "e72"
                pdir.mkdir(parents=True, exist_ok=True)
                path = pdir / f"{prefix}_{tag}_{suffix}{ext}"
                path.write_text("x\n")
                paths[prefix] = path
            hodo_root = self.data / f"{tag}_Hodo.root"
            hodo_root.write_bytes(b"hodo")
            summary = self.root / "out/calib" / tag / "hdphc_rounds"
            summary.mkdir(parents=True, exist_ok=True)
            (summary / "summary.txt").write_text("hdphc\n")
            self._write_sealed_manifest(
                exec_id="hodo_t0", family="hodo", mode="t0",
                decode_path=hodo_root, param_paths=[paths["HodoParam"]],
            )
            self._write_sealed_manifest(
                exec_id="hodo_hdprm", family="hodo", mode="hdprm",
                decode_path=hodo_root, param_paths=[paths["HodoParam"]],
            )
            self._write_sealed_manifest(
                exec_id="hodo_hdphc", family="hodo", mode="hdphc",
                decode_path=hodo_root, param_paths=[paths["HodoPHCParam"]],
            )
        if dc:
            paths = {}
            for cat, prefix, ext in (
                ("DCTDC", "DCTdcParam", ""),
                ("DCDRFT", "DCDriftParam", ".root"),
                ("DCGEO", "DCGeomParam", ""),
            ):
                pdir = self.root / "param" / cat / "e72"
                pdir.mkdir(parents=True, exist_ok=True)
                path = pdir / f"{prefix}_{tag}_{suffix}{ext}"
                path.write_bytes(b"dc")
                paths[prefix] = path
            bcout = self.data / f"{tag}_BcOut.root"
            bcout.write_bytes(b"bcout")
            # Default dispatch cfg uses dc_resi → require resi/bcout only.
            self._write_sealed_manifest(
                exec_id="dc_resi", family="dc", mode="resi",
                decode_path=bcout, param_paths=[paths["DCGeomParam"]],
                finished_at="2026-01-01T00:00:01Z",
            )
        if d5:
            pdir = self.root / "param" / "DCGEO" / "e72"
            pdir.mkdir(parents=True, exist_ok=True)
            geo = pdir / f"DCGeomParam_{tag}_{suffix}"
            geo.write_bytes(b"geo")
            (self.decode / tag).mkdir(parents=True, exist_ok=True)
            d5_root = self.decode / tag / f"{tag}_D5_{suffix}.root"
            d5_root.write_bytes(b"d5")
            (self.scratch / tag).mkdir(parents=True, exist_ok=True)
            (self.scratch / tag / f"{tag}_D5_{suffix}_bc1reso.root").write_bytes(b"ref")
            self._write_sealed_manifest(
                exec_id="d5_pos", family="d5", mode="pos_review",
                decode_path=d5_root, param_paths=[geo],
                finished_at="2026-01-01T00:00:02Z",
            )

    # Same-lane fine steps: sequential step_gate / stop / hold semantics.
    SEQ_STEPS = "  - hodo_hdprm\n  - hodo_t0\n"
    # Multi-lane coarse plan: coordinator (no intra-lane gate between families).
    PARALLEL_STEPS = "  - hodo\n  - dc_resi\n"

    def write_cfg(self, steps=None):
        if steps is None:
            steps = self.SEQ_STEPS
        self.cfg.write_text("run: 1\nsteps:\n" + steps)

    def dispatch(self, command, dry=False):
        return subprocess.run(
            [sys.executable, str(self.scripts / "calib_dispatch.py"), command, str(self.cfg)]
            + (["--dry-run"] if dry else []),
            capture_output=True, text=True, env=self.env, timeout=15)

    def _hodo_stub(self, child):
        """Body for first fine step; --only t0 touches next then exits 0."""
        next_path = str(self.root / "next")
        return (
            "import sys\n"
            "from pathlib import Path\n"
            "only = None\n"
            "if '--only' in sys.argv:\n"
            "    i = sys.argv.index('--only')\n"
            "    if i + 1 < len(sys.argv):\n"
            "        only = sys.argv[i + 1]\n"
            "if only == 't0':\n"
            f"    Path({next_path!r}).touch()\n"
            "    raise SystemExit(0)\n"
            + child
        )

    def _install_chain_scripts(self, child, steps):
        (self.scripts / "iterate_hodo_chain.py").write_text(
            self._hodo_stub(child),
        )
        (self.scripts / "iterate_dc_resi.py").write_text(
            f"from pathlib import Path\nPath({str(self.root / 'next')!r}).touch()\n")
        (self.scripts / "iterate_dc_chamber.py").write_text("pass\n")
        self.write_cfg(steps)

    def run_chain(self, child, steps=None, dry=False):
        if steps is None:
            steps = self.SEQ_STEPS
        self._install_chain_scripts(child, steps)
        return self.dispatch("start", dry=dry)

    def start_chain_async(self, child, steps=None):
        """Start dispatch in the background (for in-process step_gate wait)."""
        if steps is None:
            steps = self.SEQ_STEPS
        self._install_chain_scripts(child, steps)
        return subprocess.Popen(
            [sys.executable, str(self.scripts / "calib_dispatch.py"),
             "start", str(self.cfg)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            env=self.env,
        )

    def lane_state_path(self, lane="hodo"):
        return (
            self.root / "out/calib/run00001/lanes" / lane
            / "calib_dispatch_state.json"
        )

    def primary_state_path(self):
        """Worker lane state when present; else top-level coordinator state."""
        for lane in ("hodo", "bcout", "bcin", "d5"):
            path = self.lane_state_path(lane)
            if path.is_file():
                return path
        return self.state

    def read_primary_state(self):
        return json.loads(self.primary_state_path().read_text())

    def iter_state_paths(self):
        """Top-level and per-lane dispatch states (lane workers write under lanes/)."""
        paths = []
        if self.state.is_file():
            paths.append(self.state)
        lanes = self.root / "out/calib/run00001/lanes"
        if lanes.is_dir():
            paths.extend(sorted(lanes.glob("*/calib_dispatch_state.json")))
        return paths

    def find_held_state(self, *, source=None):
        for path in self.iter_state_paths():
            state = json.loads(path.read_text())
            hold = state.get("hold") if isinstance(state.get("hold"), dict) else {}
            if state.get("status") != "held":
                continue
            if source is not None and hold.get("source") != source:
                continue
            if state.get("gate_waiter"):
                return state, path
        return None, None

    def wait_for_step_gate(self, timeout_s=5.0):
        """Poll until held at step_gate with gate_waiter (top or lane)."""
        deadline = time.time() + timeout_s
        while time.time() < deadline:
            state, _path = self.find_held_state(source="step_gate")
            if state is not None:
                return state
            time.sleep(0.05)
        self.fail("timed out waiting for step_gate waiter")

    def dispatch_async(self, command):
        return subprocess.Popen(
            [sys.executable, str(self.scripts / "calib_dispatch.py"),
             command, str(self.cfg)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
            env=self.env,
        )

    def drive_until_done(self, proc, timeout_s=15.0):
        """While proc waits at step_gate, decide continue until it exits."""
        deadline = time.time() + timeout_s
        while time.time() < deadline:
            if proc.poll() is not None:
                return proc.returncode
            state, _path = self.find_held_state(source="step_gate")
            hold = (
                state.get("hold") if isinstance(state, dict) else {}
            ) or {}
            if (
                state is not None
                and hold.get("decision") is None
            ):
                hold_id = hold.get("id")
                args = ["decide", "--action", "continue"]
                if hold_id:
                    args.extend(["--hold-id", str(hold_id)])
                decided = self.dispatch_extra(*args)
                self.assertEqual(
                    decided.returncode, 0, decided.stderr + decided.stdout,
                )
            time.sleep(0.05)
        proc.kill()
        self.fail("timed out driving step_gate waiters")

    def run_to_awaiting(self, child="pass\n", steps=None):
        """Start chain and Proceed through all step_gates to awaiting."""
        if steps is None:
            # Multi-lane: no intra-lane gate; coordinator reaches awaiting.
            steps = self.PARALLEL_STEPS
        proc = self.start_chain_async(child, steps)
        rc = self.drive_until_done(proc)
        out, err = proc.communicate(timeout=5)
        self.assertEqual(rc, 0, err + out)
        state = self.read_primary_state()
        if state.get("status") != "awaiting_final_approval" and self.state.is_file():
            state = json.loads(self.state.read_text())
        self.assertEqual(state["status"], "awaiting_final_approval")
        return state

    def continue_through_gates(self, proc=None, max_gates=5):
        """Decide continue while in-process waiter is at step_gate."""
        if proc is not None:
            self.drive_until_done(proc)
            return
        for _ in range(max_gates):
            state, _path = self.find_held_state(source="step_gate")
            if state is None:
                return
            hold = state.get("hold") if isinstance(state.get("hold"), dict) else {}
            args = ["decide", "--action", "continue"]
            if hold.get("id"):
                args.extend(["--hold-id", str(hold["id"])])
            decided = self.dispatch_extra(*args)
            self.assertEqual(decided.returncode, 0, decided.stderr + decided.stdout)
            time.sleep(0.1)

    def test_nonzero_exit_records_stop(self):
        for rc in (1, 2, 4, 127, 137):
            with self.subTest(rc=rc):
                result = self.run_chain(f"raise SystemExit({rc})\n")
                self.assertEqual(result.returncode, rc, result.stderr)
                state = self.read_primary_state()
                self.assertEqual(state['status'], 'stopped')
                self.assertEqual(state['steps'][0]['rc'], rc)
                self.assertFalse(state['finished'])
                self.assertFalse(state.get('stop_requested'))
                self.assertFalse((self.root / 'next').exists())

    def wait_for_exit_check(self, timeout_s=5.0):
        """Poll until held at exit_check with gate_waiter."""
        deadline = time.time() + timeout_s
        while time.time() < deadline:
            state, _path = self.find_held_state(source="exit_check")
            if state is not None:
                return state
            time.sleep(0.05)
        self.fail("timed out waiting for exit_check waiter")

    def test_exit_check_records_held(self):
        proc = self.start_chain_async("raise SystemExit(3)\n")
        state = self.wait_for_exit_check()
        self.assertEqual(state["status"], "held")
        self.assertEqual(state["hold"]["source"], "exit_check")
        self.assertIsNone(state["hold"]["decision"])
        self.assertTrue(state.get("gate_waiter"))
        self.assertFalse((self.root / "next").exists())
        journal = (
            self.root / "out/calib/run00001/hold_journal.jsonl"
        ).read_text().strip().splitlines()
        events = [json.loads(line)["event"] for line in journal]
        self.assertIn("hold", events)
        hold_id = state["hold"].get("id")
        args = ["decide", "--action", "continue"]
        if hold_id:
            args.extend(["--hold-id", str(hold_id)])
        decided = self.dispatch_extra(*args)
        self.assertEqual(decided.returncode, 0, decided.stderr + decided.stdout)
        # After continue: step_gate before t0 — drive remaining gates.
        self.drive_until_done(proc)
        out, err = proc.communicate(timeout=10)
        # overall_rc keeps EXIT_CHECK from the held step; chain still advanced.
        self.assertIn(proc.returncode, (0, 3), err + out)
        self.assertTrue((self.root / "next").exists())
        final = self.read_primary_state()
        self.assertEqual(final["status"], "awaiting_final_approval")

    def test_interruption_records_stop(self):
        runner = self.root / 'common/lib/runner.py'
        original = runner.read_text()
        for exception in ('KeyboardInterrupt', 'OSError'):
            with self.subTest(exception=exception):
                runner.write_text(original + f'\ndef sh(*a, **kw): raise {exception}\n')
                result = self.run_chain('pass\n')
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertEqual(self.read_primary_state()['status'], 'stopped')
                self.assertFalse((self.root / 'next').exists())

    def test_success_and_allowed_skip(self):
        for rc in (0, 4):
            with self.subTest(rc=rc):
                next_path = self.root / "next"
                if next_path.exists():
                    next_path.unlink()
                steps = (
                    "  - hodo_hdprm: {allow_skip: true}\n"
                    "  - hodo_t0\n"
                )
                proc = self.start_chain_async(
                    f"raise SystemExit({rc})\n", steps,
                )
                state = self.wait_for_step_gate()
                self.assertEqual(state["hold"]["source"], "step_gate")
                self.assertFalse(next_path.exists())
                self.continue_through_gates(proc=proc)
                out, err = proc.communicate(timeout=10)
                self.assertEqual(proc.returncode, rc, err + out)
                self.assertTrue(next_path.exists())
                self.assertEqual(
                    self.read_primary_state()["status"],
                    "awaiting_final_approval",
                )

    def test_step_gate_holds_after_ok_before_next(self):
        proc = self.start_chain_async('pass\n')
        state = self.wait_for_step_gate()
        self.assertEqual(state['status'], 'held')
        self.assertEqual(state['hold']['source'], 'step_gate')
        self.assertEqual(state['step'], 'hodo_hdprm')
        self.assertIsNone(state['hold']['decision'])
        self.assertTrue(state.get('gate_waiter'))
        self.assertFalse((self.root / 'next').exists())
        prog = self.root / "out/calib/run00001/progress.jsonl"
        self.assertTrue(prog.is_file())
        rows = [json.loads(line) for line in prog.read_text().splitlines() if line.strip()]
        self.assertTrue(any(r.get("event") == "step_gate" for r in rows))
        # Stop button path: stay held (resume blocked while held)
        blocked = self.dispatch("resume")
        self.assertEqual(blocked.returncode, 1)
        self.assertIn('status=held', blocked.stdout)
        # Proceed: decide only → in-process waiter runs next fine step → awaiting
        hold_id = state["hold"].get("id")
        args = ["decide", "--action", "continue"]
        if hold_id:
            args.extend(["--hold-id", str(hold_id)])
        decided = self.dispatch_extra(*args)
        self.assertEqual(decided.returncode, 0, decided.stderr)
        out, err = proc.communicate(timeout=10)
        self.assertEqual(proc.returncode, 0, err + out)
        state = self.read_primary_state()
        self.assertEqual(state['status'], 'awaiting_final_approval')
        self.assertTrue((self.root / 'next').exists())

    def test_auto_skips_step_gate(self):
        (self.scripts / "iterate_hodo_chain.py").write_text(
            self._hodo_stub("pass\n"),
        )
        (self.scripts / "iterate_dc_resi.py").write_text(
            f"from pathlib import Path\nPath({str(self.root / 'next')!r}).touch()\n"
        )
        self.cfg.write_text(
            "run: 1\nauto: true\nsteps:\n"
            + self.SEQ_STEPS,
        )
        result = self.dispatch("start")
        self.assertEqual(result.returncode, 0, result.stderr)
        state = self.read_primary_state()
        self.assertEqual(state["status"], "awaiting_final_approval")
        self.assertTrue((self.root / "next").exists())
        prog = self.root / "out/calib/run00001/progress.jsonl"
        rows = [json.loads(line) for line in prog.read_text().splitlines() if line.strip()]
        self.assertFalse(any(r.get("event") == "step_gate" for r in rows))

    def test_dry_run_has_no_state(self):
        result = self.run_chain('raise SystemExit(1)\n', dry=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(self.state.exists())
        self.assertFalse(any(
            p.is_file() for p in self.iter_state_paths()
        ))
        self.assertFalse((self.root / 'next').exists())

    def test_stop_request_pauses_before_next_step(self):
        lane_state = str(self.lane_state_path("hodo"))
        child = f"""
import json
from pathlib import Path
p = Path({lane_state!r})
st = json.loads(p.read_text())
st['stop_requested'] = True
p.write_text(json.dumps(st) + '\\n')
"""
        result = self.run_chain(child)
        self.assertEqual(result.returncode, 0, result.stderr)
        state = self.read_primary_state()
        self.assertEqual(state['status'], 'stopped')
        self.assertEqual(state['step'], 'hodo_hdprm')
        self.assertFalse(state['finished'])
        self.assertFalse(state.get('stop_requested'))
        self.assertFalse((self.root / 'next').exists())

    def test_resume_retries_failed_step(self):
        first = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(first.returncode, 1, first.stderr)
        self.assertEqual(self.read_primary_state()['status'], 'stopped')
        (self.scripts / "iterate_hodo_chain.py").write_text(
            self._hodo_stub("pass\n"),
        )
        proc = self.dispatch_async("resume")
        rc = self.drive_until_done(proc)
        out, err = proc.communicate(timeout=5)
        self.assertEqual(rc, 0, err + out)
        state = self.read_primary_state()
        self.assertEqual(state['status'], 'awaiting_final_approval')
        self.assertTrue((self.root / 'next').exists())
        self.assertEqual([e['index'] for e in state['steps']], [0, 0, 1])

    def test_resume_after_cooperative_stop_continues_next(self):
        lane_state = str(self.lane_state_path("hodo"))
        child = f"""
import json
from pathlib import Path
p = Path({lane_state!r})
st = json.loads(p.read_text())
st['stop_requested'] = True
p.write_text(json.dumps(st) + '\\n')
"""
        first = self.run_chain(child)
        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertFalse((self.root / 'next').exists())
        (self.scripts / "iterate_hodo_chain.py").write_text(
            self._hodo_stub("pass\n"),
        )
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        state = self.read_primary_state()
        self.assertEqual(state['status'], 'awaiting_final_approval')
        self.assertTrue((self.root / 'next').exists())
        self.assertEqual([e['index'] for e in state['steps']], [0, 1])

    def test_resume_rejects_finished_and_not_started(self):
        self.write_cfg()
        missing = self.dispatch("resume")
        self.assertEqual(missing.returncode, 1, missing.stderr)
        self.run_to_awaiting('pass\n')
        again = self.dispatch("resume")
        self.assertEqual(again.returncode, 1, again.stderr)

    def test_stop_cli_sets_flag_and_ignores_finished(self):
        running = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(running.returncode, 1)
        flag = self.dispatch("stop")
        self.assertEqual(flag.returncode, 0, flag.stderr)
        self.assertTrue(self.read_primary_state()['stop_requested'])
        self.run_to_awaiting('pass\n')
        ignored = self.dispatch("stop")
        self.assertEqual(ignored.returncode, 0, ignored.stderr)
        state = self.read_primary_state()
        if state.get("status") != "awaiting_final_approval" and self.state.is_file():
            state = json.loads(self.state.read_text())
        self.assertEqual(state["status"], "awaiting_final_approval")

    def test_start_records_plan(self):
        result = self.run_chain(
            'raise SystemExit(1)\n', steps=self.PARALLEL_STEPS,
        )
        self.assertEqual(result.returncode, 1, result.stderr)
        # Coordinator writes top-level plan even when a lane fails.
        plan = json.loads(self.state.read_text())['plan']
        self.assertEqual([e['name'] for e in plan['steps']], ['hodo', 'dc_resi'])
        self.assertEqual(plan['calib_particle'], 'pion')
        self.assertFalse(plan['kaon'])
        self.assertFalse(plan['force'])
        self.assertEqual(plan['max_iter'], 10)

    def test_resume_rejects_steps_change(self):
        first = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(first.returncode, 1, first.stderr)
        path = self.primary_state_path()
        before = path.read_text()
        (self.scripts / "iterate_dc_chamber.py").write_text('pass\n')
        (self.scripts / "iterate_dc_resi.py").write_text(
            f"from pathlib import Path\nPath({str(self.root / 'next')!r}).touch()\n")
        self.write_cfg('  - dc_chamber\n  - dc_resi\n')
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 1, resumed.stdout + resumed.stderr)
        self.assertIn('plan mismatch', resumed.stdout)
        self.assertEqual(path.read_text(), before)
        self.assertFalse((self.root / 'next').exists())

    def test_resume_rejects_particle_change(self):
        first = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(first.returncode, 1, first.stderr)
        path = self.primary_state_path()
        before = path.read_text()
        self.cfg.write_text(
            'run: 1\nkaon: true\nsteps:\n' + self.SEQ_STEPS,
        )
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 1, resumed.stdout + resumed.stderr)
        self.assertIn('plan mismatch', resumed.stdout)
        self.assertEqual(path.read_text(), before)

    def test_resume_rejects_missing_plan(self):
        first = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(first.returncode, 1, first.stderr)
        path = self.primary_state_path()
        state = json.loads(path.read_text())
        del state['plan']
        path.write_text(json.dumps(state) + '\n')
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 1, resumed.stdout + resumed.stderr)
        self.assertIn('no plan fingerprint', resumed.stdout)
        self.assertEqual(json.loads(path.read_text())['status'], 'stopped')
        self.assertNotIn('plan', json.loads(path.read_text()))

    def dispatch_extra(self, command, *extra):
        return subprocess.run(
            [sys.executable, str(self.scripts / "calib_dispatch.py"), command, str(self.cfg)]
            + list(extra),
            capture_output=True, text=True, env=self.env, timeout=15)

    def test_held_blocks_resume_until_decide_refit(self):
        proc = self.start_chain_async('raise SystemExit(3)\n')
        state = self.wait_for_exit_check()
        blocked = self.dispatch("resume")
        self.assertEqual(blocked.returncode, 1, blocked.stdout + blocked.stderr)
        self.assertIn('status=held', blocked.stdout)
        # In-process waiter will re-run after refit — fix the child first.
        (self.scripts / "iterate_hodo_chain.py").write_text(
            self._hodo_stub("pass\n"),
        )
        hold_id = state["hold"].get("id")
        args = ["decide", "--action", "refit"]
        if hold_id:
            args.extend(["--hold-id", str(hold_id)])
        decided = self.dispatch_extra(*args)
        self.assertEqual(decided.returncode, 0, decided.stderr)
        rc = self.drive_until_done(proc)
        out, err = proc.communicate(timeout=5)
        self.assertEqual(rc, 0, err + out)
        state = self.read_primary_state()
        self.assertEqual(state['status'], 'awaiting_final_approval')
        self.assertTrue((self.root / 'next').exists())
        self.assertEqual([e['index'] for e in state['steps']], [0, 0, 1])

    def test_decide_continue_advances_past_check(self):
        proc = self.start_chain_async('raise SystemExit(3)\n')
        state = self.wait_for_exit_check()
        hold_id = state["hold"].get("id")
        args = ["decide", "--action", "continue", "--note", "accept and go"]
        if hold_id:
            args.extend(["--hold-id", str(hold_id)])
        decided = self.dispatch_extra(*args)
        self.assertEqual(decided.returncode, 0, decided.stderr)
        # In-process waiter continues; may hit step_gate before t0.
        self.drive_until_done(proc)
        out, err = proc.communicate(timeout=10)
        # CHECK stays in step history; overall worst exit remains 3.
        self.assertEqual(proc.returncode, 3, err + out)
        state = self.read_primary_state()
        self.assertEqual(state['status'], 'awaiting_final_approval')
        self.assertTrue((self.root / 'next').exists())
        self.assertEqual([e['index'] for e in state['steps']], [0, 1])
        self.assertEqual(state['rc'], 3)

    def test_manual_hold_sets_hold_requested(self):
        first = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(first.returncode, 1)
        held = self.dispatch_extra("hold", "--reason", "human flag")
        self.assertEqual(held.returncode, 0, held.stderr)
        state = self.read_primary_state()
        self.assertTrue(state.get('hold_requested'))
        self.assertEqual(state.get('hold_reason'), 'human flag')
        self.assertEqual(state['status'], 'stopped')

    def test_hold_requested_pauses_as_held(self):
        lane_state = str(self.lane_state_path("hodo"))
        child = f"""
import json
from pathlib import Path
p = Path({lane_state!r})
st = json.loads(p.read_text())
st['hold_requested'] = True
st['hold_reason'] = 'mid-chain flag'
p.write_text(json.dumps(st) + '\\n')
"""
        result = self.run_chain(child)
        self.assertEqual(result.returncode, 0, result.stderr)
        state = self.read_primary_state()
        self.assertEqual(state['status'], 'held')
        self.assertEqual(state['hold']['source'], 'manual')
        self.assertEqual(state['hold']['reason'], 'mid-chain flag')
        self.assertFalse((self.root / 'next').exists())

    def test_approve_and_revoke(self):
        state = self.run_to_awaiting('pass\n')
        self.assertTrue(state['finished'])
        blocked = self.dispatch("resume")
        self.assertEqual(blocked.returncode, 1)
        self.seed_subject_files(hodo=True, dc=True)
        approved = self.dispatch_extra("approve", "--note", "looks ok")
        self.assertEqual(approved.returncode, 0, approved.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'approved')
        self.assertIn('subject_sha', state['approval'])
        revoked = self.dispatch_extra("revoke", "--reason", "recheck")
        self.assertEqual(revoked.returncode, 0, revoked.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'awaiting_final_approval')
        self.assertNotIn('approval', state)

    def test_approve_rejects_missing_subject(self):
        self.run_to_awaiting('pass\n')
        approved = self.dispatch("approve")
        self.assertEqual(approved.returncode, 1, approved.stdout + approved.stderr)
        self.assertIn("subject incomplete", approved.stdout + approved.stderr)

    def test_approve_rejects_stale_subject_sha(self):
        self.run_to_awaiting('pass\n')
        self.seed_subject_files(hodo=True, dc=True)
        bad = self.dispatch_extra(
            "approve", "--subject-sha", "deadbeefdeadbeef",
        )
        self.assertEqual(bad.returncode, 1, bad.stdout + bad.stderr)
        self.assertIn("subject changed", bad.stdout + bad.stderr)
        # Matching sha still works (compute with the same isolated tree).
        # cwd must be the fixture root so `python -c` does not pick up the
        # real repo's common/ via the process cwd on sys.path.
        probe = subprocess.run(
            [
                sys.executable, "-c",
                "import json; from common.lib import calib_approval; "
                f"st=json.loads(open({str(self.state)!r}).read()); "
                "print(calib_approval.subject_sha("
                "calib_approval.build_subject(1, st)))",
            ],
            capture_output=True, text=True, env=self.env, timeout=15,
            cwd=str(self.root),
        )
        self.assertEqual(probe.returncode, 0, probe.stderr)
        sha = probe.stdout.strip()
        ok = self.dispatch_extra("approve", "--subject-sha", sha)
        self.assertEqual(ok.returncode, 0, ok.stderr)
        self.assertEqual(json.loads(self.state.read_text())['status'], 'approved')

    def test_legacy_finished_status_treated_as_awaiting(self):
        state = self.run_to_awaiting('pass\n')
        state['status'] = 'finished'
        self.state.write_text(json.dumps(state) + '\n')
        self.seed_subject_files(hodo=True, dc=True)
        approved = self.dispatch("approve")
        self.assertEqual(approved.returncode, 0, approved.stderr)
        self.assertEqual(json.loads(self.state.read_text())['status'], 'approved')


if __name__ == '__main__':
    unittest.main()
