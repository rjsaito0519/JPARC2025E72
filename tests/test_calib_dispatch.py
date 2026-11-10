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
        self.env = {**os.environ, "PYTHONPATH": str(self.root)}

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

    def write_cfg(self, steps="  - hodo\n  - dc_resi\n"):
        self.cfg.write_text("run: 1\nsteps:\n" + steps)

    def dispatch(self, command, dry=False):
        return subprocess.run(
            [sys.executable, str(self.scripts / "calib_dispatch.py"), command, str(self.cfg)]
            + (["--dry-run"] if dry else []),
            capture_output=True, text=True, env=self.env, timeout=15)

    def run_chain(self, child, steps="  - hodo\n  - dc_resi\n", dry=False):
        (self.scripts / "iterate_hodo_chain.py").write_text(child)
        (self.scripts / "iterate_dc_resi.py").write_text(
            f"from pathlib import Path\nPath({str(self.root / 'next')!r}).touch()\n")
        self.write_cfg(steps)
        return self.dispatch("start", dry=dry)

    def continue_through_gates(self, max_gates=5):
        """Decide continue + resume while status is held at step_gate."""
        for _ in range(max_gates):
            if not self.state.is_file():
                return
            state = json.loads(self.state.read_text())
            hold = state.get("hold") if isinstance(state.get("hold"), dict) else {}
            if state.get("status") != "held" or hold.get("source") != "step_gate":
                return
            decided = self.dispatch_extra("decide", "--action", "continue")
            self.assertEqual(decided.returncode, 0, decided.stderr + decided.stdout)
            resumed = self.dispatch("resume")
            # resume may exit with overall worst rc from history
            self.assertIn(
                resumed.returncode, (0, 2, 3, 4),
                resumed.stderr + resumed.stdout,
            )

    def test_nonzero_exit_records_stop(self):
        for rc in (1, 2, 4, 127, 137):
            with self.subTest(rc=rc):
                result = self.run_chain(f"raise SystemExit({rc})\n")
                self.assertEqual(result.returncode, rc, result.stderr)
                state = json.loads(self.state.read_text())
                self.assertEqual(state['status'], 'stopped')
                self.assertEqual(state['steps'][0]['rc'], rc)
                self.assertFalse(state['finished'])
                self.assertFalse(state.get('stop_requested'))
                self.assertFalse((self.root / 'next').exists())

    def test_exit_check_records_held(self):
        result = self.run_chain('raise SystemExit(3)\n')
        self.assertEqual(result.returncode, 3, result.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'held')
        self.assertEqual(state['hold']['source'], 'exit_check')
        self.assertIsNone(state['hold']['decision'])
        self.assertFalse((self.root / 'next').exists())
        journal = (
            self.root / "out/calib/run00001/hold_journal.jsonl"
        ).read_text().strip().splitlines()
        events = [json.loads(line)['event'] for line in journal]
        self.assertIn('hold', events)

    def test_interruption_records_stop(self):
        runner = self.root / 'common/lib/runner.py'
        original = runner.read_text()
        for exception in ('KeyboardInterrupt', 'OSError'):
            with self.subTest(exception=exception):
                runner.write_text(original + f'\ndef sh(*a, **kw): raise {exception}\n')
                result = self.run_chain('pass\n')
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertEqual(json.loads(self.state.read_text())['status'], 'stopped')
                self.assertFalse((self.root / 'next').exists())

    def test_success_and_allowed_skip(self):
        for rc in (0, 4):
            with self.subTest(rc=rc):
                result = self.run_chain(f'raise SystemExit({rc})\n',
                                        '  - hodo: {allow_skip: true}\n  - dc_resi\n')
                self.assertEqual(result.returncode, rc, result.stderr)
                if rc == 0:
                    state = json.loads(self.state.read_text())
                    self.assertEqual(state['status'], 'held')
                    self.assertEqual(state['hold']['source'], 'step_gate')
                    self.assertFalse((self.root / 'next').exists())
                    self.continue_through_gates()
                else:
                    # allow_skip EXIT_SKIP still continues; then step_gate
                    state = json.loads(self.state.read_text())
                    if state.get('status') == 'held':
                        self.continue_through_gates()
                self.assertTrue((self.root / 'next').exists())
                self.assertEqual(
                    json.loads(self.state.read_text())['status'],
                    'awaiting_final_approval',
                )

    def test_step_gate_holds_after_ok_before_next(self):
        result = self.run_chain('pass\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'held')
        self.assertEqual(state['hold']['source'], 'step_gate')
        self.assertEqual(state['step'], 'hodo')
        self.assertIsNone(state['hold']['decision'])
        self.assertFalse((self.root / 'next').exists())
        prog = self.root / "out/calib/run00001/progress.jsonl"
        self.assertTrue(prog.is_file())
        rows = [json.loads(line) for line in prog.read_text().splitlines() if line.strip()]
        self.assertTrue(any(r.get("event") == "step_gate" for r in rows))
        # Stop button path: stay held
        blocked = self.dispatch("resume")
        self.assertEqual(blocked.returncode, 1)
        self.assertIn('status=held', blocked.stdout)
        # Proceed: continue + resume → dc → awaiting
        decided = self.dispatch_extra("decide", "--action", "continue")
        self.assertEqual(decided.returncode, 0, decided.stderr)
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'awaiting_final_approval')
        self.assertTrue((self.root / 'next').exists())

    def test_auto_skips_step_gate(self):
        (self.scripts / "iterate_hodo_chain.py").write_text("pass\n")
        (self.scripts / "iterate_dc_resi.py").write_text(
            f"from pathlib import Path\nPath({str(self.root / 'next')!r}).touch()\n"
        )
        self.cfg.write_text("run: 1\nauto: true\nsteps:\n  - hodo\n  - dc_resi\n")
        result = self.dispatch("start")
        self.assertEqual(result.returncode, 0, result.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state["status"], "awaiting_final_approval")
        self.assertTrue((self.root / "next").exists())
        prog = self.root / "out/calib/run00001/progress.jsonl"
        rows = [json.loads(line) for line in prog.read_text().splitlines() if line.strip()]
        self.assertFalse(any(r.get("event") == "step_gate" for r in rows))

    def test_dry_run_has_no_state(self):
        result = self.run_chain('raise SystemExit(1)\n', dry=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(self.state.exists())
        self.assertFalse((self.root / 'next').exists())

    def test_stop_request_pauses_before_next_step(self):
        state_path = str(self.state)
        child = f"""
import json
from pathlib import Path
p = Path({state_path!r})
st = json.loads(p.read_text())
st['stop_requested'] = True
p.write_text(json.dumps(st) + '\\n')
"""
        result = self.run_chain(child)
        self.assertEqual(result.returncode, 0, result.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'stopped')
        self.assertEqual(state['step'], 'hodo')
        self.assertFalse(state['finished'])
        self.assertFalse(state.get('stop_requested'))
        self.assertFalse((self.root / 'next').exists())

    def test_resume_retries_failed_step(self):
        first = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(first.returncode, 1, first.stderr)
        self.assertEqual(json.loads(self.state.read_text())['status'], 'stopped')
        (self.scripts / "iterate_hodo_chain.py").write_text('pass\n')
        (self.scripts / "iterate_dc_resi.py").write_text(
            f"from pathlib import Path\nPath({str(self.root / 'next')!r}).touch()\n")
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        self.continue_through_gates()
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'awaiting_final_approval')
        self.assertTrue((self.root / 'next').exists())
        self.assertEqual([e['index'] for e in state['steps']], [0, 0, 1])

    def test_resume_after_cooperative_stop_continues_next(self):
        state_path = str(self.state)
        child = f"""
import json
from pathlib import Path
p = Path({state_path!r})
st = json.loads(p.read_text())
st['stop_requested'] = True
p.write_text(json.dumps(st) + '\\n')
"""
        first = self.run_chain(child)
        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertFalse((self.root / 'next').exists())
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'awaiting_final_approval')
        self.assertTrue((self.root / 'next').exists())
        self.assertEqual([e['index'] for e in state['steps']], [0, 1])

    def test_resume_rejects_finished_and_not_started(self):
        self.write_cfg()
        missing = self.dispatch("resume")
        self.assertEqual(missing.returncode, 1, missing.stderr)
        ok = self.run_chain('pass\n')
        self.assertEqual(ok.returncode, 0, ok.stderr)
        self.continue_through_gates()
        self.assertEqual(
            json.loads(self.state.read_text())['status'],
            'awaiting_final_approval',
        )
        again = self.dispatch("resume")
        self.assertEqual(again.returncode, 1, again.stderr)

    def test_stop_cli_sets_flag_and_ignores_finished(self):
        running = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(running.returncode, 1)
        flag = self.dispatch("stop")
        self.assertEqual(flag.returncode, 0, flag.stderr)
        self.assertTrue(json.loads(self.state.read_text())['stop_requested'])
        (self.scripts / "iterate_hodo_chain.py").write_text('pass\n')
        (self.scripts / "iterate_dc_resi.py").write_text('pass\n')
        done = self.dispatch("start")
        self.assertEqual(done.returncode, 0, done.stderr)
        self.continue_through_gates()
        ignored = self.dispatch("stop")
        self.assertEqual(ignored.returncode, 0, ignored.stderr)
        self.assertEqual(
            json.loads(self.state.read_text())['status'],
            'awaiting_final_approval',
        )

    def test_start_records_plan(self):
        result = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(result.returncode, 1, result.stderr)
        plan = json.loads(self.state.read_text())['plan']
        self.assertEqual([e['name'] for e in plan['steps']], ['hodo', 'dc_resi'])
        self.assertEqual(plan['calib_particle'], 'pion')
        self.assertFalse(plan['kaon'])
        self.assertFalse(plan['force'])
        self.assertEqual(plan['max_iter'], 10)

    def test_resume_rejects_steps_change(self):
        first = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(first.returncode, 1, first.stderr)
        before = self.state.read_text()
        (self.scripts / "iterate_dc_chamber.py").write_text('pass\n')
        (self.scripts / "iterate_dc_resi.py").write_text(
            f"from pathlib import Path\nPath({str(self.root / 'next')!r}).touch()\n")
        self.write_cfg('  - dc_chamber\n  - dc_resi\n')
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 1, resumed.stdout + resumed.stderr)
        self.assertIn('plan mismatch', resumed.stdout)
        self.assertEqual(self.state.read_text(), before)
        self.assertFalse((self.root / 'next').exists())

    def test_resume_rejects_particle_change(self):
        first = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(first.returncode, 1, first.stderr)
        before = self.state.read_text()
        self.cfg.write_text('run: 1\nkaon: true\nsteps:\n  - hodo\n  - dc_resi\n')
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 1, resumed.stdout + resumed.stderr)
        self.assertIn('plan mismatch', resumed.stdout)
        self.assertEqual(self.state.read_text(), before)

    def test_resume_rejects_missing_plan(self):
        first = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(first.returncode, 1, first.stderr)
        state = json.loads(self.state.read_text())
        del state['plan']
        self.state.write_text(json.dumps(state) + '\n')
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 1, resumed.stdout + resumed.stderr)
        self.assertIn('no plan fingerprint', resumed.stdout)
        self.assertEqual(json.loads(self.state.read_text())['status'], 'stopped')
        self.assertNotIn('plan', json.loads(self.state.read_text()))

    def dispatch_extra(self, command, *extra):
        return subprocess.run(
            [sys.executable, str(self.scripts / "calib_dispatch.py"), command, str(self.cfg)]
            + list(extra),
            capture_output=True, text=True, env=self.env, timeout=15)

    def test_held_blocks_resume_until_decide_refit(self):
        first = self.run_chain('raise SystemExit(3)\n')
        self.assertEqual(first.returncode, 3)
        blocked = self.dispatch("resume")
        self.assertEqual(blocked.returncode, 1, blocked.stdout + blocked.stderr)
        self.assertIn('status=held', blocked.stdout)
        decided = self.dispatch_extra("decide", "--action", "refit")
        self.assertEqual(decided.returncode, 0, decided.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'stopped')
        self.assertEqual(state['hold']['decision'], 'refit')
        (self.scripts / "iterate_hodo_chain.py").write_text('pass\n')
        (self.scripts / "iterate_dc_resi.py").write_text(
            f"from pathlib import Path\nPath({str(self.root / 'next')!r}).touch()\n")
        resumed = self.dispatch("resume")
        self.assertEqual(resumed.returncode, 0, resumed.stderr)
        self.continue_through_gates()
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'awaiting_final_approval')
        self.assertEqual([e['index'] for e in state['steps']], [0, 0, 1])

    def test_decide_continue_advances_past_check(self):
        first = self.run_chain('raise SystemExit(3)\n')
        self.assertEqual(first.returncode, 3)
        decided = self.dispatch_extra(
            "decide", "--action", "continue", "--note", "accept and go")
        self.assertEqual(decided.returncode, 0, decided.stderr)
        (self.scripts / "iterate_hodo_chain.py").write_text(
            'raise SystemExit("hodo should not re-run")\n')
        (self.scripts / "iterate_dc_resi.py").write_text(
            f"from pathlib import Path\nPath({str(self.root / 'next')!r}).touch()\n")
        resumed = self.dispatch("resume")
        # CHECK stays in step history; overall worst exit remains 3.
        self.assertEqual(resumed.returncode, 3, resumed.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'awaiting_final_approval')
        self.assertTrue((self.root / 'next').exists())
        self.assertEqual([e['index'] for e in state['steps']], [0, 1])
        self.assertEqual(state['rc'], 3)

    def test_manual_hold_on_stopped(self):
        first = self.run_chain('raise SystemExit(1)\n')
        self.assertEqual(first.returncode, 1)
        held = self.dispatch_extra("hold", "--reason", "human flag")
        self.assertEqual(held.returncode, 0, held.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'held')
        self.assertEqual(state['hold']['source'], 'manual')
        self.assertEqual(state['hold']['reason'], 'human flag')
        blocked = self.dispatch("resume")
        self.assertEqual(blocked.returncode, 1)

    def test_hold_requested_pauses_as_held(self):
        state_path = str(self.state)
        child = f"""
import json
from pathlib import Path
p = Path({state_path!r})
st = json.loads(p.read_text())
st['hold_requested'] = True
st['hold_reason'] = 'mid-chain flag'
p.write_text(json.dumps(st) + '\\n')
"""
        result = self.run_chain(child)
        self.assertEqual(result.returncode, 0, result.stderr)
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'held')
        self.assertEqual(state['hold']['source'], 'manual')
        self.assertEqual(state['hold']['reason'], 'mid-chain flag')
        self.assertFalse((self.root / 'next').exists())

    def test_approve_and_revoke(self):
        done = self.run_chain('pass\n')
        self.assertEqual(done.returncode, 0, done.stderr)
        self.continue_through_gates()
        state = json.loads(self.state.read_text())
        self.assertEqual(state['status'], 'awaiting_final_approval')
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
        done = self.run_chain('pass\n')
        self.assertEqual(done.returncode, 0, done.stderr)
        self.continue_through_gates()
        approved = self.dispatch("approve")
        self.assertEqual(approved.returncode, 1, approved.stdout + approved.stderr)
        self.assertIn("subject incomplete", approved.stdout + approved.stderr)

    def test_approve_rejects_stale_subject_sha(self):
        done = self.run_chain('pass\n')
        self.assertEqual(done.returncode, 0, done.stderr)
        self.continue_through_gates()
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
        done = self.run_chain('pass\n')
        self.assertEqual(done.returncode, 0)
        self.continue_through_gates()
        state = json.loads(self.state.read_text())
        state['status'] = 'finished'
        self.state.write_text(json.dumps(state) + '\n')
        self.seed_subject_files(hodo=True, dc=True)
        approved = self.dispatch("approve")
        self.assertEqual(approved.returncode, 0, approved.stderr)
        self.assertEqual(json.loads(self.state.read_text())['status'], 'approved')


if __name__ == '__main__':
    unittest.main()
