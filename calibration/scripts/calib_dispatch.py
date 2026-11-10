#!/usr/bin/env python3
"""Dispatch Hodo / DC / D5 iterate scripts from one YAML (see wiki Calibration-Automation).

Usage:
  calib_dispatch.py start  <config.yml> [--dry-run]
  calib_dispatch.py status <config.yml>
  calib_dispatch.py stop   <config.yml>
  calib_dispatch.py resume <config.yml> [--dry-run]

Steps: hodo | dc_resi [{bcin: true}] | dc_chamber [{bcin: true}] | d5_align
  Optional per-step: allow_skip: true  (continue after EXIT_SKIP)
  dc_chamber: t0→drift→resi→reso (see iterate_dc_chamber.py)
State: OUTPUT_DIR/calib/runXXXXX/calib_dispatch_state.json
  status=running at step start, stopped on chain halt or stop request,
  finished when all steps done.
  stop sets stop_requested; start/resume honor it between steps (no kill).
  resume continues from the stopped step (retry) or the next one after a
  continue-able stop.
Continue on exit 0; EXIT_SKIP (4) stops unless that step has allow_skip: true.
Final exit code is runner.worst of executed step codes (so a skip is not erased by a later 0).
Config: unknown top-level / step keys rejected; kaon/force/step flags are YAML bool only.
  Session YAML from calib_session.py is accepted (calib_particle, momentum_*, …).
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import config, runner  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
KNOWN_STEPS = ("hodo", "dc_resi", "dc_chamber", "d5_align")
KNOWN_TOP_KEYS = frozenset({
    "run", "kaon", "calib_particle", "force", "max_iter", "steps",
    "target", "momentum_gev_c", "momentum_source",
    "beam_particle", "beam_particle_source", "charge", "charge_source",
    "created_at", "table_updated",
})
SESSION_META_KEYS = (
    "calib_particle", "target", "momentum_gev_c", "momentum_source",
    "beam_particle", "beam_particle_source", "charge", "charge_source",
)
# Options that apply to a given step (allow_skip is common).
STEP_OPTS = {
    "hodo": frozenset({
        "skip_t0", "skip_hdphc", "ftof", "allow_skip", "hdphc_iter", "adopt",
    }),
    "dc_resi": frozenset({"bcin", "allow_skip"}),
    "dc_chamber": frozenset({
        "bcin", "allow_skip", "skip_t0", "skip_drift", "skip_resi", "skip_reso",
        "separate",
    }),
    "d5_align": frozenset({"allow_skip"}),
}
BOOL_STEP_OPTS = frozenset({
    "bcin", "skip_t0", "skip_hdphc", "skip_drift", "skip_resi", "skip_reso",
    "ftof", "allow_skip", "separate",
})
MAX_ITER_MIN = 1
MAX_ITER_MAX = 100
RUN_MIN = 1


#_____________________________________________________________________________
def state_path(run):
    return config.OUTPUT_DIR / "calib" / f"run{run:05d}" / "calib_dispatch_state.json"


#_____________________________________________________________________________
def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


#_____________________________________________________________________________
def require_bool(value, where):
    """Accept only a YAML bool (not a string or int)."""
    if not isinstance(value, bool):
        runner.fail(f"{where}: must be a YAML bool, got {value!r}")
    return value


#_____________________________________________________________________________
def require_int(value, where, lo, hi):
    """Accept an int in [lo, hi] (reject bool; YAML bool is a subclass of int)."""
    if isinstance(value, bool) or not isinstance(value, int):
        runner.fail(f"{where}: must be an int, got {value!r}")
    if value < lo or value > hi:
        runner.fail(f"{where}: must be in [{lo}, {hi}], got {value}")
    return value


#_____________________________________________________________________________
def load_config(path):
    path = Path(path).expanduser()
    if not path.is_absolute():
        path = (Path.cwd() / path).resolve()
    if not path.is_file():
        runner.fail(f"config not found: {path}")
    doc = yaml.safe_load(path.read_text())
    if not isinstance(doc, dict):
        runner.fail(f"config must be a mapping: {path}")
    unknown = sorted(set(doc) - KNOWN_TOP_KEYS)
    if unknown:
        runner.fail(f"config unknown key(s): {unknown}")
    if "run" not in doc:
        runner.fail("config missing required key: run")
    run = require_int(doc["run"], "config run", RUN_MIN, 10**9)
    steps = doc.get("steps")
    if not isinstance(steps, list) or not steps:
        runner.fail("config steps must be a non-empty list")
    kaon = False
    calib_particle = None
    if "calib_particle" in doc:
        calib_particle = str(doc["calib_particle"]).lower()
        if calib_particle not in ("pion", "kaon"):
            runner.fail(f"config calib_particle must be pion|kaon, got {doc['calib_particle']!r}")
        kaon = calib_particle == "kaon"
    if "kaon" in doc:
        kaon_flag = require_bool(doc["kaon"], "config kaon")
        if calib_particle is not None and kaon_flag != kaon:
            runner.fail("config kaon disagrees with calib_particle")
        if calib_particle is None:
            kaon = kaon_flag
            calib_particle = "kaon" if kaon else "pion"
    if calib_particle is None:
        calib_particle = "kaon" if kaon else "pion"
    force = False
    if "force" in doc:
        force = require_bool(doc["force"], "config force")
    max_iter = 10
    if "max_iter" in doc:
        max_iter = require_int(doc["max_iter"], "config max_iter", MAX_ITER_MIN, MAX_ITER_MAX)
    meta = {key: doc.get(key) for key in SESSION_META_KEYS if key in doc}
    meta["calib_particle"] = calib_particle
    return path, {
        "run": run,
        "kaon": kaon,
        "force": force,
        "max_iter": max_iter,
        "steps": steps,
        "session": meta,
    }


#_____________________________________________________________________________
def normalize_step(raw, index):
    """Return (name, opts dict) for a YAML step entry."""
    if isinstance(raw, str):
        name, opts = raw, {}
    elif isinstance(raw, dict) and len(raw) == 1:
        name, opts = next(iter(raw.items()))
        if opts is None:
            opts = {}
        if not isinstance(opts, dict):
            runner.fail(f"steps[{index}]: options must be a mapping, got {type(opts).__name__}")
    else:
        runner.fail(f"steps[{index}]: expected string or single-key mapping, got {raw!r}")
    if name not in KNOWN_STEPS:
        runner.fail(f"steps[{index}]: unknown step {name!r} (known: {', '.join(KNOWN_STEPS)})")
    allowed = STEP_OPTS[name]
    unknown = sorted(set(opts) - allowed)
    if unknown:
        runner.fail(
            f"steps[{index}] ({name}): option(s) not allowed here: {unknown} "
            f"(allowed: {sorted(allowed)})"
        )
    for key, value in opts.items():
        if key in BOOL_STEP_OPTS:
            require_bool(value, f"steps[{index}].{key}")
        if key == "hdphc_iter":
            require_int(value, f"steps[{index}].hdphc_iter", MAX_ITER_MIN, MAX_ITER_MAX)
        if key == "adopt":
            if isinstance(value, bool) or not isinstance(value, (str, int)):
                runner.fail(
                    f"steps[{index}].adopt: must be 'auto' or a positive int, got {value!r}"
                )
            if isinstance(value, str) and value != "auto":
                try:
                    value = int(value)
                except ValueError:
                    runner.fail(
                        f"steps[{index}].adopt: must be 'auto' or a positive int, got {value!r}"
                    )
            if isinstance(value, int) and value < 1:
                runner.fail(f"steps[{index}].adopt: round must be >= 1, got {value}")
            opts[key] = value
    return name, opts


#_____________________________________________________________________________
def should_continue(rc, opts):
    """True if the dispatch chain should proceed to the next step."""
    if rc == runner.EXIT_OK:
        return True
    if rc == runner.EXIT_SKIP and opts.get("allow_skip") is True:
        return True
    return False


#_____________________________________________________________________________
def step_cmd(name, opts, cfg, dry_run):
    run = cfg["run"]
    parts = ["python3", str(SCRIPT_DIR / {
        "hodo": "iterate_hodo_chain.py",
        "dc_resi": "iterate_dc_resi.py",
        "dc_chamber": "iterate_dc_chamber.py",
        "d5_align": "iterate_d5_align.py",
    }[name]), str(run)]

    if name == "dc_resi":
        if opts.get("bcin"):
            parts.append("--bcin")
        else:
            parts.append("--bcout")
        parts += ["--max-iter", str(cfg["max_iter"])]
    elif name == "dc_chamber":
        if opts.get("bcin"):
            parts.append("--bcin")
        else:
            parts.append("--bcout")
        parts += ["--max-iter", str(cfg["max_iter"])]
        if opts.get("separate"):
            parts.append("--separate")
        if opts.get("skip_t0"):
            parts.append("--skip-t0")
        if opts.get("skip_drift"):
            parts.append("--skip-drift")
        if opts.get("skip_resi"):
            parts.append("--skip-resi")
        if opts.get("skip_reso"):
            parts.append("--skip-reso")
    elif name == "d5_align":
        parts += ["--max-iter", str(cfg["max_iter"])]
    elif name == "hodo":
        if opts.get("skip_t0"):
            parts.append("--skip-t0")
        if opts.get("skip_hdphc"):
            parts.append("--skip-hdphc")
        if opts.get("ftof"):
            parts.append("--ftof")
        if "hdphc_iter" in opts:
            parts += ["--hdphc-iter", str(opts["hdphc_iter"])]
        if "adopt" in opts:
            parts += ["--adopt", str(opts["adopt"])]

    if cfg["kaon"] and name in ("hodo", "dc_resi", "dc_chamber", "d5_align"):
        parts.append("--kaon")
    if cfg["force"]:
        parts.append("--force")
    if dry_run:
        parts.append("--dry-run")
    return " ".join(parts)


#_____________________________________________________________________________
def write_state(run, payload):
    path = state_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return path


#_____________________________________________________________________________
def read_state(run):
    path = state_path(run)
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as e:
        runner.warn(f"state file unreadable: {path} ({e})")
        return None


#_____________________________________________________________________________
def plan_of(cfg):
    """Normalized execution plan used to gate resume against the YAML."""
    steps = []
    for i, raw in enumerate(cfg["steps"]):
        name, opts = normalize_step(raw, i)
        steps.append({"name": name, "opts": {k: opts[k] for k in sorted(opts)}})
    session = cfg.get("session") or {}
    calib_particle = session.get("calib_particle")
    if calib_particle is None:
        calib_particle = "kaon" if cfg["kaon"] else "pion"
    return {
        "calib_particle": str(calib_particle),
        "force": bool(cfg["force"]),
        "kaon": bool(cfg["kaon"]),
        "max_iter": int(cfg["max_iter"]),
        "steps": steps,
    }


#_____________________________________________________________________________
def plan_canonical(plan):
    return json.dumps(plan, sort_keys=True, separators=(",", ":"))


#_____________________________________________________________________________
def plans_match(saved, current):
    if not isinstance(saved, dict) or not isinstance(current, dict):
        return False
    return plan_canonical(saved) == plan_canonical(current)


#_____________________________________________________________________________
def require_resume_plan(cfg, st):
    """Fail before any step if state plan is missing or disagrees with cfg."""
    saved = st.get("plan")
    current = plan_of(cfg)
    if saved is None:
        runner.fail(
            f"run{cfg['run']:05d}: state has no plan fingerprint; "
            f"use start (not resume)"
        )
    if plans_match(saved, current):
        return
    bits = []
    if saved.get("calib_particle") != current.get("calib_particle") or (
            saved.get("kaon") != current.get("kaon")):
        bits.append(
            f"particle {saved.get('calib_particle')!r}/{saved.get('kaon')!r} -> "
            f"{current.get('calib_particle')!r}/{current.get('kaon')!r}"
        )
    if saved.get("force") != current.get("force"):
        bits.append(f"force {saved.get('force')!r} -> {current.get('force')!r}")
    if saved.get("max_iter") != current.get("max_iter"):
        bits.append(f"max_iter {saved.get('max_iter')!r} -> {current.get('max_iter')!r}")
    saved_steps = [e.get("name") for e in (saved.get("steps") or [])]
    cur_steps = [e.get("name") for e in (current.get("steps") or [])]
    if saved.get("steps") != current.get("steps"):
        bits.append(f"steps {saved_steps} -> {cur_steps}")
    detail = "; ".join(bits) if bits else "plan JSON differs"
    runner.fail(
        f"run{cfg['run']:05d}: plan mismatch on resume ({detail}); use start"
    )


#_____________________________________________________________________________
def make_payload(run, index, name, opts, rc, finished, steps_hist, status,
                 stop_requested=False, session=None, plan=None):
    payload = {
        "run": run,
        "index": index,
        "step": name,
        "opts": opts,
        "rc": rc,
        "finished": finished,
        "status": status,
        "stop_requested": bool(stop_requested),
        "steps": steps_hist,
        "updated_at": utc_now(),
    }
    if plan is not None:
        payload["plan"] = plan
    if session:
        for key in SESSION_META_KEYS:
            if key in session and session[key] is not None:
                payload[key] = session[key]
    return payload


#_____________________________________________________________________________
def stop_requested(run):
    st = read_state(run)
    return bool(st and st.get("stop_requested"))


#_____________________________________________________________________________
def overall_rc(steps_hist, fallback=runner.EXIT_OK):
    """Worst rc of the latest attempt per step index (resume may keep failed tries)."""
    if not steps_hist:
        return fallback
    by_index = {}
    for entry in steps_hist:
        by_index[entry["index"]] = entry["rc"]
    return runner.worst([by_index[i] for i in sorted(by_index)])


#_____________________________________________________________________________
def do_status(cfg):
    st = read_state(cfg["run"])
    if st is None:
        runner.info(f"run{cfg['run']:05d}: not started (no {state_path(cfg['run'])})")
        sys.exit(runner.EXIT_OK)
    mom = st.get("momentum_gev_c")
    mom_src = st.get("momentum_source")
    runner.info(
        f"run{cfg['run']:05d}: last step[{st.get('index')}]={st.get('step')} "
        f"rc={st.get('rc')} status={st.get('status')} finished={st.get('finished')} "
        f"stop_requested={st.get('stop_requested')} "
        f"calib_particle={st.get('calib_particle')} "
        f"momentum={mom} ({mom_src}) updated_at={st.get('updated_at')}"
    )
    print(json.dumps(st, indent=2, sort_keys=True))
    sys.exit(runner.EXIT_OK)


#_____________________________________________________________________________
def do_stop(cfg):
    st = read_state(cfg["run"])
    if st is None:
        runner.warn(f"run{cfg['run']:05d}: not started; nothing to stop")
        sys.exit(runner.EXIT_OK)
    if st.get("finished") or st.get("status") == "finished":
        runner.warn(f"run{cfg['run']:05d}: already finished; stop ignored")
        sys.exit(runner.EXIT_OK)
    st["stop_requested"] = True
    st["updated_at"] = utc_now()
    write_state(cfg["run"], st)
    runner.info(f"run{cfg['run']:05d}: stop requested (honored between steps)")
    sys.exit(runner.EXIT_OK)


#_____________________________________________________________________________
def resume_start_index(st):
    """Index to run next after a stopped state."""
    index = st.get("index")
    if not isinstance(index, int) or index < 0:
        runner.fail(f"state index invalid for resume: {index!r}")
    opts = st.get("opts") or {}
    if not isinstance(opts, dict):
        opts = {}
    rc = st.get("rc")
    if should_continue(rc, opts):
        return index + 1
    return index


#_____________________________________________________________________________
def run_steps(cfg, dry_run, start_index, steps_hist, label):
    """Run cfg steps from start_index; exit via sys.exit."""
    session = cfg.get("session")
    plan = plan_of(cfg)
    n_steps = len(cfg["steps"])
    if start_index > n_steps:
        runner.fail(f"{label}: start index {start_index} past end ({n_steps} steps)")
    if start_index == n_steps:
        done_rc = overall_rc(steps_hist)
        if not dry_run and steps_hist:
            last = steps_hist[-1]
            write_state(
                cfg["run"],
                make_payload(
                    cfg["run"], last["index"], last["step"], last.get("opts") or {},
                    done_rc, True, steps_hist, "finished", session=session, plan=plan,
                ),
            )
        runner.info(f"calib_dispatch {label} finished (no remaining steps)")
        sys.exit(done_rc)

    runner.info(
        f"calib_dispatch {label} run={cfg['run']:05d} "
        f"from_step={start_index + 1}/{n_steps}"
        + (" dry-run" if dry_run else "")
    )
    last_rc = runner.EXIT_OK
    last_name = None
    last_opts = {}
    last_index = start_index - 1
    for i in range(start_index, n_steps):
        raw = cfg["steps"][i]
        name, opts = normalize_step(raw, i)
        cmd = step_cmd(name, opts, cfg, dry_run)
        runner.step(f"calib_dispatch {i + 1}/{n_steps}: {name}")
        if not dry_run:
            write_state(
                cfg["run"],
                make_payload(
                    cfg["run"], i, name, opts, None, False, list(steps_hist), "running",
                    session=session, plan=plan,
                ),
            )
        if dry_run:
            print(f"[dry-run] {cmd}")
            rc = runner.EXIT_OK
        else:
            try:
                # Let the dispatcher record every child exit, including signals
                # and shell errors, before stopping the chain.
                rc, _ = runner.sh(cmd, check=False)
            except (KeyboardInterrupt, OSError) as e:
                runner.warn(f"step {name} interrupted: {type(e).__name__}")
                rc = runner.EXIT_FAIL
        last_rc = rc
        last_name = name
        last_opts = opts
        last_index = i
        steps_hist.append({"index": i, "step": name, "opts": opts, "rc": rc})
        runner.info(f"step {name} exit={rc}")
        if not should_continue(rc, opts):
            if not dry_run:
                write_state(
                    cfg["run"],
                    make_payload(
                        cfg["run"], i, name, opts, rc, False, steps_hist, "stopped",
                        session=session, plan=plan,
                    ),
                )
            runner.warn(f"stopping after step {name} (rc={rc})")
            sys.exit(rc)

        # Cooperative stop: honor stop_requested before the next step.
        if not dry_run and i + 1 < n_steps and stop_requested(cfg["run"]):
            write_state(
                cfg["run"],
                make_payload(
                    cfg["run"], i, name, opts, rc, False, steps_hist, "stopped",
                    session=session, plan=plan,
                ),
            )
            runner.warn(f"stop requested after step {name}; chain paused")
            sys.exit(rc)

    done_rc = overall_rc(steps_hist, fallback=last_rc)
    if not dry_run:
        write_state(
            cfg["run"],
            make_payload(
                cfg["run"], last_index, last_name, last_opts, done_rc, True,
                steps_hist, "finished", session=session, plan=plan,
            ),
        )
    runner.info("calib_dispatch finished")
    sys.exit(done_rc)


#_____________________________________________________________________________
def do_start(cfg, dry_run):
    # start always begins at index 0; clear any prior stop request via overwrite.
    run_steps(cfg, dry_run, 0, [], "start")


#_____________________________________________________________________________
def do_resume(cfg, dry_run):
    st = read_state(cfg["run"])
    if st is None:
        runner.fail(f"run{cfg['run']:05d}: not started; cannot resume")
    if st.get("finished") or st.get("status") == "finished":
        runner.fail(f"run{cfg['run']:05d}: already finished; cannot resume")
    if st.get("status") != "stopped":
        runner.fail(
            f"run{cfg['run']:05d}: status={st.get('status')!r} "
            f"(resume requires status=stopped)"
        )
    require_resume_plan(cfg, st)
    start_index = resume_start_index(st)
    hist = st.get("steps")
    if not isinstance(hist, list):
        hist = []
    run_steps(cfg, dry_run, start_index, list(hist), "resume")


#_____________________________________________________________________________
def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument("command", choices=("start", "status", "stop", "resume"))
    p.add_argument("config", type=Path, help="YAML config path")
    p.add_argument("--dry-run", action="store_true",
                   help="start/resume: print iterate commands only (no state write)")
    return p.parse_args()


#_____________________________________________________________________________
def main():
    a = parse_args()
    _, cfg = load_config(a.config)
    # Validate all steps up front (fail before running anything).
    for i, raw in enumerate(cfg["steps"]):
        normalize_step(raw, i)
    if a.command == "status":
        if a.dry_run:
            runner.warn("--dry-run ignored for status")
        do_status(cfg)
    if a.command == "stop":
        if a.dry_run:
            runner.warn("--dry-run ignored for stop")
        do_stop(cfg)
    if a.command == "resume":
        do_resume(cfg, a.dry_run)
    do_start(cfg, a.dry_run)


if __name__ == "__main__":
    main()
