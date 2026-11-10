#!/usr/bin/env python3
"""Dispatch Hodo / DC / D5 iterate scripts from one YAML (see wiki Calibration-Automation).

Usage:
  calib_dispatch.py start  <config.yml> [--dry-run]
  calib_dispatch.py status <config.yml>
  calib_dispatch.py stop   <config.yml>
  calib_dispatch.py hold   <config.yml> --reason TEXT
  calib_dispatch.py decide <config.yml> --action refit|continue [--note TEXT]
  calib_dispatch.py approve <config.yml> [--note TEXT] [--via cli|slack] [--user-id ID]
  calib_dispatch.py revoke  <config.yml> [--reason TEXT]
  calib_dispatch.py resume <config.yml> [--dry-run]

Steps (coarse / --auto): hodo | dc_resi | dc_chamber | d5_align
  Fine (default gates): hodo_hdprm|hodo_t0|hodo_hdphc|hodo_hdphc_adopt |
    dc_t0|dc_drift|dc_ch_resi|dc_reso | d5_pass
  Optional per-step: allow_skip: true  (continue after EXIT_SKIP)
  auto: true skips Slack step_gate between steps
State: OUTPUT_DIR/calib/runXXXXX/calib_dispatch_state.json
  status=running at step start, stopped on chain halt or stop request,
  held on EXIT_CHECK or manual hold (needs decide before resume),
  awaiting_final_approval when all steps done (not auto-approved),
  approved after human approve.
  stop sets stop_requested; hold sets hold_requested (or immediate held);
  start/resume honor them between steps (no kill).
  resume requires status=stopped (not held / awaiting / approved).
Continue on exit 0; EXIT_SKIP (4) stops unless that step has allow_skip: true.
EXIT_CHECK (3) becomes held (not a silent resume).
Final exit code is runner.worst of executed step codes (so a skip is not erased by a later 0).
Config: unknown top-level / step keys rejected; kaon/force/step flags are YAML bool only.
  Session YAML from calib_session.py is accepted (calib_particle, momentum_*, …).
"""
import json
import os
import secrets
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import yaml

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import calib_approval, config, runner  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
KNOWN_STEPS = (
    "hodo", "hodo_hdprm", "hodo_t0", "hodo_hdphc", "hodo_hdphc_adopt",
    "dc_resi", "dc_chamber", "dc_t0", "dc_drift", "dc_ch_resi", "dc_reso",
    "d5_align", "d5_pass",
)
KNOWN_TOP_KEYS = frozenset({
    "run", "kaon", "calib_particle", "force", "max_iter", "steps", "local",
    "auto",
    "target", "momentum_gev_c", "momentum_source",
    "beam_particle", "beam_particle_source", "charge", "charge_source",
    "created_at", "table_updated",
})
SESSION_META_KEYS = (
    "calib_particle", "target", "momentum_gev_c", "momentum_source",
    "beam_particle", "beam_particle_source", "charge", "charge_source",
)
HODO_FINE = frozenset({
    "hodo_hdprm", "hodo_t0", "hodo_hdphc", "hodo_hdphc_adopt",
})
DC_CHAMBER_FINE = frozenset({
    "dc_t0", "dc_drift", "dc_ch_resi", "dc_reso",
})
# Options that apply to a given step (allow_skip is common).
STEP_OPTS = {
    "hodo": frozenset({
        "skip_t0", "skip_hdphc", "ftof", "allow_skip", "hdphc_iter", "adopt",
    }),
    "hodo_hdprm": frozenset({"ftof", "allow_skip"}),
    "hodo_t0": frozenset({"ftof", "allow_skip"}),
    "hodo_hdphc": frozenset({"ftof", "allow_skip", "hdphc_round"}),
    "hodo_hdphc_adopt": frozenset({"ftof", "allow_skip", "adopt"}),
    "dc_resi": frozenset({"bcin", "allow_skip"}),
    "dc_chamber": frozenset({
        "bcin", "allow_skip", "skip_t0", "skip_drift", "skip_resi", "skip_reso",
        "separate",
    }),
    "dc_t0": frozenset({"bcin", "allow_skip", "separate"}),
    "dc_drift": frozenset({"bcin", "allow_skip", "separate"}),
    "dc_ch_resi": frozenset({"bcin", "allow_skip", "separate"}),
    "dc_reso": frozenset({"bcin", "allow_skip", "separate"}),
    "d5_align": frozenset({"allow_skip"}),
    "d5_pass": frozenset({"allow_skip", "pass"}),
}
BOOL_STEP_OPTS = frozenset({
    "bcin", "skip_t0", "skip_hdphc", "skip_drift", "skip_resi", "skip_reso",
    "ftof", "allow_skip", "separate",
})
MAX_ITER_MIN = 1
MAX_ITER_MAX = 100
RUN_MIN = 1
DECIDE_ACTIONS = ("refit", "continue")
STEP_SCRIPT = {
    "hodo": "iterate_hodo_chain.py",
    "hodo_hdprm": "iterate_hodo_chain.py",
    "hodo_t0": "iterate_hodo_chain.py",
    "hodo_hdphc": "iterate_hodo_chain.py",
    "hodo_hdphc_adopt": "iterate_hodo_chain.py",
    "dc_resi": "iterate_dc_resi.py",
    "dc_chamber": "iterate_dc_chamber.py",
    "dc_t0": "iterate_dc_chamber.py",
    "dc_drift": "iterate_dc_chamber.py",
    "dc_ch_resi": "iterate_dc_chamber.py",
    "dc_reso": "iterate_dc_chamber.py",
    "d5_align": "iterate_d5_align.py",
    "d5_pass": "iterate_d5_align.py",
}
HODO_ONLY = {
    "hodo_hdprm": "hdprm",
    "hodo_t0": "t0",
    "hodo_hdphc": "hdphc",
    "hodo_hdphc_adopt": "hdphc-adopt",
}
DC_ONLY = {
    "dc_t0": "t0",
    "dc_drift": "drift",
    "dc_ch_resi": "resi",
    "dc_reso": "reso",
}


#_____________________________________________________________________________
def state_path(run):
    return config.OUTPUT_DIR / "calib" / f"run{run:05d}" / "calib_dispatch_state.json"


#_____________________________________________________________________________
def hold_journal_path(run):
    return config.OUTPUT_DIR / "calib" / f"run{run:05d}" / "hold_journal.jsonl"


#_____________________________________________________________________________
def progress_path(run):
    from common.lib import calib_progress
    return calib_progress.progress_path(run)


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
    local = False
    if "local" in doc:
        local = require_bool(doc["local"], "config local")
    auto = False
    if "auto" in doc:
        auto = require_bool(doc["auto"], "config auto")
    max_iter = 10
    if "max_iter" in doc:
        max_iter = require_int(doc["max_iter"], "config max_iter", MAX_ITER_MIN, MAX_ITER_MAX)
    meta = {key: doc.get(key) for key in SESSION_META_KEYS if key in doc}
    meta["calib_particle"] = calib_particle
    return path, {
        "run": run,
        "kaon": kaon,
        "force": force,
        "local": local,
        "auto": auto,
        "max_iter": max_iter,
        "steps": steps,
        "session": meta,
        "momentum_gev_c": doc.get("momentum_gev_c"),
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
        if key in ("hdphc_round", "pass"):
            require_int(value, f"steps[{index}].{key}", MAX_ITER_MIN, MAX_ITER_MAX)
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
    if name == "hodo_hdphc" and "hdphc_round" not in opts:
        runner.fail(f"steps[{index}] (hodo_hdphc): hdphc_round is required")
    if name == "d5_pass" and "pass" not in opts:
        runner.fail(f"steps[{index}] (d5_pass): pass is required")
    return name, opts


# Downstream dependencies for unverified cascade (#24 Phase D).
# Keys are logical families; values are families they invalidate when changed.
CASCADE_DEPENDENTS = {
    "hodo": ("dc_chamber", "dc_resi", "d5_align"),
    "dc_chamber": ("d5_align",),
    "dc_resi": ("d5_align",),
    "d5_align": (),
}


#_____________________________________________________________________________
def step_family(name):
    if name in HODO_FINE or name == "hodo":
        return "hodo"
    if name in DC_CHAMBER_FINE or name == "dc_chamber":
        return "dc_chamber"
    if name == "dc_resi":
        return "dc_resi"
    if name in ("d5_align", "d5_pass"):
        return "d5_align"
    if name in CASCADE_DEPENDENTS:
        return name
    return name


#_____________________________________________________________________________
def dependent_indices(cfg_steps, changed_index):
    """Indices after changed_index whose step family depends on the changed step."""
    if changed_index < 0 or changed_index >= len(cfg_steps):
        return []
    changed_name, _ = normalize_step(cfg_steps[changed_index], changed_index)
    deps = set(CASCADE_DEPENDENTS.get(step_family(changed_name), ()))
    if not deps:
        return []
    out = []
    for j in range(changed_index + 1, len(cfg_steps)):
        name, _ = normalize_step(cfg_steps[j], j)
        if step_family(name) in deps:
            out.append(j)
    return out


#_____________________________________________________________________________
def mark_unverified(payload, cfg_steps, changed_index, reason):
    """Add dependent step indices to payload['unverified']."""
    payload = dict(payload)
    cur = set(int(i) for i in (payload.get("unverified") or []) if isinstance(i, int))
    added = dependent_indices(cfg_steps, changed_index)
    cur.update(added)
    payload["unverified"] = sorted(cur)
    if added:
        append_hold_journal(payload.get("run") or 0, {
            "event": "cascade_invalidate",
            "changed_index": changed_index,
            "unverified": added,
            "reason": reason,
            "at": utc_now(),
        })
    return payload


#_____________________________________________________________________________
def clear_unverified_through(payload, completed_index):
    """Drop unverified marks for indices <= completed_index."""
    payload = dict(payload)
    cur = [int(i) for i in (payload.get("unverified") or []) if isinstance(i, int)]
    payload["unverified"] = [i for i in cur if i > completed_index]
    return payload


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
    parts = ["python3", str(SCRIPT_DIR / STEP_SCRIPT[name]), str(run)]

    if name == "dc_resi":
        if opts.get("bcin"):
            parts.append("--bcin")
        else:
            parts.append("--bcout")
        parts += ["--max-iter", str(cfg["max_iter"])]
    elif name == "dc_chamber" or name in DC_CHAMBER_FINE:
        if opts.get("bcin"):
            parts.append("--bcin")
        else:
            parts.append("--bcout")
        parts += ["--max-iter", str(cfg["max_iter"])]
        if opts.get("separate"):
            parts.append("--separate")
        if name in DC_ONLY:
            parts += ["--only", DC_ONLY[name]]
        else:
            if opts.get("skip_t0"):
                parts.append("--skip-t0")
            if opts.get("skip_drift"):
                parts.append("--skip-drift")
            if opts.get("skip_resi"):
                parts.append("--skip-resi")
            if opts.get("skip_reso"):
                parts.append("--skip-reso")
    elif name in ("d5_align", "d5_pass"):
        parts += ["--max-iter", str(cfg["max_iter"])]
        mom = cfg.get("momentum_gev_c")
        if mom is None and isinstance(cfg.get("session"), dict):
            mom = cfg["session"].get("momentum_gev_c")
        if mom is not None:
            parts += ["--mom", str(mom)]
        if name == "d5_pass":
            parts += ["--only-pass", str(opts["pass"])]
    elif name == "hodo" or name in HODO_FINE:
        if name in HODO_ONLY:
            parts += ["--only", HODO_ONLY[name]]
            if name == "hodo_hdphc":
                parts += ["--hdphc-round", str(opts["hdphc_round"])]
            if name == "hodo_hdphc_adopt" and "adopt" in opts:
                parts += ["--adopt", str(opts["adopt"])]
        else:
            if opts.get("skip_t0"):
                parts.append("--skip-t0")
            if opts.get("skip_hdphc"):
                parts.append("--skip-hdphc")
            if "hdphc_iter" in opts:
                parts += ["--hdphc-iter", str(opts["hdphc_iter"])]
            if "adopt" in opts:
                parts += ["--adopt", str(opts["adopt"])]
        if opts.get("ftof"):
            parts.append("--ftof")

    kaon_steps = set(STEP_SCRIPT)
    if cfg["kaon"] and name in kaon_steps:
        parts.append("--kaon")
    if cfg["force"]:
        parts.append("--force")
    if cfg.get("local"):
        parts.append("--local")
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
def append_hold_journal(run, event):
    path = hold_journal_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    row = dict(event)
    row["at"] = row.get("at") or utc_now()
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(row, sort_keys=True) + "\n")


#_____________________________________________________________________________
def append_progress(run, event, *, step=None, index=None, rc=None,
                    status=None, hold_source=None, stage=None):
    """Append one human-readable progress line under results/calib/runXXXXX/."""
    from common.lib import calib_progress
    calib_progress.append_progress(
        run, event,
        step=step, index=index, rc=rc, status=status,
        hold_source=hold_source, stage=stage,
    )


#_____________________________________________________________________________
def note_substage(run, step, stage, *, rc=None, status="done"):
    from common.lib import calib_progress
    calib_progress.note_substage(
        run, step, stage, rc=rc, status=status,
    )


#_____________________________________________________________________________
def read_progress_tail(run, n=8):
    """Return the last n progress rows (dicts). Missing file → []."""
    from common.lib import calib_progress
    return calib_progress.read_progress_tail(run, n=n)


#_____________________________________________________________________________
def reset_progress(run):
    """Truncate progress.jsonl for a fresh start (backup previous if present)."""
    path = progress_path(run)
    if path.is_file():
        bak = path.with_suffix(path.suffix + f".bak.{utc_now().replace(':', '')}")
        try:
            path.replace(bak)
        except OSError:
            path.write_text("", encoding="utf-8")
            return path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("", encoding="utf-8")
    return path


#_____________________________________________________________________________
def clear_dispatch_state(run):
    """Remove calib_dispatch_state.json so the next start is fresh."""
    path = state_path(run)
    if path.is_file():
        path.unlink()
    return path


#_____________________________________________________________________________
def apply_decide(run, action, note=""):
    """Apply decide without exiting. Returns updated state dict."""
    if action not in DECIDE_ACTIONS:
        runner.fail(f"decide action must be one of {DECIDE_ACTIONS}, got {action!r}")
    st = read_state(run)
    if st is None:
        runner.fail(f"run{run:05d}: not started; cannot decide")
    if st.get("status") != "held":
        runner.fail(
            f"run{run:05d}: status={st.get('status')!r} "
            f"(decide requires status=held)"
        )
    hold = st.get("hold")
    if not isinstance(hold, dict):
        runner.fail(f"run{run:05d}: held state missing hold object")
    hold = dict(hold)
    hold["decision"] = action
    hold["decided_at"] = utc_now()
    if note:
        hold["note"] = note
    st["hold"] = hold
    st["status"] = "stopped"
    st["hold_requested"] = False
    st["stop_requested"] = False
    st["updated_at"] = utc_now()
    write_state(run, st)
    append_hold_journal(run, {
        "event": "decide",
        "action": action,
        "note": note,
        "source": hold.get("source"),
        "reason": hold.get("reason"),
        "step": hold.get("step"),
        "index": hold.get("index"),
        "rc": hold.get("rc"),
        "at": hold["decided_at"],
    })
    append_progress(
        run, "decide",
        step=hold.get("step"), index=hold.get("index"),
        rc=hold.get("rc"), status="stopped", hold_source=hold.get("source"),
    )
    return st


#_____________________________________________________________________________
def make_hold(source, reason, step, index, rc):
    return {
        "source": source,
        "reason": reason,
        "step": step,
        "index": index,
        "rc": rc,
        "at": utc_now(),
        "decision": None,
        "id": secrets.token_hex(8),
    }


GATE_WAIT_POLL_S = 1.0


#_____________________________________________________________________________
def mark_gate_waiter(run, *, waiting=True):
    """Set or clear in-process step_gate waiter flags on dispatch state."""
    st = read_state(run)
    if st is None:
        return None
    if waiting:
        st["gate_waiter"] = True
        st["gate_waiter_pid"] = os.getpid()
    else:
        st.pop("gate_waiter", None)
        st.pop("gate_waiter_pid", None)
    st["updated_at"] = utc_now()
    write_state(run, st)
    return st


#_____________________________________________________________________________
def wait_held_decision(run, hold_id, *, poll_s=None, label="gate"):
    """Block until Slack/CLI decide (or gate-stop) for this hold id.

    Returns 'continue' | 'refit' | 'stop'.
    Used for step_gate and exit_check in-process waits.
    """
    if poll_s is None:
        try:
            poll_s = float(os.environ.get("CALIB_GATE_WAIT_POLL_S", GATE_WAIT_POLL_S))
        except (TypeError, ValueError):
            poll_s = GATE_WAIT_POLL_S
        if poll_s <= 0:
            poll_s = GATE_WAIT_POLL_S
    from common.lib import slack_bot_ensure

    ok_bot, bot_msg = slack_bot_ensure.ensure_slack_bot(required=True)
    if not ok_bot:
        runner.fail(f"slack bot required for {label}: {bot_msg}")
    runner.info(bot_msg)
    runner.info(
        f"waiting for Slack Proceed / Refit / Stop "
        f"(hold_id={hold_id}; pid={os.getpid()}; "
        f"cluster bot tmux {slack_bot_ensure.TMUX_SESSION} "
        f"via shared heartbeat)"
    )
    mark_gate_waiter(run, waiting=True)
    try:
        while True:
            time.sleep(poll_s)
            st = read_state(run)
            if st is None:
                continue
            hold = st.get("hold") if isinstance(st.get("hold"), dict) else {}
            if hold_id and str(hold.get("id") or "") not in ("", str(hold_id)):
                # Replaced by a newer hold; keep waiting on current id only.
                continue
            decision = hold.get("decision")
            if decision in ("continue", "refit", "stop"):
                return decision
            if st.get("status") == "stopped" and decision in ("continue", "refit"):
                return decision
    finally:
        mark_gate_waiter(run, waiting=False)


#_____________________________________________________________________________
def wait_step_gate_decision(run, hold_id, *, poll_s=None):
    """Alias for wait_held_decision (step_gate label)."""
    return wait_held_decision(
        run, hold_id, poll_s=poll_s, label="step_gate",
    )


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
        "local": bool(cfg.get("local")),
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
    if isinstance(saved, dict) and "local" not in saved:
        saved = {**saved, "local": False}
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
    if saved.get("local") != current.get("local"):
        bits.append(f"local {saved.get('local')!r} -> {current.get('local')!r}")
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
                 stop_requested=False, hold_requested=False, hold=None,
                 hold_reason=None, session=None, plan=None, unverified=None):
    payload = {
        "run": run,
        "index": index,
        "step": name,
        "opts": opts,
        "rc": rc,
        "finished": finished,
        "status": status,
        "stop_requested": bool(stop_requested),
        "hold_requested": bool(hold_requested),
        "steps": steps_hist,
        "updated_at": utc_now(),
    }
    if hold is not None:
        payload["hold"] = hold
    if hold_reason is not None:
        payload["hold_reason"] = hold_reason
    if plan is not None:
        payload["plan"] = plan
    if unverified is not None:
        payload["unverified"] = list(unverified)
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
def hold_requested(run):
    st = read_state(run)
    return bool(st and st.get("hold_requested"))


#_____________________________________________________________________________
def hold_reason_pending(run):
    st = read_state(run)
    if not st:
        return "manual hold"
    reason = st.get("hold_reason")
    if isinstance(reason, str) and reason.strip():
        return reason.strip()
    return "manual hold"


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
def hold_summary(st):
    """One-line hold summary for status / Slack."""
    if not st:
        return "hold=none"
    bits = [
        f"status={st.get('status')}",
        f"hold_requested={bool(st.get('hold_requested'))}",
    ]
    hold = st.get("hold")
    if isinstance(hold, dict):
        bits.append(
            f"hold.source={hold.get('source')} decision={hold.get('decision')!r} "
            f"reason={hold.get('reason')!r}"
        )
    else:
        bits.append("hold=none")
    return " ".join(bits)


#_____________________________________________________________________________
def enter_held(run, payload, hold, event_type, cfg_steps=None):
    payload = dict(payload)
    payload["status"] = "held"
    payload["hold"] = hold
    payload["hold_requested"] = False
    payload["stop_requested"] = False
    if "hold_reason" in payload:
        del payload["hold_reason"]
    if cfg_steps is not None and isinstance(hold.get("index"), int):
        payload = mark_unverified(
            payload, cfg_steps, hold["index"], hold.get("reason") or event_type,
        )
    payload["updated_at"] = utc_now()
    write_state(run, payload)
    append_hold_journal(run, {
        "event": event_type,
        "source": hold.get("source"),
        "reason": hold.get("reason"),
        "step": hold.get("step"),
        "index": hold.get("index"),
        "rc": hold.get("rc"),
        "at": hold.get("at"),
        "unverified": payload.get("unverified") or [],
    })
    return payload


#_____________________________________________________________________________
def mark_awaiting_final(run, payload):
    """Write process-complete state and notify Slack (best-effort)."""
    payload = dict(payload)
    payload["finished"] = True
    payload["status"] = calib_approval.STATUS_AWAITING
    payload.pop("approval", None)
    payload["updated_at"] = utc_now()
    # Preserve plan/hold/steps already on payload; ensure subject can be built.
    write_state(run, payload)
    st = read_state(run) or payload
    subject = calib_approval.build_subject(run, st)
    sha = calib_approval.subject_sha(subject)
    calib_approval.append_approval_journal(run, {
        "event": "awaiting",
        "subject_sha": sha,
    })
    ok, msg = calib_approval.notify_awaiting_slack(run, sha, st)
    if not ok:
        runner.warn(f"awaiting Slack note failed: {msg}")
    elif not msg.startswith("slack skip:"):
        runner.info(msg)
    detail = calib_approval.awaiting_stage_summary(run, st)
    runner.info(
        f"run{run:05d}: awaiting final QA approval — {detail} "
        f"(subject_sha={sha[:8]}…)"
    )
    return sha


#_____________________________________________________________________________
def do_status(cfg):
    st = read_state(cfg["run"])
    if st is None:
        runner.info(f"run{cfg['run']:05d}: not started (no {state_path(cfg['run'])})")
        sys.exit(runner.EXIT_OK)
    if calib_approval.reconcile_stale_approval(cfg["run"]):
        st = read_state(cfg["run"]) or st
    mom = st.get("momentum_gev_c")
    mom_src = st.get("momentum_source")
    eff = calib_approval.effective_status(st)
    runner.info(
        f"run{cfg['run']:05d}: last step[{st.get('index')}]={st.get('step')} "
        f"rc={st.get('rc')} status={st.get('status')} effective={eff} "
        f"finished={st.get('finished')} "
        f"stop_requested={st.get('stop_requested')} "
        f"hold_requested={st.get('hold_requested')} "
        f"{hold_summary(st)} "
        f"{calib_approval.approval_summary(st)} "
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
    if calib_approval.is_process_complete(st):
        runner.warn(
            f"run{cfg['run']:05d}: process complete "
            f"({calib_approval.effective_status(st)}); stop ignored"
        )
        sys.exit(runner.EXIT_OK)
    if st.get("status") == "held":
        runner.warn(
            f"run{cfg['run']:05d}: already held; stop ignored "
            f"(use decide --action refit|continue)"
        )
        sys.exit(runner.EXIT_OK)
    st["stop_requested"] = True
    st["updated_at"] = utc_now()
    write_state(cfg["run"], st)
    runner.info(f"run{cfg['run']:05d}: stop requested (honored between steps)")
    sys.exit(runner.EXIT_OK)


#_____________________________________________________________________________
def do_hold(cfg, reason):
    reason = (reason or "").strip() or "(unspecified)"
    st = read_state(cfg["run"])
    if st is None:
        runner.fail(f"run{cfg['run']:05d}: not started; cannot hold")
    if calib_approval.is_process_complete(st):
        runner.fail(
            f"run{cfg['run']:05d}: process complete "
            f"({calib_approval.effective_status(st)}); cannot hold"
        )
    if st.get("status") == "held":
        runner.warn(f"run{cfg['run']:05d}: already held; hold ignored")
        sys.exit(runner.EXIT_OK)
    if st.get("status") == "running":
        st["hold_requested"] = True
        st["hold_reason"] = reason
        st["updated_at"] = utc_now()
        write_state(cfg["run"], st)
        runner.info(
            f"run{cfg['run']:05d}: hold requested (honored between steps): {reason}"
        )
        sys.exit(runner.EXIT_OK)
    # stopped (or other non-running): enter held immediately
    hold = make_hold(
        "manual", reason, st.get("step"), st.get("index"), st.get("rc"),
    )
    enter_held(cfg["run"], st, hold, "hold", cfg_steps=cfg["steps"])
    ok_slack, slack_msg = calib_approval.notify_held_slack(cfg["run"], hold)
    if not ok_slack:
        runner.warn(f"held Slack note failed: {slack_msg}")
    runner.info(f"run{cfg['run']:05d}: held (manual): {reason}")
    sys.exit(runner.EXIT_OK)


#_____________________________________________________________________________
def do_decide(cfg, action, note):
    st = apply_decide(cfg["run"], action, note=note or "")
    hold = st.get("hold") if isinstance(st.get("hold"), dict) else {}
    runner.info(
        f"run{cfg['run']:05d}: decided {action} "
        f"(step={hold.get('step')} index={hold.get('index')}); resume allowed"
    )
    sys.exit(runner.EXIT_OK)


#_____________________________________________________________________________
def do_flag_channel(cfg, family, mode, det, channel, reason):
    from common.lib import channel_flags
    if not mode or not det or not channel:
        runner.fail("flag-channel requires --mode --det --channel")
    if not reason or not str(reason).strip():
        runner.fail("flag-channel requires --reason")
    row = channel_flags.append_flag(
        cfg["run"],
        family=family or "hodo",
        mode=mode,
        det=det,
        channel=channel,
        source="human",
        reason=str(reason).strip(),
        status="needs_review",
    )
    st = read_state(cfg["run"])
    if st is not None and isinstance(st.get("index"), int):
        st = mark_unverified(
            st, cfg["steps"], st["index"], f"flag {mode}/{det}/{channel}",
        )
        st["updated_at"] = utc_now()
        write_state(cfg["run"], st)
    runner.info(
        f"run{cfg['run']:05d}: flagged {row['mode']}/{row['det']}/{row['channel']}"
    )
    sys.exit(runner.EXIT_OK)


#_____________________________________________________________________________
def do_clear_channel(cfg, mode, det, channel, reason):
    from common.lib import channel_flags
    if not mode or not det or not channel:
        runner.fail("clear-channel requires --mode --det --channel")
    row = channel_flags.clear_channel(
        cfg["run"],
        mode=mode,
        det=det,
        channel=channel,
        reason=str(reason or "cleared").strip(),
    )
    runner.info(
        f"run{cfg['run']:05d}: cleared {row['mode']}/{row['det']}/{row['channel']}"
    )
    sys.exit(runner.EXIT_OK)


#_____________________________________________________________________________
def do_approve(cfg, note, via, user_id, expected_subject_sha=None):
    via = via or "cli"
    st = read_state(cfg["run"])
    if st is None:
        runner.fail(f"run{cfg['run']:05d}: not started; cannot approve")
    if not calib_approval.is_process_complete(st):
        runner.fail(
            f"run{cfg['run']:05d}: status={st.get('status')!r} "
            f"(approve requires process complete)"
        )
    eff = calib_approval.effective_status(st)
    if eff == calib_approval.STATUS_APPROVED:
        runner.warn(f"run{cfg['run']:05d}: already approved")
        sys.exit(runner.EXIT_OK)
    if eff not in (
        calib_approval.STATUS_AWAITING, calib_approval.STATUS_FINISHED_LEGACY,
    ):
        runner.fail(
            f"run{cfg['run']:05d}: status={eff!r} "
            f"(approve requires awaiting_final_approval)"
        )
    unverified = [
        int(i) for i in (st.get("unverified") or []) if isinstance(i, int)
    ]
    if unverified:
        runner.fail(
            f"run{cfg['run']:05d}: unverified steps remain {unverified}; "
            f"refit/resume dependents before approve"
        )
    from common.lib import channel_flags
    open_flags = channel_flags.open_needs_review(cfg["run"])
    if open_flags:
        runner.fail(
            f"run{cfg['run']:05d}: {len(open_flags)} channel(s) still "
            f"needs_review; clear-channel before approve"
        )
    subject = calib_approval.build_subject(cfg["run"], st)
    missing = calib_approval.subject_missing(subject)
    if missing:
        runner.fail(
            f"run{cfg['run']:05d}: subject incomplete; "
            f"missing {missing} (finish decode/QA before approve)"
        )
    sha = calib_approval.subject_sha(subject)
    expected = (expected_subject_sha or "").strip() or None
    if expected is not None and expected != sha:
        runner.fail(
            f"run{cfg['run']:05d}: subject changed "
            f"(viewed={expected[:8]}… current={sha[:8]}…); "
            f"reload QA and approve again"
        )
    approval = {
        "subject_sha": sha,
        "subject": subject,
        "at": utc_now(),
        "via": via,
    }
    if note:
        approval["note"] = note
    if user_id:
        approval["user_id"] = user_id
    st["approval"] = approval
    st["status"] = calib_approval.STATUS_APPROVED
    st["finished"] = True
    st["updated_at"] = utc_now()
    write_state(cfg["run"], st)
    calib_approval.append_approval_journal(cfg["run"], {
        "event": "approve",
        "subject_sha": sha,
        "via": via,
        "user_id": user_id,
        "note": note,
    })
    runner.info(f"run{cfg['run']:05d}: approved (subject_sha={sha[:8]}… via={via})")
    sys.exit(runner.EXIT_OK)


#_____________________________________________________________________________
def do_revoke(cfg, reason):
    reason = (reason or "manual revoke").strip() or "manual revoke"
    st = read_state(cfg["run"])
    if st is None:
        runner.fail(f"run{cfg['run']:05d}: not started; cannot revoke")
    if not calib_approval.is_process_complete(st):
        runner.fail(
            f"run{cfg['run']:05d}: status={st.get('status')!r} "
            f"(revoke requires process complete)"
        )
    if calib_approval.effective_status(st) != calib_approval.STATUS_APPROVED:
        runner.warn(
            f"run{cfg['run']:05d}: not approved "
            f"({calib_approval.effective_status(st)}); nothing to revoke"
        )
        # Normalize legacy finished → awaiting
        if st.get("status") == calib_approval.STATUS_FINISHED_LEGACY:
            st["status"] = calib_approval.STATUS_AWAITING
            st["finished"] = True
            st["updated_at"] = utc_now()
            write_state(cfg["run"], st)
        sys.exit(runner.EXIT_OK)
    st["status"] = calib_approval.STATUS_AWAITING
    st["finished"] = True
    st.pop("approval", None)
    st["updated_at"] = utc_now()
    write_state(cfg["run"], st)
    calib_approval.append_approval_journal(cfg["run"], {
        "event": "revoke",
        "reason": reason,
    })
    runner.info(f"run{cfg['run']:05d}: approval revoked ({reason})")
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
    hold = st.get("hold") if isinstance(st.get("hold"), dict) else {}
    unverified = [
        int(i) for i in (st.get("unverified") or []) if isinstance(i, int)
    ]
    if hold.get("decision") == "refit" and unverified:
        return min([index] + unverified)
    if hold.get("decision") == "continue" and not should_continue(rc, opts):
        nxt = index + 1
    elif should_continue(rc, opts):
        nxt = index + 1
    else:
        nxt = index
    if unverified and hold.get("decision") == "continue":
        # Continue past the held step but re-run the earliest stale dependent.
        return min([nxt] + [u for u in unverified if u >= nxt] or [nxt])
    return nxt


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
            prev = read_state(cfg["run"])
            uv = list((prev or {}).get("unverified") or [])
            base = make_payload(
                cfg["run"], last["index"], last["step"], last.get("opts") or {},
                done_rc, True, steps_hist, calib_approval.STATUS_AWAITING,
                session=session, plan=plan, unverified=uv,
            )
            # Carry hold from prior state if present (decide continue path).
            if prev and isinstance(prev.get("hold"), dict):
                base["hold"] = prev["hold"]
            mark_awaiting_final(cfg["run"], base)
        runner.info(f"calib_dispatch {label} process complete (awaiting final approval)")
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
    prev0 = read_state(cfg["run"]) if not dry_run else None
    unverified = list((prev0 or {}).get("unverified") or []) if prev0 else []
    i = start_index
    while i < n_steps:
        raw = cfg["steps"][i]
        name, opts = normalize_step(raw, i)
        cmd = step_cmd(name, opts, cfg, dry_run)
        runner.step(f"calib_dispatch {i + 1}/{n_steps}: {name}")
        if not dry_run:
            write_state(
                cfg["run"],
                make_payload(
                    cfg["run"], i, name, opts, None, False, list(steps_hist), "running",
                    session=session, plan=plan, unverified=unverified,
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
        if should_continue(rc, opts) and not dry_run:
            # Successful step clears its own unverified mark.
            tmp = {"unverified": unverified, "run": cfg["run"]}
            tmp = clear_unverified_through(tmp, i)
            unverified = list(tmp.get("unverified") or [])
        if not should_continue(rc, opts):
            if dry_run:
                runner.warn(f"stopping after step {name} (rc={rc})")
                sys.exit(rc)
            base = make_payload(
                cfg["run"], i, name, opts, rc, False, steps_hist, "stopped",
                session=session, plan=plan, unverified=unverified,
            )
            if rc != runner.EXIT_CHECK:
                write_state(cfg["run"], base)
                append_progress(
                    cfg["run"], "step_stop",
                    step=name, index=i, rc=rc, status="stopped",
                )
                ok_slack, slack_msg = calib_approval.notify_stopped_slack(
                    cfg["run"], name, i, rc,
                )
                if not ok_slack:
                    runner.warn(f"stopped Slack note failed: {slack_msg}")
                runner.warn(f"stopping after step {name} (rc={rc})")
                sys.exit(rc)

            # EXIT_CHECK: hold and wait in-process (Slack / CLI decide).
            from common.lib import calib_anomaly
            detail = calib_anomaly.read_check_reason(cfg["run"])
            reason = (
                f"step {name} EXIT_CHECK (rc={rc})"
                + (f": {detail}" if detail else "")
            )
            hold = make_hold("exit_check", reason, name, i, rc)
            enter_held(
                cfg["run"], base, hold, "hold", cfg_steps=cfg["steps"],
            )
            append_progress(
                cfg["run"], "exit_check",
                step=name, index=i, rc=rc, status="held",
                hold_source="exit_check",
            )
            _ensure_slack_bot_for_gates(cfg, dry_run)
            ok_slack, slack_msg = calib_approval.notify_exit_check_slack(
                cfg["run"], hold,
            )
            if not ok_slack:
                runner.warn(f"held Slack note failed: {slack_msg}")
            runner.info(
                f"held after step {name} (exit_check rc={rc}); "
                f"waiting for Slack Proceed / Refit / Stop "
                f"(or ./ops.py decide {cfg['run']} continue|refit)"
            )
            decision = wait_held_decision(
                cfg["run"], hold.get("id"), label="exit_check",
            )
            if decision == "stop":
                runner.info(
                    f"exit_check stop; staying held after {name} (rc={rc})"
                )
                sys.exit(rc)
            if decision == "refit":
                runner.info(f"exit_check refit; re-running {name}")
                append_progress(
                    cfg["run"], "exit_check_refit",
                    step=name, index=i, rc=rc, status="running",
                    hold_source="exit_check",
                )
                continue
            runner.info(f"exit_check continue; next step after {name}")
            append_progress(
                cfg["run"], "exit_check_continue",
                step=name, index=i, rc=rc, status="running",
                hold_source="exit_check",
            )
            i += 1
            continue

        # Cooperative hold: honor hold_requested before the next step.
        if not dry_run and i + 1 < n_steps and hold_requested(cfg["run"]):
            reason = hold_reason_pending(cfg["run"])
            hold = make_hold("manual", reason, name, i, rc)
            base = make_payload(
                cfg["run"], i, name, opts, rc, False, steps_hist, "held",
                session=session, plan=plan, unverified=unverified,
            )
            enter_held(
                cfg["run"], base, hold, "hold", cfg_steps=cfg["steps"],
            )
            append_progress(
                cfg["run"], "hold_requested",
                step=name, index=i, rc=rc, status="held",
                hold_source="manual",
            )
            ok_slack, slack_msg = calib_approval.notify_held_slack(
                cfg["run"], hold,
            )
            if not ok_slack:
                runner.warn(f"held Slack note failed: {slack_msg}")
            runner.warn(f"hold requested after step {name}; chain held: {reason}")
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
            append_progress(
                cfg["run"], "stop_requested",
                step=name, index=i, rc=rc, status="stopped",
            )
            ok_slack, slack_msg = calib_approval.notify_stopped_slack(
                cfg["run"], name, i, rc, detail="stop_requested",
            )
            if not ok_slack:
                runner.warn(f"stopped Slack note failed: {slack_msg}")
            runner.warn(f"stop requested after step {name}; chain paused")
            sys.exit(rc)

        # Human next-step gate: wait in-process for Slack/CLI decide.
        # --auto / session auto: true skips gates (coarse or fine steps).
        if not dry_run and i + 1 < n_steps and not cfg.get("auto"):
            # Refuse the gate if Socket Mode bot cannot be ensured.
            _ensure_slack_bot_for_gates(cfg, dry_run)
            hold = make_hold(
                "step_gate",
                f"step {name} OK; next-step gate",
                name, i, rc,
            )
            base = make_payload(
                cfg["run"], i, name, opts, rc, False, steps_hist, "held",
                session=session, plan=plan, unverified=unverified,
            )
            # Do not mark_unverified: this step succeeded.
            enter_held(cfg["run"], base, hold, "hold", cfg_steps=None)
            append_progress(
                cfg["run"], "step_gate",
                step=name, index=i, rc=rc, status="held",
                hold_source="step_gate",
            )
            ok_slack, slack_msg = calib_approval.notify_step_gate_slack(
                cfg["run"], name, i, hold_id=hold.get("id"),
            )
            if not ok_slack:
                runner.warn(f"step-gate Slack note failed: {slack_msg}")
            runner.info(
                f"held after step {name} (step_gate); "
                f"waiting for Slack Proceed / Refit / Stop "
                f"(or ./ops.py decide {cfg['run']} continue|refit)"
            )
            decision = wait_held_decision(
                cfg["run"], hold.get("id"), label="step_gate",
            )
            if decision == "stop":
                runner.info(
                    f"step_gate stop; staying held after {name} (rc={rc})"
                )
                sys.exit(rc)
            if decision == "refit":
                runner.info(f"step_gate refit; re-running {name}")
                append_progress(
                    cfg["run"], "step_gate_refit",
                    step=name, index=i, rc=rc, status="running",
                    hold_source="step_gate",
                )
                continue
            runner.info(f"step_gate continue; next step after {name}")
            append_progress(
                cfg["run"], "step_gate_continue",
                step=name, index=i, rc=rc, status="running",
                hold_source="step_gate",
            )
        i += 1

    done_rc = overall_rc(steps_hist, fallback=last_rc)
    if not dry_run:
        base = make_payload(
            cfg["run"], last_index, last_name, last_opts, done_rc, True,
            steps_hist, calib_approval.STATUS_AWAITING,
            session=session, plan=plan,
        )
        prev = read_state(cfg["run"])
        if prev and isinstance(prev.get("hold"), dict):
            base["hold"] = prev["hold"]
        mark_awaiting_final(cfg["run"], base)
        append_progress(
            cfg["run"], "awaiting_final",
            step=last_name, index=last_index, rc=done_rc,
            status=calib_approval.STATUS_AWAITING,
        )
    runner.info("calib_dispatch process complete (awaiting final approval)")
    sys.exit(done_rc)


#_____________________________________________________________________________
def _ensure_slack_bot_for_gates(cfg, dry_run):
    """Start cluster Socket Mode bot when interactive gates may be needed."""
    if dry_run or cfg.get("auto"):
        return
    from common.lib import slack_bot_ensure

    ok_bot, bot_msg = slack_bot_ensure.ensure_slack_bot(required=True)
    if not ok_bot:
        runner.fail(f"slack bot required for interactive gates: {bot_msg}")
    runner.info(bot_msg)


#_____________________________________________________________________________
def do_start(cfg, dry_run):
    # start always begins at index 0; clear any prior stop/hold/approval via overwrite.
    _ensure_slack_bot_for_gates(cfg, dry_run)
    run_steps(cfg, dry_run, 0, [], "start")


#_____________________________________________________________________________
def do_resume(cfg, dry_run):
    st = read_state(cfg["run"])
    if st is None:
        runner.fail(f"run{cfg['run']:05d}: not started; cannot resume")
    if calib_approval.is_process_complete(st):
        runner.fail(
            f"run{cfg['run']:05d}: process complete "
            f"({calib_approval.effective_status(st)}); cannot resume"
        )
    if st.get("status") == "held":
        runner.fail(
            f"run{cfg['run']:05d}: status=held; "
            f"decide --action refit|continue before resume"
        )
    if st.get("status") != "stopped":
        runner.fail(
            f"run{cfg['run']:05d}: status={st.get('status')!r} "
            f"(resume requires status=stopped)"
        )
    require_resume_plan(cfg, st)
    _ensure_slack_bot_for_gates(cfg, dry_run)
    start_index = resume_start_index(st)
    hist = st.get("steps")
    if not isinstance(hist, list):
        hist = []
    run_steps(cfg, dry_run, start_index, list(hist), "resume")


#_____________________________________________________________________________
def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument(
        "command",
        choices=(
            "start", "status", "stop", "hold", "decide",
            "approve", "revoke", "resume",
            "flag-channel", "clear-channel",
        ),
    )
    p.add_argument("config", type=Path, help="YAML config path")
    p.add_argument("--dry-run", action="store_true",
                   help="start/resume: print iterate commands only (no state write)")
    p.add_argument("--reason", default=None,
                   help="hold/revoke/flag-channel: short reason")
    p.add_argument("--action", default=None, choices=DECIDE_ACTIONS,
                   help="decide: refit (retry step) or continue (advance)")
    p.add_argument("--note", default=None,
                   help="decide/approve: optional note for journal")
    p.add_argument("--via", default="cli", choices=("cli", "slack", "web"),
                   help="approve: how the approval was issued")
    p.add_argument("--user-id", default=None,
                   help="approve: Slack user id when --via slack")
    p.add_argument("--subject-sha", default=None,
                   help="approve: require this viewed subject_sha (Web/Slack pin)")
    p.add_argument("--mode", default=None,
                   help="flag-channel/clear-channel: calib mode (e.g. hdprm)")
    p.add_argument("--det", default=None,
                   help="flag-channel/clear-channel: detector name")
    p.add_argument("--channel", default=None,
                   help="flag-channel/clear-channel: channel id")
    p.add_argument("--family", default="hodo",
                   help="flag-channel: family (default hodo)")
    return p.parse_args()


#_____________________________________________________________________________
def main():
    runner.warn(
        "prefer ./calib.py <run> to launch and ./ops.py for state "
        "(calib_dispatch is compat)"
    )
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
    if a.command == "hold":
        if a.dry_run:
            runner.warn("--dry-run ignored for hold")
        do_hold(cfg, a.reason)
    if a.command == "decide":
        if a.dry_run:
            runner.warn("--dry-run ignored for decide")
        if not a.action:
            runner.fail("decide requires --action refit|continue")
        do_decide(cfg, a.action, a.note)
    if a.command == "approve":
        if a.dry_run:
            runner.warn("--dry-run ignored for approve")
        do_approve(cfg, a.note, a.via, a.user_id, a.subject_sha)
    if a.command == "revoke":
        if a.dry_run:
            runner.warn("--dry-run ignored for revoke")
        do_revoke(cfg, a.reason)
    if a.command == "flag-channel":
        if a.dry_run:
            runner.warn("--dry-run ignored for flag-channel")
        do_flag_channel(cfg, a.family, a.mode, a.det, a.channel, a.reason)
    if a.command == "clear-channel":
        if a.dry_run:
            runner.warn("--dry-run ignored for clear-channel")
        do_clear_channel(cfg, a.mode, a.det, a.channel, a.reason)
    if a.command == "resume":
        do_resume(cfg, a.dry_run)
    do_start(cfg, a.dry_run)


if __name__ == "__main__":
    main()
