#!/usr/bin/env python3
"""Dispatch Hodo / DC / D5 iterate scripts from one YAML (see wiki Calibration-Automation).

Usage:
  calib_dispatch.py start  <config.yml> [--dry-run]
  calib_dispatch.py status <config.yml>

Steps: hodo | dc_resi [{bcin: true}] | d5_align
  Optional per-step: allow_skip: true  (continue after EXIT_SKIP)
State: OUTPUT_DIR/calib/runXXXXX/calib_dispatch_state.json (no stop/resume yet).
Continue on exit 0; EXIT_SKIP (4) stops unless that step has allow_skip: true.
Final exit code is runner.worst of executed step codes (so a skip is not erased by a later 0).
"""
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import config, runner  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
KNOWN_STEPS = ("hodo", "dc_resi", "d5_align")
KNOWN_OPTS = {"bcin", "skip_t0", "skip_hdphc", "ftof", "allow_skip"}


#_____________________________________________________________________________
def state_path(run):
    return config.OUTPUT_DIR / "calib" / f"run{run:05d}" / "calib_dispatch_state.json"


#_____________________________________________________________________________
def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


#_____________________________________________________________________________
def require_bool(value, where):
    """Accept only a YAML bool (not a string)."""
    if not isinstance(value, bool):
        runner.fail(f"{where}: must be a YAML bool, got {value!r}")
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
    if "run" not in doc:
        runner.fail("config missing required key: run")
    try:
        run = int(doc["run"])
    except (TypeError, ValueError):
        runner.fail(f"config run must be int: {doc.get('run')!r}")
    steps = doc.get("steps")
    if not isinstance(steps, list) or not steps:
        runner.fail("config steps must be a non-empty list")
    return path, {
        "run": run,
        "kaon": bool(doc.get("kaon", False)),
        "force": bool(doc.get("force", False)),
        "max_iter": int(doc.get("max_iter", 10)),
        "steps": steps,
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
    if "allow_skip" in opts:
        require_bool(opts["allow_skip"], f"steps[{index}].allow_skip")
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
        "d5_align": "iterate_d5_align.py",
    }[name]), str(run)]

    if name == "dc_resi":
        if opts.get("bcin"):
            parts.append("--bcin")
        else:
            parts.append("--bcout")
        parts += ["--max-iter", str(cfg["max_iter"])]
    elif name == "d5_align":
        parts += ["--max-iter", str(cfg["max_iter"])]
    elif name == "hodo":
        if opts.get("skip_t0"):
            parts.append("--skip-t0")
        if opts.get("skip_hdphc"):
            parts.append("--skip-hdphc")
        if opts.get("ftof"):
            parts.append("--ftof")

    if cfg["kaon"] and name in ("hodo", "dc_resi", "d5_align"):
        parts.append("--kaon")
    if cfg["force"]:
        parts.append("--force")
    if dry_run:
        parts.append("--dry-run")

    unknown = set(opts) - KNOWN_OPTS
    if unknown:
        runner.fail(f"step {name}: unknown option(s) {sorted(unknown)}")
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
def make_payload(run, index, name, opts, rc, finished, steps_hist):
    return {
        "run": run,
        "index": index,
        "step": name,
        "opts": opts,
        "rc": rc,
        "finished": finished,
        "steps": steps_hist,
        "updated_at": utc_now(),
    }


#_____________________________________________________________________________
def do_status(cfg):
    st = read_state(cfg["run"])
    if st is None:
        runner.info(f"run{cfg['run']:05d}: not started (no {state_path(cfg['run'])})")
        sys.exit(runner.EXIT_OK)
    runner.info(
        f"run{cfg['run']:05d}: last step[{st.get('index')}]={st.get('step')} "
        f"rc={st.get('rc')} finished={st.get('finished')} "
        f"updated_at={st.get('updated_at')}"
    )
    print(json.dumps(st, indent=2, sort_keys=True))
    sys.exit(runner.EXIT_OK)


#_____________________________________________________________________________
def do_start(cfg, dry_run):
    runner.info(
        f"calib_dispatch start run={cfg['run']:05d} n_steps={len(cfg['steps'])}"
        + (" dry-run" if dry_run else "")
    )
    steps_hist = []
    last_rc = runner.EXIT_OK
    last_name = None
    last_opts = {}
    last_index = -1
    for i, raw in enumerate(cfg["steps"]):
        name, opts = normalize_step(raw, i)
        cmd = step_cmd(name, opts, cfg, dry_run)
        runner.step(f"calib_dispatch {i + 1}/{len(cfg['steps'])}: {name}")
        if dry_run:
            print(f"[dry-run] {cmd}")
            rc = runner.EXIT_OK
        else:
            rc, _ = runner.sh(
                cmd,
                ok=(runner.EXIT_OK, runner.EXIT_FAIL, runner.EXIT_NOT_CONVERGED,
                    runner.EXIT_CHECK, runner.EXIT_SKIP),
            )
        last_rc = rc
        last_name = name
        last_opts = opts
        last_index = i
        steps_hist.append({"index": i, "step": name, "opts": opts, "rc": rc})
        if not dry_run:
            write_state(
                cfg["run"],
                make_payload(cfg["run"], i, name, opts, rc, False, steps_hist),
            )
        runner.info(f"step {name} exit={rc}")
        if not should_continue(rc, opts):
            runner.warn(f"stopping after step {name} (rc={rc})")
            sys.exit(rc)

    overall_rc = runner.worst([e["rc"] for e in steps_hist]) if steps_hist else last_rc
    if not dry_run:
        write_state(
            cfg["run"],
            make_payload(
                cfg["run"], last_index, last_name, last_opts, overall_rc, True, steps_hist
            ),
        )
    runner.info("calib_dispatch finished")
    sys.exit(overall_rc)


#_____________________________________________________________________________
def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument("command", choices=("start", "status"))
    p.add_argument("config", type=Path, help="YAML config path")
    p.add_argument("--dry-run", action="store_true",
                   help="start: print iterate commands only (no state write)")
    return p.parse_args()


#_____________________________________________________________________________
def main():
    a = parse_args()
    _, cfg = load_config(a.config)
    # Validate all steps up front (fail before running anything).
    for i, raw in enumerate(cfg["steps"]):
        name, opts = normalize_step(raw, i)
        step_cmd(name, opts, cfg, dry_run=True)
    if a.command == "status":
        if a.dry_run:
            runner.warn("--dry-run ignored for status")
        do_status(cfg)
    do_start(cfg, a.dry_run)


if __name__ == "__main__":
    main()
