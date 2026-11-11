#!/usr/bin/env python3
"""Create / show a calibration session YAML (run, calib particle, beam momentum).

Daily entry is ./calib.py <run>. This module keeps the session helpers and a
compatibility CLI:

  calib_session.py init <run> --calib-particle pion|kaon --target hodo|chamber|both
      [--momentum GEV] [--beam-particle pion|kaon|proton] [--charge ±1]
      [--write-table] [--max-iter N] [--force] [--out PATH]
  calib_session.py show <run|session.yml>

Beam momentum comes from myanalysis beam_run_condition.yml and optional UserParam
(BeamMom etc.). Missing momentum: TTY prompt, or --momentum (no neighbor guess).
--write-table upserts myanalysis table only (not analyzer BeamRunCondition.cc).
calib_particle is the calibration window and is never used as beam particle.
"""
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import beam_run_condition as brc  # noqa: E402
from common.lib import config, runner  # noqa: E402

# Compat for old session YAML / calib_session init; new flows use stage flags.
TARGETS = ("hodo", "chamber", "both", "custom")
CALIB_PARTICLES = ("pion", "kaon")
CALIB_PARTICLE_ALIASES = {
    "pion": "pion",
    "pi": "pion",
    "kaon": "kaon",
    "k": "kaon",
}
BEAM_PARTICLE_ALIASES = {
    "pion": "pion",
    "pi": "pion",
    "kaon": "kaon",
    "k": "kaon",
    "proton": "proton",
    "p": "proton",
}


#_____________________________________________________________________________
def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


#_____________________________________________________________________________
def session_path(run):
    return config.OUTPUT_DIR / "calib" / f"run{run:05d}" / "calib_session.yml"


#_____________________________________________________________________________
def normalize_calib_particle(raw):
    """Map pion|pi / kaon|k (case-insensitive) to pion|kaon."""
    key = str(raw or "").strip().lower()
    particle = CALIB_PARTICLE_ALIASES.get(key)
    if particle is None:
        runner.fail(
            f"calib_particle must be pion|pi|kaon|k, got {raw!r}"
        )
    return particle


#_____________________________________________________________________________
def normalize_beam_particle(raw):
    """Map pion|pi / kaon|k / proton|p to canonical beam particle."""
    key = str(raw or "").strip().lower()
    particle = BEAM_PARTICLE_ALIASES.get(key)
    if particle is None:
        runner.fail(
            f"beam_particle must be pion|pi|kaon|k|proton|p, got {raw!r}"
        )
    return particle


#_____________________________________________________________________________
def default_steps(target):
    if target == "hodo":
        return ["hodo"]
    chamber = ["dc_chamber", {"dc_chamber": {"bcin": True}}, "d5_align"]
    if target == "chamber":
        return chamber
    # former "both" and any other → full chain
    return ["hodo"] + chamber


#_____________________________________________________________________________
def steps_from_flags(*, hodo=False, bcout=False, bcin=False, d5=False):
    """Build session steps from boolean stage flags.

    Multiple flags may be true. d5 forces bcout+bcin on.
    No flags → full chain (hodo + bcout + bcin + d5).
    """
    if d5:
        bcout = True
        bcin = True
    any_flag = hodo or bcout or bcin or d5
    if not any_flag:
        hodo = bcout = bcin = d5 = True
    steps = []
    if hodo:
        steps.append("hodo")
    if bcout:
        steps.append("dc_chamber")
    if bcin:
        steps.append({"dc_chamber": {"bcin": True}})
    if d5:
        steps.append("d5_align")
    if steps == default_steps("hodo"):
        target = "hodo"
    elif steps == default_steps("chamber"):
        target = "chamber"
    else:
        target = "custom"
    return steps, target


#_____________________________________________________________________________
def expand_coarse_steps(
    steps,
    *,
    auto=False,
    hdphc_iter=3,
    max_iter=10,
):
    """Expand hodo/dc_chamber/d5_align into fine steps unless auto=True."""
    if auto:
        return list(steps)
    out = []
    for raw in steps:
        if isinstance(raw, str):
            name, opts = raw, {}
        elif isinstance(raw, dict) and len(raw) == 1:
            name, opts = next(iter(raw.items()))
            opts = dict(opts or {})
        else:
            out.append(raw)
            continue

        if name == "hodo":
            ftof = bool(opts.get("ftof"))
            hodo_opts = {"ftof": True} if ftof else {}
            out.append({"hodo_hdprm": dict(hodo_opts)} if hodo_opts else "hodo_hdprm")
            if not opts.get("skip_t0"):
                out.append({"hodo_t0": dict(hodo_opts)} if hodo_opts else "hodo_t0")
            if not opts.get("skip_hdphc"):
                n = int(opts.get("hdphc_iter") or hdphc_iter)
                if n < 1:
                    n = 1
                for i in range(1, n + 1):
                    round_opts = {"hdphc_round": i}
                    round_opts.update(hodo_opts)
                    out.append({"hodo_hdphc": round_opts})
                adopt = opts.get("adopt", "auto")
                adopt_opts = {"adopt": adopt}
                adopt_opts.update(hodo_opts)
                out.append({"hodo_hdphc_adopt": adopt_opts})
            continue

        if name == "dc_chamber":
            base = {}
            if opts.get("bcin"):
                base["bcin"] = True
            if opts.get("separate"):
                base["separate"] = True
            stages = []
            if not opts.get("skip_t0"):
                stages.append("dc_t0")
            if not opts.get("skip_drift"):
                stages.append("dc_drift")
            if not opts.get("skip_resi"):
                stages.append("dc_ch_resi")
            if not opts.get("skip_reso"):
                stages.append("dc_reso")
            for stage in stages:
                out.append({stage: dict(base)} if base else stage)
            continue

        if name == "d5_align":
            n = max(1, int(max_iter))
            for i in range(1, n + 1):
                out.append({"d5_pass": {"pass": i, "allow_skip": True}})
            continue

        out.append(raw)
    return out


#_____________________________________________________________________________
def merge_session_steps(
    existing_steps,
    new_steps,
    *,
    auto=False,
    hdphc_iter=3,
    max_iter=10,
):
    """Union steps by lane: keep existing lanes, append missing lanes from new."""
    from common.lib import calib_lanes

    existing_fine = expand_coarse_steps(
        list(existing_steps or []),
        auto=auto, hdphc_iter=hdphc_iter, max_iter=max_iter,
    )
    new_fine = expand_coarse_steps(
        list(new_steps or []),
        auto=auto, hdphc_iter=hdphc_iter, max_iter=max_iter,
    )
    have = calib_lanes.split_steps_by_lane(existing_fine)
    add = calib_lanes.split_steps_by_lane(new_fine)
    out = []
    for lane in calib_lanes.LANE_ORDER:
        if lane in have:
            out.extend(have[lane])
        elif lane in add:
            out.extend(add[lane])
    return out


#_____________________________________________________________________________
def write_session(path, doc):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(doc, sort_keys=False, default_flow_style=False))
    return path


#_____________________________________________________________________________
def load_session_file(path):
    path = Path(path).expanduser()
    if not path.is_absolute():
        path = (Path.cwd() / path).resolve()
    if not path.is_file():
        runner.fail(f"session not found: {path}")
    doc = yaml.safe_load(path.read_text())
    if not isinstance(doc, dict):
        runner.fail(f"session must be a mapping: {path}")
    return path, doc


#_____________________________________________________________________________
def prompt_momentum(run):
    """Ask for beam momentum on a TTY. Fail clearly when non-interactive."""
    if not sys.stdin.isatty():
        runner.fail(
            f"run {run}: momentum unknown in table/UserParam; pass --momentum "
            f"(no neighbor guess; non-interactive)"
        )
    while True:
        try:
            raw = input(f"run {run:05d}: momentum [GeV/c]? ").strip()
        except EOFError:
            runner.fail(
                f"run {run}: momentum unknown; pass --momentum (no neighbor guess)"
            )
        try:
            momentum = float(raw)
        except ValueError:
            print("enter a positive number", file=sys.stderr)
            continue
        if momentum > 0:
            return momentum
        print("momentum must be > 0", file=sys.stderr)


#_____________________________________________________________________________
def resolve_calib_particle(*, pion=False, kaon=False, calib_particle=None):
    if pion and kaon:
        runner.fail("use only one of --pion / --kaon")
    if calib_particle is not None:
        particle = normalize_calib_particle(calib_particle)
        if pion and particle != "pion":
            runner.fail("--pion disagrees with --calib-particle")
        if kaon and particle != "kaon":
            runner.fail("--kaon disagrees with --calib-particle")
        return particle
    if kaon:
        return "kaon"
    return "pion"


#_____________________________________________________________________________
def build_session_doc(
    run,
    *,
    calib_particle="pion",
    steps=None,
    target=None,
    momentum=None,
    beam_particle=None,
    charge=None,
    write_table=False,
    max_iter=10,
    force=False,
    local=False,
    auto=False,
    prompt_if_missing=True,
    hdphc_iter=3,
):
    """Resolve beam fields and return a session mapping (does not write)."""
    if steps is None:
        if target is None:
            target = "both"
        steps = default_steps(target)
    elif target is None:
        if steps == default_steps("hodo"):
            target = "hodo"
        elif steps == default_steps("chamber"):
            target = "chamber"
        elif steps == default_steps("both"):
            target = "both"
        else:
            target = "custom"
    steps = expand_coarse_steps(
        steps, auto=auto, hdphc_iter=hdphc_iter, max_iter=max_iter,
    )

    resolved = brc.resolve(run)
    runner.info(
        f"run={run:05d} calib_particle={calib_particle} target={target}"
    )
    runner.info(
        f"beam particle={resolved['particle']} ({resolved['particle_source']}) "
        f"momentum={resolved['momentum_gev_c']} ({resolved['momentum_source']}) "
        f"charge={resolved['charge']} ({resolved['charge_source']})"
    )
    if resolved["user_param_path"]:
        runner.info(f"UserParam peeked: {resolved['user_param_path']}")

    beam_particle_out = resolved["particle"]
    beam_particle_source = resolved["particle_source"]
    if beam_particle:
        beam_particle_out = normalize_beam_particle(beam_particle)
        beam_particle_source = brc.SOURCE_USER_INPUT

    momentum_out = resolved["momentum_gev_c"]
    momentum_source = resolved["momentum_source"]
    if momentum is not None:
        momentum_out = float(momentum)
        if not (momentum_out > 0):
            runner.fail("--momentum must be > 0")
        momentum_source = brc.SOURCE_USER_INPUT

    charge_out = resolved["charge"]
    charge_source = resolved["charge_source"]
    if charge is not None:
        charge_out = int(charge)
        if charge_out not in (-1, 1):
            runner.fail("--charge must be +1 or -1")
        charge_source = brc.SOURCE_USER_INPUT

    if momentum_out is None:
        if prompt_if_missing:
            momentum_out = prompt_momentum(run)
            momentum_source = brc.SOURCE_USER_INPUT
        else:
            runner.fail(
                f"run {run}: momentum unknown in table/UserParam; pass --momentum "
                f"(no neighbor guess)"
            )

    table_updated = False
    if write_table:
        if beam_particle_out is None:
            runner.fail("--write-table needs beam particle (--beam-particle)")
        if charge_out is None:
            runner.fail("--write-table needs charge (--charge)")
        try:
            path, action = brc.upsert_run(
                run, beam_particle_out, momentum_out, charge_out,
            )
        except ValueError as e:
            runner.fail(str(e))
        table_updated = action != "unchanged"
        runner.info(f"beam table {action}: {path}")

    return {
        "run": run,
        "calib_particle": calib_particle,
        "kaon": calib_particle == "kaon",
        "target": target,
        "momentum_gev_c": momentum_out,
        "momentum_source": momentum_source,
        "beam_particle": beam_particle_out,
        "beam_particle_source": beam_particle_source,
        "charge": charge_out,
        "charge_source": charge_source,
        "force": bool(force),
        "local": bool(local),
        "auto": bool(auto),
        "max_iter": max_iter,
        "steps": steps,
        "created_at": utc_now(),
        "table_updated": table_updated,
    }


#_____________________________________________________________________________
def ensure_session(
    run,
    *,
    calib_particle="pion",
    steps=None,
    target=None,
    momentum=None,
    beam_particle=None,
    charge=None,
    write_table=False,
    max_iter=10,
    force=False,
    local=False,
    auto=False,
    out=None,
    prompt_if_missing=True,
    overwrite=False,
    hdphc_iter=3,
):
    """Create or reuse session YAML. Returns (path, doc).

    force: stored for iterate (--force on run_*).
    overwrite: rewrite the session file when it already exists.
    auto: skip Slack step gates (coarse steps kept).
    """
    out_path = Path(out).expanduser() if out else session_path(run)
    if not out_path.is_absolute():
        out_path = (Path.cwd() / out_path).resolve()
    if out_path.is_file() and not overwrite:
        path, doc = load_session_file(out_path)
        if bool(local) != bool(doc.get("local", False)):
            doc["local"] = bool(local)
            write_session(path, doc)
        return path, doc
    doc = build_session_doc(
        run,
        calib_particle=calib_particle,
        steps=steps,
        target=target,
        momentum=momentum,
        beam_particle=beam_particle,
        charge=charge,
        write_table=write_table,
        max_iter=max_iter,
        force=force,
        local=local,
        auto=auto,
        prompt_if_missing=prompt_if_missing,
        hdphc_iter=hdphc_iter,
    )
    write_session(out_path, doc)
    runner.info(f"wrote session {out_path}")
    return out_path, doc


#_____________________________________________________________________________
def set_session_local(path, local):
    path, doc = load_session_file(path)
    doc["local"] = bool(local)
    write_session(path, doc)
    return path, doc


#_____________________________________________________________________________
def do_show(target):
    text = str(target)
    if text.endswith(".yml") or text.endswith(".yaml") or "/" in text:
        path, doc = load_session_file(text)
    else:
        try:
            run = int(text)
        except ValueError:
            runner.fail(f"show: expected run number or session path, got {text!r}")
        path = session_path(run)
        if not path.is_file():
            runner.fail(f"session not found: {path}")
        doc = yaml.safe_load(path.read_text())
    runner.info(f"session {path}")
    print(yaml.safe_dump(doc, sort_keys=False, default_flow_style=False))
    sys.exit(runner.EXIT_OK)


#_____________________________________________________________________________
def do_init(a):
    steps = None
    target = a.target
    if getattr(a, "steps", None) is not None:
        steps = a.steps
    out, doc = ensure_session(
        a.run,
        calib_particle=a.calib_particle,
        steps=steps,
        target=target,
        momentum=a.momentum,
        beam_particle=a.beam_particle,
        charge=a.charge,
        write_table=a.write_table,
        max_iter=a.max_iter,
        force=a.force,
        local=getattr(a, "local", False),
        out=a.out,
        prompt_if_missing=getattr(a, "prompt_if_missing", True),
        overwrite=True,
    )
    print(yaml.safe_dump(doc, sort_keys=False, default_flow_style=False))
    runner.info(f"session ready {out}")
    sys.exit(runner.EXIT_OK)


#_____________________________________________________________________________
def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument("command", choices=("init", "show"))
    p.add_argument("target", help="run number (init/show) or session.yml (show)")
    p.add_argument("--calib-particle", choices=CALIB_PARTICLES, default=None,
                   help="calibration particle window (not beam species)")
    p.add_argument("--target", dest="calib_target", choices=TARGETS, default=None,
                   help="default steps: hodo / chamber+D5 / both")
    p.add_argument("--momentum", type=float, default=None,
                   help="beam momentum [GeV/c] when table/UserParam missing")
    p.add_argument("--beam-particle", choices=brc.PARTICLES, default=None,
                   help="beam species for session / --write-table (not calib_particle)")
    p.add_argument("--charge", type=int, default=None, help="beam charge ±1")
    p.add_argument("--write-table", action="store_true",
                   help="upsert myanalysis beam_run_condition.yml (not analyzer .cc)")
    p.add_argument("--max-iter", type=int, default=10)
    p.add_argument("--force", action="store_true")
    p.add_argument("--local", action="store_true",
                   help="myrun --local (stored on session for dispatch)")
    p.add_argument("--out", type=Path, default=None, help="session output path")
    return p.parse_args()


#_____________________________________________________________________________
def main():
    runner.warn("prefer ./calib.py <run> for daily launch (calib_session is compat)")
    a = parse_args()
    if a.command == "show":
        do_show(a.target)
    if a.calib_particle is None:
        runner.fail("init requires --calib-particle pion|kaon")
    if a.calib_target is None:
        runner.fail("init requires --target hodo|chamber|both")
    try:
        a.run = int(a.target)
    except ValueError:
        runner.fail(f"init: run must be an int, got {a.target!r}")
    if a.run < 1:
        runner.fail("run must be >= 1")
    if a.max_iter < 1 or a.max_iter > 100:
        runner.fail("--max-iter must be in [1, 100]")
    a.target = a.calib_target
    a.prompt_if_missing = a.momentum is None
    # Compatibility CLI: missing momentum still fails without TTY unless prompted;
    # keep old non-interactive fail when no --momentum (tests rely on this).
    if a.momentum is None and not sys.stdin.isatty():
        a.prompt_if_missing = False
    do_init(a)


if __name__ == "__main__":
    main()
