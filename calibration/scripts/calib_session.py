#!/usr/bin/env python3
"""Create / show a calibration session YAML (run, calib particle, beam momentum).

Usage:
  calib_session.py init <run> --calib-particle pion|kaon --target hodo|chamber|both
      [--momentum GEV] [--beam-particle pion|kaon|proton] [--charge ±1]
      [--write-table] [--max-iter N] [--force] [--out PATH]
  calib_session.py show <run|session.yml>

Beam momentum comes from myanalysis beam_run_condition.yml and optional UserParam
(BeamMom etc.). Missing momentum requires --momentum (no neighbor guess).
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

TARGETS = ("hodo", "chamber", "both")
CALIB_PARTICLES = ("pion", "kaon")


#_____________________________________________________________________________
def utc_now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


#_____________________________________________________________________________
def session_path(run):
    return config.OUTPUT_DIR / "calib" / f"run{run:05d}" / "calib_session.yml"


#_____________________________________________________________________________
def default_steps(target):
    if target == "hodo":
        return ["hodo"]
    chamber = ["dc_chamber", {"dc_chamber": {"bcin": True}}, "d5_align"]
    if target == "chamber":
        return chamber
    return ["hodo"] + chamber


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
    resolved = brc.resolve(a.run)
    runner.info(
        f"run={a.run:05d} calib_particle={a.calib_particle} target={a.target}"
    )
    runner.info(
        f"beam particle={resolved['particle']} ({resolved['particle_source']}) "
        f"momentum={resolved['momentum_gev_c']} ({resolved['momentum_source']}) "
        f"charge={resolved['charge']} ({resolved['charge_source']})"
    )
    if resolved["user_param_path"]:
        runner.info(f"UserParam peeked: {resolved['user_param_path']}")

    beam_particle = resolved["particle"]
    beam_particle_source = resolved["particle_source"]
    if a.beam_particle:
        beam_particle = a.beam_particle
        beam_particle_source = brc.SOURCE_USER_INPUT

    momentum = resolved["momentum_gev_c"]
    momentum_source = resolved["momentum_source"]
    if a.momentum is not None:
        momentum = float(a.momentum)
        if not (momentum > 0):
            runner.fail("--momentum must be > 0")
        momentum_source = brc.SOURCE_USER_INPUT

    charge = resolved["charge"]
    charge_source = resolved["charge_source"]
    if a.charge is not None:
        charge = int(a.charge)
        if charge not in (-1, 1):
            runner.fail("--charge must be +1 or -1")
        charge_source = brc.SOURCE_USER_INPUT

    if momentum is None:
        runner.fail(
            f"run {a.run}: momentum unknown in table/UserParam; pass --momentum "
            f"(no neighbor guess)"
        )

    table_updated = False
    if a.write_table:
        if beam_particle is None:
            runner.fail("--write-table needs beam particle (--beam-particle)")
        if charge is None:
            runner.fail("--write-table needs charge (--charge)")
        try:
            path, action = brc.upsert_run(
                a.run, beam_particle, momentum, charge,
            )
        except ValueError as e:
            runner.fail(str(e))
        table_updated = action != "unchanged"
        runner.info(f"beam table {action}: {path}")

    doc = {
        "run": a.run,
        "calib_particle": a.calib_particle,
        "kaon": a.calib_particle == "kaon",
        "target": a.target,
        "momentum_gev_c": momentum,
        "momentum_source": momentum_source,
        "beam_particle": beam_particle,
        "beam_particle_source": beam_particle_source,
        "charge": charge,
        "charge_source": charge_source,
        "force": bool(a.force),
        "max_iter": a.max_iter,
        "steps": default_steps(a.target),
        "created_at": utc_now(),
        "table_updated": table_updated,
    }
    out = Path(a.out).expanduser() if a.out else session_path(a.run)
    if not out.is_absolute():
        out = (Path.cwd() / out).resolve()
    write_session(out, doc)
    runner.info(f"wrote session {out}")
    print(yaml.safe_dump(doc, sort_keys=False, default_flow_style=False))
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
    p.add_argument("--out", type=Path, default=None, help="session output path")
    return p.parse_args()


#_____________________________________________________________________________
def main():
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
    do_init(a)


if __name__ == "__main__":
    main()
