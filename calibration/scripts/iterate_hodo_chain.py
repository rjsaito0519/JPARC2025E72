#!/usr/bin/env python3
"""One-pass Hodo calibration chain for one run: hdprm → decode → t0 → decode → hdphc → decode.

Usage: iterate_hodo_chain.py <run> [--kaon] [--ftof] [--force] [--dry-run]
                             [--hodo-runlist PATH] [--skip-t0] [--skip-hdphc]

Not a convergence loop: Hodo update_param writes once (empty states). This script only
runs the usual order and re-decodes Hodo after each write so the next fit sees new params.

Does not run htof (needs TPCHelixHTOF regenerate). Does not modify run_hodo.py.
Hodo runlist: --hodo-runlist or create_runlist's
runmanager/runlist/<sub>/hodo_runXXXXX.yml.
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from lib import config, runner  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
RUN_HODO = SCRIPT_DIR / "run_hodo.py"
MYRUN = config.ANALYZER_DIR / "runmanager" / "myrun.py"
RUNLIST_DIR = config.ANALYZER_DIR / "runmanager" / "runlist"
OK_DRIVER = (runner.EXIT_OK, runner.EXIT_FAIL, runner.EXIT_NOT_CONVERGED,
             runner.EXIT_CHECK, runner.EXIT_SKIP)


#_____________________________________________________________________________
def default_hodo_runlist(run):
    return RUNLIST_DIR / config.SUB_DIR / f"hodo_run{run:05d}.yml"


#_____________________________________________________________________________
def resolve_hodo_runlist(run, explicit, require_exists=True):
    if explicit is not None:
        path = Path(explicit).expanduser()
        if not path.is_absolute():
            path = (Path.cwd() / path).resolve()
        if require_exists and not path.is_file():
            runner.fail(f"hodo runlist not found: {path}")
        return path
    path = default_hodo_runlist(run)
    if path.is_file() or not require_exists:
        return path
    runner.fail(
        f"hodo runlist not found: {path} "
        f"(create_runlist.py {run} --hodo, or pass --hodo-runlist)"
    )


#_____________________________________________________________________________
def decode_cmd(runlist):
    # myrun expects cwd = analyzer root; pass an absolute runlist (same as iterate_tpc_calib).
    path = Path(runlist).expanduser().resolve()
    return f"cd {config.ANALYZER_DIR} && ./runmanager/myrun.py {path} --local"


#_____________________________________________________________________________
def hodo_cmd(run, mode, kaon, ftof, force):
    parts = ["python3", str(RUN_HODO), str(run), mode]
    if kaon:
        parts.append("--kaon")
    if ftof:
        parts.append("--ftof")
    if force:
        parts.append("--force")
    return " ".join(parts)


#_____________________________________________________________________________
def build_steps(run, runlist, kaon, ftof, force, skip_t0, skip_hdphc):
    """Ordered (label, cmd) pairs for the chain."""
    steps = [("hdprm", hodo_cmd(run, "hdprm", kaon, ftof, force)),
             ("decode", decode_cmd(runlist))]
    if not skip_t0:
        steps.append(("t0", hodo_cmd(run, "t0", kaon, ftof, force)))
        steps.append(("decode", decode_cmd(runlist)))
    if not skip_hdphc:
        steps.append(("hdphc", hodo_cmd(run, "hdphc", kaon, ftof, force)))
        steps.append(("decode", decode_cmd(runlist)))
    return steps


#_____________________________________________________________________________
def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int)
    p.add_argument("--kaon", action="store_true", help="K param suffix")
    p.add_argument("--ftof", action="store_true",
                   help="passed to run_hodo hdprm/hdphc (FTOF detector set)")
    p.add_argument("--force", action="store_true", help="passed through to run_hodo")
    p.add_argument("--dry-run", action="store_true", help="print commands only")
    p.add_argument("--hodo-runlist", type=Path, default=None,
                   help="Hodo myrun yml (default: runmanager/runlist/<sub>/hodo_runXXXXX.yml)")
    p.add_argument("--skip-t0", action="store_true", help="omit t0 and its following decode")
    p.add_argument("--skip-hdphc", action="store_true",
                   help="omit hdphc and its following decode")
    return p.parse_args()


#_____________________________________________________________________________
def main():
    a = parse_args()
    runlist = resolve_hodo_runlist(a.run, a.hodo_runlist, require_exists=not a.dry_run)
    if not MYRUN.is_file() and not a.dry_run:
        runner.fail(f"myrun not found: {MYRUN}")
    if not RUN_HODO.is_file():
        runner.fail(f"run_hodo not found: {RUN_HODO}")

    steps = build_steps(a.run, runlist, a.kaon, a.ftof, a.force, a.skip_t0, a.skip_hdphc)
    runner.info(
        f"iterate_hodo_chain run={a.run:05d} steps={len(steps)} "
        f"runlist={runlist}"
        + (" --kaon" if a.kaon else "")
        + (" --ftof" if a.ftof else "")
        + (" --skip-t0" if a.skip_t0 else "")
        + (" --skip-hdphc" if a.skip_hdphc else "")
        + (" dry-run" if a.dry_run else "")
    )

    if a.dry_run:
        runner.step("iterate_hodo_chain dry-run")
        for label, cmd in steps:
            print(f"[dry-run] ({label}) {cmd}")
        sys.exit(runner.EXIT_OK)

    codes = []
    for i, (label, cmd) in enumerate(steps, 1):
        runner.step(f"iterate_hodo_chain {i}/{len(steps)}: {label}")
        if label == "decode":
            runner.sh(cmd, ok=(runner.EXIT_OK,))
            continue
        rc, _ = runner.sh(cmd, ok=OK_DRIVER)
        runner.info(f"run_hodo {label} exit={rc}")
        if rc == runner.EXIT_FAIL:
            sys.exit(runner.EXIT_FAIL)
        codes.append(rc)

    final = runner.worst(codes) if codes else runner.EXIT_OK
    runner.info(f"done: worst exit={final}")
    sys.exit(final)


if __name__ == "__main__":
    main()
