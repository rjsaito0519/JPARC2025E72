#!/usr/bin/env python3
"""One-chamber DC chain: create_runlist → decode → fit per stage.

Usage: iterate_dc_chamber.py <run> [--bcout|--bcin] [--kaon] [--max-iter N]
                             [--separate] [--dry-run] [--force]
                             [--skip-t0] [--skip-drift] [--skip-resi] [--skip-reso]
                             [--local]

Order (wiki DC-Calibration):
  create_runlist {Pi|K}_tdc   → decode → t0
  create_runlist {Pi|K}_drift → decode → drift
  create_runlist {Pi|K}_resi  → resi loop → reso loop
"""
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parents[1]))
sys.path.insert(0, str(SCRIPT_DIR))

from common.lib import calib_runlist, config, decode_run, runner  # noqa: E402
import iterate_dc_resi as dc_fit  # noqa: E402

RUN_DCPRM = SCRIPT_DIR / "run_dcprm.py"
MYRUN = config.ANALYZER_DIR / "runmanager" / "myrun.py"
OK_DRIVER = (runner.EXIT_OK, runner.EXIT_FAIL, runner.EXIT_NOT_CONVERGED,
             runner.EXIT_CHECK, runner.EXIT_SKIP)


#_____________________________________________________________________________
def oneshot_cmd(run, mode, bcin, kaon, force):
    parts = ["python3", str(RUN_DCPRM), str(run), mode,
             "--bcin" if bcin else "--bcout"]
    if kaon:
        parts.append("--kaon")
    if force:
        parts.append("--force")
    return " ".join(parts)


#_____________________________________________________________________________
def mode_flag(bcin):
    return "--bcin" if bcin else "--bcout"


#_____________________________________________________________________________
def prepare_stage(run, bcin, kaon, stage, *, local, dry_run):
    """create_runlist + decode for tdc/drift/resi. Returns runlist path."""
    suffix = calib_runlist.dc_suffix(kaon, stage)
    calib_runlist.run_create_runlist(mode_flag(bcin), run, suffix, dry_run=dry_run)
    runlist = dc_fit.default_runlist(run, bcin)
    if dry_run:
        decode_run.run_decode(runlist, dry_run=True, local=local)
        return runlist
    if not runlist.is_file():
        runner.fail(f"runlist missing after create_runlist: {runlist}")
    decode_run.run_decode(runlist, local=local)
    return runlist


#_____________________________________________________________________________
def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int)
    g = p.add_mutually_exclusive_group()
    g.add_argument("--bcout", action="store_true", help="BcOut (default)")
    g.add_argument("--bcin", action="store_true", help="BcIn")
    p.add_argument("--kaon", action="store_true", help="K param suffix")
    p.add_argument("--max-iter", type=int, default=dc_fit.DEFAULT_MAX_ITER,
                   help=f"max resi/reso passes each (default {dc_fit.DEFAULT_MAX_ITER})")
    p.add_argument("--separate", action="store_true",
                   help="passed to reso (Res per a/b)")
    p.add_argument("--dry-run", action="store_true", help="print commands only")
    p.add_argument("--force", action="store_true", help="passed through to run_dcprm")
    p.add_argument("--skip-t0", action="store_true", help="omit t0 stage")
    p.add_argument("--skip-drift", action="store_true", help="omit drift stage")
    p.add_argument("--skip-resi", action="store_true", help="omit residual iterate loop")
    p.add_argument("--skip-reso", action="store_true", help="omit resolution iterate loop")
    p.add_argument("--local", action="store_true",
                   help="myrun --local (default: bsub)")
    a = p.parse_args()
    if a.max_iter < 1:
        runner.fail("--max-iter must be >= 1")
    a.bcin = bool(a.bcin)
    return a


#_____________________________________________________________________________
def run_oneshot(run, mode, bcin, kaon, force, dry_run):
    """Run t0 or drift only (decode already done by prepare_stage)."""
    cmd = oneshot_cmd(run, mode, bcin, kaon, force)
    runner.step(f"iterate_dc_chamber {mode}")
    if dry_run:
        print(f"[dry-run] {cmd}")
        return runner.EXIT_OK
    rc, _ = runner.sh(cmd, ok=OK_DRIVER)
    runner.info(f"run_dcprm {mode} exit={rc}")
    return rc


#_____________________________________________________________________________
def main():
    a = parse_args()

    if not a.dry_run:
        if not MYRUN.is_file():
            runner.fail(f"myrun not found: {MYRUN}")
        if not RUN_DCPRM.is_file():
            runner.fail(f"run_dcprm not found: {RUN_DCPRM}")

    runner.info(
        f"iterate_dc_chamber run={a.run:05d} {dc_fit.chamber_tag(a.bcin)} "
        f"max_iter={a.max_iter}"
        + (" dry-run" if a.dry_run else "")
        + (" --local" if a.local else "")
        + (" --skip-t0" if a.skip_t0 else "")
        + (" --skip-drift" if a.skip_drift else "")
        + (" --skip-resi" if a.skip_resi else "")
        + (" --skip-reso" if a.skip_reso else "")
    )
    runner.info("chain includes create_runlist per stage")

    runlist = None
    if not a.skip_t0:
        runlist = prepare_stage(
            a.run, a.bcin, a.kaon, "tdc", local=a.local, dry_run=a.dry_run,
        )
        rc = run_oneshot(a.run, "t0", a.bcin, a.kaon, a.force, a.dry_run)
        if rc != runner.EXIT_OK:
            sys.exit(rc)

    if not a.skip_drift:
        runlist = prepare_stage(
            a.run, a.bcin, a.kaon, "drift", local=a.local, dry_run=a.dry_run,
        )
        rc = run_oneshot(a.run, "drift", a.bcin, a.kaon, a.force, a.dry_run)
        if rc != runner.EXIT_OK:
            sys.exit(rc)

    if not a.skip_resi or not a.skip_reso:
        runlist = prepare_stage(
            a.run, a.bcin, a.kaon, "resi", local=a.local, dry_run=a.dry_run,
        )

    if not a.skip_resi:
        rc = dc_fit.run_fit_loop(
            a.run, "resi", a.bcin, a.kaon, a.force, a.max_iter, runlist, a.dry_run,
            separate=False, do_preflight=False, label="iterate_dc_chamber/resi",
            local=a.local,
        )
        if rc != runner.EXIT_OK:
            sys.exit(rc)

    if not a.skip_reso:
        rc = dc_fit.run_fit_loop(
            a.run, "reso", a.bcin, a.kaon, a.force, a.max_iter, runlist, a.dry_run,
            separate=a.separate, do_preflight=False, label="iterate_dc_chamber/reso",
            local=a.local,
        )
        if rc != runner.EXIT_OK:
            sys.exit(rc)

    runner.info("done: chamber chain EXIT_OK")
    sys.exit(runner.EXIT_OK)


if __name__ == "__main__":
    main()
