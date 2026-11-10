#!/usr/bin/env python3
"""Iterate DC residual calibration for one run: fit/update → decode → repeat.

Usage: iterate_dc_resi.py <run> [--bcout|--bcin] [--kaon] [--max-iter N]
                         [--runlist PATH] [--dry-run] [--force]

Wraps run_dcprm.py resi and, on EXIT_NOT_CONVERGED (2), re-decodes the chamber
with runmanager/myrun.py --local, then runs resi again. Stops on OK (0), FAIL (1),
CHECK (3), SKIP (4), or when --max-iter is reached (exits CHECK).

Before the first resi, if the chamber ROOT is missing or calib params are newer than
the decode (non-interactive stale check), decode once first — do not rely on --force.

Does not create runlists; pass --runlist or use create_runlist's
runmanager/runlist/<sub>/{bcout|bcin}_runXXXXX.yml.

Shared fit-loop helpers (resi/reso) live in this module for iterate_dc_reso /
iterate_dc_chamber.
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import config, decode_run, input_check, runner  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
RUN_DCPRM = SCRIPT_DIR / "run_dcprm.py"
MYRUN = config.ANALYZER_DIR / "runmanager" / "myrun.py"
DEFAULT_MAX_ITER = 10
FIT_MODES = ("resi", "reso")

# What to do after run_dcprm exits (before max-iter handling).
STOP_OK = "stop_ok"
STOP_FAIL = "stop_fail"
STOP_CHECK = "stop_check"
STOP_SKIP = "stop_skip"
DO_DECODE = "decode"


#_____________________________________________________________________________
def next_action(rc):
    """Map run_dcprm exit code to the next loop action."""
    if rc == runner.EXIT_OK:
        return STOP_OK
    if rc == runner.EXIT_FAIL:
        return STOP_FAIL
    if rc == runner.EXIT_CHECK:
        return STOP_CHECK
    if rc == runner.EXIT_SKIP:
        return STOP_SKIP
    if rc == runner.EXIT_NOT_CONVERGED:
        return DO_DECODE
    return STOP_FAIL


#_____________________________________________________________________________
def chamber_tag(bcin):
    return "BcIn" if bcin else "BcOut"


#_____________________________________________________________________________
def yml_prefix(bcin):
    return "bcin" if bcin else "bcout"


#_____________________________________________________________________________
def default_runlist(run, bcin):
    """create_runlist single-run path under runmanager/runlist/<SUB_DIR>/."""
    name = f"{yml_prefix(bcin)}_run{run:05d}.yml"
    return config.ANALYZER_DIR / "runmanager" / "runlist" / config.SUB_DIR / name


#_____________________________________________________________________________
def resolve_runlist(run, bcin, explicit, require_exists=True):
    if explicit is not None:
        path = Path(explicit).expanduser()
        if not path.is_absolute():
            path = (Path.cwd() / path).resolve()
        if require_exists and not path.is_file():
            runner.fail(f"runlist not found: {path}")
        return path
    path = default_runlist(run, bcin)
    if path.is_file() or not require_exists:
        return path
    runner.fail(
        f"runlist not found: {path} "
        f"(create_runlist.py {run} --{yml_prefix(bcin)}, or pass --runlist)"
    )


#_____________________________________________________________________________
def chamber_root(run, bcin):
    return config.DATA_DIR / f"run{run:05d}_{chamber_tag(bcin)}.root"


#_____________________________________________________________________________
def root_mtime_msg(run, bcin):
    path = chamber_root(run, bcin)
    if not path.exists():
        return f"{path.name}: missing"
    return f"{path.name} mtime_ns={path.stat().st_mtime_ns}"


#_____________________________________________________________________________
def fit_cmd(run, mode, bcin, kaon, force, separate=False):
    if mode not in FIT_MODES:
        runner.fail(f"unknown fit mode: {mode!r}")
    parts = ["python3", str(RUN_DCPRM), str(run), mode,
             "--bcin" if bcin else "--bcout"]
    if mode == "reso" and separate:
        parts.append("--separate")
    if kaon:
        parts.append("--kaon")
    if force:
        parts.append("--force")
    return " ".join(parts)


#_____________________________________________________________________________
def resi_cmd(run, bcin, kaon, force):
    """Backward-compatible alias for fit_cmd(..., mode='resi')."""
    return fit_cmd(run, "resi", bcin, kaon, force)


#_____________________________________________________________________________
def run_fit_loop(
    run, mode, bcin, kaon, force, max_iter, runlist, dry_run,
    separate=False, do_preflight=True, label=None, local=False,
):
    """Run resi/reso until terminal exit. Returns the exit code (does not sys.exit)."""
    if mode not in FIT_MODES:
        runner.fail(f"unknown fit mode: {mode!r}")
    if max_iter < 1:
        runner.fail("--max-iter must be >= 1")
    tag = label or f"iterate_dc_{mode}"
    runner.info(
        f"{tag} run={run:05d} {chamber_tag(bcin)} mode={mode} "
        f"max_iter={max_iter} runlist={runlist}"
        + (" dry-run" if dry_run else "")
        + (" --local" if local else "")
    )

    if dry_run:
        runner.step(f"{tag} dry-run")
        if do_preflight:
            print("[dry-run] preflight: decode if chamber ROOT missing or params newer than decode")
        print(f"[dry-run] {fit_cmd(run, mode, bcin, kaon, force, separate)}")
        decode_run.run_decode(runlist, dry_run=True, local=local)
        print(f"[dry-run] … repeat {mode} until exit 0/1/3/4 or max-iter={max_iter}")
        runner.info(f"chamber root: {root_mtime_msg(run, bcin)}")
        return runner.EXIT_OK

    if not MYRUN.is_file():
        runner.fail(f"myrun not found: {MYRUN}")
    if not RUN_DCPRM.is_file():
        runner.fail(f"run_dcprm not found: {RUN_DCPRM}")

    if do_preflight:
        suffix = "K" if kaon else "Pi"
        root = chamber_root(run, bcin)
        stale = input_check.stale_params(root, suffix) if root.exists() else []
        if not root.exists() or stale:
            reason = "missing" if not root.exists() else ("stale vs " + "; ".join(stale))
            runner.info(f"preflight decode ({reason})")
            decode_run.run_decode(runlist, local=local)
            runner.info(f"after preflight decode: {root_mtime_msg(run, bcin)}")

    ok_rc = (runner.EXIT_OK, runner.EXIT_FAIL, runner.EXIT_NOT_CONVERGED,
             runner.EXIT_CHECK, runner.EXIT_SKIP)
    last_rc = runner.EXIT_CHECK

    for i in range(1, max_iter + 1):
        runner.step(f"{tag} pass {i}/{max_iter}")
        runner.info(f"before {mode}: {root_mtime_msg(run, bcin)}")
        cmd = fit_cmd(run, mode, bcin, kaon, force, separate)
        rc, _ = runner.sh(cmd, ok=ok_rc)
        last_rc = rc
        action = next_action(rc)
        runner.info(f"run_dcprm {mode} exit={rc} → {action}")

        if action == STOP_OK:
            runner.info(f"{mode} converged; done")
            return runner.EXIT_OK
        if action == STOP_FAIL:
            return runner.EXIT_FAIL
        if action == STOP_CHECK:
            runner.warn("check / settled; stop")
            return runner.EXIT_CHECK
        if action == STOP_SKIP:
            runner.warn("skip; stop")
            return runner.EXIT_SKIP
        # DO_DECODE
        if i >= max_iter:
            runner.warn(f"max-iter {max_iter} reached after NOT_CONVERGED; stop as CHECK")
            return runner.EXIT_CHECK

        runner.info(f"before decode: {root_mtime_msg(run, bcin)}")
        decode_run.run_decode(runlist, local=local)
        runner.info(f"after decode: {root_mtime_msg(run, bcin)}")

    runner.warn("loop ended without terminal state; stop as CHECK")
    if last_rc in (runner.EXIT_OK, runner.EXIT_FAIL, runner.EXIT_CHECK, runner.EXIT_SKIP):
        return last_rc
    return runner.EXIT_CHECK


#_____________________________________________________________________________
def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int)
    g = p.add_mutually_exclusive_group()
    g.add_argument("--bcout", action="store_true", help="BcOut (default)")
    g.add_argument("--bcin", action="store_true", help="BcIn")
    p.add_argument("--kaon", action="store_true", help="K param suffix")
    p.add_argument("--max-iter", type=int, default=DEFAULT_MAX_ITER,
                   help=f"max resi passes (default {DEFAULT_MAX_ITER})")
    p.add_argument("--runlist", type=Path, default=None,
                   help="myrun yml (default: runmanager/runlist/<sub>/{bcout|bcin}_runXXXXX.yml)")
    p.add_argument("--dry-run", action="store_true", help="print commands only")
    p.add_argument("--force", action="store_true", help="passed through to run_dcprm")
    p.add_argument("--local", action="store_true",
                   help="myrun --local (default: bsub)")
    a = p.parse_args()
    if a.max_iter < 1:
        runner.fail("--max-iter must be >= 1")
    a.bcin = bool(a.bcin)
    return a


#_____________________________________________________________________________
def main():
    a = parse_args()
    runlist = resolve_runlist(a.run, a.bcin, a.runlist)
    rc = run_fit_loop(
        a.run, "resi", a.bcin, a.kaon, a.force, a.max_iter, runlist, a.dry_run,
        separate=False, do_preflight=True, label="iterate_dc_resi",
        local=a.local,
    )
    sys.exit(rc)


if __name__ == "__main__":
    main()
