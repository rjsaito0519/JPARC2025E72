#!/usr/bin/env python3
"""Iterate D5 BLC1 alignment for one run: ref → D5 decode → pos → BcIn → resi → …

Usage: iterate_d5_align.py <run> [--kaon] [--max-iter N] [--tol T]
                           [--bcin-runlist PATH] [--dry-run] [--force]

Wraps run_d5.py {ref,pos} and the decode steps that normally sit between them.
On pos EXIT_NOT_CONVERGED (2), re-decodes BcIn, runs run_dcprm.py resi --bcin,
re-decodes BcIn again if resi also wrote (exit 2), then starts the next outer
pass with ref. Stops on OK (0), FAIL (1), CHECK (3), SKIP (4), or when
--max-iter is reached after a written pos (exits CHECK).

Does not create runlists. BcIn defaults to create_runlist's
runmanager/runlist/<sub>/bcin_runXXXXX.yml. D5 scratch runlists are those
written by run_d5.py ref under runmanager/runlist/scratch/d5/.
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from lib import config, runner  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
RUN_D5 = SCRIPT_DIR / "run_d5.py"
RUN_DCPRM = SCRIPT_DIR / "run_dcprm.py"
MYRUN = config.ANALYZER_DIR / "runmanager" / "myrun.py"
RUNLIST_DIR = config.ANALYZER_DIR / "runmanager" / "runlist"
DEFAULT_MAX_ITER = 10
DEFAULT_TOL = 0.01
TAG = "bc1reso"

STOP_OK = "stop_ok"
STOP_FAIL = "stop_fail"
STOP_CHECK = "stop_check"
STOP_SKIP = "stop_skip"
DO_NEXT = "do_next"


#_____________________________________________________________________________
def next_action(rc):
    """Map run_d5 pos (or resi) exit code to the next loop action."""
    if rc == runner.EXIT_OK:
        return STOP_OK
    if rc == runner.EXIT_FAIL:
        return STOP_FAIL
    if rc == runner.EXIT_CHECK:
        return STOP_CHECK
    if rc == runner.EXIT_SKIP:
        return STOP_SKIP
    if rc == runner.EXIT_NOT_CONVERGED:
        return DO_NEXT
    return STOP_FAIL


#_____________________________________________________________________________
def suffix_of(kaon):
    return "K" if kaon else "Pi"


#_____________________________________________________________________________
def runs_label(run, s):
    """Same naming as run_d5.runs_label for a single run."""
    label = f"{run:05d}"
    return label if s == "K" else f"{label}_{s}"


#_____________________________________________________________________________
def d5_scratch_runlists(run, s):
    out_dir = RUNLIST_DIR / "scratch" / "d5"
    label = runs_label(run, s)
    return (
        out_dir / f"d5_{label}.yml",
        out_dir / f"d5{TAG}_{label}.yml",
    )


#_____________________________________________________________________________
def default_bcin_runlist(run):
    return RUNLIST_DIR / config.SUB_DIR / f"bcin_run{run:05d}.yml"


#_____________________________________________________________________________
def resolve_bcin_runlist(run, explicit, require_exists=True):
    if explicit is not None:
        path = Path(explicit).expanduser()
        if not path.is_absolute():
            path = (Path.cwd() / path).resolve()
        if require_exists and not path.is_file():
            runner.fail(f"bcin runlist not found: {path}")
        return path
    path = default_bcin_runlist(run)
    if path.is_file() or not require_exists:
        return path
    runner.fail(
        f"bcin runlist not found: {path} "
        f"(create_runlist.py {run} --bcin, or pass --bcin-runlist)"
    )


#_____________________________________________________________________________
def decode_cmd(runlist):
    # myrun expects cwd = analyzer root (same as iterate_dc_resi / iterate_tpc_calib).
    try:
        rel = runlist.relative_to(config.ANALYZER_DIR / "runmanager")
    except ValueError:
        rel = runlist
    return f"cd {config.ANALYZER_DIR} && ./runmanager/myrun.py {rel} --local"


#_____________________________________________________________________________
def ref_cmd(run, kaon):
    parts = ["python3", str(RUN_D5), str(run), "ref"]
    if kaon:
        parts.append("--kaon")
    return " ".join(parts)


#_____________________________________________________________________________
def pos_cmd(run, kaon, tol, force):
    parts = ["python3", str(RUN_D5), str(run), "pos", f"--tol={tol}"]
    if kaon:
        parts.append("--kaon")
    if force:
        parts.append("--force")
    return " ".join(parts)


#_____________________________________________________________________________
def resi_cmd(run, kaon, force):
    parts = ["python3", str(RUN_DCPRM), str(run), "resi", "--bcin"]
    if kaon:
        parts.append("--kaon")
    if force:
        parts.append("--force")
    return " ".join(parts)


#_____________________________________________________________________________
def require_scratch_runlists(run, s, dry_run):
    paths = d5_scratch_runlists(run, s)
    if dry_run:
        return paths
    missing = [p for p in paths if not p.is_file()]
    if missing:
        runner.fail(
            "D5 scratch runlist(s) missing after ref: "
            + ", ".join(str(p) for p in missing)
        )
    return paths


#_____________________________________________________________________________
def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int)
    p.add_argument("--kaon", action="store_true", help="K param suffix")
    p.add_argument("--max-iter", type=int, default=DEFAULT_MAX_ITER,
                   help=f"max outer ref→pos passes (default {DEFAULT_MAX_ITER})")
    p.add_argument("--tol", type=float, default=DEFAULT_TOL,
                   help=f"passed to run_d5 pos (default {DEFAULT_TOL})")
    p.add_argument("--bcin-runlist", type=Path, default=None,
                   help="BcIn myrun yml (default: runmanager/runlist/<sub>/bcin_runXXXXX.yml)")
    p.add_argument("--dry-run", action="store_true", help="print commands only")
    p.add_argument("--force", action="store_true",
                   help="passed through to run_d5 pos and run_dcprm resi")
    a = p.parse_args()
    if a.max_iter < 1:
        runner.fail("--max-iter must be >= 1")
    if a.tol <= 0:
        runner.fail("--tol must be > 0")
    return a


#_____________________________________________________________________________
def stop_with(action, rc_fallback=None):
    if action == STOP_OK:
        runner.info("converged; done")
        sys.exit(runner.EXIT_OK)
    if action == STOP_FAIL:
        sys.exit(runner.EXIT_FAIL)
    if action == STOP_CHECK:
        runner.warn("check / settled; stop")
        sys.exit(runner.EXIT_CHECK)
    if action == STOP_SKIP:
        runner.warn("skip; stop")
        sys.exit(runner.EXIT_SKIP)
    if rc_fallback is not None:
        sys.exit(rc_fallback)
    sys.exit(runner.EXIT_CHECK)


#_____________________________________________________________________________
def main():
    a = parse_args()
    s = suffix_of(a.kaon)
    bcin_runlist = resolve_bcin_runlist(a.run, a.bcin_runlist, require_exists=not a.dry_run)
    if not MYRUN.is_file() and not a.dry_run:
        runner.fail(f"myrun not found: {MYRUN}")
    if not RUN_D5.is_file():
        runner.fail(f"run_d5 not found: {RUN_D5}")
    if not RUN_DCPRM.is_file():
        runner.fail(f"run_dcprm not found: {RUN_DCPRM}")

    p_norm, p_tag = d5_scratch_runlists(a.run, s)
    runner.info(
        f"iterate_d5_align run={a.run:05d} suffix={s} "
        f"max_iter={a.max_iter} tol={a.tol} bcin_runlist={bcin_runlist}"
        + (" dry-run" if a.dry_run else "")
    )

    if a.dry_run:
        runner.step("iterate_d5_align dry-run")
        print(f"[dry-run] {ref_cmd(a.run, a.kaon)}")
        print(f"[dry-run] {decode_cmd(p_norm)}")
        print(f"[dry-run] {decode_cmd(p_tag)}")
        print(f"[dry-run] {pos_cmd(a.run, a.kaon, a.tol, a.force)}")
        print(f"[dry-run] {decode_cmd(bcin_runlist)}   # if pos exit 2")
        print(f"[dry-run] {resi_cmd(a.run, a.kaon, a.force)}")
        print(f"[dry-run] {decode_cmd(bcin_runlist)}   # if resi exit 2")
        print(f"[dry-run] … repeat until pos exit 0/1/3/4 or max-iter={a.max_iter}")
        sys.exit(runner.EXIT_OK)

    ok_driver = (runner.EXIT_OK, runner.EXIT_FAIL, runner.EXIT_NOT_CONVERGED,
                 runner.EXIT_CHECK, runner.EXIT_SKIP)
    last_rc = runner.EXIT_CHECK

    for i in range(1, a.max_iter + 1):
        runner.step(f"iterate_d5_align pass {i}/{a.max_iter}")

        rc, _ = runner.sh(ref_cmd(a.run, a.kaon), ok=(runner.EXIT_OK,))
        runner.info(f"run_d5 ref exit={rc}")
        p_norm, p_tag = require_scratch_runlists(a.run, s, dry_run=False)
        runner.sh(decode_cmd(p_norm), ok=(runner.EXIT_OK,))
        runner.sh(decode_cmd(p_tag), ok=(runner.EXIT_OK,))

        rc, _ = runner.sh(pos_cmd(a.run, a.kaon, a.tol, a.force), ok=ok_driver)
        last_rc = rc
        action = next_action(rc)
        runner.info(f"run_d5 pos exit={rc} → {action}")

        if action != DO_NEXT:
            stop_with(action)

        if i >= a.max_iter:
            runner.warn(f"max-iter {a.max_iter} reached after NOT_CONVERGED; stop as CHECK")
            sys.exit(runner.EXIT_CHECK)

        runner.sh(decode_cmd(bcin_runlist), ok=(runner.EXIT_OK,))
        rc_resi, _ = runner.sh(resi_cmd(a.run, a.kaon, a.force), ok=ok_driver)
        action_resi = next_action(rc_resi)
        runner.info(f"run_dcprm resi --bcin exit={rc_resi} → {action_resi}")
        if action_resi in (STOP_FAIL, STOP_CHECK, STOP_SKIP):
            stop_with(action_resi)
        if action_resi == DO_NEXT:
            runner.sh(decode_cmd(bcin_runlist), ok=(runner.EXIT_OK,))
        # EXIT_OK or after re-decode: continue outer loop with ref

    runner.warn("loop ended without terminal state; stop as CHECK")
    sys.exit(last_rc if last_rc in (runner.EXIT_OK, runner.EXIT_FAIL,
                                    runner.EXIT_CHECK, runner.EXIT_SKIP)
             else runner.EXIT_CHECK)


if __name__ == "__main__":
    main()
