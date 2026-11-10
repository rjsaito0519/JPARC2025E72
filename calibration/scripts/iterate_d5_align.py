#!/usr/bin/env python3
"""Iterate D5 BLC1 alignment for one run: create_runlist → ref → decode → pos → …

Usage: iterate_d5_align.py <run> [--kaon] [--max-iter N] [--tol T]
                           [--bcin-runlist PATH] [--mom P] [--dry-run] [--force]
                           [--local]

Order:
  create_runlist --d5 RUN {Pi|K} --mom …
  create_runlist --bcin RUN {Pi|K}_resi   (for BcIn re-decode in the loop)
  then: ref → D5×2 decode → pos → (exit 2) BcIn decode → resi --bcin → …

D5 scratch runlists are those written by run_d5.py ref under
runmanager/runlist/scratch/d5/.
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import calib_runlist, config, decode_run, runner  # noqa: E402

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
def prepare_runlists(run, kaon, mom, *, dry_run):
    """create_runlist --d5 and --bcin (resi) before the align loop."""
    d5_sfx = calib_runlist.d5_suffix(kaon)
    calib_runlist.run_create_runlist(
        "--d5", run, d5_sfx, mom=mom, dry_run=dry_run,
    )
    bcin_sfx = calib_runlist.dc_suffix(kaon, "resi")
    calib_runlist.run_create_runlist(
        "--bcin", run, bcin_sfx, dry_run=dry_run,
    )


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
    p.add_argument("--mom", type=float, default=None,
                   help="D5 central momentum [GeV/c] for create_runlist "
                        "(default: table/UserParam)")
    p.add_argument("--dry-run", action="store_true", help="print commands only")
    p.add_argument("--force", action="store_true",
                   help="passed through to run_d5 pos and run_dcprm resi")
    p.add_argument("--local", action="store_true",
                   help="myrun --local (default: bsub)")
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
    mom = calib_runlist.resolve_d5_mom(a.run, a.mom)
    if not MYRUN.is_file() and not a.dry_run:
        runner.fail(f"myrun not found: {MYRUN}")
    if not RUN_D5.is_file():
        runner.fail(f"run_d5 not found: {RUN_D5}")
    if not RUN_DCPRM.is_file():
        runner.fail(f"run_dcprm not found: {RUN_DCPRM}")

    p_norm, p_tag = d5_scratch_runlists(a.run, s)
    runner.info(
        f"iterate_d5_align run={a.run:05d} suffix={s} "
        f"max_iter={a.max_iter} tol={a.tol} mom={mom}"
        + (" dry-run" if a.dry_run else "")
        + (" --local" if a.local else "")
    )
    runner.info("chain includes create_runlist --d5 and --bcin before loop")

    prepare_runlists(a.run, a.kaon, mom, dry_run=a.dry_run)
    bcin_runlist = resolve_bcin_runlist(
        a.run, a.bcin_runlist, require_exists=not a.dry_run,
    )

    if a.dry_run:
        runner.step("iterate_d5_align dry-run")
        print(f"[dry-run] {ref_cmd(a.run, a.kaon)}")
        decode_run.run_decode(p_norm, dry_run=True, local=a.local)
        decode_run.run_decode(p_tag, dry_run=True, local=a.local)
        print(f"[dry-run] {pos_cmd(a.run, a.kaon, a.tol, a.force)}")
        print("[dry-run] (if pos exit 2)")
        decode_run.run_decode(bcin_runlist, dry_run=True, local=a.local)
        print(f"[dry-run] {resi_cmd(a.run, a.kaon, a.force)}")
        print("[dry-run] (if resi exit 2)")
        decode_run.run_decode(bcin_runlist, dry_run=True, local=a.local)
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
        decode_run.run_decode(p_norm, local=a.local)
        decode_run.run_decode(p_tag, local=a.local)

        rc, _ = runner.sh(pos_cmd(a.run, a.kaon, a.tol, a.force), ok=ok_driver)
        last_rc = rc
        action = next_action(rc)
        runner.info(f"run_d5 pos exit={rc} → {action}")

        if action != DO_NEXT:
            stop_with(action)

        if i >= a.max_iter:
            runner.warn(f"max-iter {a.max_iter} reached after NOT_CONVERGED; stop as CHECK")
            sys.exit(runner.EXIT_CHECK)

        decode_run.run_decode(bcin_runlist, local=a.local)
        rc_resi, _ = runner.sh(resi_cmd(a.run, a.kaon, a.force), ok=ok_driver)
        action_resi = next_action(rc_resi)
        runner.info(f"run_dcprm resi --bcin exit={rc_resi} → {action_resi}")
        if action_resi in (STOP_FAIL, STOP_CHECK, STOP_SKIP):
            stop_with(action_resi)
        if action_resi == DO_NEXT:
            decode_run.run_decode(bcin_runlist, local=a.local)
        # EXIT_OK or after re-decode: continue outer loop with ref

    runner.warn("loop ended without terminal state; stop as CHECK")
    sys.exit(last_rc if last_rc in (runner.EXIT_OK, runner.EXIT_FAIL,
                                    runner.EXIT_CHECK, runner.EXIT_SKIP)
             else runner.EXIT_CHECK)


if __name__ == "__main__":
    main()
