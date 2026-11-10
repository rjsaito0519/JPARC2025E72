#!/usr/bin/env python3
"""Iterate DC resolution calibration for one run: fit/update → decode → repeat.

Usage: iterate_dc_reso.py <run> [--bcout|--bcin] [--kaon] [--max-iter N]
                         [--separate] [--runlist PATH] [--dry-run] [--force]

Same loop as iterate_dc_resi, but calls run_dcprm.py reso. --separate is passed
through to run_dcprm (a/b Res groups). Does not use --all.
"""
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.append(str(SCRIPT_DIR.parents[1]))
sys.path.insert(0, str(SCRIPT_DIR))

from common.lib import runner  # noqa: E402
import iterate_dc_resi as dc_fit  # noqa: E402


#_____________________________________________________________________________
def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int)
    g = p.add_mutually_exclusive_group()
    g.add_argument("--bcout", action="store_true", help="BcOut (default)")
    g.add_argument("--bcin", action="store_true", help="BcIn")
    p.add_argument("--kaon", action="store_true", help="K param suffix")
    p.add_argument("--max-iter", type=int, default=dc_fit.DEFAULT_MAX_ITER,
                   help=f"max reso passes (default {dc_fit.DEFAULT_MAX_ITER})")
    p.add_argument("--separate", action="store_true",
                   help="passed to run_dcprm reso (Res per a/b)")
    p.add_argument("--runlist", type=Path, default=None,
                   help="myrun yml (default: runmanager/runlist/<sub>/{bcout|bcin}_runXXXXX.yml)")
    p.add_argument("--dry-run", action="store_true", help="print commands only")
    p.add_argument("--force", action="store_true", help="passed through to run_dcprm")
    a = p.parse_args()
    if a.max_iter < 1:
        runner.fail("--max-iter must be >= 1")
    a.bcin = bool(a.bcin)
    return a


#_____________________________________________________________________________
def main():
    a = parse_args()
    runlist = dc_fit.resolve_runlist(a.run, a.bcin, a.runlist)
    rc = dc_fit.run_fit_loop(
        a.run, "reso", a.bcin, a.kaon, a.force, a.max_iter, runlist, a.dry_run,
        separate=a.separate, do_preflight=True, label="iterate_dc_reso",
    )
    sys.exit(rc)


if __name__ == "__main__":
    main()
