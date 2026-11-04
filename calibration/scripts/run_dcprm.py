#!/usr/bin/env python3
"""BC (BLC) drift chamber calibration: fit, then write the params.

Usage: run_dcprm.py <run...> {t0,drift,resi,reso} [--bcout|--bcin|--all] [--separate] [--kaon] [--debug] [--force]
  t0     BLC_TDC -> DCTDC (update_param.py dctdc)
  drift  BLC_DRFT -> DCDRFT (BLC_DRFT writes param/DCDRFT/<sub>/DCDriftParam_runXXXXX_<S>.root itself)
  resi   BLC_residual -> DCGEO Ofs (update_param.py residual)
  reso   BLC_pull: exclusive pull width p -> common Res_new = p * mean(Res in group) in DCGEO
         (default one Res for a+b; --separate one each for a and b; --all one for BLC1+BLC2).
         Iterate: decode BcIn/BcOut -> reso -> decode.
Chamber: BcOut by default, --bcin for BcIn. Input: DATA_DIR/runXXXXX_<BcIn|BcOut>.root (create_runlist.py).
Runs are processed one by one. resi and reso judge convergence (calibration/scripts/calib_convergence.py);
the exit code is the worst of the runs: 3 no longer improving (look at it) over 2 written (decode again
and repeat) over 0 converged; 1 failed.
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from lib import config, runner  # noqa: E402


def chamber(args):
    return "BcIn" if args.bcin else "BcOut"


def chamber_flag(args):
    return "--bcin" if args.bcin else "--bcout"


def fit_and_update(run, suffix, args, exe_name, ptype):
    """One chamber: run exe_name on its decode, then update_param.py ptype for that chamber."""
    (exe,) = runner.require_bins(exe_name)
    root = runner.input_root(run, chamber(args), suffix, args.force)
    runner.step(f"run{run:05d}: {exe_name} ({chamber(args)})")
    runner.sh(f"{exe} {root} {suffix}")
    return runner.update_param(run, suffix, ptype, chamber_flag(args), debug=args.debug)


def do_t0(run, suffix, args):
    return fit_and_update(run, suffix, args, "BLC_TDC", "dctdc")


def do_resi(run, suffix, args):
    return fit_and_update(run, suffix, args, "BLC_residual", "residual")


def do_drift(run, suffix, args):
    """BLC_DRFT writes the DCDRFT ROOT param itself; back it up first."""
    (exe,) = runner.require_bins("BLC_DRFT")
    root = runner.input_root(run, chamber(args), suffix, args.force)
    if not args.debug:
        runner.backup_file(config.PARAM_DIR / "DCDRFT" / config.SUB_DIR / f"DCDriftParam_run{run:05d}_{suffix}.root")
    runner.step(f"run{run:05d}: BLC_DRFT ({chamber(args)})")
    runner.sh(f"{exe} {root} {suffix}" + (" --debug" if args.debug else ""))
    return runner.EXIT_OK


def do_reso(run, suffix, args):
    (exe,) = runner.require_bins("BLC_pull")
    dets = ["BcIn", "BcOut"] if args.all else [chamber(args)]
    roots = [runner.input_root(run, d, suffix, args.force) for d in dets]
    for d, root in zip(dets, roots):
        runner.step(f"run{run:05d}: BLC_pull ({d})")
        runner.sh(f"{exe} {root} {suffix}")
    opts = ["--all" if args.all else chamber_flag(args)] + (["--separate"] if args.separate else [])
    return runner.update_param(run, suffix, "resolution", *opts, debug=args.debug)


MODES = {"t0": do_t0, "drift": do_drift, "resi": do_resi, "reso": do_reso}


def main():
    parser = runner.base_parser(__doc__, list(MODES))
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--bcout", action="store_true", help="BcOut (default)")
    group.add_argument("--bcin", action="store_true", help="BcIn")
    group.add_argument("--all", action="store_true",
                       help="reso only: one common Res for BLC1+BLC2 (p = mean of both); with --separate: 4 groups")
    parser.add_argument("--separate", action="store_true",
                        help="reso only: common Res per a/b (default one Res for a+b)")
    args = parser.parse_args()
    if (args.all or args.separate) and args.mode != "reso":
        runner.fail("--all and --separate are for mode reso only")
    suffix = runner.suffix_of(args)
    runner.info(f"mode {args.mode}, suffix {suffix}, " + ("BcIn+BcOut" if args.all else chamber(args))
                + (", --debug (no param written)" if args.debug else ""))
    codes = {run: MODES[args.mode](run, suffix, args) for run in args.runs}
    runner.info(f"done: {args.mode}; exit code per run {codes}")
    sys.exit(runner.worst(codes.values()))


if __name__ == "__main__":
    main()
