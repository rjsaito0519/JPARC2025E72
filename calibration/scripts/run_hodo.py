#!/usr/bin/env python3
"""Hodoscope calibration: fit, then write the params (update_param.py).

Usage: run_hodo.py <run...> {hdprm,t0,hdphc,htof} [--kaon] [--ftof] [--debug] [--force]
  hdprm  <det>_HDPRM for the detector set of lib/config.py HODO_DETECTORS -> HDPRM
  t0     T0_Offset -> HDPRM (BH2 T0 offsets)
  hdphc  <det>_PHC for the detector set -> HDPHC
  htof   HTOF_Calib (TDC, ADC, PHC in one pass) -> HDPRM and HDPHC
Input: DATA_DIR/runXXXXX_Hodo.root (create_runlist.py --hodo); htof also DATA_DIR/runXXXXX_TPCHelixHTOF.root.
Runs are processed one by one, except htof, which merges the runs and labels the result with the
first run. For htof, --force also allows a TPCHelixHTOF older than that run's Hodo (logged as a warning).
Exit code: worst over runs (and htofprm/htofphc); 0 OK, 4 skip/debug, 3 check, 2 not converged, 1 fail.
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import config, runner  # noqa: E402


#_____________________________________________________________________________
def fit_and_update(run, suffix, args, mode, bin_suffix):
    """Fit every detector of the set of `mode`, then write that mode's params."""
    dets = config.hodo_detectors(mode, args.ftof)
    bins = runner.require_bins(*[f"{d}_{bin_suffix}" for d in dets])
    root = runner.input_root(run, "Hodo", suffix, args.force)
    runner.step(f"run{run:05d}: {bin_suffix} fit for {dets}")
    for b in bins:
        runner.sh(f"{b} {root} {suffix}")
    return runner.update_param(
        run, suffix, mode, *(["--ftof"] if args.ftof else []), debug=args.debug
    )


#_____________________________________________________________________________
def do_hdprm(run, suffix, args):
    return fit_and_update(run, suffix, args, "hdprm", "HDPRM")


#_____________________________________________________________________________
def do_hdphc(run, suffix, args):
    return fit_and_update(run, suffix, args, "hdphc", "PHC")


#_____________________________________________________________________________
def do_t0(run, suffix, args):
    (exe,) = runner.require_bins("T0_Offset")
    root = runner.input_root(run, "Hodo", suffix, args.force)
    runner.step(f"run{run:05d}: T0_Offset")
    runner.sh(f"{exe} {root} {suffix}")
    return runner.update_param(run, suffix, "t0", debug=args.debug)


#_____________________________________________________________________________
def htof_dst(run, force=False):
    """DstTPCHelixHTOF output of the run (dE/dx-tagged helix x HTOF match); it must follow the Hodo decode."""
    dst = config.DATA_DIR / f"run{run:05d}_TPCHelixHTOF.root"
    if not dst.exists():
        runner.fail(f"{dst} not found (dst_create_runlist.py {run} --tpchelixhtof, then run the job)")
    hodo = config.DATA_DIR / f"run{run:05d}_Hodo.root"
    if hodo.exists() and dst.resolve().stat().st_mtime < hodo.resolve().stat().st_mtime:
        msg = (f"{dst.name} is older than {hodo.name}; regenerate it after the Hodo re-decode")
        if not force:
            runner.fail(f"{msg} (use --force to run anyway)")
        runner.warn(f"--force: allowing stale DST ({msg})")
    runner.info(f"input: {dst}")
    return dst


#_____________________________________________________________________________
def do_htof(runs, suffix, args):
    """HTOF_Calib on all runs merged; params labeled with the first run."""
    if args.kaon:
        runner.fail("htof supports Pi only (TPC dE/dx pion tag); --kaon is not implemented")
    if args.ftof:
        runner.warn("--ftof does not apply to htof; ignored")
    (exe,) = runner.require_bins("HTOF_Calib")
    pairs = [(runner.input_root(r, "Hodo", suffix, args.force), htof_dst(r, args.force)) for r in runs]
    if len(runs) > 1:
        runner.info(f"merging {len(runs)} runs; results/params labeled with run {runs[0]:05d}")
    runner.step("HTOF_Calib")
    runner.sh(f"{exe} " + " ".join(f"{h} {d}" for h, d in pairs))
    rc_prm = runner.update_param(runs[0], suffix, "htofprm", debug=args.debug)
    rc_phc = runner.update_param(runs[0], suffix, "htofphc", debug=args.debug)
    return runner.worst([rc_prm, rc_phc])


PER_RUN = {"hdprm": do_hdprm, "t0": do_t0, "hdphc": do_hdphc}


#_____________________________________________________________________________
def main():
    parser = runner.base_parser(__doc__, ["hdprm", "t0", "hdphc", "htof"])
    parser.add_argument("--ftof", action="store_true",
                        help="hdprm/hdphc: the FTOF detector set instead of the default one "
                             f"(lib/config.py HODO_DETECTORS: {config.HODO_DETECTORS}); t0 is unchanged")
    args = parser.parse_args()
    suffix = runner.suffix_of(args)
    runner.info(f"mode {args.mode}, suffix {suffix}" + (", FTOF set" if args.ftof else "")
                + (", --debug (no param written)" if args.debug else ""))
    codes = []
    if args.mode == "htof":
        codes.append(do_htof(args.runs, suffix, args))
    else:
        for run in args.runs:
            codes.append(PER_RUN[args.mode](run, suffix, args))
    final = runner.worst(codes) if codes else runner.EXIT_OK
    runner.info(f"done: {args.mode} for {' '.join(map(str, args.runs))} exit={final}")
    sys.exit(final)


if __name__ == "__main__":
    main()
