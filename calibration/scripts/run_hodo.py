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
from common.lib import config, hodo_qa_pdf, qa_session, runner, slack_notify  # noqa: E402


#_____________________________________________________________________________
def hodo_param_path(run, mode, suffix):
    """Per-run Hodo param path written by update_param for this mode."""
    if mode in ("hdprm", "t0"):
        return (
            config.PARAM_DIR / "HDPRM" / config.SUB_DIR
            / f"HodoParam_run{run:05d}_{suffix}"
        )
    if mode == "hdphc":
        return (
            config.PARAM_DIR / "HDPHC" / config.SUB_DIR
            / f"HodoPHCParam_run{run:05d}_{suffix}"
        )
    raise ValueError(f"unsupported hodo mode for param path: {mode!r}")


#_____________________________________________________________________________
def after_qa_merge(run, mode, args, rc, session=None):
    """Merge detector QA PDFs; incomplete book raises the run to EXIT_CHECK."""
    if getattr(args, "skip_qa", False):
        return rc
    if mode not in hodo_qa_pdf.MERGE_MODES:
        return rc
    ok, message, out = hodo_qa_pdf.merge(
        run, mode, kaon=args.kaon, ftof=args.ftof, dry_run=False, session=session,
    )
    if ok:
        runner.info(message)
        manifest = None
        post_pdf = out
        if session is not None and out is not None:
            _ok, _msg, hashes = qa_session.verify_sources(session)
            sealed = qa_session.seal_outputs(session, [out])
            suffix = "K" if args.kaon else "Pi"
            qa_session.attach_params(
                session, [hodo_param_path(run, mode, suffix)], role="written",
            )
            manifest_path = qa_session.write_manifest(
                session, source_hashes=hashes, output_pdfs=sealed,
            )
            manifest = qa_session.load_manifest(manifest_path)
            post_pdf = sealed[0]
            runner.info(f"QA manifest: {manifest_path}")
        posted, post_msg = slack_notify.try_post_hodo_after_merge(
            run, mode, kaon=args.kaon, ftof=args.ftof, pdf_path=post_pdf,
            manifest=manifest,
        )
        if not posted:
            runner.warn(f"Slack QA post failed: {post_msg}")
        elif not post_msg.startswith("slack skip:"):
            runner.info(post_msg)
        return rc
    runner.warn(f"QA PDF merge failed: {message}")
    suffix = "K" if args.kaon else "Pi"
    ok_n, nmsg = slack_notify.try_notify_failure(
        run,
        f"hodo_merge_{mode}_{suffix}",
        f"[CHECK] run{run:05d} hodo {mode} QA PDF merge failed: {message}",
    )
    if not ok_n:
        runner.warn(f"Slack failure notify failed: {nmsg}")
    elif nmsg.startswith("notified"):
        runner.info(nmsg)
    if rc in (runner.EXIT_OK, runner.EXIT_SKIP, runner.EXIT_NOT_CONVERGED):
        return runner.EXIT_CHECK
    return rc


#_____________________________________________________________________________
def fit_and_update(run, suffix, args, mode, bin_suffix):
    """Fit every detector of the set of `mode`, then write that mode's params."""
    dets = config.hodo_detectors(mode, args.ftof)
    bins = runner.require_bins(*[f"{d}_{bin_suffix}" for d in dets])
    root = runner.input_root(run, "Hodo", suffix, args.force)
    session = qa_session.begin(
        run,
        f"hodo_{mode}",
        hodo_qa_pdf.source_paths(run, mode, args.kaon, args.ftof),
        input_path=root,
        meta={
            "family": "hodo",
            "mode": mode,
            "particle": "kaon" if args.kaon else "pion",
            "ftof": bool(args.ftof),
        },
    )
    qa_session.attach_params(
        session, [hodo_param_path(run, mode, suffix)], role="fit_input",
    )
    runner.step(f"run{run:05d}: {bin_suffix} fit for {dets}")
    import os
    from common.lib import fit_overrides
    flat = fit_overrides.write_flat_mip_left(run)
    os.environ["E72_FIT_OVERRIDE"] = str(flat)
    for b in bins:
        runner.sh(f"{b} {root} {suffix}")
    rc = runner.update_param(
        run, suffix, mode, *(["--ftof"] if args.ftof else []), debug=args.debug
    )
    return after_qa_merge(run, mode, args, rc, session=session)


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
    session = qa_session.begin(
        run,
        "hodo_t0",
        hodo_qa_pdf.source_paths(run, "t0", args.kaon, False),
        input_path=root,
        meta={
            "family": "hodo",
            "mode": "t0",
            "particle": "kaon" if args.kaon else "pion",
            "ftof": False,
        },
    )
    qa_session.attach_params(
        session, [hodo_param_path(run, "t0", suffix)], role="fit_input",
    )
    runner.step(f"run{run:05d}: T0_Offset")
    runner.sh(f"{exe} {root} {suffix}")
    rc = runner.update_param(run, suffix, "t0", debug=args.debug)
    return after_qa_merge(run, "t0", args, rc, session=session)


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
    parser.add_argument("--skip-qa", action="store_true",
                        help="skip QA PDF merge/Slack (for intermediate hdphc iterations)")
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
