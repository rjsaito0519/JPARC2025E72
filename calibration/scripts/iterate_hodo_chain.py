#!/usr/bin/env python3
"""Hodo calibration chain: create_runlist → decode → fit → recreate runlist …

Usage: iterate_hodo_chain.py <run> [--kaon] [--ftof] [--force] [--dry-run]
                             [--skip-t0] [--skip-hdphc] [--hdphc-iter N]
                             [--adopt auto|K] [--skip-hdphc-adopt] [--local]
                             [--assume-window] [--skip-window-bootstrap]
                             [--only hdprm|t0|hdphc|hdphc-adopt]
                             [--hdphc-round I]

Full order (wiki Hodoscope-Calibration / TDC-Window):
  0) if TDC window unset: create_runlist --hodo RUN 0 → decode → EXIT_CHECK
  1) create_runlist {Pi|K}_hdprm → decode → run_hodo hdprm
  2) create_runlist {Pi|K}_t0    → decode → run_hodo t0
  3) create_runlist {Pi|K}_hdphc → hdphc×N (fit→decode) → adopt → decode → qa-only

Does not run htof. Does not roll back hdprm/t0.
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import (  # noqa: E402
    calib_runlist,
    config,
    decode_run,
    hodo_btof_monitor,
    hodo_hdphc_rounds,
    hodo_qa_pdf,
    runner,
)

SCRIPT_DIR = Path(__file__).resolve().parent
RUN_HODO = SCRIPT_DIR / "run_hodo.py"
MYRUN = config.ANALYZER_DIR / "runmanager" / "myrun.py"
OK_DRIVER = (runner.EXIT_OK, runner.EXIT_FAIL, runner.EXIT_NOT_CONVERGED,
             runner.EXIT_CHECK, runner.EXIT_SKIP)


#_____________________________________________________________________________
def hodo_cmd(
    run, mode, kaon, ftof, force, skip_qa=False, qa_only=False, adopted_round=None,
):
    parts = ["python3", str(RUN_HODO), str(run), mode]
    if kaon:
        parts.append("--kaon")
    if ftof:
        parts.append("--ftof")
    if force:
        parts.append("--force")
    if skip_qa:
        parts.append("--skip-qa")
    if qa_only:
        parts.append("--qa-only")
        if adopted_round is not None:
            parts.append(f"--adopted-round {int(adopted_round)}")
    return " ".join(parts)


#_____________________________________________________________________________
def note_hodo_substage(run, stage, *, rc=None, status="done"):
    from common.lib import calib_progress
    calib_progress.note_substage(
        run, "hodo", stage, rc=rc, status=status,
    )


#_____________________________________________________________________________
def run_hodo_ok(run, cmd, label):
    rc, _ = runner.sh(cmd, ok=OK_DRIVER)
    runner.info(f"run_hodo {label} exit={rc}")
    if rc != runner.EXIT_OK:
        note_hodo_substage(run, label, rc=rc, status="stopped")
        runner.warn(f"stopping after {label} (rc={rc})")
        sys.exit(rc)
    note_hodo_substage(run, label, rc=rc, status="done")


#_____________________________________________________________________________
def stage_prepare_and_decode(run, kaon, stage, *, local, dry_run):
    """create_runlist for hodo stage, then myrun decode. Returns runlist path."""
    suffix = calib_runlist.hodo_suffix(kaon, stage)
    calib_runlist.run_create_runlist("--hodo", run, suffix, dry_run=dry_run)
    runlist = (
        config.ANALYZER_DIR / "runmanager" / "runlist" / config.SUB_DIR
        / f"hodo_run{run:05d}.yml"
    )
    if dry_run:
        decode_run.run_decode(runlist, dry_run=True, local=local)
        return runlist
    if not runlist.is_file():
        runner.fail(f"hodo runlist missing after create_runlist: {runlist}")
    decode_run.run_decode(runlist, local=local)
    return runlist


#_____________________________________________________________________________
def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int)
    p.add_argument("--kaon", action="store_true", help="K param suffix")
    p.add_argument("--ftof", action="store_true",
                   help="passed to run_hodo hdprm/hdphc (FTOF detector set)")
    p.add_argument("--force", action="store_true", help="passed through to run_hodo")
    p.add_argument("--dry-run", action="store_true", help="print commands only")
    p.add_argument("--skip-t0", action="store_true", help="omit t0 stage")
    p.add_argument("--skip-hdphc", action="store_true",
                   help="omit hdphc loop and adopt")
    p.add_argument("--hdphc-iter", type=int, default=3,
                   help="number of hdphc→decode passes (default 3)")
    p.add_argument("--adopt", default="auto",
                   help="auto (min CBTOF width) or 1-based round index to keep")
    p.add_argument("--skip-hdphc-adopt", action="store_true",
                   help="leave the last hdphc write as-is (no restore/final decode)")
    p.add_argument("--local", action="store_true",
                   help="myrun --local (default: bsub)")
    p.add_argument("--assume-window", action="store_true",
                   help="skip TDC window gate (windows already set)")
    p.add_argument("--skip-window-bootstrap", action="store_true",
                   help="do not run create_runlist 0 / decode for window setup")
    p.add_argument(
        "--only", default=None,
        choices=("hdprm", "t0", "hdphc", "hdphc-adopt"),
        help="run a single stage (for fine-grained dispatch step_gate)",
    )
    p.add_argument(
        "--hdphc-round", type=int, default=None,
        help="with --only hdphc: 1-based round index (required)",
    )
    return p.parse_args()


#_____________________________________________________________________________
def ensure_window_or_exit(a):
    ready, upath, detail = calib_runlist.tdc_windows_ready(a.run)
    if a.assume_window:
        ready = True
        detail = "assumed via --assume-window"
    runner.info(f"TDC window: ready={ready} ({detail}); UserParam={upath}")
    if ready:
        return
    if a.dry_run:
        runner.step("iterate_hodo_chain dry-run (window unset)")
        if not a.skip_window_bootstrap:
            calib_runlist.run_create_runlist(
                "--hodo", a.run, "0", dry_run=True,
            )
            print("[dry-run] myrun decode (suffix 0)")
        calib_runlist.announce_tdc_window(a.run)
        print("[dry-run] would stop with EXIT_CHECK until windows are set")
        sys.exit(runner.EXIT_OK)
    if not a.skip_window_bootstrap:
        stage_prepare_and_decode(
            a.run, a.kaon, "0", local=a.local, dry_run=False,
        )
    calib_runlist.announce_tdc_window(a.run)
    sys.exit(runner.EXIT_CHECK)


#_____________________________________________________________________________
def run_hdprm_stage(a):
    stage_prepare_and_decode(a.run, a.kaon, "hdprm", local=a.local, dry_run=False)
    run_hodo_ok(
        a.run, hodo_cmd(a.run, "hdprm", a.kaon, a.ftof, a.force), "hdprm",
    )


#_____________________________________________________________________________
def run_t0_stage(a):
    stage_prepare_and_decode(a.run, a.kaon, "t0", local=a.local, dry_run=False)
    run_hodo_ok(
        a.run, hodo_cmd(a.run, "t0", a.kaon, a.ftof, a.force), "t0",
    )


#_____________________________________________________________________________
def run_hdphc_one_round(a, round_i, *, clear_journal):
    suffix = "K" if a.kaon else "Pi"
    runlist = stage_prepare_and_decode(
        a.run, a.kaon, "hdphc", local=a.local, dry_run=False,
    )
    if clear_journal:
        hodo_hdphc_rounds.clear_rounds(a.run)
    hodo_path = hodo_btof_monitor.hodo_root(a.run)
    runner.step(f"iterate_hodo_chain hdphc {round_i}")
    run_hodo_ok(
        a.run,
        hodo_cmd(a.run, "hdphc", a.kaon, a.ftof, a.force, skip_qa=True),
        f"hdphc[{round_i}]",
    )
    decode_run.run_decode(runlist, local=a.local)
    score = hodo_btof_monitor.score_btof(
        a.run, particle_suffix=suffix, ftof=a.ftof, root_path=hodo_path,
    )
    entry = hodo_hdphc_rounds.snapshot_round(
        a.run, round_i, suffix, score, input_path=hodo_path,
        qa_paths=hodo_qa_pdf.source_paths(a.run, "hdphc", a.kaon, a.ftof),
    )
    wtxt = (
        f"{float(entry['width']):.4f}"
        if hodo_hdphc_rounds.width_ok(entry.get("width"))
        else "nan"
    )
    runner.info(
        f"hdphc round {round_i}: width={wtxt} n_entries={entry.get('n_entries')} "
        f"ok={entry.get('ok')}"
    )


#_____________________________________________________________________________
def run_hdphc_adopt_stage(a):
    suffix = "K" if a.kaon else "Pi"
    runlist = (
        config.ANALYZER_DIR / "runmanager" / "runlist" / config.SUB_DIR
        / f"hodo_run{a.run:05d}.yml"
    )
    journal = hodo_hdphc_rounds.load_journal(a.run)
    try:
        chosen = hodo_hdphc_rounds.pick_round(journal, a.adopt)
    except ValueError as exc:
        runner.fail(str(exc))
    runner.info(f"adopting hdphc round {chosen}")
    hodo_hdphc_rounds.restore_round(a.run, chosen, suffix)
    if not runlist.is_file():
        runlist = stage_prepare_and_decode(
            a.run, a.kaon, "hdphc", local=a.local, dry_run=False,
        )
    else:
        decode_run.run_decode(runlist, local=a.local)
    hodo_hdphc_rounds.write_summary(a.run, adopted=chosen)
    run_hodo_ok(
        a.run,
        hodo_cmd(
            a.run, "hdphc", a.kaon, a.ftof, a.force,
            qa_only=True, adopted_round=chosen,
        ),
        f"hdphc[qa-only adopt={chosen}]",
    )
    runner.info(f"done: hdphc adopted round {chosen}")


#_____________________________________________________________________________
def main():
    a = parse_args()
    if a.hdphc_iter < 1 and not a.skip_hdphc and a.only is None:
        runner.fail("--hdphc-iter must be >= 1")
    if a.adopt != "auto":
        try:
            k = int(a.adopt)
        except ValueError:
            runner.fail("--adopt must be 'auto' or a positive integer")
        if k < 1:
            runner.fail("--adopt round must be >= 1")
    if a.only == "hdphc":
        if a.hdphc_round is None or a.hdphc_round < 1:
            runner.fail("--only hdphc requires --hdphc-round >= 1")

    if not MYRUN.is_file() and not a.dry_run:
        runner.fail(f"myrun not found: {MYRUN}")
    if not RUN_HODO.is_file():
        runner.fail(f"run_hodo not found: {RUN_HODO}")

    suffix = "K" if a.kaon else "Pi"
    n_hdphc = 0 if a.skip_hdphc else a.hdphc_iter
    only_s = f" only={a.only}" if a.only else ""
    runner.info(
        f"iterate_hodo_chain run={a.run:05d} hdphc_iter={n_hdphc} adopt={a.adopt}"
        + only_s
        + (" --kaon" if a.kaon else "")
        + (" --ftof" if a.ftof else "")
        + (" --skip-t0" if a.skip_t0 else "")
        + (" --skip-hdphc" if a.skip_hdphc else "")
        + (" --local" if a.local else "")
        + (" --assume-window" if a.assume_window else "")
        + (" dry-run" if a.dry_run else "")
    )
    runner.info(
        "chain includes create_runlist per stage "
        "(UserParam / conf / runlist → decode → fit)"
    )

    ensure_window_or_exit(a)

    if a.only is not None:
        if a.dry_run:
            runner.info(f"[dry-run] --only {a.only}")
            sys.exit(runner.EXIT_OK)
        if a.only == "hdprm":
            run_hdprm_stage(a)
        elif a.only == "t0":
            run_t0_stage(a)
        elif a.only == "hdphc":
            run_hdphc_one_round(
                a, int(a.hdphc_round),
                clear_journal=(int(a.hdphc_round) == 1),
            )
        elif a.only == "hdphc-adopt":
            run_hdphc_adopt_stage(a)
        sys.exit(runner.EXIT_OK)

    if a.dry_run:
        runner.step("iterate_hodo_chain dry-run")
        stage_prepare_and_decode(a.run, a.kaon, "hdprm", local=a.local, dry_run=True)
        print(f"[dry-run] {hodo_cmd(a.run, 'hdprm', a.kaon, a.ftof, a.force)}")
        if not a.skip_t0:
            stage_prepare_and_decode(a.run, a.kaon, "t0", local=a.local, dry_run=True)
            print(f"[dry-run] {hodo_cmd(a.run, 't0', a.kaon, a.ftof, a.force)}")
        if n_hdphc:
            stage_prepare_and_decode(
                a.run, a.kaon, "hdphc", local=a.local, dry_run=True,
            )
            for i in range(1, n_hdphc + 1):
                print(
                    f"[dry-run] ({i}/{n_hdphc}) "
                    f"{hodo_cmd(a.run, 'hdphc', a.kaon, a.ftof, a.force, skip_qa=True)}"
                )
                print("[dry-run] myrun decode + score_btof + snapshot")
            if not a.skip_hdphc_adopt:
                print(
                    f"[dry-run] adopt {a.adopt} → restore HDPHC + decode + "
                    f"run_hodo hdphc --qa-only"
                )
        sys.exit(runner.EXIT_OK)

    run_hdprm_stage(a)
    if not a.skip_t0:
        run_t0_stage(a)
    if a.skip_hdphc:
        runner.info("done: skipped hdphc")
        sys.exit(runner.EXIT_OK)

    for i in range(1, n_hdphc + 1):
        run_hdphc_one_round(a, i, clear_journal=(i == 1))

    if a.skip_hdphc_adopt:
        hodo_hdphc_rounds.write_summary(a.run)
        runner.info("done: hdphc rounds recorded; adopt skipped")
        sys.exit(runner.EXIT_OK)

    run_hdphc_adopt_stage(a)
    sys.exit(runner.EXIT_OK)


if __name__ == "__main__":
    main()
