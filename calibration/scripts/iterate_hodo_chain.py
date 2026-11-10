#!/usr/bin/env python3
"""Hodo calibration chain for one run: hdprm → decode → t0 → decode → hdphc×N → adopt.

Usage: iterate_hodo_chain.py <run> [--kaon] [--ftof] [--force] [--dry-run]
                             [--hodo-runlist PATH] [--skip-t0] [--skip-hdphc]
                             [--hdphc-iter N] [--adopt auto|K] [--skip-hdphc-adopt]

hdprm/t0 run once. hdphc repeats N times (default 3): each pass fits, writes HDPHC,
re-decodes, scores CBTOF robust width, and snapshots the param plus fit-time QA PDFs.
After the loop, the round with the smallest width is restored (or --adopt K), then a
final decode runs and `run_hodo hdphc --qa-only` regenerates post-adopt QA (no param write)
with a manifest tied to the adopted round.

Continue only on run_hodo EXIT_OK (0) for each step.

Does not run htof. Does not roll back hdprm/t0.
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import (  # noqa: E402
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
def run_hodo_ok(cmd, label):
    rc, _ = runner.sh(cmd, ok=OK_DRIVER)
    runner.info(f"run_hodo {label} exit={rc}")
    if rc != runner.EXIT_OK:
        runner.warn(f"stopping after {label} (rc={rc})")
        sys.exit(rc)


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
                   help="omit hdphc loop and adopt")
    p.add_argument("--hdphc-iter", type=int, default=3,
                   help="number of hdphc→decode passes (default 3)")
    p.add_argument("--adopt", default="auto",
                   help="auto (min CBTOF width) or 1-based round index to keep")
    p.add_argument("--skip-hdphc-adopt", action="store_true",
                   help="leave the last hdphc write as-is (no restore/final decode)")
    return p.parse_args()


#_____________________________________________________________________________
def main():
    a = parse_args()
    if a.hdphc_iter < 1 and not a.skip_hdphc:
        runner.fail("--hdphc-iter must be >= 1")
    if a.adopt != "auto":
        try:
            k = int(a.adopt)
        except ValueError:
            runner.fail("--adopt must be 'auto' or a positive integer")
        if k < 1:
            runner.fail("--adopt round must be >= 1")

    runlist = resolve_hodo_runlist(a.run, a.hodo_runlist, require_exists=not a.dry_run)
    if not MYRUN.is_file() and not a.dry_run:
        runner.fail(f"myrun not found: {MYRUN}")
    if not RUN_HODO.is_file():
        runner.fail(f"run_hodo not found: {RUN_HODO}")

    suffix = "K" if a.kaon else "Pi"
    n_hdphc = 0 if a.skip_hdphc else a.hdphc_iter
    runner.info(
        f"iterate_hodo_chain run={a.run:05d} hdphc_iter={n_hdphc} adopt={a.adopt} "
        f"runlist={runlist}"
        + (" --kaon" if a.kaon else "")
        + (" --ftof" if a.ftof else "")
        + (" --skip-t0" if a.skip_t0 else "")
        + (" --skip-hdphc" if a.skip_hdphc else "")
        + (" --skip-hdphc-adopt" if a.skip_hdphc_adopt else "")
        + (" dry-run" if a.dry_run else "")
    )

    if a.dry_run:
        runner.step("iterate_hodo_chain dry-run")
        print(f"[dry-run] {hodo_cmd(a.run, 'hdprm', a.kaon, a.ftof, a.force)}")
        print(f"[dry-run] decode {runlist}")
        if not a.skip_t0:
            print(f"[dry-run] {hodo_cmd(a.run, 't0', a.kaon, a.ftof, a.force)}")
            print(f"[dry-run] decode {runlist}")
        for i in range(1, n_hdphc + 1):
            print(
                f"[dry-run] ({i}/{n_hdphc}) "
                f"{hodo_cmd(a.run, 'hdphc', a.kaon, a.ftof, a.force, skip_qa=True)}"
            )
            print(f"[dry-run] decode {runlist} + score_btof + snapshot round_{i:02d}")
        if n_hdphc and not a.skip_hdphc_adopt:
            print(
                f"[dry-run] adopt {a.adopt} -> restore HDPHC + decode + "
                f"run_hodo hdphc --qa-only"
            )
        sys.exit(runner.EXIT_OK)

    # hdprm
    run_hodo_ok(hodo_cmd(a.run, "hdprm", a.kaon, a.ftof, a.force), "hdprm")
    decode_run.run_decode(runlist)

    if not a.skip_t0:
        run_hodo_ok(hodo_cmd(a.run, "t0", a.kaon, a.ftof, a.force), "t0")
        decode_run.run_decode(runlist)

    if a.skip_hdphc:
        runner.info("done: skipped hdphc")
        sys.exit(runner.EXIT_OK)

    hodo_hdphc_rounds.clear_rounds(a.run)
    hodo_path = hodo_btof_monitor.hodo_root(a.run)
    for i in range(1, n_hdphc + 1):
        runner.step(f"iterate_hodo_chain hdphc {i}/{n_hdphc}")
        run_hodo_ok(
            hodo_cmd(a.run, "hdphc", a.kaon, a.ftof, a.force, skip_qa=True),
            f"hdphc[{i}]",
        )
        decode_run.run_decode(runlist)
        score = hodo_btof_monitor.score_btof(
            a.run, particle_suffix=suffix, ftof=a.ftof, root_path=hodo_path,
        )
        entry = hodo_hdphc_rounds.snapshot_round(
            a.run, i, suffix, score, input_path=hodo_path,
            qa_paths=hodo_qa_pdf.source_paths(a.run, "hdphc", a.kaon, a.ftof),
        )
        wtxt = (
            f"{float(entry['width']):.4f}"
            if hodo_hdphc_rounds.width_ok(entry.get("width"))
            else "nan"
        )
        runner.info(
            f"hdphc round {i}: width={wtxt} n_entries={entry.get('n_entries')} "
            f"ok={entry.get('ok')}"
        )

    if a.skip_hdphc_adopt:
        hodo_hdphc_rounds.write_summary(a.run)
        runner.info("done: hdphc rounds recorded; adopt skipped")
        sys.exit(runner.EXIT_OK)

    journal = hodo_hdphc_rounds.load_journal(a.run)
    try:
        chosen = hodo_hdphc_rounds.pick_round(journal, a.adopt)
    except ValueError as exc:
        runner.fail(str(exc))
    runner.info(f"adopting hdphc round {chosen}")
    hodo_hdphc_rounds.restore_round(a.run, chosen, suffix)
    decode_run.run_decode(runlist)
    hodo_hdphc_rounds.write_summary(a.run, adopted=chosen)

    # Post-adopt QA: re-fit plots on the restored param / final decode; no param write.
    run_hodo_ok(
        hodo_cmd(
            a.run, "hdphc", a.kaon, a.ftof, a.force,
            qa_only=True, adopted_round=chosen,
        ),
        f"hdphc[qa-only adopt={chosen}]",
    )

    runner.info(f"done: hdphc adopted round {chosen}")
    sys.exit(runner.EXIT_OK)


if __name__ == "__main__":
    main()
