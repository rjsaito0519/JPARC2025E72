#!/usr/bin/env python3
"""Run one additional hdphc→decode pass and append an hdphc_rounds journal entry.

Usage: extra_hdphc_iter.py <run> [--kaon] [--ftof] [--force] [--dry-run]
                           [--hodo-runlist PATH]

Does not clear prior rounds. Does not auto-adopt (live HDPHC stays the new write).
Invalidates final approval when the dispatch chain was already complete.
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import (  # noqa: E402
    calib_approval,
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
OK_DRIVER = (
    runner.EXIT_OK, runner.EXIT_FAIL, runner.EXIT_NOT_CONVERGED,
    runner.EXIT_CHECK, runner.EXIT_SKIP,
)


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
def next_round_index(run):
    journal = hodo_hdphc_rounds.load_journal(run)
    if not journal:
        return 1
    return max(int(entry["round"]) for entry in journal) + 1


#_____________________________________________________________________________
def hodo_cmd(run, kaon, ftof, force):
    parts = ["python3", str(RUN_HODO), str(run), "hdphc", "--skip-qa"]
    if kaon:
        parts.append("--kaon")
    if ftof:
        parts.append("--ftof")
    if force:
        parts.append("--force")
    return " ".join(parts)


#_____________________________________________________________________________
def main():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int)
    p.add_argument("--kaon", action="store_true")
    p.add_argument("--ftof", action="store_true")
    p.add_argument("--force", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--hodo-runlist", type=Path, default=None)
    a = p.parse_args()
    if a.run < 1:
        runner.fail("run must be >= 1")

    runlist = resolve_hodo_runlist(a.run, a.hodo_runlist, require_exists=not a.dry_run)
    if not MYRUN.is_file() and not a.dry_run:
        runner.fail(f"myrun not found: {MYRUN}")
    if not RUN_HODO.is_file():
        runner.fail(f"run_hodo not found: {RUN_HODO}")

    suffix = "K" if a.kaon else "Pi"
    round_i = next_round_index(a.run)
    cmd = hodo_cmd(a.run, a.kaon, a.ftof, a.force)
    runner.info(
        f"extra_hdphc_iter run={a.run:05d} round={round_i} "
        f"runlist={runlist}"
        + (" --kaon" if a.kaon else "")
        + (" --ftof" if a.ftof else "")
        + (" dry-run" if a.dry_run else "")
    )
    if a.dry_run:
        print(f"[dry-run] {cmd}")
        print(f"[dry-run] decode {runlist} + score + snapshot round_{round_i:02d}")
        sys.exit(runner.EXIT_OK)

    hodo_hdphc_rounds.rounds_dir(a.run).mkdir(parents=True, exist_ok=True)
    rc, _ = runner.sh(cmd, ok=OK_DRIVER)
    runner.info(f"run_hodo hdphc[extra] exit={rc}")
    if rc != runner.EXIT_OK:
        runner.warn(f"stopping after hdphc extra (rc={rc})")
        sys.exit(rc)

    decode_run.run_decode(runlist)
    hodo_path = hodo_btof_monitor.hodo_root(a.run)
    score = hodo_btof_monitor.score_btof(
        a.run, particle_suffix=suffix, ftof=a.ftof, root_path=hodo_path,
    )
    entry = hodo_hdphc_rounds.snapshot_round(
        a.run, round_i, suffix, score, input_path=hodo_path,
        qa_paths=hodo_qa_pdf.source_paths(a.run, "hdphc", a.kaon, a.ftof),
    )
    hodo_hdphc_rounds.write_summary(a.run)
    wtxt = (
        f"{float(entry['width']):.4f}"
        if hodo_hdphc_rounds.width_ok(entry.get("width"))
        else "nan"
    )
    runner.info(
        f"hdphc extra round {round_i}: width={wtxt} "
        f"n_entries={entry.get('n_entries')} ok={entry.get('ok')}"
    )
    if calib_approval.invalidate_approval(
        a.run, f"extra_hdphc_iter round {round_i}",
    ):
        runner.info(
            f"run{a.run:05d}: final approval cleared after extra hdphc "
            f"(awaiting re-approval)"
        )
    sys.exit(runner.EXIT_OK)


if __name__ == "__main__":
    main()
