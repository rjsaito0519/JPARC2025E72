#!/usr/bin/env python3
"""Restore an hdphc round snapshot to the live HDPHC param (manual override).

Usage: adopt_hdphc_round.py <run> <round> [--kaon] [--decode] [--hodo-runlist PATH]

Round directories come from iterate_hodo_chain (OUTPUT_DIR/calib/runXXXXX/hdphc_rounds/).
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import decode_run, hodo_hdphc_rounds, runner  # noqa: E402
from common.lib import config  # noqa: E402

RUNLIST_DIR = config.ANALYZER_DIR / "runmanager" / "runlist"


#_____________________________________________________________________________
def main():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int)
    p.add_argument("round", type=int, help="1-based round index from journal")
    p.add_argument("--kaon", action="store_true")
    p.add_argument("--decode", action="store_true",
                   help="re-decode Hodo after restoring the param")
    p.add_argument("--hodo-runlist", type=Path, default=None)
    a = p.parse_args()
    if a.run < 1 or a.round < 1:
        runner.fail("run and round must be >= 1")
    suffix = "K" if a.kaon else "Pi"
    journal = hodo_hdphc_rounds.load_journal(a.run)
    rounds = {int(e["round"]) for e in journal}
    if a.round not in rounds:
        runner.fail(f"round {a.round} not in journal {sorted(rounds)}")
    dest = hodo_hdphc_rounds.restore_round(a.run, a.round, suffix)
    hodo_hdphc_rounds.write_summary(a.run, adopted=a.round)
    runner.info(f"restored {dest.name} from round_{a.round:02d}")
    if a.decode:
        if a.hodo_runlist is not None:
            runlist = Path(a.hodo_runlist).expanduser()
            if not runlist.is_absolute():
                runlist = (Path.cwd() / runlist).resolve()
        else:
            runlist = RUNLIST_DIR / config.SUB_DIR / f"hodo_run{a.run:05d}.yml"
        if not runlist.is_file():
            runner.fail(f"hodo runlist not found: {runlist}")
        decode_run.run_decode(runlist)
    sys.exit(runner.EXIT_OK)


if __name__ == "__main__":
    main()
