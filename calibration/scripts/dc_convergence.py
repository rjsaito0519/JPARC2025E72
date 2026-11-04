"""Pass history and the per-group decision of iterated calibrations (DC resi/reso, D5 pos).

One text file per run, chamber (or ALL), mode and particle in OUTPUT_DIR/calib/runXXXXX/:
  dc_<BLC1|BLC2|ALL>_<resi|reso>_<K|Pi>.txt (update_param.py), d5_pos_<K|Pi>.txt (run_d5.py)
one line per pass and group: time, group, measured value, state[, details].

States (lib/runner.py exit codes):
  converged  within the criterion; nothing written (0)
  settled    outside the criterion but not getting better over SETTLE_PASSES passes; nothing
             written, a person looks (3)
  written    params written; decode again and run once more (2)
"""
import datetime
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from lib import config, runner  # noqa: E402

CONVERGED, SETTLED, WRITTEN = "converged", "settled", "written"
SETTLE_PASSES = 3
EXIT_CODE = {CONVERGED: runner.EXIT_OK, WRITTEN: runner.EXIT_NOT_CONVERGED, SETTLED: runner.EXIT_CHECK}


def history_path(run, group, mode, suffix):
    chamber = "ALL" if group == "ALL" else group[:4]
    return config.OUTPUT_DIR / "calib" / f"run{run:05d}" / f"dc_{chamber}_{mode}_{suffix}.txt"


def previous_values(path, group):
    """Values of `group` measured since its last converged pass, oldest first."""
    values = []
    if not path.exists():
        return values
    for line in path.read_text().splitlines():
        f = line.split()
        if len(f) >= 4 and f[1] == group:
            values = [] if f[3] == CONVERGED else values + [float(f[2])]
    return values


def decide(value, within, previous, not_better):
    """State of this pass; not_better(a, b) says pass b did not improve on pass a."""
    if within:
        return CONVERGED
    recent = previous[-(SETTLE_PASSES - 1):] + [value]
    if len(recent) == SETTLE_PASSES and all(not_better(a, b) for a, b in zip(recent, recent[1:])):
        return SETTLED
    return WRITTEN


def record(path, group, value, state, details=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as f:
        f.write(f"{datetime.datetime.now():%Y-%m-%dT%H:%M:%S}\t{group}\t{value:.6g}\t{state}"
                + (f"\t{details}" if details else "") + "\n")


def judge_at(path, group, value, within, not_better, details=""):
    """decide() with the group's history in path, record the pass, and return the state."""
    state = decide(value, within, previous_values(path, group), not_better)
    record(path, group, value, state, details)
    return state


def judge(run, group, mode, suffix, value, within, not_better):
    """judge_at() on the DC history file of this run, group, mode and particle."""
    return judge_at(history_path(run, group, mode, suffix), group, value, within, not_better)


def exit_code(states):
    """Worst state of the groups: written > settled > converged."""
    if WRITTEN in states:
        return EXIT_CODE[WRITTEN]
    if SETTLED in states:
        return EXIT_CODE[SETTLED]
    return EXIT_CODE[CONVERGED]
