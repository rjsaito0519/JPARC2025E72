"""Pass history and the per-group decision of the iterated DC calibrations (resi, reso).

One text file per run, chamber (or ALL), mode and particle:
  OUTPUT_DIR/calib/runXXXXX/dc_<BLC1|BLC2|ALL>_<resi|reso>_<K|Pi>.txt
one line per pass and group: time, group, measured value, state.

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
        if len(f) == 4 and f[1] == group:
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


def record(path, group, value, state):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as f:
        f.write(f"{datetime.datetime.now():%Y-%m-%dT%H:%M:%S}\t{group}\t{value:.6g}\t{state}\n")


def judge(run, group, mode, suffix, value, within, not_better):
    """decide() with this group's history, record the pass, and return the state."""
    path = history_path(run, group, mode, suffix)
    state = decide(value, within, previous_values(path, group), not_better)
    record(path, group, value, state)
    return state


def exit_code(states):
    """Worst state of the groups: written > settled > converged."""
    if WRITTEN in states:
        return EXIT_CODE[WRITTEN]
    if SETTLED in states:
        return EXIT_CODE[SETTLED]
    return EXIT_CODE[CONVERGED]
