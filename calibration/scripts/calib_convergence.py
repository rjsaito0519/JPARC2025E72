"""Pass history and the per-group decision of iterated calibrations (DC resi/reso, D5 pos).

One text file per run, chamber (or ALL), mode and particle in OUTPUT_DIR/calib/runXXXXX/:
  dc_<BLC1|BLC2|ALL>_<resi|reso>_<K|Pi>.txt (update_param.py), d5_pos_<K|Pi>.txt (run_d5.py)
one line per pass and group: time, group, measured value, state[, details].

States (lib/runner.py exit codes; several groups or runs give the worst,
  settled/check > written > skip > converged):
  converged  within the criterion; nothing written (0)
  settled    outside the criterion but not getting better over SETTLE_PASSES passes; nothing
             written, a person looks (3)
  written    params written; decode again and run once more (2)
  skipped    same fit already applied; nothing written this pass (4)

`written` must be recorded only after a successful param write. Use decide_* then record().
"""
import datetime
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from lib import config, runner  # noqa: E402

CONVERGED, SETTLED, WRITTEN, SKIPPED = "converged", "settled", "written", "skipped"
SETTLE_PASSES = 3
EXIT_CODE = {
    CONVERGED: runner.EXIT_OK,
    WRITTEN: runner.EXIT_NOT_CONVERGED,
    SETTLED: runner.EXIT_CHECK,
    SKIPPED: runner.EXIT_SKIP,
}


def history_path(run, group, mode, suffix):
    chamber = "ALL" if group == "ALL" else group[:4]
    return config.OUTPUT_DIR / "calib" / f"run{run:05d}" / f"dc_{chamber}_{mode}_{suffix}.txt"


def previous_values(path, group):
    """Values of `group` measured since its last converged pass, oldest first."""
    values = []
    if not path.exists():
        return values
    for line in path.read_text().splitlines():
        f = line.split("\t") if "\t" in line else line.split()
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


def decide_at(path, group, value, within, not_better):
    """Return the state for this pass without appending to the history file."""
    return decide(value, within, previous_values(path, group), not_better)


def decide_pass(run, group, mode, suffix, value, within, not_better):
    """decide_at() on the DC history file of this run, group, mode and particle."""
    return decide_at(history_path(run, group, mode, suffix), group, value, within, not_better)


def fit_token(path):
    """Stable id of a fit ROOT (name + mtime_ns) for idempotent apply checks."""
    path = Path(path)
    return f"{path.name}:{path.stat().st_mtime_ns}"


def details_with_fit(fit_tok, param_path=None, extra=""):
    """Build a details string: optional extra, fit=..., optional param_mtime=..."""
    parts = []
    if extra:
        parts.append(extra.strip())
    parts.append(f"fit={fit_tok}")
    if param_path is not None and Path(param_path).exists():
        parts.append(f"param_mtime={Path(param_path).stat().st_mtime_ns}")
    return " ".join(parts)


def already_applied(hist_path, group, fit_tok, param_path):
    """True when a successful written pass for this fit is already reflected in param_path.

    Looks for the newest `written` line of `group` whose details contain fit=<fit_tok>.
    If that line has param_mtime=N and the param file's mtime is >= N, the apply succeeded.
    If param_mtime is missing (legacy line), treat as applied to avoid double correction.
    """
    hist_path = Path(hist_path)
    param_path = Path(param_path)
    if not hist_path.exists() or not param_path.exists():
        return False
    param_now = param_path.stat().st_mtime_ns
    needle = f"fit={fit_tok}"
    for line in reversed(hist_path.read_text().splitlines()):
        if not line.strip():
            continue
        fields = line.split("\t") if "\t" in line else line.split()
        if len(fields) < 4 or fields[1] != group or fields[3] != WRITTEN:
            continue
        details = fields[4] if len(fields) > 4 else ""
        if needle not in details:
            continue
        pm = None
        for part in details.split():
            if part.startswith("param_mtime="):
                try:
                    pm = int(part.split("=", 1)[1])
                except ValueError:
                    pm = None
                break
        if pm is None:
            return True
        return param_now >= pm
    return False


def record(path, group, value, state, details=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as f:
        f.write(f"{datetime.datetime.now():%Y-%m-%dT%H:%M:%S}\t{group}\t{value:.6g}\t{state}"
                + (f"\t{details}" if details else "") + "\n")


def judge_at(path, group, value, within, not_better, details=""):
    """Decide the state without recording (details ignored; kept for call-site compatibility)."""
    del details
    return decide_at(path, group, value, within, not_better)


def judge(run, group, mode, suffix, value, within, not_better):
    """Decide on the DC history path without recording."""
    return decide_pass(run, group, mode, suffix, value, within, not_better)


def exit_code(states):
    """Exit code of several groups (runner.worst; 0 when there are none)."""
    return runner.worst(EXIT_CODE[s] for s in states)
