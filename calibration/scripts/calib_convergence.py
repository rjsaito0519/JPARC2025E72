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
  skipped    same apply already reflected in param; nothing written this pass (4)
  pending    apply planned (before write); not an exit state; not counted in settle history

`written` must be recorded only after a successful param write. Record `pending` first with
before=/after= content hashes, then write, then `written`. Idempotency uses apply_key
(input decode token + correction content hash), not fit ROOT mtime.
`previous_values` ignores pending so one apply is one measurement for settle.
"""
import datetime
import hashlib
import json
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import config, runner  # noqa: E402

CONVERGED, SETTLED, WRITTEN, SKIPPED, PENDING = (
    "converged", "settled", "written", "skipped", "pending",
)
SETTLE_PASSES = 3
EXIT_CODE = {
    CONVERGED: runner.EXIT_OK,
    WRITTEN: runner.EXIT_NOT_CONVERGED,
    SETTLED: runner.EXIT_CHECK,
    SKIPPED: runner.EXIT_SKIP,
}


#_____________________________________________________________________________
def history_path(run, group, mode, suffix):
    chamber = "ALL" if group == "ALL" else group[:4]
    return config.OUTPUT_DIR / "calib" / f"run{run:05d}" / f"dc_{chamber}_{mode}_{suffix}.txt"


#_____________________________________________________________________________
def previous_values(path, group):
    """Values of `group` measured since its last converged pass, oldest first.

    Journal `pending` lines are ignored; only completed pass states count toward settle.
    """
    values = []
    if not path.exists():
        return values
    for line in path.read_text().splitlines():
        f = line.split("\t") if "\t" in line else line.split()
        if len(f) >= 4 and f[1] == group:
            if f[3] == PENDING:
                continue
            values = [] if f[3] == CONVERGED else values + [float(f[2])]
    return values


#_____________________________________________________________________________
def decide(value, within, previous, not_better):
    """State of this pass; not_better(a, b) says pass b did not improve on pass a."""
    if within:
        return CONVERGED
    recent = previous[-(SETTLE_PASSES - 1):] + [value]
    if len(recent) == SETTLE_PASSES and all(not_better(a, b) for a, b in zip(recent, recent[1:])):
        return SETTLED
    return WRITTEN


#_____________________________________________________________________________
def decide_at(path, group, value, within, not_better):
    """Return the state for this pass without appending to the history file."""
    return decide(value, within, previous_values(path, group), not_better)


#_____________________________________________________________________________
def decide_pass(run, group, mode, suffix, value, within, not_better):
    """decide_at() on the DC history file of this run, group, mode and particle."""
    return decide_at(history_path(run, group, mode, suffix), group, value, within, not_better)


#_____________________________________________________________________________
def sha256_hex(data, n=16):
    """Short hex digest of bytes or str."""
    if isinstance(data, str):
        data = data.encode()
    return hashlib.sha256(data).hexdigest()[:n]


#_____________________________________________________________________________
def file_sha256(path, n=16):
    path = Path(path)
    return sha256_hex(path.read_bytes(), n=n)


#_____________________________________________________________________________
def input_token(path):
    """Id of an input ROOT (name + size + mtime_ns); stable across refit of the same decode."""
    path = Path(path)
    st = path.stat()
    return f"{path.name}:{st.st_size}:{st.st_mtime_ns}"


#_____________________________________________________________________________
def corr_token(mapping):
    """Hash of the numeric correction payload about to be applied."""
    payload = json.dumps(mapping, sort_keys=True, separators=(",", ":"), default=str)
    return sha256_hex(payload)


#_____________________________________________________________________________
def make_apply_key(input_tok, corr_tok):
    return f"input={input_tok}+corr={corr_tok}"


#_____________________________________________________________________________
def fit_token(path):
    """Legacy id of a fit ROOT (name + mtime_ns); for logs only, not apply checks."""
    path = Path(path)
    return f"{path.name}:{path.stat().st_mtime_ns}"


#_____________________________________________________________________________
def _detail_map(details):
    out = {}
    for part in details.split():
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


#_____________________________________________________________________________
def details_with_apply(apply_key, before=None, after=None, extra=""):
    """Build details: optional extra, apply=..., optional before=/after=."""
    parts = []
    if extra:
        parts.append(extra.strip())
    parts.append(f"apply={apply_key}")
    if before is not None:
        parts.append(f"before={before}")
    if after is not None:
        parts.append(f"after={after}")
    return " ".join(parts)


#_____________________________________________________________________________
def details_with_fit(fit_tok, param_path=None, extra=""):
    """Legacy details helper (fit= / param_mtime=). Prefer details_with_apply for new writes."""
    parts = []
    if extra:
        parts.append(extra.strip())
    parts.append(f"fit={fit_tok}")
    if param_path is not None and Path(param_path).exists():
        parts.append(f"param_mtime={Path(param_path).stat().st_mtime_ns}")
    return " ".join(parts)


#_____________________________________________________________________________
def already_applied(hist_path, group, apply_key, param_path):
    """True when param_path already reflects a successful apply for apply_key.

    New journal: look for apply=<apply_key> with pending/written and content hashes.
    Legacy lines with fit=<token> and param_mtime= keep the old mtime rule.
    """
    hist_path = Path(hist_path)
    param_path = Path(param_path)
    if not hist_path.exists() or not param_path.exists():
        return False
    current = file_sha256(param_path)
    for line in reversed(hist_path.read_text().splitlines()):
        if not line.strip():
            continue
        fields = line.split("\t") if "\t" in line else line.split()
        if len(fields) < 4 or fields[1] != group:
            continue
        state = fields[3]
        details = fields[4] if len(fields) > 4 else ""
        dm = _detail_map(details)
        if dm.get("apply") == apply_key:
            before = dm.get("before")
            after = dm.get("after")
            if state == WRITTEN and after is not None:
                return after == current
            if state == PENDING:
                if before is not None and before == current:
                    return False
                if before is not None and before != current:
                    return True
                if after is not None and after == current:
                    return True
                return False
            if state == WRITTEN and after is None and before is not None:
                return before != current
            continue
        # Legacy: fit=<token> + optional param_mtime=
        if dm.get("fit") != apply_key or state != WRITTEN:
            continue
        if "param_mtime" not in dm:
            return True
        try:
            pm = int(dm["param_mtime"])
        except ValueError:
            return True
        return param_path.stat().st_mtime_ns >= pm
    return False


#_____________________________________________________________________________
def record(path, group, value, state, details=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as f:
        f.write(f"{datetime.datetime.now():%Y-%m-%dT%H:%M:%S}\t{group}\t{value:.6g}\t{state}"
                + (f"\t{details}" if details else "") + "\n")


#_____________________________________________________________________________
def judge_at(path, group, value, within, not_better, details=""):
    """Decide the state without recording (details ignored; kept for call-site compatibility)."""
    del details
    return decide_at(path, group, value, within, not_better)


#_____________________________________________________________________________
def judge(run, group, mode, suffix, value, within, not_better):
    """Decide on the DC history path without recording."""
    return decide_pass(run, group, mode, suffix, value, within, not_better)


#_____________________________________________________________________________
def exit_code(states):
    """Exit code of several groups (runner.worst; 0 when there are none)."""
    return runner.worst(EXIT_CODE[s] for s in states if s in EXIT_CODE)
