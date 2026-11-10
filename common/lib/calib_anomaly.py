"""Provisional auto-QA anomaly checks (#24).

Escalate existing skip/flag/bad-segment signals to EXIT_CHECK so dispatch
enters held. Thresholds are provisional — confirm against data before
tightening. Does not invent new fit-width algorithms.
"""
from __future__ import annotations

import json
from pathlib import Path

from . import config, runner

# Provisional ratios (confirm with data; wiki: 要データ確認).
HDPRM_SKIP_RATIO = 0.25
PHC_BAD_RATIO = 0.25
HTOF_FLAG_RATIO = 0.25
DRIFT_MIN_GRAPHS = 1


#_____________________________________________________________________________
def last_check_reason_path(run: int) -> Path:
    return (
        Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}"
        / "last_check_reason.txt"
    )


#_____________________________________________________________________________
def remember_check_reason(run: int, reason: str) -> None:
    path = last_check_reason_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(reason.strip() + "\n", encoding="utf-8")
    runner.warn(f"[ANOMALY] {reason}")


#_____________________________________________________________________________
def read_check_reason(run: int) -> str:
    path = last_check_reason_path(run)
    if not path.is_file():
        return ""
    try:
        return path.read_text(encoding="utf-8").strip()
    except OSError:
        return ""


#_____________________________________________________________________________
def ratio_too_high(bad: int, total: int, limit: float) -> bool:
    if total <= 0:
        return False
    return (bad / total) >= limit


#_____________________________________________________________________________
def bac_skip_counts(tree) -> tuple[int, int, list[str]]:
    """Return (bad, total, channel labels) for BAC adc_ok/tdc_ok skips."""
    labels = []
    total = 0
    if "tdc_ok" in tree:
        total += 1
        if not tree["tdc_ok"][0]:
            labels.append("BAC:tdc_sum")
    if "adc_ok" in tree:
        for i, ch in enumerate(tree["ch"]):
            total += 1
            if not tree["adc_ok"][i]:
                labels.append(f"BAC:seg{int(ch)}:adc")
    return len(labels), total, labels


#_____________________________________________________________________________
def kvc_skip_counts(tree) -> tuple[int, int, list[str]]:
    labels = []
    total = 0
    if "tdc_ok" not in tree:
        return 0, 0, []
    for i, ch in enumerate(tree["ch"]):
        total += 1
        if not tree["tdc_ok"][i]:
            labels.append(f"KVC:seg{int(ch)}:tdc")
    return len(labels), total, labels


#_____________________________________________________________________________
def htof_adc_flag_counts(flag_rows: list) -> tuple[int, int, list[str]]:
    """flag_rows from update_hdprm.htof_adc_quality_report; total = len(rows)*3 sides
    is unknown here — caller passes total_sides."""
    labels = [
        f"HTOF:seg{seg}{side}:flag{flag}" for seg, side, flag, _mip in flag_rows
    ]
    return len(labels), labels


#_____________________________________________________________________________
def phc_bad_counts(bad_segs: list, n_seg: int) -> tuple[int, int, list[str]]:
    labels = [f"PHC:seg{ch}" for ch in bad_segs]
    return len(labels), int(n_seg), labels


#_____________________________________________________________________________
def decide_ratio(
    family: str,
    bad: int,
    total: int,
    limit: float,
    labels: list[str],
) -> str | None:
    if not ratio_too_high(bad, total, limit):
        return None
    sample = ",".join(labels[:8])
    more = "" if len(labels) <= 8 else f",…(+{len(labels) - 8})"
    return (
        f"{family} anomaly: {bad}/{total} "
        f"(>= {limit:.0%}) skipped/flagged [{sample}{more}]"
    )


#_____________________________________________________________________________
def count_drift_graphs(param_root: Path) -> int:
    """Count TGraph objects in a DCDriftParam ROOT file."""
    if not param_root.is_file():
        return 0
    try:
        import uproot
    except ImportError:
        return 0
    try:
        with uproot.open(param_root) as f:
            n = 0
            for key in f.keys():
                try:
                    obj = f[key]
                except Exception:
                    continue
                cls = getattr(obj, "classname", lambda: "")()
                if "TGraph" in str(cls):
                    n += 1
            return n
    except Exception:
        return 0


#_____________________________________________________________________________
def drift_graph_anomaly(param_root: Path, min_graphs: int = DRIFT_MIN_GRAPHS) -> str | None:
    n = count_drift_graphs(param_root)
    if n >= min_graphs:
        return None
    return f"drift anomaly: only {n} TGraph in {param_root.name} (need >= {min_graphs})"


#_____________________________________________________________________________
def candidates_payload(
    family: str,
    mode: str,
    det: str,
    labels: list[str],
    reason: str,
) -> list[dict]:
    """Channel-flag style dicts for Phase C (source=auto)."""
    out = []
    for lab in labels:
        out.append({
            "family": family,
            "mode": mode,
            "det": det,
            "channel": lab,
            "source": "auto",
            "reason": reason,
            "status": "needs_review",
        })
    return out


#_____________________________________________________________________________
def dump_candidates(run: int, rows: list[dict]) -> None:
    """Append auto candidates for later channel_flags ingestion."""
    if not rows:
        return
    path = (
        Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}"
        / "anomaly_candidates.jsonl"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row, ensure_ascii=False) + "\n")
