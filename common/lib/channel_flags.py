"""Per-channel review flags for calib (#24 Phase C).

Persist under OUTPUT_DIR/calib/runXXXXX/channel_flags.jsonl.
Open needs_review rows block param writes for that mode.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from . import config

STATUSES = ("needs_review", "cleared", "locked_calib")
SOURCES = ("auto", "human")


#_____________________________________________________________________________
def flags_path(run: int) -> Path:
    return (
        Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}"
        / "channel_flags.jsonl"
    )


#_____________________________________________________________________________
def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


#_____________________________________________________________________________
def append_flag(
    run: int,
    *,
    family: str,
    mode: str,
    det: str,
    channel: str,
    source: str = "human",
    reason: str = "",
    status: str = "needs_review",
) -> dict[str, Any]:
    if source not in SOURCES:
        raise ValueError(f"source must be {SOURCES}, got {source!r}")
    if status not in STATUSES:
        raise ValueError(f"status must be {STATUSES}, got {status!r}")
    row = {
        "at": _utc_now(),
        "family": str(family),
        "mode": str(mode),
        "det": str(det),
        "channel": str(channel),
        "source": source,
        "reason": str(reason or ""),
        "status": status,
    }
    path = flags_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(row, ensure_ascii=False) + "\n")
    return row


#_____________________________________________________________________________
def read_all(run: int) -> list[dict[str, Any]]:
    path = flags_path(run)
    if not path.is_file():
        return []
    rows = []
    try:
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                obj = json.loads(line)
            except json.JSONDecodeError:
                continue
            if isinstance(obj, dict):
                rows.append(obj)
    except OSError:
        return []
    return rows


#_____________________________________________________________________________
def latest_by_channel(run: int) -> dict[tuple[str, str, str], dict[str, Any]]:
    """Latest row per (mode, det, channel)."""
    out: dict[tuple[str, str, str], dict[str, Any]] = {}
    for row in read_all(run):
        key = (
            str(row.get("mode") or ""),
            str(row.get("det") or ""),
            str(row.get("channel") or ""),
        )
        out[key] = row
    return out


#_____________________________________________________________________________
def open_needs_review(
    run: int, mode: Optional[str] = None,
) -> list[dict[str, Any]]:
    rows = []
    for row in latest_by_channel(run).values():
        if row.get("status") != "needs_review":
            continue
        if mode is not None and str(row.get("mode") or "") != str(mode):
            continue
        rows.append(row)
    return rows


#_____________________________________________________________________________
def write_blocked_by_flags(
    run: int, mode: Optional[str] = None,
) -> tuple[bool, str]:
    open_rows = open_needs_review(run, mode=mode)
    if not open_rows:
        return False, ""
    sample = open_rows[0]
    return True, (
        f"channel needs_review "
        f"({sample.get('mode')}/{sample.get('det')}/{sample.get('channel')}"
        f"; +{len(open_rows) - 1} more)"
        if len(open_rows) > 1
        else f"channel needs_review "
        f"({sample.get('mode')}/{sample.get('det')}/{sample.get('channel')})"
    )


#_____________________________________________________________________________
def clear_channel(
    run: int, *, mode: str, det: str, channel: str, reason: str = "cleared",
) -> dict[str, Any]:
    return append_flag(
        run,
        family="hodo",
        mode=mode,
        det=det,
        channel=channel,
        source="human",
        reason=reason,
        status="cleared",
    )


#_____________________________________________________________________________
def ingest_anomaly_candidates(run: int) -> int:
    """Copy anomaly_candidates.jsonl into channel_flags (once each line)."""
    path = (
        Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}"
        / "anomaly_candidates.jsonl"
    )
    if not path.is_file():
        return 0
    n = 0
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return 0
    for line in lines:
        line = line.strip()
        if not line:
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        if not isinstance(row, dict):
            continue
        append_flag(
            run,
            family=str(row.get("family") or "hodo"),
            mode=str(row.get("mode") or ""),
            det=str(row.get("det") or ""),
            channel=str(row.get("channel") or ""),
            source="auto",
            reason=str(row.get("reason") or ""),
            status="needs_review",
        )
        n += 1
    return n
