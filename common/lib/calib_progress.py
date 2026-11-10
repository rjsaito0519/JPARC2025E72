"""Calibration progress.jsonl helpers (dispatch + iterate substages)."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from . import config


#_____________________________________________________________________________
def progress_path(run: int) -> Path:
    return Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}" / "progress.jsonl"


#_____________________________________________________________________________
def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


#_____________________________________________________________________________
def append_progress(
    run: int,
    event: str,
    *,
    step: Any = None,
    index: Any = None,
    rc: Any = None,
    status: Any = None,
    hold_source: Any = None,
    stage: Any = None,
) -> None:
    """Append one progress line under OUTPUT_DIR/calib/runXXXXX/."""
    path = progress_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "at": _utc_now(),
        "event": event,
        "step": step,
        "index": index,
        "rc": rc,
        "status": status,
        "hold_source": hold_source,
        "stage": stage,
    }
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(row, sort_keys=True) + "\n")


#_____________________________________________________________________________
def note_substage(
    run: int,
    step: str,
    stage: str,
    *,
    rc: Optional[int] = None,
    status: str = "done",
) -> None:
    """Record a fine-grained substage inside a dispatch step."""
    append_progress(
        int(run),
        "substage",
        step=str(step),
        stage=str(stage),
        rc=rc,
        status=str(status),
    )


#_____________________________________________________________________________
def read_progress_tail(run: int, n: int = 8) -> list[dict[str, Any]]:
    """Return the last n progress rows (dicts). Missing file → []."""
    path = progress_path(run)
    if not path.is_file():
        return []
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return []
    out: list[dict[str, Any]] = []
    for line in lines[-max(1, int(n)):]:
        line = line.strip()
        if not line:
            continue
        try:
            row = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(row, dict):
            out.append(row)
    return out


#_____________________________________________________________________________
def latest_substage(run: int, *, scan: int = 80) -> Optional[dict[str, Any]]:
    """Most recent event=substage row, or None."""
    for row in reversed(read_progress_tail(run, n=scan)):
        if row.get("event") == "substage" and row.get("stage"):
            return row
    return None
