"""Explicit unattended calibration policy and a separate automatic-decision journal."""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path

from . import config
from .calib_state_io import locked_path


def enabled():
    return os.environ.get("E72_CALIB_DISPATCH") == "1" and os.environ.get("E72_CALIB_AUTO") == "1"


def note(run, kind, **details):
    path = Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}" / "auto_decisions.jsonl"
    row = {"at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "run": int(run), "source": "auto", "kind": kind, **details}
    with locked_path(path):
        with path.open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + "\n")
    return row


def retire_previous_auto_flags(run):
    """Recompute automatic findings on restart; retain human flags and calibration locks."""
    from . import channel_flags
    retired = 0
    for row in channel_flags.open_needs_review(int(run)):
        if row.get("source") != "auto":
            continue
        channel_flags.append_flag(
            int(run), family=row.get("family") or "hodo", mode=row["mode"],
            det=row["det"], channel=row["channel"], source="auto", status="cleared",
            reason="explicit auto restart; findings will be recomputed from fresh fits",
        )
        retired += 1
    if retired:
        note(run, "retired_previous_flags", count=retired)
    return retired
