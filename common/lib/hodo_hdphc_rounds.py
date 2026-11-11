"""Snapshot / journal / adopt helpers for multi-pass hdphc."""
from __future__ import annotations

import json
import math
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from . import config, qa_session


#_____________________________________________________________________________
def rounds_dir(run: int) -> Path:
    return Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}" / "hdphc_rounds"


#_____________________________________________________________________________
def journal_path(run: int) -> Path:
    return rounds_dir(run) / "journal.jsonl"


#_____________________________________________________________________________
def summary_path(run: int) -> Path:
    return rounds_dir(run) / "summary.txt"


#_____________________________________________________________________________
def round_dir(run: int, round_i: int) -> Path:
    return rounds_dir(run) / f"round_{int(round_i):02d}"


#_____________________________________________________________________________
def hdphc_param_path(run: int, suffix: str) -> Path:
    return (
        Path(config.PARAM_DIR)
        / "HDPHC"
        / config.SUB_DIR
        / f"HodoPHCParam_run{int(run):05d}_{suffix}"
    )


#_____________________________________________________________________________
def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


#_____________________________________________________________________________
def clear_rounds(run: int) -> None:
    path = rounds_dir(run)
    if path.is_dir():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)
    # Starting hdphc from round 1 again: re-confirm ranges next collect_phc.
    from . import phc_range_gate

    phc_range_gate.clear_range_confirmed(int(run))


#_____________________________________________________________________________
def snapshot_round(
    run: int,
    round_i: int,
    suffix: str,
    score: dict[str, Any],
    *,
    input_path: Path | None = None,
    qa_paths=None,
) -> dict[str, Any]:
    """Copy current HDPHC param (and optional fit QA PDFs) into round_XX/."""
    dest = round_dir(run, round_i)
    dest.mkdir(parents=True, exist_ok=True)
    src = hdphc_param_path(run, suffix)
    if not src.is_file():
        raise FileNotFoundError(f"HDPHC param missing: {src}")
    shutil.copy2(src, dest / src.name)
    param_sha = qa_session.file_sha256(src)
    tok = None
    if input_path is not None and Path(input_path).is_file():
        tok = qa_session.input_token(input_path)
    qa_files: list[dict[str, Any]] = []
    if qa_paths:
        qa_dir = dest / "qa"
        qa_dir.mkdir(parents=True, exist_ok=True)
        for path in qa_paths:
            path = Path(path)
            if not path.is_file():
                qa_files.append({
                    "name": path.name,
                    "sha": None,
                    "missing": True,
                })
                continue
            sealed = qa_dir / path.name
            shutil.copy2(path, sealed)
            qa_files.append({
                "name": path.name,
                "sha": qa_session.file_sha256(sealed),
                "path": str(sealed),
                "missing": False,
            })
    entry = {
        "round": int(round_i),
        "width": score.get("width"),
        "n_entries": score.get("n_entries"),
        "ok": bool(score.get("ok")),
        "input_token": tok,
        "param_sha": param_sha,
        "param_file": src.name,
        "qa_kind": "fit",
        "qa_files": qa_files,
        "at": _utc_now(),
        "per_det": score.get("per_det") or {},
    }
    with journal_path(run).open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry, sort_keys=True) + "\n")
    return entry


#_____________________________________________________________________________
def load_journal(run: int) -> list[dict[str, Any]]:
    path = journal_path(run)
    if not path.is_file():
        return []
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        rows.append(json.loads(line))
    return rows


#_____________________________________________________________________________
def pick_round(journal: list[dict[str, Any]], adopt: str) -> int:
    """Return 1-based round index. adopt is 'auto' or an integer string."""
    if not journal:
        raise ValueError("empty hdphc journal")
    if adopt != "auto":
        k = int(adopt)
        rounds = {int(e["round"]) for e in journal}
        if k not in rounds:
            raise ValueError(f"adopt round {k} not in journal {sorted(rounds)}")
        return k
    best = None
    best_w = None
    for entry in journal:
        if not entry.get("ok"):
            continue
        width = entry.get("width")
        if width is None:
            continue
        try:
            w = float(width)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(w):
            continue
        if best_w is None or w < best_w:
            best_w = w
            best = int(entry["round"])
    if best is None:
        raise ValueError("no finite width in journal for auto adopt")
    return best


#_____________________________________________________________________________
def width_ok(value) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


#_____________________________________________________________________________
def restore_round(run: int, round_i: int, suffix: str) -> Path:
    """Copy round snapshot back to the live HDPHC param path."""
    src = round_dir(run, round_i) / f"HodoPHCParam_run{int(run):05d}_{suffix}"
    if not src.is_file():
        raise FileNotFoundError(f"round snapshot missing: {src}")
    dest = hdphc_param_path(run, suffix)
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    return dest


#_____________________________________________________________________________
def write_summary(run: int, adopted: Optional[int] = None) -> Path:
    rows = load_journal(run)
    lines = [
        f"hdphc rounds run{int(run):05d}",
        f"{'round':>5}  {'width':>10}  {'n_entries':>12}  {'ok':>3}  param_sha",
    ]
    for entry in rows:
        mark = " *" if adopted is not None and int(entry["round"]) == int(adopted) else ""
        width = entry.get("width")
        wtxt = f"{float(width):.4f}" if width_ok(width) else "nan"
        lines.append(
            f"{int(entry['round']):5d}  {wtxt:>10}  {float(entry.get('n_entries') or 0):12.0f}  "
            f"{'Y' if entry.get('ok') else 'N':>3}  {entry.get('param_sha')}{mark}"
        )
    if adopted is not None:
        lines.append(f"adopted: round {int(adopted)}")
    path = summary_path(run)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
