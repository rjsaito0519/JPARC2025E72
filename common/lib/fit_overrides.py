"""Run-local Hodo fit overrides and calib locks (#23)."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional

from . import config


#_____________________________________________________________________________
def overrides_path(run: int) -> Path:
    return (
        Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}"
        / "fit_overrides.json"
    )


#_____________________________________________________________________________
def flat_override_path(run: int) -> Path:
    """Text file consumed by C++ via E72_FIT_OVERRIDE."""
    return (
        Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}"
        / "fit_override_mip_left.txt"
    )


#_____________________________________________________________________________
def empty_doc() -> dict[str, Any]:
    return {"fit": {}, "calib_locks": {}}


#_____________________________________________________________________________
def load(run: int) -> dict[str, Any]:
    path = overrides_path(run)
    if not path.is_file():
        return empty_doc()
    try:
        doc = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return empty_doc()
    if not isinstance(doc, dict):
        return empty_doc()
    doc.setdefault("fit", {})
    doc.setdefault("calib_locks", {})
    return doc


#_____________________________________________________________________________
def save(run: int, doc: dict[str, Any]) -> Path:
    path = overrides_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(doc, indent=2, sort_keys=True) + "\n", encoding="utf-8",
    )
    write_flat_mip_left(run, doc)
    return path


#_____________________________________________________________________________
def channel_key(det: str, seg: int, side: str) -> str:
    """C++ / params.h style key: bh2-6-u."""
    return f"{det.lower()}-{int(seg)}-{side.lower()}"


#_____________________________________________________________________________
def set_mip_left(
    run: int, det: str, seg: int, side: str, mip_left: float, mode: str = "hdprm",
) -> dict[str, Any]:
    doc = load(run)
    fit = doc.setdefault("fit", {})
    det_map = fit.setdefault(det.upper(), {})
    mode_map = det_map.setdefault(mode, {})
    key = channel_key(det, seg, side)
    entry = mode_map.get(key) if isinstance(mode_map.get(key), dict) else {}
    entry["mip_left"] = float(mip_left)
    mode_map[key] = entry
    save(run, doc)
    return doc


#_____________________________________________________________________________
def write_flat_mip_left(run: int, doc: Optional[dict[str, Any]] = None) -> Path:
    doc = load(run) if doc is None else doc
    lines = ["# key mip_left  (E72_FIT_OVERRIDE)"]
    fit = doc.get("fit") if isinstance(doc.get("fit"), dict) else {}
    for det, modes in sorted(fit.items()):
        if not isinstance(modes, dict):
            continue
        hdprm = modes.get("hdprm") if isinstance(modes.get("hdprm"), dict) else {}
        for key, entry in sorted(hdprm.items()):
            if not isinstance(entry, dict):
                continue
            if "mip_left" not in entry:
                continue
            lines.append(f"{key} {float(entry['mip_left'])}")
    path = flat_override_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


#_____________________________________________________________________________
def calib_lock_keys(run: int, det: str, mode: str = "hdprm") -> set[str]:
    doc = load(run)
    locks = doc.get("calib_locks") if isinstance(doc.get("calib_locks"), dict) else {}
    det_map = locks.get(det.upper()) if isinstance(locks.get(det.upper()), dict) else {}
    rows = det_map.get(mode) if isinstance(det_map.get(mode), list) else []
    return {str(x) for x in rows}


#_____________________________________________________________________________
def set_calib_lock(
    run: int, det: str, param_key: str, locked: bool, mode: str = "hdprm",
) -> dict[str, Any]:
    doc = load(run)
    locks = doc.setdefault("calib_locks", {})
    det_map = locks.setdefault(det.upper(), {})
    rows = list(det_map.get(mode) or [])
    key = str(param_key)
    if locked and key not in rows:
        rows.append(key)
    if not locked:
        rows = [r for r in rows if r != key]
    det_map[mode] = rows
    save(run, doc)
    return doc
