"""Run-local Hodo fit overrides and calib locks (#23)."""
from __future__ import annotations

import json
import math
import os
import tempfile
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
    from .calib_state_io import write_json
    write_json(path, doc)
    write_flat_mip_left(run, doc)
    write_flat_mip_fit(run, doc)
    write_flat_t0_fit(run, doc)
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


#_____________________________________________________________________________
def validate_mip_fit(settings):
    """Validate the four free MIP parameters and fit range, independently of sliders."""
    if not isinstance(settings, dict):
        raise ValueError("MIP settings must be an object")
    parameters = settings.get("parameters")
    fit_range = settings.get("range")
    if not isinstance(parameters, list) or len(parameters) != 4:
        raise ValueError("four MIP parameters are required")
    if not isinstance(fit_range, list) or len(fit_range) != 2:
        raise ValueError("two fit bounds are required")
    numbers = parameters + fit_range
    if any(isinstance(value, bool) or not isinstance(value, (float, int))
           or not math.isfinite(value) for value in numbers):
        raise ValueError("MIP settings must contain finite numbers")
    if parameters[0] <= 0 or parameters[2] <= 0 or parameters[3] < 0:
        raise ValueError("MIP normalization/width must be positive; background must be nonnegative")
    if fit_range[0] >= fit_range[1]:
        raise ValueError("fit lower bound must be below upper bound")
    if settings.get("function") not in ("auto", "gaussian", "landau", "mean"):
        raise ValueError("unknown MIP function")
    if settings.get("range_mode") not in ("auto", "fixed"):
        raise ValueError("unknown MIP range mode")
    return {"parameters": [float(value) for value in parameters],
            "range": [float(value) for value in fit_range],
            "function": settings["function"], "range_mode": settings["range_mode"]}


#_____________________________________________________________________________
def flat_mip_fit_path(run):
    return overrides_path(run).with_name("fit_override_mip.txt")


#_____________________________________________________________________________
def mip_fit_line(key, particle, settings):
    settings = validate_mip_fit(settings)
    numbers = settings["range"] + settings["parameters"]
    return " ".join([f"{key}:{particle}", settings["function"], settings["range_mode"],
                     *[format(value, ".17g") for value in numbers]])


#_____________________________________________________________________________
def write_flat_mip_fit(run, doc=None):
    doc = load(run) if doc is None else doc
    fit = doc.get("fit") if isinstance(doc.get("fit"), dict) else {}
    bh2 = fit.get("BH2") if isinstance(fit.get("BH2"), dict) else {}
    channels = bh2.get("hdprm") if isinstance(bh2.get("hdprm"), dict) else {}
    lines = ["# key:particle function range_mode lower upper p0 p1 p2 p3"]
    for key, entry in sorted(channels.items()):
        if not isinstance(entry, dict):
            continue
        for particle, settings in sorted(entry.get("mip_fit", {}).items()):
            if particle not in ("Pi", "K"):
                raise ValueError("unknown MIP particle")
            lines.append(mip_fit_line(key, particle, settings))
    path = flat_mip_fit_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w") as stream:
            stream.write("\n".join(lines) + "\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return path


#_____________________________________________________________________________
def validate_t0_fit(settings):
    if not isinstance(settings, dict):
        raise ValueError("T0 settings must be an object")
    parameters, bounds = settings.get("parameters"), settings.get("range")
    if not isinstance(parameters, list) or len(parameters) != 3:
        raise ValueError("three T0 Gaussian parameters are required")
    if not isinstance(bounds, list) or len(bounds) != 2:
        raise ValueError("two T0 fit bounds are required")
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
           for value in parameters+bounds):
        raise ValueError("T0 parameters and bounds must be finite")
    if parameters[0] <= 0 or parameters[2] <= 0 or bounds[0] >= bounds[1]:
        raise ValueError("T0 normalization and width must be positive, with ordered bounds")
    if settings.get("function") != "gaussian" or settings.get("range_mode") not in ("auto", "fixed"):
        raise ValueError("T0 supports Gaussian with auto or fixed range")
    return {"parameters": [float(value) for value in parameters], "range": [float(value) for value in bounds],
            "function": "gaussian", "range_mode": settings["range_mode"]}


#_____________________________________________________________________________
def t0_fit_line(key, particle, settings):
    settings = validate_t0_fit(settings)
    return " ".join([f"{key}:{particle}", "gaussian", settings["range_mode"],
                     *[format(value, ".17g") for value in settings["range"]+settings["parameters"]+[0.]]])


#_____________________________________________________________________________
def flat_t0_fit_path(run):
    return overrides_path(run).with_name("fit_override_t0.txt")


#_____________________________________________________________________________
def write_flat_t0_fit(run, doc=None):
    doc = load(run) if doc is None else doc
    channels = doc.get("fit", {}).get("BH2", {}).get("t0", {})
    lines = ["# key function range_mode lower upper p0 p1 p2 0 (E72_T0_FIT_OVERRIDE)"]
    for key, entry in sorted(channels.items()):
        for particle, settings in sorted(entry.get("gaussian_fit", {}).items()):
            if particle not in ("Pi", "K"):
                raise ValueError("unknown T0 calibration particle")
            lines.append(t0_fit_line(key, particle, settings))
    path = flat_t0_fit_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines)+"\n", encoding="utf-8")
    return path
