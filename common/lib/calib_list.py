"""Load calibration/caliblist/*.yml (batch DEFAULT + RUN)."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from . import runner


#_____________________________________________________________________________
def _as_bool(value, *, default=False):
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    text = str(value).strip().lower()
    if text in ("1", "true", "yes", "on"):
        return True
    if text in ("0", "false", "no", "off", ""):
        return False
    raise ValueError(f"expected bool, got {value!r}")


#_____________________________________________________________________________
def _parse_run_key(key) -> int:
    text = str(key).strip()
    if text.lower().startswith("run"):
        text = text[3:]
    try:
        run = int(text)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"bad RUN key {key!r}") from exc
    if run < 1:
        raise ValueError(f"run must be >= 1, got {run}")
    return run


#_____________________________________________________________________________
def _merge_opts(default: dict[str, Any], override: dict[str, Any] | None) -> dict[str, Any]:
    out = dict(default)
    if isinstance(override, dict):
        for key, value in override.items():
            if value is not None:
                out[key] = value
    return out


#_____________________________________________________________________________
def normalize_run_opts(raw: dict[str, Any]) -> dict[str, Any]:
    """Normalize one run's options (after DEFAULT merge)."""
    opts = dict(raw or {})
    hodo = _as_bool(opts.get("hodo"), default=False)
    bcout = _as_bool(opts.get("bcout"), default=False)
    bcin = _as_bool(opts.get("bcin"), default=False)
    d5 = _as_bool(opts.get("d5"), default=False)
    if d5:
        bcout = True
        bcin = True
    if not (hodo or bcout or bcin or d5):
        # Same default as calib_session.steps_from_flags with no flags.
        hodo = bcout = bcin = d5 = True
    particle = opts.get("calib_particle", opts.get("particle", "pion"))
    max_iter = int(opts.get("max_iter", 10))
    hdphc_iter = int(opts.get("hdphc_iter", 3))
    if max_iter < 1 or max_iter > 100:
        raise ValueError(f"max_iter must be in [1, 100], got {max_iter}")
    if hdphc_iter < 1:
        raise ValueError(f"hdphc_iter must be >= 1, got {hdphc_iter}")
    out = {
        "calib_particle": particle,
        "hodo": hodo,
        "bcout": bcout,
        "bcin": bcin,
        "d5": d5,
        "local": _as_bool(opts.get("local"), default=False),
        "force": _as_bool(opts.get("force"), default=False),
        "auto": _as_bool(opts.get("auto"), default=False),
        "dry_run": _as_bool(
            opts.get("dry_run", opts.get("dry")), default=False,
        ),
        "max_iter": max_iter,
        "hdphc_iter": hdphc_iter,
        "write_table": _as_bool(opts.get("write_table"), default=False),
    }
    if opts.get("momentum") is not None:
        out["momentum"] = float(opts["momentum"])
    if opts.get("beam_particle") is not None:
        out["beam_particle"] = opts["beam_particle"]
    if opts.get("charge") is not None:
        out["charge"] = int(opts["charge"])
    return out


#_____________________________________________________________________________
def load_caliblist(path) -> list[tuple[int, dict[str, Any]]]:
    """Return ordered [(run, opts), ...] from a caliblist YAML."""
    path = Path(path).expanduser()
    if not path.is_file():
        runner.fail(f"caliblist not found: {path}")
    try:
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        runner.fail(f"caliblist read failed: {path}: {exc}")
    if not isinstance(doc, dict):
        runner.fail(f"caliblist root must be a mapping: {path}")
    default = doc.get("DEFAULT") or doc.get("default") or {}
    if default is None:
        default = {}
    if not isinstance(default, dict):
        runner.fail(f"DEFAULT must be a mapping: {path}")
    runs = doc.get("RUN") or doc.get("run") or {}
    if not isinstance(runs, dict) or not runs:
        runner.fail(f"RUN must be a non-empty mapping: {path}")
    out: list[tuple[int, dict[str, Any]]] = []
    for key, override in runs.items():
        try:
            run = _parse_run_key(key)
            if override is None:
                override = {}
            if not isinstance(override, dict):
                raise ValueError(f"RUN.{key} must be a mapping or null")
            merged = _merge_opts(default, override)
            opts = normalize_run_opts(merged)
        except (TypeError, ValueError) as exc:
            runner.fail(f"caliblist {path}: RUN {key!r}: {exc}")
        out.append((run, opts))
    return out
