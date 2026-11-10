"""Score corrected BTOF / TOF width from a Hodo decode ROOT.

Uses CBTOF (BHT/BH2), CTOF (T1), CFTOF (CVC) vs DeltaE TH2s. The scalar score
is the median over segments of the robust half-width (q84-q16)/2 of the TOF
projection (summed over dE). No Gaussian auto-fit.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import uproot

from . import config

# Detector name in HODO_DETECTORS -> (hist name prefix in Hodo.root, TOF tag)
_HIST_PREFIX = {
    "BHT": ("BHT", "CBTOF"),
    "BH2": ("T0", "CBTOF"),
    "T1": ("T1", "CTOF"),
    "CVC": ("CVC", "CFTOF"),
}


#_____________________________________________________________________________
def hodo_root(run: int) -> Path:
    return Path(config.DATA_DIR) / f"run{int(run):05d}_Hodo.root"


#_____________________________________________________________________________
def _weighted_quantile(centers: np.ndarray, weights: np.ndarray, q: float) -> float:
    mask = weights > 0
    if not np.any(mask):
        return float("nan")
    x = centers[mask]
    w = weights[mask].astype(float)
    order = np.argsort(x)
    x = x[order]
    w = w[order]
    cdf = np.cumsum(w)
    total = float(cdf[-1])
    if total <= 0:
        return float("nan")
    cdf = cdf / total
    return float(np.interp(q, cdf, x))


#_____________________________________________________________________________
def robust_half_width(centers: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    """Return ((q84-q16)/2, total_weight)."""
    n = float(np.sum(weights))
    if n <= 0:
        return float("nan"), 0.0
    q16 = _weighted_quantile(centers, weights, 0.16)
    q84 = _weighted_quantile(centers, weights, 0.84)
    if not (np.isfinite(q16) and np.isfinite(q84)):
        return float("nan"), n
    return 0.5 * (q84 - q16), n


#_____________________________________________________________________________
def _pick_keys(all_keys: list[str], prefix: str, tof_tag: str, particle: str) -> list[str]:
    """Prefer particle-suffixed seg TH2s; fall back to unsuffixed."""
    base = []
    suffixed = []
    needle = f"_{tof_tag}_vs_DeltaE"
    for key in all_keys:
        name = key.split(";")[0]
        if not name.startswith(f"{prefix}_seg"):
            continue
        if needle not in name:
            continue
        if name.endswith(f"_{particle}"):
            suffixed.append(name)
        elif name.endswith(needle) or name.endswith(f"{needle}"):
            # unsuffixed: ..._CBTOF_vs_DeltaE
            if not name.endswith("_Pi") and not name.endswith("_K"):
                base.append(name)
    return suffixed if suffixed else base


#_____________________________________________________________________________
def project_tof_width(hist) -> tuple[float, float]:
    """TH2: axis0=dE, axis1=TOF. Return (half_width, entries)."""
    values = hist.values(flow=False)
    tof_edges = hist.axes[1].edges()
    weights = np.sum(values, axis=0).astype(float)
    centers = 0.5 * (tof_edges[:-1] + tof_edges[1:])
    return robust_half_width(centers, weights)


#_____________________________________________________________________________
def score_btof(
    run: int,
    *,
    particle_suffix: str = "Pi",
    ftof: bool = False,
    root_path: Path | None = None,
) -> dict[str, Any]:
    """Median CBTOF/CTOF robust half-width over hdphc detector segments."""
    path = Path(root_path) if root_path is not None else hodo_root(run)
    if not path.is_file():
        return {
            "ok": False,
            "width": float("nan"),
            "n_entries": 0.0,
            "message": f"missing {path.name}",
            "per_det": {},
        }
    particle = "K" if particle_suffix in ("K", "kaon") else "Pi"
    dets = config.hodo_detectors("hdphc", ftof)
    try:
        root = uproot.open(path)
    except Exception as exc:  # noqa: BLE001
        return {
            "ok": False,
            "width": float("nan"),
            "n_entries": 0.0,
            "message": f"open failed: {type(exc).__name__}",
            "per_det": {},
        }

    all_keys = list(root.keys())
    widths: list[float] = []
    total_n = 0.0
    per_det: dict[str, Any] = {}
    for det in dets:
        spec = _HIST_PREFIX.get(det)
        if spec is None:
            per_det[det] = {"width": float("nan"), "n_seg": 0, "n_entries": 0.0}
            continue
        prefix, tof_tag = spec
        names = _pick_keys(all_keys, prefix, tof_tag, particle)
        seg_w = []
        seg_n = 0.0
        for name in names:
            try:
                hist = root[name]
            except Exception:  # noqa: BLE001
                continue
            width, n_ent = project_tof_width(hist)
            if np.isfinite(width) and n_ent > 0:
                seg_w.append(width)
                seg_n += n_ent
        det_width = float(np.median(seg_w)) if seg_w else float("nan")
        per_det[det] = {
            "width": det_width,
            "n_seg": len(seg_w),
            "n_entries": seg_n,
        }
        if np.isfinite(det_width):
            widths.append(det_width)
        total_n += seg_n

    if not widths:
        return {
            "ok": False,
            "width": float("nan"),
            "n_entries": total_n,
            "message": "no CBTOF/CTOF segments scored",
            "per_det": per_det,
        }
    return {
        "ok": True,
        "width": float(np.median(widths)),
        "n_entries": total_n,
        "message": "ok",
        "per_det": per_det,
    }
