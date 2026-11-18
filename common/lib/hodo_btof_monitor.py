"""Score corrected BTOF / TOF width from a Hodo decode ROOT.

Uses CBTOF (BHT/BH2), CTOF (T1), CFTOF (CVC) vs DeltaE TH2s. The scalar score
is the median over segments of the robust half-width (q84-q16)/2 of the TOF
projection (summed over dE). No Gaussian auto-fit.
"""
from __future__ import annotations

import hashlib
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
    return list(dict.fromkeys(suffixed if suffixed else base))


#_____________________________________________________________________________
def project_tof_width(hist) -> tuple[float, float]:
    """TH2: axis0=dE, axis1=TOF. Return (half_width, entries)."""
    values = hist.values(flow=False)
    tof_edges = hist.axes[1].edges()
    weights = np.sum(values, axis=0).astype(float)
    centers = 0.5 * (tof_edges[:-1] + tof_edges[1:])
    return robust_half_width(centers, weights)


#_____________________________________________________________________________
def _score_btof_impl(
    run: int,
    *,
    particle_suffix: str = "Pi",
    ftof: bool = False,
    root_path: Path | None = None,
    scoring_version: int = 2,
    quality_policy: dict | None = None,
) -> dict[str, Any]:
    """Score global BTOF/FTOF width; version 1 preserves historical segment medians."""
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

    from . import phc_quality
    if scoring_version not in (1, 2):
        root.close()
        raise ValueError("unknown PHC scoring version")
    validity = _time_validity(root, dets, particle) if scoring_version == 2 else {}
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
        segments = []
        tof_edges = None
        profiles_consistent = True
        for name in names:
            try:
                hist = root[name]
            except Exception:  # noqa: BLE001
                continue
            width, n_ent = project_tof_width(hist)
            flow_values = hist.values(flow=True)
            axes_sha = hashlib.sha256(b"".join(
                np.asarray(axis.edges(), dtype="<f8").tobytes() for axis in hist.axes
            )).hexdigest()
            segments.append({
                "histogram": name, "width": float(width) if np.isfinite(width) else None,
                "n_entries": n_ent, "n_total": float(np.sum(flow_values)),
                "n_tof_flow": float(np.sum(flow_values[1:-1, [0, -1]])),
                "n_de_flow": float(np.sum(flow_values[[0, -1], 1:-1])),
                "n_both_flow": float(np.sum(flow_values[np.ix_([0, -1], [0, -1])])),
                "axes_sha": axes_sha,
            })
            if scoring_version == 2:
                edges = hist.axes[1].edges().tolist()
                if tof_edges is not None and edges != tof_edges:
                    profiles_consistent = False
                tof_edges = edges
                segments[-1]["tof_counts"] = flow_values.sum(axis=0).astype(float).tolist()
            if np.isfinite(width) and n_ent > 0:
                seg_w.append(width)
                seg_n += n_ent
        det_width = float(np.median(seg_w)) if seg_w else float("nan")
        per_det[det] = {
            "width": det_width,
            "n_seg": len(seg_w),
            "n_entries": seg_n,
            "segments": segments,
        }
        if scoring_version == 2:
            per_det[det]["tof_edges"] = tof_edges if profiles_consistent else None
            per_det[det]["time_validity"] = validity.get(det)
            try:
                per_det[det]["global_tof"] = phc_quality.aggregate(per_det[det])
            except (ValueError, TypeError, KeyError):
                per_det[det]["global_tof"] = None
        if np.isfinite(det_width):
            widths.append(det_width)
        total_n += seg_n

    root.close()
    context = {"version": scoring_version, "particle": particle, "ftof": bool(ftof),
               "required_detectors": list(dets)}
    if scoring_version == 2:
        primary = ["BHT", "CVC"] if ftof else ["BHT"]
        context.update(primary_detectors=primary, quality_policy=phc_quality.policy(quality_policy))
        widths = [(per_det.get(det, {}).get("global_tof") or {}).get("half68") for det in primary]
        if any(value is None or not np.isfinite(value) for value in widths):
            widths = []
    if not widths:
        return {
            "ok": False,
            "width": float("nan"),
            "n_entries": total_n,
            "message": "no CBTOF/CTOF segments scored",
            "scoring_context": context,
            "per_det": per_det,
        }
    return {
        "ok": True,
        "width": float(np.mean(widths) if scoring_version == 2 else np.median(widths)),
        "n_entries": total_n,
        "message": "ok",
        "scoring_context": context,
        "per_det": per_det,
    }


def _time_validity(root, detectors, particle):
    """Track non-finite corrected times by event and ordinal, without a TOF cut."""
    import awkward as ak
    result = {}
    try:
        tree = root["hodo"]
        fields = ["event_number", "beam_flag"] + [
            det.lower() + suffix for det in detectors for suffix in ("_cmt", "_mt", "_hit_seg")]
        arrays = tree.arrays(fields, library="ak")
        arrays = arrays[arrays.beam_flag == (2 if particle == "K" else 1)]
        for det in detectors:
            values = ak.flatten(arrays[det.lower() + "_cmt"], axis=2)
            numbers = ak.broadcast_arrays(arrays.event_number, values)[0]
            ordinal = ak.local_index(values, axis=1)
            flat = ak.to_numpy(ak.flatten(values, axis=None))
            invalid = ~np.isfinite(flat)
            event_numbers = ak.to_numpy(ak.flatten(numbers, axis=None))[invalid]
            positions = ak.to_numpy(ak.flatten(ordinal, axis=None))[invalid]
            raw_signature = hashlib.sha256()
            raw_mt = arrays[det.lower() + "_mt"]
            raw_parts = [
                np.asarray(ak.to_numpy(arrays.event_number), dtype="<i8"),
                np.asarray(ak.to_numpy(ak.num(raw_mt, axis=1)), dtype="<i8"),
                np.asarray(ak.to_numpy(ak.flatten(arrays[det.lower() + "_hit_seg"], axis=None)), dtype="<i8"),
                np.asarray(ak.to_numpy(ak.flatten(ak.num(raw_mt, axis=2), axis=None)), dtype="<i8"),
            ]
            raw_times = np.asarray(ak.to_numpy(ak.flatten(raw_mt, axis=None)), dtype="<f8").copy()
            raw_times[np.isnan(raw_times)] = np.nan
            for part in [*raw_parts, raw_times]:
                raw_signature.update(part.tobytes())
            result[det] = {"n_cmt": int(len(flat)), "raw_hit_signature": raw_signature.hexdigest(),
                           "n_events": int(ak.sum(ak.num(values, axis=1) > 0)),
                           "nonfinite_keys": [[int(event), int(index)] for event, index in zip(event_numbers, positions)]}
    except (KeyError, TypeError, ValueError, IndexError, AttributeError):
        return {}
    return result


def score_btof(run, *, particle_suffix="Pi", ftof=False, root_path=None,
               scoring_version=2, quality_policy=None):
    return _score_btof_impl(run, particle_suffix=particle_suffix, ftof=ftof, root_path=root_path,
                           scoring_version=scoring_version, quality_policy=quality_policy)


def initial_integrity(run, *, particle_suffix="Pi", ftof=False, root_path=None):
    """Snapshot pre-update hit integrity independently of later quality ranking."""
    from . import phc_quality
    score = _score_btof_impl(run, particle_suffix=particle_suffix, ftof=ftof, root_path=root_path)
    if score.get("scoring_context", {}).get("version") != 2 or not score.get("per_det"):
        return None
    return phc_quality.integrity_from_score(score)
