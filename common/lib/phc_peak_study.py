"""Offline PHC comparisons on identical BHT/BH2 hit pairs; never writes parameters."""
from __future__ import annotations

from pathlib import Path

import awkward as ak
import numpy as np
import uproot
from scipy.optimize import least_squares
from scipy.special import ndtr, xlogy

PAIR_KEY = ("event", "bht_seg", "bht_entry", "bh2_seg", "bh2_entry")


#_____________________________________________________________________________
def bh2_offsets(parameter):
    offsets = {}
    for line in Path(parameter).read_text().splitlines():
        fields = line.split("#", 1)[0].split()
        if len(fields) >= 7 and tuple(map(int, fields[:2])) == (3, 0):
            if tuple(map(int, fields[3:5])) == (1, 2):
                offsets[int(fields[2])] = -float(fields[5]) * float(fields[6])
    if not offsets:
        raise ValueError("No BH2 TimeOffset rows in HDPRM")
    return offsets


#_____________________________________________________________________________
def skim_pairs(root_path, hdprm, particle):
    """Match native BHT_CBTOF: all hit pairs, with no corrected-time selection."""
    if particle not in ("Pi", "K"):
        raise ValueError("particle must be Pi or K")
    offsets = bh2_offsets(hdprm)
    columns = {name: [] for name in (*PAIR_KEY, "bht_cmt", "bh2_cmt", "offset", "raw_btof")}
    seen_events = set()
    branches = ["run_number", "event_number", "beam_flag"]
    for detector in ("bht", "bh2"):
        branches.extend(f"{detector}_{suffix}" for suffix in ("hit_seg", "mt", "cmt"))
    with uproot.open(root_path) as source:
        tree = source["hodo"]
        for chunk in tree.iterate(branches, step_size=10000, library="ak"):
            selected = chunk[chunk.beam_flag == (2 if particle == "K" else 1)]
            for event in ak.to_list(selected):
                event_number = int(event["event_number"])
                if event_number in seen_events:
                    raise ValueError(f"Duplicate event number: {event_number}")
                seen_events.add(event_number)
                for detector in ("bht", "bh2"):
                    segments = event[f"{detector}_hit_seg"]
                    if len(set(segments)) != len(segments):
                        raise ValueError("Repeated segment in one event")
                for bht_index, bht_segment in enumerate(event["bht_hit_seg"]):
                    for bht_entry, bht_time in enumerate(event["bht_cmt"][bht_index]):
                        for bh2_index, bh2_segment in enumerate(event["bh2_hit_seg"]):
                            for bh2_entry, bh2_time in enumerate(event["bh2_cmt"][bh2_index]):
                                offset = offsets[bh2_segment]
                                raw_btof = (event["bh2_mt"][bh2_index][bh2_entry] + offset
                                            - event["bht_mt"][bht_index][bht_entry])
                                values = (event_number, bht_segment, bht_entry, bh2_segment, bh2_entry,
                                          bht_time, bh2_time, offset, raw_btof)
                                for name, value in zip(columns, values):
                                    columns[name].append(value)
        histogram = source[f"BHT_CBTOF_vs_DeltaE_{particle}"]
        native_counts = np.sum(histogram.values(flow=True), axis=0)
        edges = histogram.axes[1].edges()
    pairs = {name: np.asarray(values, dtype=np.int64 if name in PAIR_KEY else float)
             for name, values in columns.items()}
    btof = pairs["bh2_cmt"] + pairs["offset"] - pairs["bht_cmt"]
    finite = btof[np.isfinite(btof)]
    counts = np.r_[np.sum(finite < edges[0]), np.histogram(finite, edges)[0],
                   np.sum(finite >= edges[-1])]
    # ROOT treats the upper axis edge as overflow, unlike numpy.histogram.
    counts[-2] -= np.sum(finite == edges[-1])
    # ROOT places NaN in the upper overflow bin. Keep it visible as invalid data.
    counts[0] += np.sum(np.isneginf(btof))
    counts[-1] += np.sum(np.isnan(btof) | np.isposinf(btof))
    if not np.array_equal(counts, native_counts):
        raise ValueError(f"Hit reconstruction differs from native histogram: "
                         f"{np.sum(np.abs(counts-native_counts))} counts; "
                         f"nonfinite={np.sum(~np.isfinite(btof))}; native={native_counts.sum()} vs finite={counts.sum()}; "
                         f"flow_native={native_counts[[0,-1]]} flow_rebuilt={counts[[0,-1]]}")
    return pairs, {"native_histogram": histogram.name, "native_match": True,
                   "n_pairs": len(btof), "n_events": len(seen_events),
                   "nonfinite_pairs": int(np.sum(~np.isfinite(btof)))}


#_____________________________________________________________________________
def pair_keys(pairs):
    keys = np.empty(len(pairs["event"]), dtype=[(name, "i8") for name in PAIR_KEY])
    for name in PAIR_KEY:
        keys[name] = pairs[name]
    if len(np.unique(keys)) != len(keys):
        raise ValueError("Duplicate hit-pair identity")
    return keys


#_____________________________________________________________________________
def align_pairs(candidates):
    keys = [pair_keys(pairs) for pairs in candidates]
    common = keys[0]
    for current in keys[1:]:
        common = np.intersect1d(common, current, assume_unique=True)
    if len(common) != len(keys[0]) or any(len(current) != len(common) for current in keys):
        raise ValueError("Candidate event/hit-pair membership changed; no silent intersection")
    aligned = []
    for pairs, current in zip(candidates, keys):
        order = np.argsort(current)
        indices = order[np.searchsorted(current[order], common)]
        aligned.append({name: values[indices] for name, values in pairs.items()})
    for candidate in aligned[1:]:
        if not np.array_equal(candidate["raw_btof"], aligned[0]["raw_btof"], equal_nan=True):
            raise ValueError("Uncorrected times changed between candidates")
        if not np.array_equal(candidate["offset"], aligned[0]["offset"]):
            raise ValueError("BH2 TimeOffset changed between candidates")
    return aligned


#_____________________________________________________________________________
def event_partition(event_numbers):
    """A deterministic event split keeps every pair of one event in the same set."""
    hashed = np.asarray(event_numbers, dtype=np.uint64).copy()
    with np.errstate(over="ignore"):
        hashed = (hashed ^ (hashed >> 30)) * np.uint64(0xbf58476d1ce4e5b9)
        hashed = (hashed ^ (hashed >> 27)) * np.uint64(0x94d049bb133111eb)
    return ((hashed ^ (hashed >> 31)) & np.uint64(1)) == 0


#_____________________________________________________________________________
def gaussian_core(values, half_range=1.0, bin_width=0.025, min_entries=300, event_ids=None):
    """Binned Poisson fit of an integrated Gaussian plus constant in a fixed window.

    This is a study diagnostic. No production acceptance or zero-position cut is set.
    The returned covariance ignores correlations between hit pairs in one event.
    """
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    edges = np.linspace(-half_range, half_range, round(2*half_range/bin_width)+1)
    counts = np.histogram(finite, edges)[0].astype(float)
    statistics = {"n_window": int(counts.sum()), "n_nonfinite": int(len(values)-len(finite))}
    if event_ids is not None:
        event_ids = np.asarray(event_ids)
        if event_ids.shape != values.shape:
            raise ValueError("Event identifiers must match the sample shape")
        statistics["n_events_window"] = len(np.unique(event_ids[np.isfinite(values) & (np.abs(values) <= half_range)]))
    if counts.sum() < min_entries or statistics.get("n_events_window", min_entries) < min_entries:
        return dict(statistics, ok=False, reason="low_statistics")
    centers = (edges[:-1]+edges[1:])/2
    search = np.abs(centers) < min(0.5, half_range*0.8)
    peak = float(centers[search][np.argmax(counts[search])])
    initial = [counts.sum()*0.8, peak, min(0.2, half_range/3), np.percentile(counts, 20)]
    lower = [0, -half_range*0.9, bin_width/2, 0]
    upper = [counts.sum()*100, half_range*0.9, half_range*2, counts.max()*10+1]
    def prediction(parameters):
        area, mean, sigma, background = parameters
        return area*np.diff(ndtr((edges-mean)/sigma)) + background
    def residual(parameters):
        expected = np.maximum(prediction(parameters), 1e-12)
        deviance = 2*(expected-counts+xlogy(counts, counts/expected))
        return np.sign(counts-expected)*np.sqrt(np.maximum(deviance, 0))
    result = least_squares(residual, initial, bounds=(lower, upper), max_nfev=500)
    area, mean, sigma, background = result.x
    covariance = np.linalg.pinv(result.jac.T @ result.jac)
    bounded = bool(np.any(result.active_mask[:3]))
    return {"ok": bool(result.success and not bounded), "reason": "bound" if bounded else result.message,
            "sigma": float(sigma), "mean": float(mean), "sigma_error_pair_naive": float(np.sqrt(covariance[2, 2])),
            "area": float(area), "background_per_bin": float(background),
            "deviance_ndf": float(np.sum(result.fun**2)/(len(counts)-4)),
            "n_finite": len(finite), "half_range": half_range, "bin_width": bin_width, **statistics}


#_____________________________________________________________________________
def btof_values(pairs):
    return pairs["bh2_cmt"] + pairs["offset"] - pairs["bht_cmt"]


#_____________________________________________________________________________
def channel_mixture(candidates, train, half_range=1.0, min_entries=300):
    """Study a per-segment choice with the other detector fixed to saved round 1.

    All readout sides of a segment come from the same round. Low statistics keep
    the saved round-1 reference; this is not the pre-calibration parameter set.
    """
    reference = candidates[0]
    chosen_times = {}
    choices = {}
    for detector in ("bht", "bh2"):
        chosen_times[detector] = reference[f"{detector}_cmt"].copy()
        for segment in np.unique(reference[f"{detector}_seg"]):
            selected = train & (reference[f"{detector}_seg"] == segment)
            fits = []
            for candidate in candidates:
                bht = candidate["bht_cmt"] if detector == "bht" else reference["bht_cmt"]
                bh2 = candidate["bh2_cmt"] if detector == "bh2" else reference["bh2_cmt"]
                fits.append(gaussian_core((bh2+reference["offset"]-bht)[selected], half_range,
                                         min_entries=min_entries, event_ids=reference["event"][selected]))
            eligible = [index for index, fit in enumerate(fits)
                        if fit["ok"] and fit["n_nonfinite"] <= fits[0]["n_nonfinite"]]
            chosen = min(eligible, key=lambda index: fits[index]["sigma"]) if eligible else 0
            mask = reference[f"{detector}_seg"] == segment
            chosen_times[detector][mask] = candidates[chosen][f"{detector}_cmt"][mask]
            choices[f"{detector}:{segment}"] = {"candidate_index": chosen,
                    "reason": "train_minimum" if eligible else "keep_reference_unresolved", "fits": fits}
    values = chosen_times["bh2"] + reference["offset"] - chosen_times["bht"]
    return values, choices
