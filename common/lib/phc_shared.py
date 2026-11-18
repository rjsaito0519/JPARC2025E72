"""Shared Type1 PHC shape from native raw-time TH2 counts, separated by side.

The model is A*(-1/sqrt(abs(E-B)) + 1/sqrt(abs(1-B))) + c_segment.
Channel intercepts are fitted together with the shape; coefficients are never
averaged. Existing valid individual fits remain the preferred correction.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

from . import config

SUPPORTED = {"BHT": ("BHT", "BTOF", (-4., 2.)),
             "BH2": ("T0", "BTOF", (-2., 1.)),
             "CVC": ("CVC", "FTOF", (-4., 2.))}
ENERGY_RANGE = (0.6, 1.4)
MIN_SOURCE_COUNTS = 300
HUBER_NS = 0.3


def basis(energy, shift):
    return -1. / np.sqrt(np.maximum(np.abs(energy - shift), 1.e-12)) + 1. / np.sqrt(abs(1. - shift))


def histogram_name(det, segment, side, suffix):
    prefix, kind, _ = SUPPORTED[det]
    return f"{prefix}_seg{segment}{'UD'[side]}_{kind}_vs_DeltaE_{suffix}"


def cells(histogram, det):
    """Native time-window bins and the studied 0.6--1.4 MIP fit interval."""
    counts, x_edges, y_edges = histogram.to_numpy(flow=False)
    x = (x_edges[:-1] + x_edges[1:]) * .5
    y = (y_edges[:-1] + y_edges[1:]) * .5
    lo, hi = SUPPORTED[det][2]
    # ROOT FindBin endpoints are inclusive in the native ProfileX.
    first = np.searchsorted(y_edges, lo, side="right") - 1
    last = np.searchsorted(y_edges, hi, side="right") - 1
    selected = ((x[:, None] >= ENERGY_RANGE[0]) & (x[:, None] <= ENERGY_RANGE[1])
                & (np.arange(len(y))[None, :] >= first)
                & (np.arange(len(y))[None, :] <= last)
                & np.isfinite(counts) & (counts > 0))
    xi, yi = np.nonzero(selected)
    return x[xi], y[yi], counts[xi, yi].astype(float)


def weighted_median(values, weights):
    order = np.argsort(values)
    return float(np.interp(.5 * weights.sum(),
                          np.cumsum(weights[order]) - .5 * weights[order], values[order]))


def location(values, weights):
    center = weighted_median(values, weights)
    for _ in range(12):
        robust = weights * np.minimum(1., HUBER_NS / np.maximum(abs(values - center), 1.e-12))
        updated = float(np.average(values, weights=robust))
        if abs(updated - center) < 1.e-7:
            break
        center = updated
    return center


def fit_shape(samples):
    """Variable projection with Huber weights and a free intercept per source."""
    channels = sorted(samples)
    if len(channels) < 2:
        return {"ok": False, "reason": "fewer_than_two_source_segments", "sources": channels}
    energy = np.concatenate([samples[ch][0] for ch in channels])
    times = np.concatenate([samples[ch][1] for ch in channels])
    weights = np.concatenate([samples[ch][2] * min(1., 5000. / samples[ch][2].sum()) for ch in channels])
    groups = np.concatenate([np.full(len(samples[ch][0]), i) for i, ch in enumerate(channels)])
    number = len(channels)

    def solve(shift, details=False):
        shape = basis(energy, shift)
        robust = weights.copy()
        amplitude = 0.
        constants = np.zeros(number)
        for _ in range(12):
            total = np.bincount(groups, weights=robust, minlength=number)
            mx = np.bincount(groups, weights=robust * shape, minlength=number) / total
            my = np.bincount(groups, weights=robust * times, minlength=number) / total
            dx = shape - mx[groups]
            information = float(np.dot(robust, dx * dx))
            if information <= 1.e-9:
                return None if details else 1.e100
            amplitude = float(np.clip(np.dot(robust * dx, times - my[groups]) / information, .001, 15.))
            constants = my - amplitude * mx
            residual = times - amplitude * shape - constants[groups]
            robust = weights * np.minimum(1., HUBER_NS / np.maximum(abs(residual), 1.e-12))
        loss = float(np.average(np.where(abs(residual) <= HUBER_NS, residual**2 / 2.,
                                        HUBER_NS * (abs(residual) - HUBER_NS / 2.)), weights=weights))
        if details:
            return amplitude, constants, loss
        return loss

    grid = np.array([-5., -3., -2., -1., -.5, -.2, 0., .2, .4])
    losses = [solve(value) for value in grid]
    index = int(np.argmin(losses))
    left, right = grid[max(index - 1, 0)], grid[min(index + 1, len(grid) - 1)]
    result = minimize_scalar(solve, bounds=(left, right), method="bounded",
                             options={"xatol": 1.e-5})
    shift = float(result.x) if result.success and result.fun < losses[index] else float(grid[index])
    solution = solve(shift, details=True)
    if solution is None:
        return {"ok": False, "reason": "unidentifiable_shape", "sources": channels}
    amplitude, constants, loss = solution
    boundary = amplitude <= .00101 or amplitude >= 14.999 or shift <= -4.999 or shift >= .3999
    return {"ok": not boundary, "reason": "shape_at_limit" if boundary else "fitted",
            "sources": channels, "amplitude": amplitude, "shift": shift, "loss": loss,
            "constants": {str(ch): float(c) for ch, c in zip(channels, constants)},
            "energy_range": list(ENERGY_RANGE), "huber_ns": HUBER_NS}


def fallback_rows(root, det, suffix, channels, sides, donor_sides, good_ch_range=None):
    """Return fallback parameters and an auditable per-side decision report."""
    output = {}
    report = {"detector": det, "range_origin": "manual" if good_ch_range is not None else "automatic",
              "source_range": list(good_ch_range) if good_ch_range is not None else None,
              "models": {}, "channels": []}
    for side in range(sides):
        samples = {}
        for ch in channels:
            name = histogram_name(det, ch, side, suffix)
            if root is not None and name in root:
                samples[ch] = cells(root[name], det)
        donors = {}
        for ch in channels:
            if (ch, side) not in donor_sides or ch not in samples:
                continue
            x, _, w = samples[ch]
            if w.sum() >= MIN_SOURCE_COUNTS and len(np.unique(x)) >= 6 and np.ptp(x) >= .2:
                donors[ch] = samples[ch]
        model = fit_shape(donors)
        report["models"]["UD"[side]] = model
        centers = np.array(list(model.get("constants", {}).values()))
        group_center = float(np.median(centers)) if len(centers) else None
        spread = max(.05, float(1.4826 * np.median(abs(centers - group_center)))) if len(centers) else None
        for ch in channels:
            if (ch, side) in donor_sides:
                report["channels"].append({"segment": ch, "side": "UD"[side], "mode": "individual",
                                           "counts": float(samples[ch][2].sum()) if ch in samples else 0.})
                continue
            row = {"segment": ch, "side": "UD"[side], "mode": "unresolved",
                   "reason": model["reason"], "counts": float(samples[ch][2].sum()) if ch in samples else 0.}
            if model["ok"]:
                x, y, w = samples.get(ch, (np.array([]), np.array([]), np.array([])))
                n = float(w.sum())
                if n:
                    residual = y - model["amplitude"] * basis(x, model["shift"])
                    target_center = location(residual, w)
                    mad = weighted_median(abs(residual - target_center), w)
                    # Counts are histogram fills, not independent-event statistics.
                    # A floor avoids treating the coarse native y bins as exact.
                    error = max(.1, 1.4826 * mad) / np.sqrt(n)
                    shrink = float(spread**2 / (spread**2 + error**2))
                    center = shrink * target_center + (1. - shrink) * group_center
                else:
                    target_center, error, shrink, center = None, None, 0., group_center
                pars = [model["amplitude"], model["shift"],
                        model["amplitude"] / np.sqrt(abs(1. - model["shift"])) + center]
                output[(ch, side)] = [float(value) for value in pars]
                row.update(mode="shared" if n else "shared_prediction", reason="fitted" if n else "no_target_counts",
                           counts=n, parameters=output[(ch, side)], constant=float(center), target_constant=target_center,
                           constant_error=error, shrink_weight=shrink)
            report["channels"].append(row)
    report["counts"] = {mode: sum(row["mode"] == mode for row in report["channels"])
                        for mode in ("individual", "shared", "shared_prediction", "unresolved")}
    return output, report


def preserve_bh2_reference(rows, previous, reports):
    """Fix the common clock when BHT and BH2 are fitted against each other.

    Relative BH2 shapes/offsets are retained. The most populated usable BH2
    segment keeps its previous mean correction at 1 MIP, removing the global
    offset degree of freedom that otherwise updates twice in one round.
    """
    if not any(key.startswith("1-0-") for key in rows) or not any(key.startswith("3-0-") for key in rows):
        return None
    detector = next((item for item in reports if item["detector"] == "BH2"), None)
    if detector is None:
        return None
    counts = {(item["segment"], item["side"]): item.get("counts", 0.) for item in detector["channels"]}
    sources = set(detector["models"]["U"]["sources"]) & set(detector["models"]["D"]["sources"])
    candidates = []
    for segment in sources:
        keys = [f"3-0-{segment}-{side}" for side in (0, 1)]
        if all(key in rows and key in previous for key in keys):
            values = [float(v) for mapping in (previous, rows) for key in keys for v in mapping[key][:5]]
            if (len(values) == 20 and np.isfinite(values).all()
                    and all(float(mapping[key][0]) == 1 and float(mapping[key][1]) == 3
                            and abs(1. - float(mapping[key][3])) > 1.e-8
                            for mapping in (previous, rows) for key in keys)):
                candidates.append(segment)
    if not candidates:
        raise RuntimeError("no usable BH2 reference for simultaneous BHT/BH2 update")
    anchor = max(candidates, key=lambda segment: (min(counts.get((segment, side), 0.) for side in "UD"), -segment))
    def pivot(value):
        a, b, c = map(float, value[2:5])
        return -a / np.sqrt(abs(1. - b)) + c
    keys = [f"3-0-{anchor}-{side}" for side in (0, 1)]
    old = float(np.mean([pivot(previous[key]) for key in keys]))
    fitted = float(np.mean([pivot(rows[key]) for key in keys]))
    shift = old - fitted
    if not np.isfinite(shift):
        raise RuntimeError("non-finite BH2 reference adjustment")
    for key, values in rows.items():
        if key.startswith("3-0-"):
            values[4] += shift
    for channel in detector["channels"]:
        if "parameters" in channel:
            channel["fit_parameters"] = list(channel["parameters"])
            channel["parameters"][2] += shift
    record = {"detector": "BH2", "segment": anchor, "energy_mip": 1.,
              "previous_mean_correction_ns": old, "fitted_mean_correction_ns": fitted,
              "common_shift_ns": shift}
    detector["reference_anchor"] = record
    return record


def report_path(run, suffix):
    return Path(config.OUTPUT_DIR) / "calib" / f"run{run:05d}" / f"phc_shared_{suffix}.json"


def write_report(run, suffix, report):
    path = report_path(run, suffix)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    try:
        temporary.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return path


def mark_written(run, suffix, parameter):
    path = report_path(run, suffix)
    if not path.is_file():
        return
    report = json.loads(path.read_text())
    if report.get("decision") != "rows_prepared":
        return
    from . import qa_session
    report["decision"] = "written"
    report["param_sha"] = qa_session.file_sha256(parameter)
    write_report(run, suffix, report)


def summary_line(report):
    counts = report.get("counts", {})
    origin = {"manual": "手動設定", "automatic": "自動選択",
              "session_automatic": "今回のみ採用"}.get(report.get("range_origin"), "不明")
    reference = report.get("reference_anchor")
    reference_note = f"; 基準 seg{reference['segment']} 固定" if reference else ""
    return (f"{report.get('detector', '?')}: 個別 {counts.get('individual', 0)} / "
            f"共有 {counts.get('shared', 0)} / 予測のみ {counts.get('shared_prediction', 0)} / "
            f"未解決 {counts.get('unresolved', 0)} "
            f"({origin}; U/D別{reference_note})")

def render_pdf(report, source_root, output, execution_id=""):
    """Show the actual fallback curves; the original individual-fit pages follow."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LogNorm
    import uproot

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(output) as pdf, uproot.open(source_root) as root:
        fig = plt.figure(figsize=(11.7, 8.3))
        lines = [f"Shared PHC fallback -- run{report['run']:05d} ({report['suffix']})",
                 f"Execution: {execution_id}", f"Decision: {report['decision']}",
                 "", "Good individual fits are retained. Counts below refer to readout sides.",
                 "Shared = raw histogram fit with common shape and channel offsets.",
                 "Prediction only = no target counts; not an individually measured calibration.",
                 "Unresolved = no update; previous parameter retained.", ""]
        for detector in report["detectors"]:
            counts = detector["counts"]
            lines.append(f"{detector['detector']}: individual {counts['individual']}, shared {counts['shared']}, "
                         f"prediction {counts['shared_prediction']}, unresolved {counts['unresolved']}")
            lines.append(f"  Source range: {detector['range_origin']} {detector['source_range']}")
            if detector.get("reference_anchor"):
                anchor = detector["reference_anchor"]
                lines.append(f"  Clock reference: BH2 seg{anchor['segment']} at 1 MIP; common shift {anchor['common_shift_ns']:+.4f} ns")
            for side, model in detector["models"].items():
                lines.append(f"  {side}: {model['reason']}; sources {model['sources']}")
        lines += ["", "Native raw-time histograms; fit region 0.6--1.4 MIP; U/D and detectors separate.",
                  "Sparse offsets use shrinkage. Histogram fills may include correlated pairs.",
                  "Curve quality must also be checked after native re-decode.",
                  "Following detector pages show the original individual fits."]
        fig.text(.06, .94, "\n".join(lines), va="top", family="monospace", fontsize=10)
        pdf.savefig(fig); plt.close(fig)
        for detector in report["detectors"]:
            rows = [r for r in detector["channels"] if r["mode"] != "individual" and r.get("counts", 0) > 0]
            for first in range(0, len(rows), 8):
                fig, axes = plt.subplots(4, 2, figsize=(11.7, 12))
                for ax, row in zip(axes.flat, rows[first:first + 8]):
                    name = histogram_name(detector["detector"], row["segment"], "UD".index(row["side"]), report["suffix"])
                    hist = root[name]
                    counts, xedges, yedges = hist.to_numpy(flow=False)
                    ax.pcolormesh(xedges, yedges, np.ma.masked_less_equal(counts.T, 0), norm=LogNorm(), rasterized=True)
                    if "parameters" in row:
                        a, b, c = row["parameters"]
                        x = np.linspace(*ENERGY_RANGE, 200)
                        ax.plot(x, -a / np.sqrt(abs(x - b)) + c, color="red", lw=1.5, label="applied")
                        if "fit_parameters" in row:
                            fa, fb, fc = row["fit_parameters"]
                            ax.plot(x, -fa / np.sqrt(abs(x - fb)) + fc, color="gray", ls="--", label="raw fit")
                            ax.legend(fontsize=7)
                    ax.set(xlim=(0, 3), ylim=SUPPORTED[detector["detector"]][2],
                           title=f"{detector['detector']} {row['segment']}{row['side']} -- {row['mode']} (n={row['counts']:g})",
                           xlabel="DeltaE [MIP]", ylabel="reference - raw time [ns]")
                for ax in list(axes.flat)[len(rows[first:first + 8]):]:
                    ax.set_visible(False)
                fig.tight_layout()
                pdf.savefig(fig); plt.close(fig)
    return output
