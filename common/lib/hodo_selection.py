"""Explicit Hodo input selection, independent of the parameter particle."""
from pathlib import Path
import math
import os
import yaml

SELECTIONS = Path(__file__).resolve().parents[2] / "calibration/config/hodo_input_selection.yml"
ENVIRONMENT_KEYS = ("E72_CALIB_HIST_SELECTION", "E72_CALIB_PHC_TIME_WINDOW", "E72_CALIB_T0_TIME_WINDOW")


def selection(run, particle="Pi"):
    document = (yaml.safe_load(SELECTIONS.read_text()) or {}) if SELECTIONS.is_file() else {}
    if not isinstance(document, dict):
        raise ValueError("Hodo input selections must be a mapping")
    run_selection = document.get(int(run), {})
    if run_selection is None:
        run_selection = {}
    if not isinstance(run_selection, dict):
        raise ValueError("Per-run Hodo input selection must be a mapping")
    entry = run_selection.get(particle)
    if entry is None:
        return None
    if (not isinstance(entry, dict) or not {"histograms", "peak_window_ns"} <= set(entry)
            or set(entry) - {"histograms", "peak_window_ns", "fit_time_windows_ns", "t0_fit_window_ns"} or entry["histograms"] != "All"):
        raise ValueError("Invalid Hodo input selection")
    bounds = entry["peak_window_ns"]
    if (not isinstance(bounds, list) or len(bounds) != 2
            or any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) for x in bounds)
            or bounds[0] >= bounds[1] or not bounds[0] < 0 < bounds[1]):
        raise ValueError("Invalid zero-centered Hodo peak window")
    result = {"histograms": "All", "peak_window_ns": list(map(float, bounds))}
    if "fit_time_windows_ns" in entry:
        windows = entry["fit_time_windows_ns"]
        if not isinstance(windows, dict) or set(windows) - {"BHT", "BH2", "T1", "CVC"}:
            raise ValueError("Invalid detector PHC fit windows")
        for detector, window in windows.items():
            if (not isinstance(window, list) or len(window) != 2
                    or any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) for x in window)
                    or window[0] >= window[1]):
                raise ValueError("Invalid detector PHC fit window")
        result["fit_time_windows_ns"] = {detector: list(map(float, window)) for detector, window in windows.items()}
    if "t0_fit_window_ns" in entry:
        window = entry["t0_fit_window_ns"]
        if (not isinstance(window, list) or len(window) != 2
                or any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) for x in window)
                or window[0] >= window[1]):
            raise ValueError("Invalid raw-time T0 fit window")
        result["t0_fit_window_ns"] = list(map(float, window))
    return result


def histogram_particle(run, particle="Pi"):
    return "All" if selection(run, particle) else particle


def histogram_name(base, run, particle="Pi"):
    return base if selection(run, particle) else f"{base}_{particle}"


def environment(run, particle="Pi", base=None, *, detector=None):
    result = dict(os.environ if base is None else base)
    for key in ENVIRONMENT_KEYS:
        result.pop(key, None)
    entry = selection(run, particle)
    if entry:
        result["E72_CALIB_HIST_SELECTION"] = "All"
        window = entry.get("fit_time_windows_ns", {}).get(detector, entry["peak_window_ns"])
        result["E72_CALIB_PHC_TIME_WINDOW"] = " ".join(map(str, window))
        if "t0_fit_window_ns" in entry:
            result["E72_CALIB_T0_TIME_WINDOW"] = " ".join(map(str, entry["t0_fit_window_ns"]))
    return result


def metadata(run, particle="Pi"):
    entry = selection(run, particle)
    if not entry:
        return {}
    return {"input_selection": "All", "input_peak_window_ns": entry["peak_window_ns"],
            "output_parameter_particle": particle,
            "input_fit_time_windows_ns": entry.get("fit_time_windows_ns", {}),
            "input_t0_fit_window_ns": entry.get("t0_fit_window_ns")}


def cover_pdf(run, particle, mode, root_path):
    from . import config, qa_session
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    entry = selection(run, particle)
    if not entry:
        return None
    output = Path(config.OUTPUT_DIR) / "img" / f"run{run:05d}" / f"run{run:05d}_input_All_{particle}_{mode}.pdf"
    output.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"run{run:05d} / {mode}", "Input: All histograms (no beam PID selection)",
             f"Output parameter set: {particle}",
             f"Corrected BTOF peak evaluation window (BHT/BH2): {entry['peak_window_ns']} ns",
             f"Raw T0 fit window: {entry.get('t0_fit_window_ns')} ns",
             f"Raw-time PHC fit windows: {entry.get('fit_time_windows_ns', {})}",
             "Full profiles, total fills, outside-peak counts and nonfinite times remain checked.",
             f"Input: {Path(root_path).name}", qa_session.input_token(root_path)]
    with PdfPages(output) as pdf:
        figure = plt.figure(figsize=(11.7, 8.3))
        figure.text(.08, .9, "\n\n".join(lines), va="top", fontsize=13)
        pdf.savefig(figure)
        plt.close(figure)
    return output


def require_fit(path, run, particle):
    if not selection(run, particle):
        return
    import uproot
    with uproot.open(path) as output:
        if ("calib_input_selection" not in output or str(output["calib_input_selection"]) != "All"
                or "calib_parameter_particle" not in output or str(output["calib_parameter_particle"]) != particle):
            raise ValueError(f"{Path(path).name}: Fit input selection differs; rerun the stage with All input")


def require_sources(run, particle, mode, root_path, ftof=False):
    if not selection(run, particle):
        return
    import uproot
    from . import config
    names=[]
    if mode == "t0":
        names=[f"T0_seg{segment}_TimeOffset" for segment in range(config.NUM_OF_CH["BH2"])]
    else:
        detectors = config.hodo_detectors(mode, ftof)
        entry = selection(run, particle)
        if mode == "hdphc" and any(det not in entry.get("fit_time_windows_ns", {}) for det in detectors):
            raise ValueError("All PHC requires declared raw-time fit windows for each detector")
        for detector in detectors:
            if detector == "BAC":
                names += [f"BAC_ADC_seg{segment}U" for segment in range(5)] + ["BAC_TDC_seg4U"]
                continue
            if detector == "KVC":
                names += [f"KVC_ADC_seg{segment}{side}" for segment in range(8) for side in "abcdS"]
                names += [f"KVC_TDC_seg{segment}S" for segment in range(8)]
                continue
            sides="U" if detector=="T1" else "UD"
            for segment in range(config.NUM_OF_CH[detector]):
                for side in sides:
                    if mode=="hdprm":
                        names.append(f"{detector}_TDC_seg{segment}{side}")
                        names.append(f"BHT_Hit_TOT_seg{segment}{side}" if detector=="BHT"
                                     else f"{detector}_ADC_seg{segment}{side}")
                    else:
                        prefix,raw,corrected={"BHT":("BHT","BTOF","CBTOF"),"BH2":("T0","BTOF","CBTOF"),
                                              "T1":("T1","TOF","CTOF"),"CVC":("CVC","FTOF","CFTOF")}[detector]
                        names += [f"{prefix}_seg{segment}{side}_{kind}_vs_DeltaE" for kind in (raw,corrected)]
    with uproot.open(root_path) as root:
        missing=[name for name in names if name not in root]
    if missing:
        raise ValueError("Required All histograms missing: "+", ".join(missing[:5]))
