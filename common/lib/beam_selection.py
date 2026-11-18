"""Check explicitly declared calibration selections against native beam tags."""
from pathlib import Path

import awkward as ak
import numpy as np
import uproot
import yaml

REQUIREMENTS = Path(__file__).resolve().parents[2] / "calibration/config/beam_pid_requirements.yml"


def requirement(run, particle):
    document = (yaml.safe_load(REQUIREMENTS.read_text()) or {}) if REQUIREMENTS.is_file() else {}
    if not isinstance(document, dict):
        raise ValueError("beam PID requirements must be a mapping")
    record = document.get(int(run), {})
    if not isinstance(record, dict):
        raise ValueError("per-run beam PID requirements must be a mapping")
    entry = record.get(particle)
    if entry is not None and (not isinstance(entry, dict) or entry.get("bac") not in ("require", "veto")):
        raise ValueError("invalid explicit BAC calibration requirement")
    return entry


def check(run, particle, root_path, *, parameter_files=None):
    expected = requirement(run, particle)
    if expected is None:
        return None
    result = {"ok": False, "bac": expected["bac"], "particle": particle, "selected_events": 0,
              "violations": None, "reason": "beam_selection_unavailable"}
    try:
        with uproot.open(root_path) as root:
            if "param/USER" in root:
                item = root["param/USER"]
                text = "\n".join(str(line) for line in item.member("fLines")) if item.classname == "TMacro" else str(item)
            else:
                text = Path((parameter_files or {})["USER"]).read_text()
            windows = []
            for line in text.splitlines():
                fields = line.split("#", 1)[0].split()
                if fields and fields[0] == "BAC_TDC":
                    windows.append(tuple(map(float, fields[1:])))
            if (len(windows) != 1 or len(windows[0]) != 2 or not np.isfinite(windows[0]).all()
                    or windows[0][0] >= windows[0][1]):
                raise ValueError("BAC_TDC window unavailable or ambiguous")
            lo, hi = windows[0]
            selected = violations = 0
            for values in root["hodo"].iterate(["beam_flag", "bac_tdc_u"], step_size=10000, library="ak"):
                mask = values.beam_flag == (2 if particle == "K" else 1)
                tdcs = values.bac_tdc_u[mask]
                hit = ak.any(ak.any((tdcs > lo) & (tdcs < hi), axis=-1), axis=-1)
                selected += int(ak.sum(mask))
                violations += int(ak.sum(hit if expected["bac"] == "veto" else ~hit))
            result.update(ok=selected > 0 and violations == 0, selected_events=selected,
                          violations=violations, reason="ok" if selected > 0 and violations == 0 else
                          "no_selected_particle" if selected == 0 else "beam_selection_mismatch",
                          bac_window=[lo, hi])
    except (OSError, KeyError, ValueError, TypeError, IndexError) as error:
        result["detail"] = str(error)
    return result


def require(run, particle, root_path, *, parameter_files=None):
    result = check(run, particle, root_path, parameter_files=parameter_files)
    if result is not None and not result["ok"]:
        from . import calib_anomaly, runner
        reason = (f"beam_selection_mismatch: run{int(run):05d} {particle}, BAC {result['bac']}; "
                  f"selected={result['selected_events']}, violations={result['violations']}. "
                  "粒子選別が指定条件と不一致です。analyzerの選別設定を確認し、新しくdecodeしてください。")
        calib_anomaly.remember_check_reason(run, reason)
        runner.warn(reason)
        raise SystemExit(runner.EXIT_CHECK)
    return result
