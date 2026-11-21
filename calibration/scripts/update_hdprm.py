"""HDPRM rows {"CId-PlId-SegId-AorT-UorD": [values]} from the Hodo fit outputs (results/root).

AorT: 0 = ADC (values pedestal, gain), 1 = TDC (values offset, TDC_GAIN).
The detector is given by the caller (update_param.py); CId and the number of UorD rows come from
lib/config.py.
"""
import sys
from pathlib import Path

import uproot

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import config, fit_overrides, fit_validity, runner  # noqa: E402


def _drop_locked(run, det, data, mode="hdprm"):
    """Remove calib-locked param keys so auto update does not overwrite them."""
    if run is None:
        return data
    locked = fit_overrides.calib_lock_keys(int(run), det, mode=mode)
    if not locked:
        return data
    out = {}
    for key, val in data.items():
        if key in locked:
            runner.warn(f"{det}: calib lock keeps current value for {key}")
            continue
        out[key] = val
    return out

TDC_GAIN = -0.0009765625  # ns per TDC channel (-1/1024), written as the TDC gain column


def read_tree(root_file_path):
    return uproot.open(root_file_path)["tree"].arrays(library="np")


#_____________________________________________________________________________
def _finite_row(values, *, label):
    """Return values if all finite; else warn and return None (keep current param)."""
    if fit_validity.all_finite(values):
        return list(values)
    runner.warn(f"{label}: non-finite fit value {values}; not written (current value kept)")
    return None


def _branch_value(tree, name, row_index, value_index):
    value = tree[name][row_index]
    if value_index is None:
        return value
    return value[value_index]


def _require_fit_fields(tree, prefixes, *, label):
    missing = [
        f"{prefix}_fit_{suffix}"
        for prefix in prefixes
        for suffix in ("status", "valid", "method")
        if f"{prefix}_fit_{suffix}" not in tree
    ]
    if missing:
        runner.fail(
            f"{label}: fit-quality branches {missing} are missing; "
            "rerun the corresponding Hodo calibration executable"
        )


def _fit_is_usable(tree, prefix, row_index, value_index, *, label, branch_names=None):
    """Require the fit-quality branches added by the current Hodo calibration binaries."""
    names = branch_names or [f"{prefix}_fit_{suffix}" for suffix in ("status", "valid", "method")]
    missing = [name for name in names if name not in tree]
    if missing:
        runner.fail(
            f"{label}: fit-quality branches {missing} are missing; "
            "rerun the corresponding Hodo calibration executable"
        )
    status_name, valid_name, method_name = names
    valid = bool(_branch_value(tree, valid_name, row_index, value_index))
    method_value = _branch_value(tree, method_name, row_index, value_index)
    method = method_value.decode() if isinstance(method_value, bytes) else str(method_value)
    status = int(_branch_value(tree, status_name, row_index, value_index))
    if valid and method and method != "invalid":
        return True
    runner.warn(
        f"{label}: fit unusable (status={status}, method={method}); "
        "not written (current value kept)"
    )
    return False


def ud_counter_rows(tree, cid, n_ud, adc_kind="mip"):
    """Counters with U/D (or single) ADC and TDC per segment: ADC (pedestal, gain), TDC offset."""
    data = {}
    for i, ch in enumerate(tree["ch"]):
        for ud in range(n_ud):
            if "adc_p0_val" in tree:
                label = f"cid{cid} seg{int(ch)} ud{ud} ADC"
                if adc_kind == "tot":
                    adc_usable = _fit_is_usable(tree, "tot", i, ud, label=label)
                else:
                    adc_usable = (
                        _fit_is_usable(tree, "adc_ped", i, ud, label=label)
                        and _fit_is_usable(tree, "adc_mip", i, ud, label=label)
                    )
                if adc_usable:
                    row = _finite_row(
                        [tree["adc_p0_val"][i][ud], tree["adc_p1_val"][i][ud]],
                        label=label,
                    )
                    if row is not None:
                        data[f"{cid}-0-{ch:.0f}-0-{ud}"] = row
            if "tdc_p0_val" in tree:
                label = f"cid{cid} seg{int(ch)} ud{ud} TDC"
                if _fit_is_usable(tree, "tdc", i, ud, label=label):
                    row = _finite_row(
                        [tree["tdc_p0_val"][i][ud], TDC_GAIN],
                        label=label,
                    )
                    if row is not None:
                        data[f"{cid}-0-{ch:.0f}-1-{ud}"] = row
    return data


def bac_rows(tree, cid):
    """BAC: ADC pedestal per segment (gain kept); one TDC offset on the sum channel (seg 4).
    Rows with adc_ok / tdc_ok == 0 (not determined, see params.h bac_*) are not written; the current value stays."""
    if "adc_ok" not in tree or "tdc_ok" not in tree:
        runner.fail("BAC fit output has no adc_ok/tdc_ok (made by an older BAC_HDPRM); run BAC_HDPRM again")
    data = {}
    if tree["tdc_ok"][0]:
        if _fit_is_usable(tree, "tdc", 0, 0, label="BAC TDC"):
            row = _finite_row([tree["tdc_p0_val"][0][0], TDC_GAIN], label="BAC TDC")
            if row is not None:
                data[f"{cid}-0-4-1-0"] = row
    else:
        runner.warn(f"BAC TDC: offset {tree['tdc_p0_val'][0][0]:.0f} outside the window, not written (current value kept)")
    for i, ch in enumerate(tree["ch"]):
        if tree["adc_ok"][i]:
            label = f"BAC seg{ch:.0f} ADC"
            if _fit_is_usable(tree, "adc_ped", i, 0, label=label):
                row = _finite_row([tree["adc_p0_val"][i][0]], label=label)
                if row is not None:
                    data[f"{cid}-0-{ch:.0f}-0-0"] = row
        else:
            runner.warn(f"BAC seg{ch:.0f} ADC: pedestal not determined (signal-dominated), not written (current value kept)")
    return data


def kvc_rows(tree, cid, n_ud):
    """KVC: ADC pedestal of a-d and SUM (one-photon gain kept); TDC offset on SUM (UorD 4).
    A TDC row with tdc_ok == 0 (too few hits around the peak) is not written; the current value stays."""
    if "tdc_ok" not in tree:
        runner.fail("KVC fit output has no tdc_ok (made by an older KVC_HDPRM); run KVC_HDPRM again")
    data = {}
    for i, ch in enumerate(tree["ch"]):
        if tree["tdc_ok"][i]:
            label = f"KVC seg{ch:.0f}S TDC"
            if _fit_is_usable(tree, "tdc", i, 0, label=label):
                row = _finite_row([tree["tdc_p0_val"][i][0], TDC_GAIN], label=label)
                if row is not None:
                    data[f"{cid}-0-{ch:.0f}-1-4"] = row
        else:
            runner.warn(f"KVC seg{ch:.0f}S TDC: too few hits around the peak, not written (current value kept)")
        for ud in range(n_ud):
            if "adc_p0_val" in tree:
                label = f"KVC seg{ch:.0f} ud{ud} ADC"
                if _fit_is_usable(tree, "adc_p0", i, ud, label=label):
                    row = _finite_row([tree["adc_p0_val"][i][ud]], label=label)
                    if row is not None:
                        data[f"{cid}-0-{ch:.0f}-0-{ud}"] = row
    return data


def htof_tdc_rows(tree, cid):
    """HTOF TDC offsets of U, D, S (ADC is make_htofprm_dictdata)."""
    data = {}
    for i, ch in enumerate(tree["ch"]):
        for ud in range(3):
            if "tdc_p0_val" in tree:
                label = f"HTOF seg{int(ch)} ud{ud} TDC"
                if _fit_is_usable(tree, "tdc", i, ud, label=label):
                    row = _finite_row([tree["tdc_p0_val"][i][ud], TDC_GAIN], label=label)
                    if row is not None:
                        data[f"{cid}-0-{ch:.0f}-1-{ud}"] = row
    return data


def make_dictdata(root_file_path, det, run=None):
    """HDPRM rows of one detector's *_HDPRM fit output."""
    tree, cid = read_tree(root_file_path), config.DETECTOR_I_LIST[det]
    required_fits = {
        "HTOF": ("tdc",),
        "BAC": ("tdc", "adc_ped"),
        "KVC": ("tdc", "adc_p0"),
        "BHT": ("tdc", "tot"),
        "SFV": ("tdc",),
    }.get(det, ("tdc", "adc_ped", "adc_mip"))
    _require_fit_fields(tree, required_fits, label=f"{det} HDPRM")
    if det == "HTOF":
        data = htof_tdc_rows(tree, cid)
    elif det == "BAC":
        data = bac_rows(tree, cid)
    elif det == "KVC":
        data = kvc_rows(tree, cid, config.DETECTOR_N_UD_LIST[det])
    else:
        adc_kind = "tot" if det == "BHT" else "mip"
        data = ud_counter_rows(tree, cid, config.DETECTOR_N_UD_LIST[det], adc_kind=adc_kind)
    return _drop_locked(run, det, data, mode="hdprm")


def make_t0_offset_dictdata(root_file_path):
    """BH2 time-zero offset per segment (UorD 2) from the T0_Offset output."""
    tree, cid = read_tree(root_file_path), config.DETECTOR_I_LIST["BH2"]
    root_file = uproot.open(root_file_path)
    if "t0_fit" not in root_file:
        runner.fail("BH2 T0 output has no t0_fit status tree; rerun T0_Offset")
    fit_tree = root_file["t0_fit"].arrays(library="np")
    if "ch" not in fit_tree:
        runner.fail("BH2 T0 fit-status tree has no ch branch; rerun T0_Offset")
    fit_rows = {int(ch): i for i, ch in enumerate(fit_tree["ch"])}
    data = {}
    for i, ch in enumerate(tree["ch"]):
        fit_row = fit_rows.get(int(ch))
        if fit_row is None:
            runner.fail(f"BH2 seg{int(ch)} has no T0 fit-status row; rerun T0_Offset")
        if not _fit_is_usable(
                fit_tree, "t0", fit_row, None,
                label=f"BH2 seg{int(ch)} T0",
                branch_names=["status", "valid", "method"]):
            continue
        row = _finite_row(
            [tree["offset_p0_val"][i][0], 1.0],
            label=f"BH2 seg{int(ch)} T0 offset",
        )
        if row is not None:
            data[f"{cid}-0-{ch:.0f}-1-2"] = row
    return data


# -- prepare HTOF ADC data (param_type "htofprm") -----------------------------
def make_htofprm_dictdata(root_file_path):
    """HTOF ADC pedestal / MIP (HDPRM "gain" column) from the HTOF_Calib ADC output
    (run{run}_HTOF_ADC_Pi.root; MIP from dE/dx pion-tagged DstTPCHelixHTOF tracks).
    UorD: 0=U, 1=D, 2=S."""
    detector_id = config.DETECTOR_I_LIST["HTOF"]
    file = uproot.open(root_file_path)
    tree = file["tree"].arrays(library="np")

    _require_fit_fields(tree, ("adc_ped", "adc_mip"), label="HTOF ADC")
    data = dict()
    for i in range(len(tree["ch"])):
        ch = tree["ch"][i]
        for UorD in range(3):
            label = f"HTOF seg{int(ch)} ud{UorD} ADC"
            if (_fit_is_usable(tree, "adc_ped", i, UorD, label=label)
                    and _fit_is_usable(tree, "adc_mip", i, UorD, label=label)):
                row = _finite_row(
                    [tree["adc_p0_val"][i][UorD], tree["adc_p1_val"][i][UorD]],
                    label=label,
                )
                if row is not None:
                    key = f"{detector_id}-0-{ch:.0f}-0-{UorD:.0f}"
                    data[key] = row
    return data
# ---------------------------------------------------------------------------

# -- HTOF fit-quality reports -------------------------------------------------
_HTOF_SIDE = "UDS"

def htof_adc_quality_report(root_file_path):
    """Channels of the HTOF_Calib ADC output whose adc_flag has bits other than 32 set.
    Returns a list of (seg, side, flag, mip). Flag bits (ana_helper::htof_adc_fit):
    1 raw seed unusable, 2 low statistics, 4 ndf < 10, 8 at fit-range edge / no fit,
    16 candidate two-peak structure; htof_adc_fit_weak (beam-window weak side): 32 two-component
    refit used (informational only, not reported alone), 64 hump not separated (MIP uncertain);
    128 MIP range / model from params.h htof_adc_fit_hint (informational only, not reported alone)."""
    tree = uproot.open(root_file_path)["tree"].arrays(library="np")
    out = []
    if "adc_flag" not in tree:
        return out
    for i in range(len(tree["ch"])):
        for s in range(3):
            flag = int(tree["adc_flag"][i][s])
            if flag & ~(32 | 128):
                out.append((int(tree["ch"][i]), _HTOF_SIDE[s], flag, float(tree["adc_p1_val"][i][s])))
    return out

def htof_phc_quality_report(root_file_path, min_offset_n=100):
    """Channels of the HTOF_Calib PHC output whose p0/p1 sit at the fit limits and segments
    whose absolute TOF offset was not applied (tof_offset_n < min_offset_n). Limits: walk from
    pion-bit-only tracks (tof_sub = 1, fit_walk_shared in HTOF_Calib.cpp): p0 at 0, p1 at the
    scan ends -3 or 0.25; two-stage fit (tof_sub = 0, ana_helper::htof_phc_fit): p0 at 0.001
    or 15, p1 at -5.
    Returns (limit_list[(seg, side, p0, p1)], offset_list[(seg, n)])."""
    tree = uproot.open(root_file_path)["tree"].arrays(library="np")
    lim, off = [], []
    for i in range(len(tree["ch"])):
        ch = int(tree["ch"][i])
        for s in range(2):
            p0, p1 = float(tree["p0_val"][i][s]), float(tree["p1_val"][i][s])
            if int(tree["tof_sub"][i]) == 1:
                at_limit = p0 <= 0.0 or p1 < -2.995 or p1 > 0.245
            else:
                at_limit = p0 < 0.0011 or p0 > 14.99 or p1 < -4.999
            if at_limit:
                lim.append((ch, _HTOF_SIDE[s], p0, p1))
        if "tof_offset_n" in tree and int(tree["tof_offset_n"][i]) < min_offset_n:
            off.append((ch, int(tree["tof_offset_n"][i])))
    return lim, off
# ---------------------------------------------------------------------------
