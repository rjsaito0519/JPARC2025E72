"""HDPRM rows {"CId-PlId-SegId-AorT-UorD": [values]} from the Hodo fit outputs (results/root).

AorT: 0 = ADC (values pedestal, gain), 1 = TDC (values offset, TDC_GAIN).
The detector is given by the caller (update_param.py); CId and the number of UorD rows come from
lib/config.py.
"""
import sys
from pathlib import Path

import uproot

sys.path.append(str(Path(__file__).resolve().parents[2]))
from lib import config, runner  # noqa: E402

TDC_GAIN = -0.0009765625  # ns per TDC channel (-1/1024), written as the TDC gain column


def read_tree(root_file_path):
    return uproot.open(root_file_path)["tree"].arrays(library="np")


def ud_counter_rows(tree, cid, n_ud):
    """Counters with U/D (or single) ADC and TDC per segment: ADC (pedestal, gain), TDC offset."""
    data = {}
    for i, ch in enumerate(tree["ch"]):
        for ud in range(n_ud):
            if "adc_p0_val" in tree:
                data[f"{cid}-0-{ch:.0f}-0-{ud}"] = [tree["adc_p0_val"][i][ud], tree["adc_p1_val"][i][ud]]
            if "tdc_p0_val" in tree:
                data[f"{cid}-0-{ch:.0f}-1-{ud}"] = [tree["tdc_p0_val"][i][ud], TDC_GAIN]
    return data


def bac_rows(tree, cid):
    """BAC: ADC pedestal per segment (gain kept); one TDC offset on the sum channel (seg 4).
    Rows with adc_ok / tdc_ok == 0 (not determined, see params.h bac_*) are not written; the current value stays."""
    if "adc_ok" not in tree or "tdc_ok" not in tree:
        runner.fail("BAC fit output has no adc_ok/tdc_ok (made by an older BAC_HDPRM); run BAC_HDPRM again")
    data = {}
    if tree["tdc_ok"][0]:
        data[f"{cid}-0-4-1-0"] = [tree["tdc_p0_val"][0][0], TDC_GAIN]
    else:
        runner.warn(f"BAC TDC: offset {tree['tdc_p0_val'][0][0]:.0f} outside the window, not written (current value kept)")
    for i, ch in enumerate(tree["ch"]):
        if tree["adc_ok"][i]:
            data[f"{cid}-0-{ch:.0f}-0-0"] = [tree["adc_p0_val"][i][0]]
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
            data[f"{cid}-0-{ch:.0f}-1-4"] = [tree["tdc_p0_val"][i][0], TDC_GAIN]
        else:
            runner.warn(f"KVC seg{ch:.0f}S TDC: too few hits around the peak, not written (current value kept)")
        for ud in range(n_ud):
            if "adc_p0_val" in tree:
                data[f"{cid}-0-{ch:.0f}-0-{ud}"] = [tree["adc_p0_val"][i][ud]]
    return data


def htof_tdc_rows(tree, cid):
    """HTOF TDC offsets of U, D, S (ADC is make_htofprm_dictdata)."""
    data = {}
    for i, ch in enumerate(tree["ch"]):
        for ud in range(3):
            if "tdc_p0_val" in tree:
                data[f"{cid}-0-{ch:.0f}-1-{ud}"] = [tree["tdc_p0_val"][i][ud], TDC_GAIN]
    return data


def make_dictdata(root_file_path, det):
    """HDPRM rows of one detector's *_HDPRM fit output."""
    tree, cid = read_tree(root_file_path), config.DETECTOR_I_LIST[det]
    if det == "HTOF":
        return htof_tdc_rows(tree, cid)
    if det == "BAC":
        return bac_rows(tree, cid)
    if det == "KVC":
        return kvc_rows(tree, cid, config.DETECTOR_N_UD_LIST[det])
    return ud_counter_rows(tree, cid, config.DETECTOR_N_UD_LIST[det])


def make_t0_offset_dictdata(root_file_path):
    """BH2 time-zero offset per segment (UorD 2) from the T0_Offset output."""
    tree, cid = read_tree(root_file_path), config.DETECTOR_I_LIST["BH2"]
    return {f"{cid}-0-{ch:.0f}-1-2": [tree["offset_p0_val"][i][0], 1.0] for i, ch in enumerate(tree["ch"])}


# -- prepare HTOF ADC data (param_type "htofprm") -----------------------------
def make_htofprm_dictdata(root_file_path):
    """HTOF ADC pedestal / MIP (HDPRM "gain" column) from the HTOF_Calib ADC output
    (run{run}_HTOF_ADC_Pi.root; MIP from dE/dx pion-tagged DstTPCHelixHTOF tracks).
    UorD: 0=U, 1=D, 2=S."""
    detector_id = config.DETECTOR_I_LIST["HTOF"]
    file = uproot.open(root_file_path)
    tree = file["tree"].arrays(library="np")

    data = dict()
    for i in range(len(tree["ch"])):
        ch = tree["ch"][i]
        for UorD in range(3):
            key = f"{detector_id}-0-{ch:.0f}-0-{UorD:.0f}"
            data[key] = [ tree["adc_p0_val"][i][UorD], tree["adc_p1_val"][i][UorD] ]
    return data
# ---------------------------------------------------------------------------

# -- HTOF fit-quality reports (warnings only; all channels are still written) --
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
