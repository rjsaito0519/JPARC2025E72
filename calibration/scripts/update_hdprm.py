from re import U
import uproot
import sys
import numpy as np

detector_id_list = {
    "BHT":  1,
    "T0":   2,
    "BH2":  3,
    "BAC":  4,
    "HTOF": 5,
    "KVC":  6,
    "T1":   7, 
    "CVC":  8,
    "SAC3": 9,
    "SFV": 10,
}

detector_n_ud_list = {
    "BHT":  2,
    "T0":   2,
    "BH2":  2,
    "BAC":  1,
    "HTOF": -1,
    "KVC":  5,
    "T1":   1, 
    "CVC":  2,
    "SAC3": 1,
    "SFV":  1,
}

# -- prepare HDPRM data  -----------------------------------
def make_dictdata(root_file_path, good_ch_range = [-np.inf, np.inf], is_t0_offset = False):

    file = uproot.open(root_file_path)
    tree = file["tree"].arrays(library="np")
    detector_id = -1
    n_ud = -1
    for key, det_id in detector_id_list.items():
        if key in root_file_path:
            detector_id = det_id
            n_ud = detector_n_ud_list[key]
            break

    if detector_id == -1:
        print("something wrong")
        sys.exit()

    data = dict()
    if is_t0_offset:
        detector_id = detector_id_list["BH2"]
        for i in range(len(tree["offset_p0_val"])):
            # CId - PlId - SegId - AorT(0:adc, 1:tdc) - UorD(0:u, 1:d)
            ch = tree["ch"][i]
            key = f"{detector_id}-0-{ch:.0f}-1-2"
            data[key] = [ tree["offset_p0_val"][i][0], 1.0 ]
    else:
        for i in range(len(tree["ch"])):
            # CId - PlId - SegId - AorT(0:adc, 1:tdc) - UorD(0:u, 1:d)
            ch = tree["ch"][i]
            if n_ud != -1:
                if detector_id == detector_id_list["BAC"]:
                    if "tdc_p0_val" in tree.keys():
                        key = f"{detector_id}-0-4-1-0"
                        data[key] = [ tree["tdc_p0_val"][0][0], -0.0009765625 ]
                elif detector_id == detector_id_list["KVC"]:
                    if "tdc_p0_val" in tree.keys():
                        key = f"{detector_id}-0-{ch:.0f}-1-4"
                        data[key] = [ tree["tdc_p0_val"][i][0], -0.0009765625 ]

                for UorD in range(n_ud):
                    if detector_id in [detector_id_list["BAC"], detector_id_list["KVC"]] :
                        if "adc_p0_val" in tree.keys():
                            key = f"{detector_id}-0-{ch:.0f}-0-{UorD:.0f}"
                            data[key] = [ tree["adc_p0_val"][i][UorD] ]
                    else:
                        # -- ADC -----
                        if "adc_p0_val" in tree.keys():
                            key = f"{detector_id}-0-{ch:.0f}-0-{UorD:.0f}"
                            data[key] = [ tree["adc_p0_val"][i][UorD], tree["adc_p1_val"][i][UorD] ]

                        # -- TDC -----
                        if "tdc_p0_val" in tree.keys():
                            key = f"{detector_id}-0-{ch:.0f}-1-{UorD:.0f}"
                            data[key] = [ tree["tdc_p0_val"][i][UorD], -0.0009765625 ]
            else:
                if detector_id == detector_id_list["HTOF"]:
                    # HTOF HDPRM (this function/param_type) is TDC-only (All-event input,
                    # 3 UorD entries: U, D, S). HTOF ADC is handled separately by
                    # make_htofprm_dictdata() below (param_type "htofprm"), which uses the
                    # dE/dx-tagged DstTPCHelixHTOF output instead of the Hodo-level tree.
                    for UorD in range(3):
                        if "tdc_p0_val" in tree.keys():
                            key = f"{detector_id}-0-{ch:.0f}-1-{UorD:.0f}"
                            data[key] = [ tree["tdc_p0_val"][i][UorD], -0.0009765625 ]

    return data
# ---------------------------------------------------------------------------

# -- prepare HTOF ADC data (param_type "htofprm") -----------------------------
def make_htofprm_dictdata(root_file_path):
    """HTOF ADC pedestal / MIP (HDPRM "gain" column) from the HTOF_Calib ADC output
    (run{run}_HTOF_ADC_Pi.root; MIP from dE/dx pion-tagged DstTPCHelixHTOF tracks).
    UorD: 0=U, 1=D, 2=S (see update_hdprm.py detector_n_ud_list / AorT=0 rows)."""
    detector_id = detector_id_list["HTOF"]
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
    """Channels of the HTOF_Calib PHC output whose p0/p1 sit at the fit limits
    (p0 at 0.001 or 15, p1 at -5; see ana_helper::htof_phc_fit) and segments whose
    absolute TOF offset was not applied (tof_offset_n < min_offset_n).
    Returns (limit_list[(seg, side, p0, p1)], offset_list[(seg, n)])."""
    tree = uproot.open(root_file_path)["tree"].arrays(library="np")
    lim, off = [], []
    for i in range(len(tree["ch"])):
        ch = int(tree["ch"][i])
        for s in range(2):
            p0, p1 = float(tree["p0_val"][i][s]), float(tree["p1_val"][i][s])
            if p0 < 0.0011 or p0 > 14.99 or p1 < -4.999:
                lim.append((ch, _HTOF_SIDE[s], p0, p1))
        if "tof_offset_n" in tree and int(tree["tof_offset_n"][i]) < min_offset_n:
            off.append((ch, int(tree["tof_offset_n"][i])))
    return lim, off
# ---------------------------------------------------------------------------

# -- write HDPRM file  -----------------------------------
def update_file(target_file, data):
    buf = []
    n_update = 0
    with open(target_file) as f:
        for line in f:
            s_list = line.split()
            if s_list[0][0] == "#":
                buf.append(s_list)
                continue
            # key structure
            # CId - PlId - SegId - AorT(0:adc, 1:tdc) - UorD(0:u, 1:d)
            key_length = 5
            key = s_list[0]
            for i in range(1, key_length):
                key += "-"+s_list[i]
            if key in data.keys():
                for i in range(len(data[key])):
                    s_list[i+key_length] = data[key][i]
                n_update += 1               
            buf.append(s_list)

    with open(target_file, mode='w') as f:
        for l in buf:
            f.write('\t'.join(str(item) for item in l))
            f.write("\n")

    return len(data) == n_update
# ---------------------------------------------------------------------------
