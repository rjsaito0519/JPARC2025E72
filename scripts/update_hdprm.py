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
    "BAC":  -1,
    "HTOF": -1,
    "KVC":  -1,
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

    if detector_id == detector_id_list["HTOF"] and "adc_p0_val" in tree.keys():
        import os
        import conf
        calib_root_file_path = os.path.join(
            conf.data_dir,
            f"root/run02603_HTOF_HDPRM_Pi.root"
        )
        calib_tree = uproot.open(calib_root_file_path)["tree"].arrays(library="np")
        factor = [[], [], []]
        for UorD in range(3):
            for i in range(len(tree["ch"])):
                ch = tree["ch"][i]
                if good_ch_range[0] <= ch <= good_ch_range[1]:
                    mip = tree["adc_p1_val"][i][UorD] - tree["adc_p0_val"][i][UorD]
                    calib_mip = calib_tree["adc_p1_val"][i][UorD] - calib_tree["adc_p0_val"][i][UorD]
                    factor[UorD].append(mip/calib_mip)

        factor = [ np.mean(x) if len(x)>0 else 1 for x in factor ]
        print(factor)

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
                for UorD in range(n_ud):
                    # -- ADC -----
                    if "adc_p0_val" in tree.keys():
                        key = f"{detector_id}-0-{ch:.0f}-0-{UorD:.0f}"
                        data[key] = [ tree["adc_p0_val"][i][UorD], tree["adc_p1_val"][i][UorD] ]

                    # -- TDC -----
                    if "tdc_p0_val" in tree.keys():
                        key = f"{detector_id}-0-{ch:.0f}-1-{UorD:.0f}"
                        data[key] = [ tree["tdc_p0_val"][i][UorD], -0.0009765625 ]
            else:
                if detector_id == detector_id_list["BAC"]:
                    if "tdc_p0_val" in tree.keys():
                        key = f"{detector_id}-0-4-1-0"
                        data[key] = [ tree["tdc_p0_val"][0][0], -0.0009765625 ]
                elif detector_id == detector_id_list["KVC"]:
                    if "tdc_p0_val" in tree.keys():
                        key = f"{detector_id}-0-{ch:.0f}-1-4"
                        data[key] = [ tree["tdc_p0_val"][i][0], -0.0009765625 ]
                elif detector_id == detector_id_list["HTOF"]:
                    for UorD in range(3):
                        # -- ADC -----
                        if good_ch_range[0] <= ch <= good_ch_range[1]:
                            key = f"{detector_id}-0-{ch:.0f}-0-{UorD:.0f}"
                            data[key] = [ tree["adc_p0_val"][i][UorD], tree["adc_p1_val"][i][UorD] ]
                        else:
                            key = f"{detector_id}-0-{ch:.0f}-0-{UorD:.0f}"
                            calib_mip = calib_tree["adc_p1_val"][i][UorD] - calib_tree["adc_p0_val"][i][UorD]
                            data[key] = [ tree["adc_p0_val"][i][UorD], tree["adc_p0_val"][i][UorD] + calib_mip*factor[UorD] ]
                                
                        # -- TDC -----
                        if "tdc_p0_val" in tree.keys():
                            key = f"{detector_id}-0-{ch:.0f}-1-{UorD:.0f}"
                            data[key] = [ tree["tdc_p0_val"][i][UorD], -0.0009765625 ]

    return data
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
