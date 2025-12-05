import uproot
import sys

detector_id_list = {
    "BHT":  1,
    "T0":   2,
    "BH2":  3,
    "BH2":  4,
    "HTOF": 5,
    "KVC":  6,
    "T1":   7, 
    "CVC":  8,
    "SAC3": 9,
    "SFV": 10,
}

# -- prepare HDPRM data  -----------------------------------
def make_dictdata(root_file_path, is_t0_offset = False):

    file = uproot.open(root_file_path)
    tree = file["tree"].arrays(library="np")
    detector_id = -1
    for key, det_id in detector_id_list.items():
        if key in root_file_path:
            detector_id = det_id
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
        print("adc_p0_val" in tree.keys())
        for i in range(len(tree["ch"])):
            # CId - PlId - SegId - AorT(0:adc, 1:tdc) - UorD(0:u, 1:d)
            ch = tree["ch"][i]
            for UorD in [0, 1]:
                # -- ADC -----
                key = f"{detector_id}-0-{ch:.0f}-0-{UorD:.0f}"
                data[key] = [ tree["adc_p0_val"][i][UorD], tree["adc_p1_val"][i][UorD] ]

                # -- TDC -----
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
