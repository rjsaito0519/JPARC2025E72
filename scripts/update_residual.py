import uproot
import sys

detector_id_list = {
    "BLC1": None,
    "BLC2": 1001,
}

# -- prepare HDPRM data  -----------------------------------
def make_dictdata(root_file_path):

    file = uproot.open(root_file_path)
    tree = file["tree"].arrays(library="np")
    detector_id = -1
    if "BLC1" in root_file_path:
        detector_id = detector_id_list["BLC1"]
    elif "BLC2" in root_file_path:
        detector_id = detector_id_list["BLC2"]
        # mu_value = -1.0 * (Mean - ref_off)

    if detector_id == -1:
        print("something wrong")
        sys.exit()

    data = dict()
    ref_offset = tree["tdc_p0_val"][0][0]
    for i in range(len(tree["tdc_p0_val"])):
        # id
        id = detector_id + tree["ch"][i]
        for WireId in range(32):
            # -- TDC BLCa -----
            key = f"{detector_id}"
            data[key] = [ -1.0*(tree["tdc_p0_val"][i][0] - ref_offset) ]
            # -- TDC BLCb -----
            key = f"{detector_id+8}"
            data[key] = [ -1.0*(tree["tdc_p0_val"][i][0] - ref_offset) ]

    return data
# ---------------------------------------------------------------------------

# -- write HDPRM file  -----------------------------------
def update_file(target_file, data):
    buf = []
    n_update = 0
    with open(target_file) as f:
        for line in f:
            s_list = line.split()
            if len(s_list) != 0 and s_list[0] in data.keys():
                s_list[12] = data[s_list[0]][0]
                n_update += 1               
            buf.append(s_list)

    with open(target_file, mode='w') as f:
        for l in buf:
            f.write('\t'.join(str(item) for item in l))
            f.write("\n")

    return len(data) == n_update
# ---------------------------------------------------------------------------
