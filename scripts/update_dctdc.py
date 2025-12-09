import uproot
import sys

detector_id_list = {
    "BLC1": 101, # (1a: 101, 1b: 102)
    "BLC2": 103, # (1a: 103, 1b: 104)
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

    if detector_id == -1:
        print("something wrong")
        sys.exit()

    data = dict()
    for i in range(len(tree["tdc_p0_val"])):
        # Cid - Plid - WireId
        ch = tree["ch"][i]
        for WireId in range(32):
            # -- TDC BLCa -----
            key = f"{detector_id}-{ch:.0f}-{WireId:.0f}"
            data[key] = [ tree["tdc_p0_val"][i][0], -0.8333333333 ]
            # -- TDC BLCb -----
            key = f"{detector_id+1}-{ch:.0f}-{WireId:.0f}"
            data[key] = [ tree["tdc_p0_val"][i][1], -0.8333333333 ]

    return data
# ---------------------------------------------------------------------------

# -- write HDPRM file  -----------------------------------
def update_file(target_file, data):
    buf = []
    n_update = 0
    with open(target_file) as f:
        for line in f:
            s_list = line.split()
            if s_list[0] == "#":
                buf.append(s_list)
                continue
            # key structure
            # Cid - Plid - WireId
            key_length = 3
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
