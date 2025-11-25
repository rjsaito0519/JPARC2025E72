import uproot
import sys
import numpy as np
import statistics

detector_id_list = {
    "BHT":  1,
    "T0":   2,
    "BH2":  3,
    "HTOF": 5,
    "CVC":  8
}

# -- prepare HDPRM data  -----------------------------------
def make_dictdata(root_file_path, good_ch_range = [-np.inf, np.inf]):

    file = uproot.open(root_file_path)
    tree = file["tree"].arrays(library="np")
    detector_id = -1
    if "BHT" in root_file_path:
        detector_id   = detector_id_list["BHT"]
    elif "T0" in root_file_path:
        detector_id = detector_id_list["T0"]
    elif "BH2" in root_file_path:
        detector_id = detector_id_list["BH2"]
    elif "HTOF" in root_file_path:
        detector_id = detector_id_list["HTOF"]
    elif "CVC" in root_file_path:
        detector_id = detector_id_list["CVC"]

    if detector_id == -1:
        print("something wrong")
        sys.exit()

    data = dict()
    buf = [
        [[], [], []], # up
        [[], [], []]  # down
    ]
    bad_ch_index = []
    for i in range(len(tree["p0_val"])):
        ch = tree["ch"][i]
        if good_ch_range[0] <= ch and ch <= good_ch_range[1]:
            for UorD in [0, 1]:
                # CId - PlId - SegId - UorD(0:u, 1:d)
                key = f"{detector_id}-0-{ch:.0f}-{UorD:.0f}"
                data[key] = [ 1, 3, tree["p0_val"][i][UorD], tree["p1_val"][i][UorD], tree["p2_val"][i][UorD] ]
                for j in range(3):
                    buf[UorD][j].append(tree[f"p{j}_val"][i][UorD])
        else:
            bad_ch_index.append(i)

    for i in bad_ch_index:
        ch = tree["ch"][i]
        for UorD in [0, 1]:
            # CId - PlId - SegId - UorD(0:u, 1:d)
            key = f"{detector_id}-0-{ch:.0f}-{UorD:.0f}"
            data[key] = [ 1, 3, statistics.mean(buf[UorD][0]), statistics.mean(buf[UorD][1]), statistics.mean(buf[UorD][2]) ]

    return data
# ---------------------------------------------------------------------------


# -- write HDPHC file  -----------------------------------
def update_file(target_file, data):
    buf = []
    n_update = 0
    with open(target_file) as f:
        for line in f:
            s_list = line.split()
            # key structure
            # CId - PlId - SegId - UorD(0:u, 1:d)
            key_length = 4
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
