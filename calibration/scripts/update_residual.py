import uproot
import sys
import numpy as np

detector_id_list = {
    "BLC1": 1,
    "BLC2": 1001,
}

# -- prepare HDPRM data  -----------------------------------
def make_dictdata(root_file_path):
    file = uproot.open(root_file_path)
    tree = file["tree"].arrays(library="np")
    detector_id = -1
    
    # Identify detector. 
    # BLC1/2 files often named specifically, but common DC files need checking too.
    if "BLC1" in root_file_path:
        detector_id = detector_id_list["BLC1"]
    elif "BLC2" in root_file_path:
        detector_id = detector_id_list["BLC2"]
    elif "DC" in root_file_path:
        # Fallback for combined DC file: we might need to handle both trees or identify which one.
        # But usually update_param.py iterates over them if separate.
        # For combined, we default to BLC1 and assume caller handles the offset.
        print(colored("Warning: Identifying combined DC file as BLC1 (ID 1)", "yellow"))
        detector_id = detector_id_list["BLC1"]

    if detector_id == -1:
        print(colored(f"Error: Could not determine detector from path {root_file_path}", "red"))
        sys.exit(1)

    data = dict()
    # residual_p0_val is usually [plane_idx][in_or_out_idx] (0:a, 1:b)
    # Ref offset is the first plane's first value in some cases, or fixed 0.
    # In E72, values are often absolute offsets to be ADDED or replaced.
    # The C++ code writes blca_residual[i].par[1] (the Mean of the residual fit).
    # If the residual is +0.1mm, the offset should be Adjusted by -0.1mm.
    # DCGeom Ofs = OldOfs - Mean? Or purely overwrite?
    # Usually, we overwrite with the calculated center.
    
    for i in range(len(tree["residual_p0_val"])):
        # CID for BLC planes: 
        # BLC1: 1..8 (a), 9..16 (b)
        # BLC2: 1001..1008 (a), 1009..1016 (b)
        
        # BLCa
        key_a = f"{detector_id + i}"
        data[key_a] = [ tree["residual_p0_val"][i][0] ]
        
        # BLCb
        key_b = f"{detector_id + i + 8}"
        data[key_b] = [ tree["residual_p0_val"][i][1] ]

    return data
# ---------------------------------------------------------------------------
