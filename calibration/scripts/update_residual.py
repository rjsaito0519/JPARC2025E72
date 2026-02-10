import uproot
import sys
import numpy as np
from termcolor import colored

detector_id_list = {
    "BLC1": 1,
    "BLC2": 1001,
}

# -- prepare HDPRM data  -----------------------------------
def make_dictdata(root_file_path):
    try:
        file = uproot.open(root_file_path)
        if "tree" not in file:
            print(colored(f"  [Warning] 'tree' not found in {root_file_path}. Skipping.", "yellow"))
            return {}
        tree = file["tree"].arrays(library="np")
    except Exception as e:
        print(colored(f"  [Error] Opening {root_file_path}: {e}", "red"))
        return {}

    detector_id = -1
    if "BLC1" in root_file_path:
        detector_id = detector_id_list["BLC1"]
    elif "BLC2" in root_file_path:
        detector_id = detector_id_list["BLC2"]
    elif "DC" in root_file_path:
        print(colored("  [Warning] Identifying combined DC file as BLC1 (ID 1)", "yellow"))
        detector_id = detector_id_list["BLC1"]

    if detector_id == -1:
        print(colored(f"  [Error] Could not determine detector from path {root_file_path}", "red"))
        sys.exit(1)

    # --- Macro-like relative correction logic ---
    # We reference the first plane (U1) to avoid overall detector shift/rotation.
    # mu_value = -1.0 * (current_mu - ref_mu)
    
    ref_mu = 0.0
    if len(tree["residual_p0_val"]) > 0:
        # Reference is side 0 (a-side), plane 0 (U1)
        ref_mu = tree["residual_p0_val"][0][0]
        print(colored(f"  [INFO] Reference Mean (Plane 0): {ref_mu:.4f} mm", "cyan"))

    data = dict()
    has_sigma = "residual_p1_val" in tree
    
    for i in range(len(tree["residual_p0_val"])):
        # BLCa (Planes 1-8 or 1001-1008)
        mu_a = tree["residual_p0_val"][i][0]
        sig_a = tree["residual_p1_val"][i][0] if has_sigma else 0.0
        key_a = f"{detector_id + i}"
        data[key_a] = [ -1.0 * (mu_a - ref_mu), sig_a ]
        
        # BLCb (Planes 9-16 or 1009-1016)
        mu_b = tree["residual_p0_val"][i][1]
        sig_b = tree["residual_p1_val"][i][1] if has_sigma else 0.0
        key_b = f"{detector_id + i + 8}"
        data[key_b] = [ -1.0 * (mu_b - ref_mu), sig_b ]

    return data
# ---------------------------------------------------------------------------
