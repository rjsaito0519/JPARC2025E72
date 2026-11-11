import math
import sys
from pathlib import Path

import numpy as np
import uproot
from termcolor import colored

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import fit_validity  # noqa: E402

detector_id_list = {
    "BLC1": 1,
    "BLC2": 1001,
}


class ResidualFitInvalid(ValueError):
    """Residual means/sigmas are not all finite or the reference plane is bad."""


# -- prepare residual correction data  -----------------------------------
def make_dictdata(root_file_path):
    try:
        file = uproot.open(root_file_path)
        if "tree" not in file:
            print(colored(f"  [Warning] 'tree' not found in {root_file_path}. Skipping.", "yellow"))
            return {}
        tree = file["tree"].arrays(library="np")
    except ResidualFitInvalid:
        raise
    except Exception as e:
        print(colored(f"  [Error] Opening {root_file_path}: {e}", "red"))
        return {}
    detector_id = -1

    if "BLC1" in root_file_path:
        detector_id = detector_id_list["BLC1"]
    elif "BLC2" in root_file_path:
        detector_id = detector_id_list["BLC2"]

    if detector_id == -1:
        print(colored(f"Error: Could not determine detector from path {root_file_path}", "red"))
        sys.exit(1)

    # Relative to first plane (U1a) to avoid overall detector shift/rotation.
    if len(tree["residual_p0_val"]) < 1:
        raise ResidualFitInvalid(f"{root_file_path}: empty residual_p0_val")
    ref_mu = float(tree["residual_p0_val"][0][0])
    if not fit_validity.is_finite_number(ref_mu):
        raise ResidualFitInvalid(
            f"{root_file_path}: reference plane mean is not finite ({ref_mu!r})"
        )
    print(colored(f"  [INFO] Reference Mean (Plane 0): {ref_mu:.4f} mm", "cyan"))

    data = dict()
    has_sigma = "residual_p1_val" in tree
    invalid = []

    for i in range(len(tree["residual_p0_val"])):
        mu_a = float(tree["residual_p0_val"][i][0])
        sig_a = float(tree["residual_p1_val"][i][0]) if has_sigma else 0.0
        mu_b = float(tree["residual_p0_val"][i][1])
        sig_b = float(tree["residual_p1_val"][i][1]) if has_sigma else 0.0
        for label, mu, sig in (
            (f"{detector_id + i}", mu_a, sig_a),
            (f"{detector_id + i + 8}", mu_b, sig_b),
        ):
            if not fit_validity.all_finite([mu, sig]):
                invalid.append(label)
                continue
            corr = -1.0 * (mu - ref_mu)
            if not math.isfinite(corr):
                invalid.append(label)
                continue
            data[label] = [corr, sig]

    if invalid:
        raise ResidualFitInvalid(
            f"{root_file_path}: non-finite residual mean/sigma on plane(s) {invalid}"
        )
    return data
# ---------------------------------------------------------------------------
