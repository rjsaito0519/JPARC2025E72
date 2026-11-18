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
def make_dictdata(root_file_path, input_path=None):
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

    # New native outputs carry primary ROOT fit status. Legacy outputs retain
    # the finite/width/bounds checks below.
    if "residual_fit_status" in tree:
        statuses = tree["residual_fit_status"]
        if len(statuses) != len(tree["residual_p0_val"]):
            raise ResidualFitInvalid(f"{root_file_path}: residual fit status length mismatch")
        for plane, sides in enumerate(statuses):
            try:
                usable = len(sides) == 2 and all(fit_validity.root_fit_status_ok(v) for v in sides)
            except TypeError:
                usable = False
            if not usable:
                raise ResidualFitInvalid(
                    f"{root_file_path}: unusable residual fit status at plane {plane}: {sides}"
                )

    # Relative to first plane (U1a) to avoid overall detector shift/rotation.
    if len(tree["residual_p0_val"]) < 1:
        raise ResidualFitInvalid(f"{root_file_path}: empty residual_p0_val")
    ref_mu = float(tree["residual_p0_val"][0][0])
    if not fit_validity.is_finite_number(ref_mu):
        raise ResidualFitInvalid(
            f"{root_file_path}: reference plane mean is not finite ({ref_mu!r})"
        )
    print(colored(f"  [INFO] Reference Mean (Plane 0): {ref_mu:.4f} mm", "cyan"))

    ranges = {}
    if input_path is not None:
        try:
            with uproot.open(input_path) as source:
                for plane in range(len(tree["residual_p0_val"])):
                    for side in ("a", "b"):
                        suffix = "K" if Path(root_file_path).stem.endswith("_K") else "Pi"
                        name = f"BLC{1 if detector_id == 1 else 2}{side}_Track_Residual_plane{plane}_{suffix}"
                        edges = source[name].axis().edges()
                        ranges[(plane, side)] = (float(edges[0]), float(edges[-1]))
        except Exception as error:
            raise ResidualFitInvalid(f"cannot verify residual histogram bounds: {input_path}: {error}") from error

    data = dict()
    has_sigma = "residual_p1_val" in tree
    invalid = []

    for i in range(len(tree["residual_p0_val"])):
        mu_a = float(tree["residual_p0_val"][i][0])
        sig_a = float(tree["residual_p1_val"][i][0]) if has_sigma else 0.0
        mu_b = float(tree["residual_p0_val"][i][1])
        sig_b = float(tree["residual_p1_val"][i][1]) if has_sigma else 0.0
        for side, label, mu, sig in (
            ("a", f"{detector_id + i}", mu_a, sig_a),
            ("b", f"{detector_id + i + 8}", mu_b, sig_b),
        ):
            if not fit_validity.all_finite([mu, sig]):
                invalid.append(label)
                continue
            bounds = ranges.get((i, side))
            if ((has_sigma and sig <= 0.) or
                    (bounds is not None and (not bounds[0] <= mu <= bounds[1]
                                             or (has_sigma and sig > bounds[1] - bounds[0])))):
                invalid.append(label)
                continue
            corr = -1.0 * (mu - ref_mu)
            if not math.isfinite(corr):
                invalid.append(label)
                continue
            data[label] = [corr, sig]

    if invalid:
        raise ResidualFitInvalid(
            f"{root_file_path}: invalid residual mean/sigma (non-finite, nonpositive width, or outside histogram bounds) on plane(s) {invalid}"
        )
    return data
# ---------------------------------------------------------------------------
