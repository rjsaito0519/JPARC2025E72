"""HDPHC rows {"CId-PlId-SegId-UorD": [type, npar, p0, p1, p2]} from a *_PHC fit output.

The detector is given by the caller (update_param.py); its CId comes from lib/config.py.
Which segments are usable:
  - hand range (phc_conf.py limits_dict "RUN_SUFFIX_det": [first, last]) if the caller gives one;
  - otherwise automatic: a segment is not used if a side has fewer than MIN_FIT_ENTRIES entries in
    the fit or no fit at all (ndf 0). chi2/ndf is not used as a cut: it grows with the statistics
    because the PHC function does not follow the profile exactly.
  - In both modes, fit_status (when present), finite parameters, and chi2_ndf >= 0 are required.
    Failed / non-finite fits keep the previous parameter (not written) and do not enter the mean.
Segments not used (out of range / low stats) get the mean parameters of the used segments of the
same side when at least one valid fit exists; otherwise nothing is written for them.
"""
import statistics
import sys
from pathlib import Path

import numpy as np
import uproot

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import config, fit_validity  # noqa: E402

PHC_TYPE, PHC_NPAR = 1, 3
MIN_FIT_ENTRIES = 300


#_____________________________________________________________________________
def _side_arrays(tree, index):
    n = np.asarray(tree["n_fit"][index], dtype=float)
    chi = np.asarray(tree["chi2_ndf"][index], dtype=float)
    n_ud = len(tree["p0_val"][index])
    statuses = None
    if "fit_status" in tree:
        statuses = np.asarray(tree["fit_status"][index])
    pars_by_ud = []
    for ud in range(n_ud):
        pars_by_ud.append([
            float(tree[f"p{j}_val"][index][ud]) for j in range(3)
        ])
    return n, chi, statuses, pars_by_ud


#_____________________________________________________________________________
def side_fit_usable(n_entries, chi2_ndf, status, pars) -> bool:
    """True when one U/D side may be written or contribute to the mean."""
    if not fit_validity.is_finite_number(n_entries):
        return False
    if float(n_entries) < MIN_FIT_ENTRIES:
        return False
    if not fit_validity.is_finite_number(chi2_ndf) or float(chi2_ndf) < 0:
        return False
    if not fit_validity.phc_fit_status_ok(status):
        return False
    if not fit_validity.all_finite(pars):
        return False
    return True


#_____________________________________________________________________________
def segment_fit_usable(tree, index) -> bool:
    """True when every side of the segment has a usable fit."""
    n, chi, statuses, pars_by_ud = _side_arrays(tree, index)
    for ud, pars in enumerate(pars_by_ud):
        status = None if statuses is None else statuses[ud]
        n_ud = float(n[ud]) if ud < len(n) else float("nan")
        chi_ud = float(chi[ud]) if ud < len(chi) else float("nan")
        if not side_fit_usable(n_ud, chi_ud, status, pars):
            return False
    return True


#_____________________________________________________________________________
def auto_bad_segments(tree):
    """Segments whose fit is unusable on any side (stats / status / finite)."""
    if "n_fit" not in tree:
        raise RuntimeError(
            "fit output without n_fit/chi2_ndf: fit again with the current *_PHC binary"
        )
    bad = []
    for i, ch in enumerate(tree["ch"]):
        if not segment_fit_usable(tree, i):
            bad.append(int(ch))
    return bad


#_____________________________________________________________________________
def make_dictdata(root_file_path, det, good_ch_range=None):
    """(rows, not-used segments). good_ch_range = hand range, or None for automatic.

    Hand range selects which segments are intended; fit validity still applies
    inside that range (failed fits are not written and do not seed the mean).
    """
    tree = uproot.open(root_file_path)["tree"].arrays(library="np")
    cid = config.DETECTOR_I_LIST[det]
    if "n_fit" not in tree:
        raise RuntimeError(
            "fit output without n_fit/chi2_ndf: fit again with the current *_PHC binary"
        )

    out_of_range = set()
    if good_ch_range is not None:
        out_of_range = {
            int(ch) for ch in tree["ch"]
            if not good_ch_range[0] <= ch <= good_ch_range[1]
        }

    data, good = {}, {0: [[], [], []], 1: [[], [], []]}
    unused = set(out_of_range)
    for i, ch in enumerate(tree["ch"]):
        ch_i = int(ch)
        if ch_i in out_of_range:
            continue
        if not segment_fit_usable(tree, i):
            unused.add(ch_i)
            continue
        for ud in range(len(tree["p0_val"][i])):
            pars = [float(tree[f"p{j}_val"][i][ud]) for j in range(3)]
            data[f"{cid}-0-{ch_i}-{ud}"] = [PHC_TYPE, PHC_NPAR] + pars
            if ud in good:
                for j in range(3):
                    good[ud][j].append(pars[j])

    for ch in sorted(unused):
        for ud, lists in good.items():
            if not lists[0]:
                continue
            if not fit_validity.all_finite([statistics.mean(v) for v in lists]):
                continue
            data[f"{cid}-0-{ch}-{ud}"] = (
                [PHC_TYPE, PHC_NPAR] + [statistics.mean(v) for v in lists]
            )
    return data, sorted(unused)
