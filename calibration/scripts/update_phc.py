"""HDPHC rows {"CId-PlId-SegId-UorD": [type, npar, p0, p1, p2]} from a *_PHC fit output.

The detector is given by the caller (update_param.py); its CId comes from lib/config.py.
Which segments are usable:
  - hand range (phc_conf.py limits_dict "RUN_SUFFIX_det": [first, last]) if the caller gives one;
  - otherwise automatic: a segment is not used if a side has fewer than MIN_FIT_ENTRIES entries in
    the fit or no fit at all (ndf 0). chi2/ndf is not used: it grows with the statistics because the
    PHC function does not follow the profile exactly.
Segments not used get the mean parameters of the used segments of the same side.
"""
import statistics
import sys
from pathlib import Path

import numpy as np
import uproot

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import config  # noqa: E402

PHC_TYPE, PHC_NPAR = 1, 3
MIN_FIT_ENTRIES = 300


def auto_bad_segments(tree):
    """Segments whose fit has too few entries, or no fit, on any side."""
    if "n_fit" not in tree:
        raise RuntimeError("fit output without n_fit/chi2_ndf: fit again with the current *_PHC binary")
    bad = []
    for i, ch in enumerate(tree["ch"]):
        n, chi = np.asarray(tree["n_fit"][i]), np.asarray(tree["chi2_ndf"][i])
        if np.any(n < MIN_FIT_ENTRIES) or np.any(chi < 0):
            bad.append(int(ch))
    return bad


def make_dictdata(root_file_path, det, good_ch_range=None):
    """(rows, not-used segments). good_ch_range = hand range, or None for the automatic choice."""
    tree = uproot.open(root_file_path)["tree"].arrays(library="np")
    cid = config.DETECTOR_I_LIST[det]
    if good_ch_range is None:
        bad = auto_bad_segments(tree)
    else:
        bad = [int(ch) for ch in tree["ch"] if not good_ch_range[0] <= ch <= good_ch_range[1]]
    data, good = {}, {0: [[], [], []], 1: [[], [], []]}
    for i, ch in enumerate(tree["ch"]):
        if int(ch) in bad:
            continue
        for ud in range(len(tree["p0_val"][i])):
            pars = [tree[f"p{j}_val"][i][ud] for j in range(3)]
            data[f"{cid}-0-{ch:.0f}-{ud}"] = [PHC_TYPE, PHC_NPAR] + pars
            if ud in good:
                for j in range(3):
                    good[ud][j].append(pars[j])
    for ch in bad:
        for ud, lists in good.items():
            if lists[0]:
                data[f"{cid}-0-{ch}-{ud}"] = [PHC_TYPE, PHC_NPAR] + [statistics.mean(v) for v in lists]
    return data, bad
