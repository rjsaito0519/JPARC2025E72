"""HDPHC rows {"CId-PlId-SegId-UorD": [type, npar, p0, p1, p2]} from a *_PHC fit output.

The detector is given by the caller (update_param.py); its CId comes from lib/config.py.
Segments outside good_ch_range get the mean parameters of the segments inside it.
"""
import statistics
import sys
from pathlib import Path

import numpy as np
import uproot

sys.path.append(str(Path(__file__).resolve().parents[2]))
from lib import config  # noqa: E402

PHC_TYPE, PHC_NPAR = 1, 3


def make_dictdata(root_file_path, det, good_ch_range=(-np.inf, np.inf)):
    tree = uproot.open(root_file_path)["tree"].arrays(library="np")
    cid = config.DETECTOR_I_LIST[det]
    data, good, bad = {}, {0: [[], [], []], 1: [[], [], []]}, []
    for i, ch in enumerate(tree["ch"]):
        if not good_ch_range[0] <= ch <= good_ch_range[1]:
            bad.append(ch)
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
                data[f"{cid}-0-{ch:.0f}-{ud}"] = [PHC_TYPE, PHC_NPAR] + [statistics.mean(v) for v in lists]
    return data
