#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Offline prediction of a TPCPOS map change on beam-through momenta (no re-reconstruction).

Hits of a TPCHelix file reconstructed with map P are shifted by (M - P) at their
cluster centre pad (Y-flat maps), the bending-plane circle is refitted with free
centre and radius, and mom0 is rescaled by R'/R. Track finding, clustering and the
full helix fit are not redone, so this is only a proxy for a real reconstruction.

Usage:
  tpcpos_offline_refit.py <bt_helix.root> --applied P --cand LABEL=M [--cand ...]
                          [--max-entries N]
Prints per-run median(mom0')/p_beam and the charge-pair difference per label.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tpc_position_update import SPEC, Z_TARGET, read_map  # noqa: E402


def map_lookup(path: Path) -> Dict[int, np.ndarray]:
    _, _, _, grid = read_map(path)
    return {k[0] * 1000 + k[1]: np.array(v[0][1:]) for k, v in grid.items()}


def free_circle(hX, hY, cx, cy, R, n_iter=6):
    """Gauss-Newton fit of (cx, cy, R) per track (jagged hits)."""
    for _ in range(n_iter):
        dx = hX - cx
        dy = hY - cy
        rr = np.sqrt(dx * dx + dy * dy)
        d = rr - R
        jx, jy, jr = -dx / rr, -dy / rr, -1.0 + 0.0 * rr
        J = [jx, jy, jr]
        A = [[ak.to_numpy(ak.sum(J[i] * J[j], axis=-1)) for j in range(3)] for i in range(3)]
        b = [ak.to_numpy(ak.sum(J[i] * d, axis=-1)) for i in range(3)]
        A = np.moveaxis(np.array(A), -1, 0)
        b = np.moveaxis(np.array(b), -1, 0)
        step = np.linalg.solve(A + 1e-9 * np.eye(3), b[..., None])[..., 0]
        cx = cx - step[:, 0]
        cy = cy - step[:, 1]
        R = R - step[:, 2]
    return R


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("helix_root", type=Path)
    ap.add_argument("--applied", type=Path, required=True)
    ap.add_argument("--cand", action="append", default=[], help="LABEL=map")
    ap.add_argument("--max-entries", type=int, default=None)
    ap.add_argument("--stride-chunks", type=int, default=0,
                    help="if >0, read this many evenly spaced chunks of 20k entries instead")
    args = ap.parse_args()

    applied = map_lookup(args.applied)
    cands = [("refit_same", applied)] + [(c.split("=", 1)[0], map_lookup(Path(c.split("=", 1)[1])))
                                          for c in args.cand]
    t = uproot.open(args.helix_root)["tpc"]
    br = ["run_number", "nhtrack", "helix_dz", "mom0", "charge", "helix_cx", "helix_cy",
          "helix_r", "hitlayer", "hitpos_x", "hitpos_z", "track_cluster_row_center"]
    if args.stride_chunks > 0:
        starts = np.linspace(0, t.num_entries - 20000, args.stride_chunks).astype(int)
        ranges = [(int(s), int(s) + 20000) for s in starts]
    else:
        ranges = [(0, args.max_entries or t.num_entries)]
    res: Dict[str, Dict[int, List[np.ndarray]]] = {lab: {} for lab, _ in cands}
    res["helix"] = {}
    for start, stop in ranges:
        for a in t.iterate(br, step_size="300 MB", entry_start=start, entry_stop=stop):
            runs = ak.to_numpy(a["run_number"])
            pb = np.array([SPEC[int(r)].p_beam if int(r) in SPEC else np.nan for r in runs])
            qe = np.array([SPEC[int(r)].expect_charge if int(r) in SPEC else 0 for r in runs])
            pb_t = ak.broadcast_arrays(pb, a["mom0"])[0]
            q_t = ak.broadcast_arrays(qe, a["mom0"])[0]
            sel = ((a["nhtrack"] >= 15) & (abs(a["helix_dz"]) < 0.05) & (a["charge"] == q_t)
                   & (abs(a["mom0"] - pb_t) < np.maximum(0.12, 0.2 * pb_t)))
            run_t = ak.to_numpy(ak.flatten(ak.broadcast_arrays(runs, a["mom0"])[0][sel]))
            mom = ak.to_numpy(ak.flatten(a["mom0"][sel]))
            cx = ak.to_numpy(ak.flatten(a["helix_cx"][sel]))
            cy = ak.to_numpy(ak.flatten(a["helix_cy"][sel]))
            R = ak.to_numpy(ak.flatten(a["helix_r"][sel]))
            hx = ak.flatten(a["hitpos_x"][sel], axis=1)
            hz = ak.flatten(a["hitpos_z"][sel], axis=1)
            key = (np.rint(ak.flatten(a["hitlayer"][sel], axis=1)) * 1000
                   + np.rint(ak.flatten(a["track_cluster_row_center"][sel], axis=1)))
            counts = ak.num(hx)
            kf = ak.to_numpy(ak.flatten(key)).astype(np.int64)
            for r in np.unique(run_t):
                res["helix"].setdefault(int(r), []).append(mom[run_t == r])
            for lab, M in cands:
                dvx = np.array([(M.get(k, np.zeros(2)) - applied.get(k, np.zeros(2)))[0] for k in kf])
                dvz = np.array([(M.get(k, np.zeros(2)) - applied.get(k, np.zeros(2)))[1] for k in kf])
                hx2 = hx + ak.unflatten(dvx, counts)
                hz2 = hz + ak.unflatten(dvz, counts)
                R2 = free_circle(-hx2, hz2 - Z_TARGET, cx, cy, R)
                m2 = mom * R2 / R
                for r in np.unique(run_t):
                    res[lab].setdefault(int(r), []).append(m2[run_t == r])
    labels = ["helix"] + [lab for lab, _ in cands]
    med = {lab: {r: float(np.median(np.concatenate(v))) for r, v in res[lab].items()} for lab in labels}
    print(f"{'run':>6} {'q':>3} {'p':>6} " + " ".join(f"{lab:>12}" for lab in labels))
    for r in sorted(med["helix"]):
        s = SPEC[r]
        print(f"{r:6d} {s.expect_charge:+3d} {s.p_beam:6.3f} "
              + " ".join(f"{med[lab].get(r, np.nan) / s.p_beam:12.4f}" for lab in labels))
    pairs: Dict[tuple, Dict[int, List[int]]] = {}
    for r in med["helix"]:
        s = SPEC[r]
        pairs.setdefault((s.species, s.p_beam), {}).setdefault(s.expect_charge, []).append(r)
    print("charge-pair difference ratio(+) - ratio(-)")
    for (sp, p), byq in sorted(pairs.items()):
        if +1 not in byq or -1 not in byq:
            continue
        cells = []
        for lab in labels:
            rp = np.mean([med[lab][r] for r in byq[+1]]) / p
            rm = np.mean([med[lab][r] for r in byq[-1]]) / p
            cells.append(f"{rp - rm:+12.4f}")
        print(f"{sp:>4} {p:6.3f}      " + " ".join(cells))
    return 0


if __name__ == "__main__":
    sys.exit(main())
