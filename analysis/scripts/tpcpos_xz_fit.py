#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TPCPOS update with both in-plane components (x, z) per pad from inclusive helix residuals.

The DST residual (residual_x, residual_z) of a hit lies along the track normal n in the
bending plane (centre -> hit). A pad offset d shifts the residual by r_n = n . d, so per pad
  M = sum n n^T,  b = sum n r_n,  d = M^-1 b.
The component of d along the track direction is not seen by a single track; it is
determined only where the tracks crossing the pad have different directions. Pads whose
smallest eigenvalue of M per hit is below --cond keep the parent value.
Note: per pad, the radial part of d is degenerate with a cluster-position bias linear in the
track crossing angle.

Pad basis: e_t = (cos th, -sin th), e_r = (sin th, cos th) in (x, z), th = pad angle.
Update (additive): node += -alpha * d for every y node of the pad (--y-bins 1), or for the
nodes with y < 0 / y >= 0 from the hits in that half (--y-bins 2).
Center-frame pads and --zero-layers are set to (0, 0).

Usage:
  tpcpos_xz_fit.py HELIX.root [HELIX.root ...] --parent MAP --out MAP
                   [--alpha A] [--cond C] [--threshold N] [--y-bins 1|2] [--zero-layers L0-L1]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tpc_position_update import Z_TARGET, read_map, write_map  # noqa: E402

Y_OFFSET = 10_000_000  # key offset of the y >= 0 half


def flat(v):
    return ak.to_numpy(ak.flatten(v, axis=None))


def accumulate(path: Path, pad_phi: Dict[int, float], acc: Dict[int, np.ndarray], y_bins: int, step: str):
    br = ["nhtrack", "helix_dz", "helix_cx", "helix_cy", "hitlayer", "track_cluster_row_center",
          "hitpos_x", "hitpos_y", "hitpos_z", "residual_x", "residual_z"]
    known = np.array(sorted(pad_phi))
    for a in uproot.open(path)["tpc"].iterate(br, step_size=step):
        sel = (a["nhtrack"] >= 15) & (abs(a["helix_dz"]) < 0.05)
        cx = flat(ak.broadcast_arrays(a["helix_cx"][sel], a["hitpos_x"][sel])[0])
        cy = flat(ak.broadcast_arrays(a["helix_cy"][sel], a["hitpos_x"][sel])[0])
        x, y, z = flat(a["hitpos_x"][sel]), flat(a["hitpos_y"][sel]), flat(a["hitpos_z"][sel])
        rx, rz = flat(a["residual_x"][sel]), flat(a["residual_z"][sel])
        key = (np.rint(flat(a["hitlayer"][sel])).astype(np.int64) * 1000
               + np.rint(flat(a["track_cluster_row_center"][sel])).astype(np.int64))
        # Helix frame X = -x, Y = z - Z_TARGET; the normal (X - cx, Y - cy) in (x, z).
        nx, nz = x + cx, (z - Z_TARGET) - cy
        nn = np.hypot(nx, nz)
        nx, nz = nx / nn, nz / nn
        rn = rx * nx + rz * nz
        ok = np.isfinite(rn) & (np.abs(rn) < 5.0) & np.isin(key, known)
        key, y, nx, nz, rn = key[ok], y[ok], nx[ok], nz[ok], rn[ok]
        th = np.array([pad_phi[int(k)] for k in key])
        nt = nx * np.cos(th) - nz * np.sin(th)
        nr = nx * np.sin(th) + nz * np.cos(th)
        if y_bins == 2:
            key = key + np.where(y >= 0.0, Y_OFFSET, 0)
        uk, inv = np.unique(key, return_inverse=True)
        cols = [np.bincount(inv, weights=c) for c in
                (nt * nt, nt * nr, nr * nr, nt * rn, nr * rn, np.ones_like(rn))]
        for i, k in enumerate(uk):
            v = acc.setdefault(int(k), np.zeros(6))
            v += [c[i] for c in cols]
        print(f"[xz] {path.name}: {len(rn)} hits", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("helix_root", type=Path, nargs="+")
    ap.add_argument("--parent", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--pad-dir", type=Path, default=Path("/group/had/sks/Users/sryuta/tmp/debug/tpcpos/cfit"),
                    help="dir with center_frame.txt and pad_phi.txt (dump_pads.C)")
    ap.add_argument("--alpha", type=float, default=0.7)
    ap.add_argument("--cond", type=float, default=0.03,
                    help="minimum smallest eigenvalue of M per hit (0: radial unseen, 0.5: isotropic)")
    ap.add_argument("--threshold", type=int, default=100)
    ap.add_argument("--y-bins", type=int, choices=[1, 2], default=1)
    ap.add_argument("--zero-layers", default="")
    ap.add_argument("--step", default="300 MB")
    args = ap.parse_args()

    pad_phi = {}
    with open(args.pad_dir / "pad_phi.txt") as f:
        for line in f:
            s = line.split()
            pad_phi[int(s[0]) * 1000 + int(s[1])] = float(s[2])
    cf = set()
    with open(args.pad_dir / "center_frame.txt") as f:
        for line in f:
            s = line.split()
            if len(s) >= 2:
                cf.add((int(s[0]), int(s[1])))
    zl = tuple(int(v) for v in args.zero_layers.split("-")) if args.zero_layers else None

    acc: Dict[int, np.ndarray] = {}
    for p in args.helix_root:
        accumulate(p, pad_phi, acc, args.y_bins, args.step)

    n_y, y0, dy, grid = read_map(args.parent)
    n_upd, n_low, dts, drs = 0, 0, [], []
    for k, v in acc.items():
        half = k // Y_OFFSET
        pk = k % Y_OFFSET
        key = (pk // 1000, pk % 1000)
        if v[5] < args.threshold or key in cf or key not in grid or (zl and zl[0] <= key[0] <= zl[1]):
            continue
        M = np.array([[v[0], v[1]], [v[1], v[2]]])
        if np.linalg.eigvalsh(M)[0] / v[5] < args.cond:
            n_low += 1
            continue
        d_t, d_r = np.linalg.solve(M, v[3:5])
        th = pad_phi[pk]
        vx = -args.alpha * (d_t * np.cos(th) + d_r * np.sin(th))
        vz = -args.alpha * (-d_t * np.sin(th) + d_r * np.cos(th))
        for node in grid[key]:
            if args.y_bins == 2 and (node[0] >= 0.0) != (half == 1):
                continue
            node[1] += vx
            node[2] += vz
        n_upd += 1
        dts.append(d_t)
        drs.append(d_r)
    zero = {key for key in grid if zl and zl[0] <= key[0] <= zl[1]}
    for key in cf | zero:
        for node in grid.get(key, []):
            node[1] = 0.0
            node[2] = 0.0
    write_map(args.out, n_y, y0, dy, grid)
    dts, drs = np.array(dts), np.array(drs)
    print(f"y-bins {args.y_bins}: updated {n_upd} pad(-halves), below cond {n_low}")
    print(f"median|d_t|={np.median(np.abs(dts)):.4f} median|d_r|={np.median(np.abs(drs)):.4f} "
          f"rms d_t={dts.std():.4f} rms d_r={drs.std():.4f}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
