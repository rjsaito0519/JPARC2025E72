#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""One (x, z) update step of the TPC hit-position correction map (TPCPOS).

Per pad, the in-plane offset d = d_t e_t + d_r e_r is fitted from the residuals of two
kinds of tracks:

  --bt FILE    beam-through helix DST (runs in compare_tpcpos_beamthrough.RUNS); each
               track is refitted in the bending plane with a Gaussian prior on its radius
               around the run median momentum (relative width --momentum-width), as
               tpc_position_update.py --momentum-constraint soft.
  --incl FILE  any helix DST (multi-pion, K beam, ...); the DST residual (residual_x/z) of
               every well-measured track (--min-nhit, --max-dip, --min-mom, --max-chisqr).

A residual r of a hit lies along the track normal n in the bending plane, and a pad
offset d shifts it by n . d, so per pad (t, r basis)
  M = sum w n n^T,  b = sum w n r_n,  (M + ridge * N * e_r e_r^T) d = b
with w = --bt-weight for beam-through hits and 1 otherwise, N = sum w. The ridge pulls
the radial part, which a pad sees only from the spread of track crossing angles (and
which is degenerate with a cluster bias linear in that angle), towards 0; ridge -> inf
gives the tangential-only update. Pads whose smallest eigenvalue of M / N stays below
--cond are tangential-only.

Pads with fewer than --threshold hits get the median update of the fitted pads within
--fill-rows rows in the same layer (if at least two exist); others keep the parent.
Update (additive, y-flat): node += -alpha * d. Centre-frame pads and --zero-layers are
set to (0, 0).

Usage:
  tpcpos_xz_update.py --parent MAP --out MAP [--bt BT.root] [--incl HELIX.root ...]
                      [--alpha A] [--ridge R] [--momentum-width S] [--fill-rows K]
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
from tpc_position_update import (SPEC, Z_TARGET, mom_window, read_map, run_medians,  # noqa: E402
                                 soft_residual, write_map)

# per-pad sums: [tt, tr, rr, t*rn, r*rn, w, n_hits]
N_ACC = 7


def flat(v):
    return ak.to_numpy(ak.flatten(v, axis=None))


def add_hits(acc, key, nx, nz, rn, w, pad_phi):
    ok = np.isfinite(rn) & (np.abs(rn) < 5.0) & np.isin(key, pad_phi["keys"])
    key, nx, nz, rn = key[ok], nx[ok], nz[ok], rn[ok]
    th = pad_phi["th"][np.searchsorted(pad_phi["keys"], key)]
    nt = nx * np.cos(th) - nz * np.sin(th)
    nr = nx * np.sin(th) + nz * np.cos(th)
    ww = np.full(len(rn), w)
    uk, inv = np.unique(key, return_inverse=True)
    cols = [np.bincount(inv, weights=c) for c in
            (ww * nt * nt, ww * nt * nr, ww * nr * nr, ww * nt * rn, ww * nr * rn, ww, np.ones_like(rn))]
    for i, k in enumerate(uk):
        v = acc.setdefault(int(k), np.zeros(N_ACC))
        v += [c[i] for c in cols]
    return len(rn)


def pad_key(a, sel):
    return (np.rint(flat(a["hitlayer"][sel])).astype(np.int64) * 1000
            + np.rint(flat(a["track_cluster_row_center"][sel])).astype(np.int64))


def accumulate_bt(path, acc, pad_phi, width, hit_res, weight, step):
    tree = uproot.open(path)["tpc"]
    med = run_medians(tree, step)
    br = ["run_number", "nhtrack", "helix_dz", "mom0", "charge", "helix_cx", "helix_cy",
          "helix_r", "hitlayer", "hitpos_x", "hitpos_z", "track_cluster_row_center"]
    n_trk = n_hit = 0
    for a in tree.iterate(br, step_size=step):
        runs = ak.to_numpy(a["run_number"])
        pr = np.array([med.get(int(r), np.nan) for r in runs])
        pb = np.array([SPEC[int(r)].p_beam if int(r) in SPEC else np.nan for r in runs])
        qe = np.array([SPEC[int(r)].expect_charge if int(r) in SPEC else 0 for r in runs])
        pr_t = ak.broadcast_arrays(pr, a["mom0"])[0]
        pb_t = ak.broadcast_arrays(pb, a["mom0"])[0]
        q_t = ak.broadcast_arrays(qe, a["mom0"])[0]
        sel = ((a["nhtrack"] >= 15) & (abs(a["helix_dz"]) < 0.05) & (a["charge"] == q_t)
               & np.isfinite(pr_t) & (abs(a["mom0"] - pb_t) < np.maximum(0.12, 0.2 * pb_t)))
        cx = ak.flatten(a["helix_cx"][sel])
        cy = ak.flatten(a["helix_cy"][sel])
        R = ak.flatten(a["helix_r"][sel] * pr_t[sel] / a["mom0"][sel])
        hx = ak.flatten(a["hitpos_x"][sel], axis=1)
        hz = ak.flatten(a["hitpos_z"][sel], axis=1)
        d, nX, nY = soft_residual(-hx, hz - Z_TARGET, cx, cy, R, width, hit_res)
        # residual vector in (x, z) is d * (-nX, nY)
        n_hit += add_hits(acc, pad_key(a, sel), flat(-nX), flat(nY), flat(d), weight, pad_phi)
        n_trk += len(cx)
    print(f"[bt] {path.name}: runs {len(med)} tracks {n_trk} hits {n_hit}", flush=True)


def accumulate_incl(path, acc, pad_phi, args):
    br = ["nhtrack", "helix_dz", "helix_cx", "helix_cy", "mom0", "chisqr", "hitlayer",
          "track_cluster_row_center", "hitpos_x", "hitpos_z", "residual_x", "residual_z"]
    n_trk = n_hit = 0
    for a in uproot.open(path)["tpc"].iterate(br, step_size=args.step):
        sel = ((a["nhtrack"] >= args.min_nhit) & (abs(a["helix_dz"]) < args.max_dip)
               & (a["mom0"] > args.min_mom) & (a["chisqr"] < args.max_chisqr))
        cx = flat(ak.broadcast_arrays(a["helix_cx"][sel], a["hitpos_x"][sel])[0])
        cy = flat(ak.broadcast_arrays(a["helix_cy"][sel], a["hitpos_x"][sel])[0])
        x, z = flat(a["hitpos_x"][sel]), flat(a["hitpos_z"][sel])
        # helix frame X = -x, Y = z - Z_TARGET; the normal (X - cx, Y - cy) in (x, z)
        nx, nz = x + cx, (z - Z_TARGET) - cy
        nn = np.hypot(nx, nz)
        nx, nz = nx / nn, nz / nn
        rn = flat(a["residual_x"][sel]) * nx + flat(a["residual_z"][sel]) * nz
        n_hit += add_hits(acc, pad_key(a, sel), nx, nz, rn, 1.0, pad_phi)
        n_trk += int(ak.sum(sel))
    print(f"[incl] {path.name}: tracks {n_trk} hits {n_hit}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--parent", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--bt", type=Path, action="append", default=[])
    ap.add_argument("--incl", type=Path, action="append", default=[])
    ap.add_argument("--pad-dir", type=Path, default=Path("/group/had/sks/Users/sryuta/tmp/debug/tpcpos/cfit"),
                    help="dir with center_frame.txt and pad_phi.txt (dump_pads.C)")
    ap.add_argument("--alpha", type=float, default=0.7)
    ap.add_argument("--ridge", type=float, default=0.02,
                    help="radial ridge per unit weight (0: free, large: tangential only)")
    ap.add_argument("--cond", type=float, default=0.01,
                    help="minimum smallest eigenvalue of M/N for the radial part")
    ap.add_argument("--threshold", type=int, default=100)
    ap.add_argument("--fill-rows", type=int, default=3)
    ap.add_argument("--momentum-width", type=float, default=0.01)
    ap.add_argument("--hit-resolution", type=float, default=0.6)
    ap.add_argument("--bt-weight", type=float, default=1.0)
    ap.add_argument("--min-nhit", type=int, default=15)
    ap.add_argument("--max-dip", type=float, default=0.05)
    ap.add_argument("--min-mom", type=float, default=0.15)
    ap.add_argument("--max-chisqr", type=float, default=1e9)
    ap.add_argument("--zero-layers", default="")
    ap.add_argument("--npz", type=Path, default=None, help="per-pad fit results")
    ap.add_argument("--step", default="300 MB")
    args = ap.parse_args()

    phi = {}
    with open(args.pad_dir / "pad_phi.txt") as f:
        for line in f:
            s = line.split()
            phi[int(s[0]) * 1000 + int(s[1])] = float(s[2])
    keys = np.array(sorted(phi), dtype=np.int64)
    pad_phi = {"keys": keys, "th": np.array([phi[int(k)] for k in keys])}
    cf = set()
    with open(args.pad_dir / "center_frame.txt") as f:
        for line in f:
            s = line.split()
            if len(s) >= 2:
                cf.add((int(s[0]), int(s[1])))
    zl = tuple(int(v) for v in args.zero_layers.split("-")) if args.zero_layers else None

    acc: Dict[int, np.ndarray] = {}
    for p in args.bt:
        accumulate_bt(p, acc, pad_phi, args.momentum_width, args.hit_resolution, args.bt_weight, args.step)
    for p in args.incl:
        accumulate_incl(p, acc, pad_phi, args)

    n_y, y0, dy, grid = read_map(args.parent)
    upd: Dict[tuple, np.ndarray] = {}
    rec = []
    n_tan = 0
    for k, v in acc.items():
        key = (k // 1000, k % 1000)
        if v[6] < args.threshold or key in cf or key not in grid or (zl and zl[0] <= key[0] <= zl[1]):
            continue
        M = np.array([[v[0], v[1]], [v[1], v[2]]])
        ev = np.linalg.eigvalsh(M / v[5])
        if ev[0] < args.cond:
            d_t, d_r = v[3] / v[0], 0.0
            n_tan += 1
        else:
            d_t, d_r = np.linalg.solve(M + np.diag([0.0, args.ridge * v[5]]), v[3:5])
        upd[key] = np.array([d_t, d_r])
        rec.append((k, v[6], ev[0], d_t, d_r))
    # fill thin pads from fitted neighbours in the same layer
    fill = {}
    for key in grid:
        if key in upd or key in cf or (zl and zl[0] <= key[0] <= zl[1]):
            continue
        nb = [upd[(key[0], key[1] + j)] for j in range(-args.fill_rows, args.fill_rows + 1)
              if (key[0], key[1] + j) in upd]
        if len(nb) >= 2:
            fill[key] = np.median(np.array(nb), axis=0)
    for key, (d_t, d_r) in list(upd.items()) + list(fill.items()):
        th = phi[key[0] * 1000 + key[1]]
        vx = -args.alpha * (d_t * np.cos(th) + d_r * np.sin(th))
        vz = -args.alpha * (-d_t * np.sin(th) + d_r * np.cos(th))
        for node in grid[key]:
            node[1] += vx
            node[2] += vz
    zero = {key for key in grid if zl and zl[0] <= key[0] <= zl[1]}
    for key in cf | zero:
        for node in grid.get(key, []):
            node[1] = 0.0
            node[2] = 0.0
    write_map(args.out, n_y, y0, dy, grid)
    r = np.array(rec)
    if args.npz:
        np.savez(args.npz, key=r[:, 0].astype(np.int64), n=r[:, 1], ev=r[:, 2], d_t=r[:, 3], d_r=r[:, 4])
    print(f"fitted pads {len(upd)} (tangential-only {n_tan}), filled {len(fill)}, "
          f"untouched {len(grid) - len(upd) - len(fill) - len(cf | zero)}")
    print(f"median|d_t|={np.median(np.abs(r[:, 3])):.4f} median|d_r|={np.median(np.abs(r[:, 4])):.4f} "
          f"rms d_t={r[:, 3].std():.4f} rms d_r={r[:, 4].std():.4f}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
