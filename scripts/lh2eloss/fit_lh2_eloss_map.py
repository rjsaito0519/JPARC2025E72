#!/usr/bin/env python3
"""Fit an LH2 eloss map (Lh2ElossMan 5-column format) to per-sample G4 truth.

Input: the TTree "samples" written by `kpsc_g4_lh2_eloss_map --dump-samples`
(p1 = |p| at the first TPC hit [MeV/c], L = LH2 path length [mm],
loss = p_vertex - p1 [MeV/c], cap = end-cap exit flag, holdout flag).

For each exit surface (side, cap) the grid node values are chosen so that the
clamped bilinear interpolation used by Lh2ElossMan::Loss matches the samples
in the L1 (median) sense, with a weak curvature penalty along p and L that
keeps nodes with few samples smooth. Solved by iteratively reweighted sparse
least squares.

Usage:
  fit_lh2_eloss_map.py samples_proton.root --particle proton \
      --pedges 80,90,... --ledges 0,2,5,... -o Lh2ElossMap_proton [--lam 0.3]
      [--exclude-holdout] [--min-cap N]
"""
import argparse
import sys

import numpy as np
import scipy.sparse as sp
import uproot
from scipy.sparse.linalg import lsqr

PID = {"proton": 0, "pion": 1, "kaon": 2}
TARGET = (("TargetRadius", 40.0), ("TargetHalfHeight", 50.0), ("TargetCenterX", 0.0),
          ("TargetCenterY", 0.0), ("TargetCenterZ", -143.0))


def design(pe, le, p, L):
    """Sparse matrix of the clamped bilinear weights (rows: samples, cols: nodes)."""
    p = np.clip(p, pe[0], pe[-1])
    L = np.clip(L, le[0], le[-1])
    ip = np.clip(np.searchsorted(pe, p, side="right") - 1, 0, len(pe) - 2)
    il = np.clip(np.searchsorted(le, L, side="right") - 1, 0, len(le) - 2)
    tp = (p - pe[ip]) / (pe[ip + 1] - pe[ip])
    tl = (L - le[il]) / (le[il + 1] - le[il])
    nl, n = len(le), len(p)
    rows = np.repeat(np.arange(n), 4)
    cols = np.stack([ip * nl + il, (ip + 1) * nl + il,
                     ip * nl + il + 1, (ip + 1) * nl + il + 1], 1).ravel()
    vals = np.stack([(1 - tp) * (1 - tl), tp * (1 - tl),
                     (1 - tp) * tl, tp * tl], 1).ravel()
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, len(pe) * nl))


def curvature(pe, le, lam):
    """Second differences along p and L on the non-uniform grid, each scaled by
    the squared median spacing of its axis times lam."""
    npn, nl = len(pe), len(le)
    r, c, v = [], [], []
    k = 0
    sp_p = lam * np.median(np.diff(pe)) ** 2
    sp_l = lam * np.median(np.diff(le)) ** 2
    for i in range(1, npn - 1):
        h1, h2 = pe[i] - pe[i - 1], pe[i + 1] - pe[i]
        for j in range(nl):
            for cc, vv in (((i - 1) * nl + j, 2 / (h1 * (h1 + h2))), (i * nl + j, -2 / (h1 * h2)),
                           ((i + 1) * nl + j, 2 / (h2 * (h1 + h2)))):
                r.append(k); c.append(cc); v.append(vv * sp_p)
            k += 1
    for i in range(npn):
        for j in range(1, nl - 1):
            h1, h2 = le[j] - le[j - 1], le[j + 1] - le[j]
            for cc, vv in ((i * nl + j - 1, 2 / (h1 * (h1 + h2))), (i * nl + j, -2 / (h1 * h2)),
                           (i * nl + j + 1, 2 / (h2 * (h1 + h2)))):
                r.append(k); c.append(cc); v.append(vv * sp_l)
            k += 1
    return sp.csr_matrix((v, (r, c)), shape=(k, npn * nl))


def fit_grid(pe, le, p, L, y, lam, n_iter=40, eps=0.05):
    A = design(pe, le, p, L)
    D = curvature(pe, le, lam)
    zeros = np.zeros(D.shape[0])
    x = lsqr(sp.vstack([A, D]), np.r_[y, zeros], atol=1e-10, btol=1e-10)[0]
    for _ in range(n_iter):
        w = 1.0 / np.sqrt(np.maximum(np.abs(y - A @ x), eps))
        x_new = lsqr(sp.vstack([sp.diags(w) @ A, D]), np.r_[w * y, zeros],
                     atol=1e-10, btol=1e-10, x0=x)[0]
        done = np.max(np.abs(x_new - x)) < 1e-4
        x = x_new
        if done:
            break
    return np.maximum(x.reshape(len(pe), len(le)), 0.0)


def edges(text):
    e = np.array([float(x) for x in text.split(",")])
    if len(e) < 2 or np.any(np.diff(e) <= 0):
        sys.exit(f"edges must be increasing: {text}")
    return e


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("samples")
    ap.add_argument("--particle", required=True, choices=sorted(PID))
    ap.add_argument("--pedges", required=True)
    ap.add_argument("--ledges", required=True)
    ap.add_argument("--lam", type=float, default=0.3, help="curvature penalty weight (default 0.3)")
    ap.add_argument("--exclude-holdout", action="store_true", help="fit only samples with holdout==0")
    ap.add_argument("--min-cap", type=int, default=2000,
                    help="fit a separate cap grid only with at least this many cap samples (else cap = side grid)")
    ap.add_argument("--vessel-radius", type=float, default=48.5)
    ap.add_argument("-o", "--output", required=True)
    a = ap.parse_args()

    pe, le = edges(a.pedges), edges(a.ledges)
    d = uproot.open(a.samples)["samples"].arrays(["p1", "L", "loss", "cap", "holdout"], library="np")
    use = np.ones(len(d["p1"]), bool)
    if a.exclude_holdout:
        use &= d["holdout"] == 0
    grids, nsamp = {}, {}
    for ex in (0, 1):
        s = use & (d["cap"] == ex)
        nsamp[ex] = int(s.sum())
        if ex == 1 and nsamp[ex] < a.min_cap:
            grids[1] = grids[0]
            continue
        grids[ex] = fit_grid(pe, le, d["p1"][s].astype(float), d["L"][s].astype(float),
                             d["loss"][s].astype(float), a.lam)

    pid = PID[a.particle]
    with open(a.output, "w") as f:
        f.write("#\n")
        f.write(f"#  Lh2ElossMap_{a.particle}: LH2 momentum-loss map (Lh2ElossMan 5-column format)\n")
        f.write("#  Loss = p_vertex - p_TPC1st [MeV/c]; apply |p|_vtx = |p|_TPC + Loss\n")
        f.write("#  p axis = p_TPC1st; L = path length in the LH2 cylinder from the vertex [mm]\n")
        f.write("#  node values: L1 (median) fit of clamped bilinear interpolation to G4 truth,\n")
        f.write(f"#    curvature penalty lam={a.lam} (fit_lh2_eloss_map.py)\n")
        f.write(f"#  samples: {a.samples}{' (holdout excluded)' if a.exclude_holdout else ''}\n")
        f.write(f"#  N side={nsamp[0]}  N cap={nsamp[1]}"
                f"{'' if nsamp[1] >= a.min_cap else '  (cap rows = side grid)'}\n#\n")
        for k, v in TARGET:
            f.write(f"{k:<18} {v}    # mm\n")
        f.write(f"{'VesselRadius':<18} {a.vessel_radius}    # mm\n\n")
        f.write("# particle  exit  p[MeV/c]  L[mm]  Loss[MeV/c]  (particle 0=proton 1=pion 2=kaon, exit 0=side 1=cap)\n")
        for ex in (0, 1):
            for i, p in enumerate(pe):
                for j, L in enumerate(le):
                    f.write(f"{pid}  {ex}  {p:.6f}  {L:.6f}  {grids[ex][i, j]:.6f}\n")
    print(f"wrote {a.output}: {a.particle} np={len(pe)} nL={len(le)} N side={nsamp[0]} cap={nsamp[1]}")


if __name__ == "__main__":
    main()
