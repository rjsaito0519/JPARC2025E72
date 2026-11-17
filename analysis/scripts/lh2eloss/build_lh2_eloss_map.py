#!/usr/bin/env python3
"""Build the LH2 momentum-loss map of the analyzer's Lh2ElossMan from G4 single-particle gun samples.

Input per particle: two sample files written by `lh2eloss_g4_samples` (TTree "samples"):
  INNER  gun vertices r < 48.4 mm (liquid and vacuum inside the vessel)
  OUTER  gun vertices r < 53.5 mm (adds the gap between vessel and holder)
The vessel grids use INNER only; the gap and gap_cross grids use INNER + OUTER.
Every sample is classified with lh2eloss_map.classify (region, L, exit, holder, cos); samples with
holdout == 0 are fitted, holdout == 1 are only used for the residual table printed at the end.

Fit per grid: clamped multilinear L1 (median) fit of loss = p0 - p1 against (p1, L, cos) with a curvature
penalty LAMBDA along each axis (IRLS, normal equations). The gap grids have no L axis (single node L = 0).
A cap grid with fewer than MIN_SAMPLES samples is written as a copy of the side grid of the same region
and holder (cap x window); any other grid below MIN_SAMPLES is an error.

Usage:
  build_lh2_eloss_map.py -o Lh2ElossMap_e72 --samples proton INNER OUTER --samples pion INNER OUTER
"""
import argparse
import datetime
import os
import sys

import numpy as np
import scipy.sparse as sp
import uproot
from scipy.sparse.linalg import spsolve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lh2eloss_map as M  # noqa: E402

# G4 target geometry (hyptpc E72) [mm, deg]
GEO_E72 = dict(TargetRadius=40.0, TargetHalfHeight=50.0, TargetCenterX=0.0, TargetCenterY=0.0, TargetCenterZ=-143.0,
               VesselRadius=48.5, VesselHalfHeight=50.0, HolderInnerRadius=53.5, HolderRadius=55.0,
               HolderWindowHalfPhi=58.0, HolderWindowHalfHeight=50.0, HolderBottomY=-119.5)

P_AXIS = {
    M.PROTON: np.r_[np.arange(80, 301, 10), np.arange(320, 381, 20), 400, 450, 500, 550, 600, 700, 800, 900, 1000, 1200],
    M.PION: np.r_[np.arange(30, 121, 10), 135, 150, 175, 200, 250, 300, 350, 400, 500, 600],
    M.KAON: np.r_[120, 150, 175, 200, 225, 250, 275, np.arange(300, 701, 50), 800, 900, 1000, 1200],
}
L_AXIS = np.array([0, 2, 5, 7.5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 115, 130.])
COS_AXIS = np.array([0, .3, .5, .65, .8, .9, .97, 1.0])
LAMBDA = 0.3
MIN_SAMPLES = 2000
LOW_P = {M.PROTON: 300, M.PION: 150, M.KAON: 250}     # residual table: "low p" = p1 below this [MeV/c]
BRANCHES = ["p1", "loss", "vx", "vy", "vz", "ux", "uy", "uz", "holdout"]


# ---------------------------------------------------------------- fit

def _locate(e, x):
    x = np.clip(x, e[0], e[-1])
    i = np.clip(np.searchsorted(e, x, side="right") - 1, 0, len(e) - 2)
    return i, (x - e[i]) / (e[i + 1] - e[i])


def design(axes, coords):
    """Sparse matrix of the clamped multilinear weights (rows: samples, columns: grid nodes)."""
    shape = [len(e) for e in axes]
    strides = np.cumprod([1] + shape[::-1])[:-1][::-1]
    n = len(coords[0])
    loc = [_locate(e, np.asarray(c, float)) for e, c in zip(axes, coords)]
    rows, cols, vals = [], [], []
    for corner in np.ndindex(*([2] * len(axes))):
        idx = np.zeros(n, np.int64)
        w = np.ones(n)
        for d, (i, t) in enumerate(loc):
            idx += (i + corner[d]) * strides[d]
            w *= t if corner[d] else (1 - t)
        rows.append(np.arange(n)); cols.append(idx); vals.append(w)
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(n, int(np.prod(shape))))


def curvature(axes, lam):
    """Second differences along each axis on its non-uniform nodes, scaled by lam * (median spacing)^2."""
    shape = [len(e) for e in axes]
    strides = np.cumprod([1] + shape[::-1])[:-1][::-1]
    r, c, v = [], [], []
    k = 0
    for d, e in enumerate(axes):
        if len(e) < 3:
            continue
        scale = lam * np.median(np.diff(e)) ** 2
        for base in np.ndindex(*[s if dd != d else 1 for dd, s in enumerate(shape)]):
            b = sum(bi * strides[dd] for dd, bi in enumerate(base))
            for i in range(1, len(e) - 1):
                h1, h2 = e[i] - e[i - 1], e[i + 1] - e[i]
                for off, w in ((-1, 2 / (h1 * (h1 + h2))), (0, -2 / (h1 * h2)), (1, 2 / (h2 * (h1 + h2)))):
                    r.append(k); c.append(b + (i + off) * strides[d]); v.append(w * scale)
                k += 1
    return sp.csr_matrix((v, (r, c)), shape=(k, int(np.prod(shape))))


def fit(axes, coords, y, lam=LAMBDA, n_iter=40, eps=0.05, ridge=1e-9):
    """Node values minimising sum |y - interp| + curvature penalty (IRLS), clipped to >= 0."""
    A = design(axes, coords).tocsr()
    D = curvature(axes, lam)
    DtD = (D.T @ D).tocsr()
    reg = ridge * sp.identity(A.shape[1], format="csr")
    w = np.ones(len(y))
    x = None
    for _ in range(n_iter + 1):
        Aw = sp.diags(w) @ A
        xn = spsolve((A.T @ Aw + DtD + reg).tocsc(), Aw.T @ y)
        done = x is not None and np.max(np.abs(xn - x)) < 1e-4
        x = xn
        if done:
            break
        w = 1.0 / np.maximum(np.abs(y - A @ x), eps)
    return np.maximum(x.reshape([len(e) for e in axes]), 0.0)


# ---------------------------------------------------------------- samples

def load(path, geo, source):
    t = uproot.open(path)["samples"]
    s = {k: np.asarray(v, float) for k, v in t.arrays(BRANCHES, library="np").items()}
    s.update(M.classify(geo, np.c_[s["vx"], s["vy"], s["vz"]], np.c_[s["ux"], s["uy"], s["uz"]]))
    s["source"] = np.full(len(s["p1"]), source)
    return s


def join(a, b):
    return {k: np.r_[a[k], b[k]] for k in a}


def build_particle(pa, s):
    """Grids {(pa, region, exit, holder): grid} from the classified samples s."""
    pe = P_AXIS[pa]
    train = s["holdout"] == 0
    grids, copies = {}, []
    for rg, ex, ho in M.REQUIRED:
        key = (pa, rg, ex, ho)
        m = train & (s["region"] == rg) & (s["exit"] == ex) & (s["holder"] == ho)
        if rg == M.VESSEL:
            m &= s["source"] == 0
        n = int(m.sum())
        print(f"  {M.grid_name(key):32s} n {n:8d}")
        if n < MIN_SAMPLES:
            if ex == M.CAP:
                copies.append(key)
                continue
            sys.exit(f"error: {M.grid_name(key)} has {n} samples (< {MIN_SAMPLES})")
        if rg == M.GAP:
            T = fit([pe, COS_AXIS], [s["p1"][m], s["cos_alpha"][m]], s["loss"][m])
            grids[key] = dict(p=pe, L=np.array([0.0]), c=COS_AXIS, loss=T.reshape(len(pe), 1, len(COS_AXIS)))
        else:
            T = fit([pe, L_AXIS, COS_AXIS], [s["p1"][m], s["L"][m], s["cos_alpha"][m]], s["loss"][m])
            grids[key] = dict(p=pe, L=L_AXIS, c=COS_AXIS, loss=T)
    for key in copies:
        side = (key[0], key[1], M.SIDE, key[3])
        print(f"  {M.grid_name(key)}: copy of {M.grid_name(side)}")
        grids[key] = grids[side]
    return grids, copies


def residual_table(pa, grids, geo, s):
    ho = s["holdout"] == 1
    trk = {k: v[ho] for k, v in s.items()}
    res = trk["loss"] - M.loss(grids, trk, trk["p1"], pa)
    low = trk["p1"] < LOW_P[pa]
    print(f"  holdout: true - map [MeV/c] (median, sig68); low p = p1 < {LOW_P[pa]} MeV/c")
    for rg in (M.VESSEL, M.GAP, M.GAP_CROSS, M.OUTSIDE):
        for h in (M.WINDOW, M.G10, -1):
            m = (trk["region"] == rg) & (trk["holder"] == h)
            for tag, mm in (("all ", m), ("lowp", m & low)):
                if mm.sum() < 20:
                    continue
                r = res[mm]
                rname = M.REGION_NAMES[rg] if rg != M.OUTSIDE else "outside"
                hname = M.HOLDER_NAMES[h] if h >= 0 else "-"
                print(f"    {rname:9s} {hname:6s} {tag} N {mm.sum():8d}  med {np.median(r):+6.2f}"
                      f"  sig68 {np.diff(np.percentile(r, [16, 84]))[0] / 2:6.2f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-o", "--output", required=True, help="map file to write")
    ap.add_argument("--samples", nargs=3, action="append", required=True, metavar=("PARTICLE", "INNER", "OUTER"),
                    help="particle (proton|pion|kaon) and its two sample files")
    ap.add_argument("--note", action="append", default=[], help="extra header line (repeatable)")
    a = ap.parse_args()

    geo = GEO_E72
    grids, header_src, copied = {}, [], []
    for name, inner, outer in a.samples:
        pa = M.PARTICLE_NAMES.index(name)
        print(f"{name}: {inner} + {outer}")
        s = join(load(inner, geo, 0), load(outer, geo, 1))
        g, c = build_particle(pa, s)
        residual_table(pa, g, geo, s)
        grids.update(g)
        copied += c
        header_src.append(f"  {name}: {os.path.abspath(inner)}")
        header_src.append(f"  {' ' * len(name)}  {os.path.abspath(outer)}")

    names = ", ".join(M.PARTICLE_NAMES[p] for p in sorted({k[0] for k in grids}))
    header = [
        f"{os.path.basename(a.output)}: LH2 momentum-loss map for Lh2ElossMan",
        "Loss [MeV/c] = |p| at the vertex - |p| at the first TPC hit; the analyzer adds it to the TPC |p|.",
        f"Source: G4 single-particle gun ({names}), built {datetime.date.today().isoformat()} by",
        f"  {os.path.abspath(__file__)}",
        "Samples (inner: vertex r < 48.4 mm, outer: r < 53.5 mm):",
        *header_src,
        f"Fit: L1 median per node, curvature penalty {LAMBDA}; vessel grids from the inner samples only.",
        "region: vessel = vertex in the vacuum vessel; gap = between vessel and holder, ray misses the vessel (L = 0);",
        "        gap_cross = between vessel and holder, ray enters the vessel",
    ]
    if copied:
        header.append("Copies of the side grid (too few samples): "
                      + ", ".join(M.grid_name(k) for k in copied))
    header += a.note
    M.write_map(a.output, grids, geo, header)
    g2, geo2 = M.read_map(a.output)
    assert set(g2) == set(grids) and geo2 == geo
    print(f"wrote {a.output}: {len(grids)} grids")


if __name__ == "__main__":
    main()
