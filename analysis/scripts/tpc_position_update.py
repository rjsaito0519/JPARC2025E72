#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""One update step of the TPC hit-position correction map (TPCPOS) from beam-through residuals.

Each beam-through track is refitted in the bending plane; the residual of its hits to that
circle, projected on the pad tangent, is averaged per pad and a fraction of it is added to
the parent map. The inclusive helix residual alone cannot see curvature-type distortions (a
fixed sagitta shift is absorbed into the fitted radius), so the radius may be constrained to
a reference momentum (--momentum-constraint):

  none         no radius constraint: pads are updated from the inclusive helix residual
               (DST residual_x/z) of all tracks, like the multi-pion fill.
  soft         the radius is free but has a Gaussian prior around the run median momentum:
               chi2 = sum(d^2)/hit_resolution^2 + ((R - R_ref)/(momentum_width R_ref))^2.
               momentum_width is the allowed relative momentum spread (beam momentum bite);
               momentum_width -> 0 is fixed, large is none.
  fixed        the radius is fixed to the run median momentum (p_ref = median(mom0) of the
               run): only the track-to-track curvature spread is constrained.
  charge-pair  the radius is fixed to a common |p| of the q=+ and q=- runs at the same
               (species, |p_beam|):
                 p_ref = p_beam * (w R+ + (1-w) R-),  R+- = median(mom0)/p_beam of the q=+-1 run(s),
               w = --pair-weight (0.5: pair mean, 1: q=+ beam). Runs without a partner of the
               opposite polarity are not used.

Helix frame (DstTPCHelixTracking): X = -x, Y = z - Z_TARGET, centre (helix_cx,
helix_cy), radius helix_r; mom0 is proportional to helix_r, so
R_ref = helix_r * p_ref / mom0.

Residual: r = hit - closest point on the circle, projected on the pad tangent
  r_t = r_x cos(phi) - r_z sin(phi), phi = atan2(x, z - Z_TARGET).
Update (tangential, y-independent):  new = parent + (-f * <r_t>) * e_t,
f = --update-fraction, e_t = (cos phi_pad, -sin phi_pad); centre-frame pads set to (0, 0).
--overwrite replaces the touched pads by (-f * <r_t>) * e_t instead.
--zero-layers L0-L1 keeps those layers at (0, 0) (e.g. the hit-starved inner layers under
the target holder) and excludes them from the update.

Multi-pion calibration runs (HTOF trigger, 2601-2607, pi- 1.0 GeV/c beam on target):
  --multipion FILE      : pads without a beam-through update get -f*<r_t> of the
                          inclusive helix residual (DST residual_x/z) of all tracks;
  --multipion-beam FILE : beam pi- tracks added to the constrained residual with
                          p_ref = median(run) * ref(pi 1.0)/R-(pi 1.0) of the beam-through pair,
                          i.e. the same relative bias as the beam-through pi- 1.0 beam is removed.

Usage:
  tpc_position_update.py <beamthrough_helix_hadd.root> --parent MAP --out MAP
                         [--momentum-constraint none|soft|fixed|charge-pair] [--momentum-width S]
                         [--update-fraction F] [--min-hits N] [--multipion FILE] [--zero-layers L0-L1]
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "macro"))
from compare_tpcpos_beamthrough import RUNS  # noqa: E402

Z_TARGET = -143.0
SPEC = {r.run: r for r in RUNS}
HTOF_RUNS = (2601, 2602, 2603, 2604, 2606, 2607)
HTOF_SPEC = {r: type(RUNS[0])(r, -1, 1.000, "pi") for r in HTOF_RUNS}


def mom_window(p_beam: float) -> float:
    return max(0.12, 0.20 * p_beam)


def read_map(path: Path):
    with open(path) as f:
        # lines starting with '#' are comments (allowed before the header and in the body)
        lines = (ln for ln in f if ln.strip() and not ln.lstrip().startswith("#"))
        head = next(lines).split()
        n_y, y0, dy = int(head[0]), float(head[1]), float(head[2])
        grid: Dict[Tuple[int, int], List[List[float]]] = {}
        for line in lines:
            s = line.split()
            if len(s) != 5:
                continue
            grid.setdefault((int(s[0]), int(s[1])), []).append(
                [float(s[2]), float(s[3]), float(s[4])])
    return n_y, y0, dy, grid


def write_map(path: Path, n_y, y0, dy, grid) -> None:
    with open(path, "w") as f:
        f.write(f"{n_y}  {y0:g}  {dy:g}\n")
        for key in sorted(grid):
            for y, vx, vz in grid[key]:
                f.write(f"{key[0]}   {key[1]}   {y:g}   {vx:.10g}   {vz:.10g}\n")


def run_medians(tree, step, spec_map=None) -> Dict[int, float]:
    spec_map = SPEC if spec_map is None else spec_map
    vals: Dict[int, List[np.ndarray]] = {}
    for a in tree.iterate(["run_number", "nhtrack", "helix_dz", "mom0", "charge"], step_size=step):
        runs = ak.to_numpy(a["run_number"])
        for run in np.unique(runs):
            spec = spec_map.get(int(run))
            if spec is None:
                continue
            m = runs == run
            sel = ((a["nhtrack"][m] >= 15) & (abs(a["helix_dz"][m]) < 0.05)
                   & (a["charge"][m] == spec.expect_charge)
                   & (abs(a["mom0"][m] - spec.p_beam) < mom_window(spec.p_beam)))
            vals.setdefault(int(run), []).append(ak.to_numpy(ak.flatten(a["mom0"][m][sel])))
    return {r: float(np.median(np.concatenate(v))) for r, v in vals.items()}


def reference_p(med: Dict[int, float], w: float) -> Dict[int, float]:
    pairs: Dict[Tuple[str, float], Dict[int, List[int]]] = {}
    for r in med:
        s = SPEC[r]
        pairs.setdefault((s.species, s.p_beam), {}).setdefault(s.expect_charge, []).append(r)
    pref: Dict[int, float] = {}
    for (sp, p), byq in pairs.items():
        if +1 not in byq or -1 not in byq:
            continue
        rp = np.mean([med[r] for r in byq[+1]]) / p
        rm = np.mean([med[r] for r in byq[-1]]) / p
        for runs in byq.values():
            for r in runs:
                pref[r] = p * (w * rp + (1.0 - w) * rm)
    return pref


def constrained_residual(hX, hY, cx, cy, R, n_iter=6):
    """Gauss-Newton fit of circle centre with fixed radius R (jagged hits per track).
    Returns per-hit signed distance d (|h-c|-R) and unit normal components."""
    for _ in range(n_iter):
        dx = hX - cx
        dy = hY - cy
        rr = np.sqrt(dx * dx + dy * dy)
        d = rr - R
        jx = -dx / rr
        jy = -dy / rr
        a11 = ak.sum(jx * jx, axis=-1)
        a12 = ak.sum(jx * jy, axis=-1)
        a22 = ak.sum(jy * jy, axis=-1)
        b1 = ak.sum(jx * d, axis=-1)
        b2 = ak.sum(jy * d, axis=-1)
        det = a11 * a22 - a12 * a12
        det = ak.where(abs(det) > 1e-12, det, 1e-12)
        cx = cx - (a22 * b1 - a12 * b2) / det
        cy = cy - (a11 * b2 - a12 * b1) / det
    dx = hX - cx
    dy = hY - cy
    rr = np.sqrt(dx * dx + dy * dy)
    return rr - R, dx / rr, dy / rr


def soft_residual(hX, hY, cx, cy, R_ref, sigma_p, sigma_hit, n_iter=8):
    """Gauss-Newton circle fit (cx, cy, R) with a Gaussian prior on R around R_ref
    (relative width sigma_p). Returns per-hit signed distance and unit normal components."""
    cx = ak.to_numpy(cx).astype(float)
    cy = ak.to_numpy(cy).astype(float)
    r_ref = ak.to_numpy(R_ref).astype(float)
    r = r_ref.copy()
    w_hit = 1.0 / sigma_hit ** 2
    w_r = 1.0 / (sigma_p * r_ref) ** 2
    for _ in range(n_iter):
        dx = hX - cx
        dy = hY - cy
        rr = np.sqrt(dx * dx + dy * dy)
        d = rr - r
        jx = -dx / rr
        jy = -dy / rr
        f = lambda x: ak.to_numpy(ak.sum(x, axis=-1)).astype(float)
        a = np.zeros((len(cx), 3, 3))
        a[:, 0, 0] = w_hit * f(jx * jx)
        a[:, 0, 1] = a[:, 1, 0] = w_hit * f(jx * jy)
        a[:, 1, 1] = w_hit * f(jy * jy)
        a[:, 0, 2] = a[:, 2, 0] = -w_hit * f(jx)
        a[:, 1, 2] = a[:, 2, 1] = -w_hit * f(jy)
        a[:, 2, 2] = w_hit * f(d * 0 + 1) + w_r
        b = np.stack([w_hit * f(jx * d), w_hit * f(jy * d), -w_hit * f(d) + w_r * (r - r_ref)], axis=1)
        step = np.linalg.solve(a, b[..., None])[..., 0]
        cx = cx - step[:, 0]
        cy = cy - step[:, 1]
        r = r - step[:, 2]
    dx = hX - cx
    dy = hY - cy
    rr = np.sqrt(dx * dx + dy * dy)
    return rr - r, dx / rr, dy / rr


def pad_keys(a, sel):
    ly = ak.to_numpy(ak.flatten(ak.flatten(a["hitlayer"][sel], axis=1)))
    rc = ak.to_numpy(ak.flatten(ak.flatten(a["track_cluster_row_center"][sel], axis=1)))
    return np.rint(ly).astype(np.int64) * 1000 + np.rint(rc).astype(np.int64)


def add_pad_sums(key, rt, sums, cnts):
    ok = np.isfinite(rt) & (np.abs(rt) < 5.0)
    uk, inv = np.unique(key[ok], return_inverse=True)
    s = np.bincount(inv, weights=rt[ok])
    c = np.bincount(inv)
    for k, sv, cv in zip(uk, s, c):
        sums[int(k)] = sums.get(int(k), 0.0) + float(sv)
        cnts[int(k)] = cnts.get(int(k), 0) + int(cv)


def accumulate_constrained(tree, pref, spec_map, sums, cnts, step, max_entries, soft=None) -> int:
    """Per-pad sums of the tangential residual to the radius-constrained circle."""
    br = ["run_number", "nhtrack", "helix_dz", "mom0", "charge", "helix_cx", "helix_cy",
          "helix_r", "hitlayer", "hitpos_x", "hitpos_z", "track_cluster_row_center"]
    n_trk = 0
    n_done = 0
    for a in tree.iterate(br, step_size=step, entry_stop=max_entries):
        runs = ak.to_numpy(a["run_number"])
        pr = np.array([pref.get(int(r), np.nan) for r in runs])
        pb = np.array([spec_map[int(r)].p_beam if int(r) in spec_map else np.nan for r in runs])
        qexp = np.array([spec_map[int(r)].expect_charge if int(r) in spec_map else 0 for r in runs])
        pr_t = ak.broadcast_arrays(pr, a["mom0"])[0]
        pb_t = ak.broadcast_arrays(pb, a["mom0"])[0]
        q_t = ak.broadcast_arrays(qexp, a["mom0"])[0]
        sel = ((a["nhtrack"] >= 15) & (abs(a["helix_dz"]) < 0.05) & (a["charge"] == q_t)
               & np.isfinite(pr_t) & (abs(a["mom0"] - pb_t) < np.maximum(0.12, 0.2 * pb_t)))
        cx = ak.flatten(a["helix_cx"][sel])
        cy = ak.flatten(a["helix_cy"][sel])
        R = ak.flatten(a["helix_r"][sel] * pr_t[sel] / a["mom0"][sel])
        hx = ak.flatten(a["hitpos_x"][sel], axis=1)
        hz = ak.flatten(a["hitpos_z"][sel], axis=1)
        n_trk += len(cx)
        if soft is None:
            d, nX, nY = constrained_residual(-hx, hz - Z_TARGET, cx, cy, R)
        else:
            d, nX, nY = soft_residual(-hx, hz - Z_TARGET, cx, cy, R, *soft)
        rx = ak.to_numpy(ak.flatten(-d * nX))
        rz = ak.to_numpy(ak.flatten(d * nY))
        x = ak.to_numpy(ak.flatten(hx))
        z = ak.to_numpy(ak.flatten(hz))
        phi = np.arctan2(x, z - Z_TARGET)
        add_pad_sums(pad_keys(a, sel), rx * np.cos(phi) - rz * np.sin(phi), sums, cnts)
        n_done += len(runs)
        print(f"[fit] entries {n_done}  tracks {n_trk}", flush=True)
    return n_trk


def accumulate_inclusive(tree, sums, cnts, step) -> None:
    """Per-pad sums of the inclusive helix residual (DST residual_x/z) of all tracks."""
    br = ["nhtrack", "helix_dz", "hitlayer", "residual_x", "residual_z", "hitpos_x", "hitpos_z",
          "track_cluster_row_center"]
    for a in tree.iterate(br, step_size=step):
        sel = (a["nhtrack"] >= 15) & (abs(a["helix_dz"]) < 0.05)
        rx = ak.to_numpy(ak.flatten(ak.flatten(a["residual_x"][sel], axis=1)))
        rz = ak.to_numpy(ak.flatten(ak.flatten(a["residual_z"][sel], axis=1)))
        x = ak.to_numpy(ak.flatten(ak.flatten(a["hitpos_x"][sel], axis=1)))
        z = ak.to_numpy(ak.flatten(ak.flatten(a["hitpos_z"][sel], axis=1)))
        phi = np.arctan2(x, z - Z_TARGET)
        add_pad_sums(pad_keys(a, sel), rx * np.cos(phi) - rz * np.sin(phi), sums, cnts)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("helix_root", type=Path)
    ap.add_argument("--parent", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--pad-dir", type=Path, default=Path("/group/had/sks/Users/sryuta/tmp/debug/tpcpos/cfit"),
                    help="dir with center_frame.txt ('layer row') and pad_phi.txt ('layer row theta[rad]'),"
                    " dumped from TPCPadHelper by dump_pads.C")
    ap.add_argument("--min-hits", type=int, default=100,
                    help="minimum number of hits in a pad to update it")
    ap.add_argument("--update-fraction", type=float, default=1.0,
                    help="fraction of the mean residual added to the map in this step")
    ap.add_argument("--momentum-constraint", choices=["none", "soft", "fixed", "charge-pair"],
                    default="charge-pair",
                    help="none: inclusive residual, no radius constraint; soft: run median as a Gaussian "
                    "prior (--momentum-width); fixed: radius fixed to the run median; "
                    "charge-pair: common |p| of the +/- runs (--pair-weight)")
    ap.add_argument("--momentum-width", type=float, default=0.01,
                    help="soft: allowed relative momentum spread around the run median")
    ap.add_argument("--hit-resolution", type=float, default=0.6,
                    help="soft: hit resolution [mm] used to weight the prior")
    ap.add_argument("--pair-weight", type=float, default=0.5,
                    help="charge-pair: weight w of the q=+ ratio in the reference (0.5: mean, 1: q=+)")
    ap.add_argument("--overwrite", action="store_true",
                    help="replace touched pads instead of adding to the parent map")
    ap.add_argument("--zero-layers", default="",
                    help="layer range kept at (0,0), e.g. 0-4")
    ap.add_argument("--max-entries", type=int, default=None)
    ap.add_argument("--multipion-beam", type=Path, action="append", default=[])
    ap.add_argument("--multipion", type=Path, action="append", default=[])
    ap.add_argument("--step", default="300 MB")
    args = ap.parse_args()

    tree = uproot.open(args.helix_root)["tpc"]
    med = run_medians(tree, args.step)
    pref = reference_p(med, args.pair_weight) if args.momentum_constraint == "charge-pair" else dict(med)
    soft = (args.momentum_width, args.hit_resolution) if args.momentum_constraint == "soft" else None
    if soft is not None:
        print(f"soft prior: momentum width {args.momentum_width}  hit resolution {args.hit_resolution} mm")
    print("run  q   p_beam  med/p   p_ref/p")
    for r in sorted(med):
        s = SPEC[r]
        print(f"{r:5d} {s.expect_charge:+d} {s.p_beam:6.3f} {med[r] / s.p_beam:7.4f} "
              f"{pref.get(r, float('nan')) / s.p_beam:7.4f}")

    sums: Dict[int, float] = {}
    cnts: Dict[int, int] = {}
    if args.momentum_constraint == "none":
        accumulate_inclusive(tree, sums, cnts, args.step)
        n_trk = 0
    else:
        n_trk = accumulate_constrained(tree, pref, SPEC, sums, cnts, args.step, args.max_entries, soft)
    for hf in args.multipion_beam:
        ht = uproot.open(hf)["tpc"]
        hmed = run_medians(ht, args.step, HTOF_SPEC)
        pi_bt = [r for r in pref if SPEC[r].species == "pi" and SPEC[r].p_beam == 1.0
                 and SPEC[r].expect_charge == -1]
        k = np.mean([pref[r] / med[r] for r in pi_bt]) if pi_bt else 1.0
        href = {r: m * k for r, m in hmed.items()}
        for r in sorted(hmed):
            print(f"[multipion-beam] run {r} med/p {hmed[r]:.4f} p_ref/p {href[r]:.4f} (k={k:.4f})")
        n_trk += accumulate_constrained(ht, href, HTOF_SPEC, sums, cnts, args.step, None)
    isums: Dict[int, float] = {}
    icnts: Dict[int, int] = {}
    for hf in args.multipion:
        accumulate_inclusive(uproot.open(hf)["tpc"], isums, icnts, args.step)

    n_y, y0, dy, grid = read_map(args.parent)
    cf = set()
    with open(args.pad_dir / "center_frame.txt") as f:
        for line in f:
            s = line.split()
            if len(s) >= 2:
                cf.add((int(s[0]), int(s[1])))
    n_fit = 0
    n_incl = 0
    upd = []
    for k, n in cnts.items():
        if n < args.min_hits:
            continue
        key = (k // 1000, k % 1000)
        if key in cf or key not in grid:
            continue
        mean_rt = sums[k] / n
        ct = -args.update_fraction * mean_rt
        upd.append((key, ct))
    done = {key for key, _ in upd}
    for k, n in icnts.items():
        key = (k // 1000, k % 1000)
        if n < args.min_hits or key in done or key in cf or key not in grid:
            continue
        upd.append((key, -args.update_fraction * isums[k] / n))
        n_incl += 1
    # e_t = (cos th, -sin th) with the pad angle th, as in tpc_pos_calib.
    pad_phi = {}
    with open(args.pad_dir / "pad_phi.txt") as f:
        for line in f:
            s = line.split()
            pad_phi[(int(s[0]), int(s[1]))] = float(s[2])
    zl = None
    if args.zero_layers:
        zl = tuple(int(x) for x in args.zero_layers.split("-"))
    upd = [(key, ct) for key, ct in upd if not (zl and zl[0] <= key[0] <= zl[1])]
    for key, ct in upd:
        ph = pad_phi[key]
        ex, ez = math.cos(ph), -math.sin(ph)
        for node in grid[key]:
            if args.overwrite:
                node[1] = ct * ex
                node[2] = ct * ez
            else:
                node[1] += ct * ex
                node[2] += ct * ez
        n_fit += 1
    zero = set()
    if args.zero_layers:
        lo, hi = (int(x) for x in args.zero_layers.split("-"))
        zero = {key for key in grid if lo <= key[0] <= hi}
    for key in cf | zero:
        for node in grid.get(key, []):
            node[1] = 0.0
            node[2] = 0.0
    print(f"zero-layers {args.zero_layers or '-'}: {len(zero)} pads set to (0,0)")
    write_map(args.out, n_y, y0, dy, grid)
    cts = np.array([c for _, c in upd])
    print(f"multipion-filled pads={n_incl}")
    print(f"tracks={n_trk} fitted_pads={n_fit} median|ct|={np.median(np.abs(cts)):.4f} "
          f"rms ct={np.std(cts):.4f}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
