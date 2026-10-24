#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TPCPOS iteration QA on beam-through TPCHelix hadd files (one file per map generation).

Usage:
  tpcpos_iter_eval.py LABEL=path/to/bt17_iterN.root [LABEL=...] [--tag NAME]
                      [--htof LABEL=path/to/htof6.root ...]
                      [--res-chunks N] [--res-chunk-size N]

  --htof: multi-pion (HTOF 2601-2607) helix, all tracks (nhtrack>=15, |dz|<0.05,
          no momentum window): inclusive residual per pad and per layer.
Momentum width per run: (q84-q16)/2 of mom0/p_beam (track-to-track curvature spread).

Per label (= map generation):
  * momentum: per run median(mom0)/p_beam for tracks with the run polarity.
    Charge-pair difference (q=+ minus q=-) at the same |p| isolates a fixed
    curvature offset (position distortion) from a common scale (B x D5).
    Fit over runs with p_beam >= 0.4:  1/p_meas - 1/p_beam = eps/p_beam + q*b
      eps : common scale-type term,  b : curvature offset [1/(GeV/c)]
  * residual (sampled entries): per-layer mean of q*r_t and per-pad |<r_t>|
      r_t = residual_x*cos(phi) - residual_z*sin(phi),
      phi = atan2(hitpos_x, hitpos_z - Z_TARGET)  (helix_inclusive_tangent_phi_v1)
Selection = compare_charge_beamthrough (|dz|<0.05, nhtrack>=15,
|mom0-p_beam|<max(0.12,0.2 p_beam); no is_beam cut).

Output: /group/had/sks/Users/sryuta/tmp/tpcpos_iter/tpcpos_iter_eval[_TAG].pdf and .txt
        (with --img the PDF goes to RESULTS_IMG/run00000/tpcpos_iter/; the .txt always stays in tmp)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List

import awkward as ak
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "macro"))
from compare_tpcpos_beamthrough import RUNS  # noqa: E402

RESULTS_IMG = Path("/home/had/sryuta/JPARC2025E72/results/img")
SCRATCH_OUT = Path("/group/had/sks/Users/sryuta/tmp/tpcpos_iter")
Z_TARGET = -143.0
N_LAYER = 32
MIN_P_FIT = 0.4
SPEC = {r.run: r for r in RUNS}


def mom_window(p_beam: float) -> float:
    return max(0.12, 0.20 * p_beam)


def load_momentum(path: Path) -> Dict[int, np.ndarray]:
    """run -> selected mom0 array (tracks with the run polarity)."""
    out: Dict[int, List[np.ndarray]] = {}
    br = ["run_number", "nhtrack", "helix_dz", "mom0", "charge"]
    for arr in uproot.open(path)["tpc"].iterate(br, step_size="400 MB"):
        runs = ak.to_numpy(arr["run_number"])
        nh, dz, mom, q = arr["nhtrack"], arr["helix_dz"], arr["mom0"], arr["charge"]
        for run in np.unique(runs):
            spec = SPEC.get(int(run))
            if spec is None:
                continue
            m = runs == run
            sel = (
                (nh[m] >= 15)
                & (abs(dz[m]) < 0.05)
                & (q[m] == spec.expect_charge)
                & (abs(mom[m] - spec.p_beam) < mom_window(spec.p_beam))
            )
            out.setdefault(int(run), []).append(ak.to_numpy(ak.flatten(mom[m][sel])))
    return {r: np.concatenate(v) for r, v in out.items()}


def load_residual(path: Path, n_chunks: int, chunk: int, all_tracks: bool = False):
    """Sampled per-hit (layer, row-free pad key not needed) q*r_t and per-pad mean r_t."""
    t = uproot.open(path)["tpc"]
    n = t.num_entries
    starts = np.linspace(0, max(0, n - chunk), n_chunks).astype(np.int64)
    br = ["run_number", "nhtrack", "helix_dz", "mom0", "charge", "hitlayer",
          "residual_x", "residual_z", "hitpos_x", "hitpos_z", "track_cluster_row_center"]
    lay_sum = np.zeros(N_LAYER)
    lay_n = np.zeros(N_LAYER)
    pad_sum: Dict[int, float] = {}
    pad_n: Dict[int, int] = {}
    rt_all = []
    for s in starts:
        a = t.arrays(br, entry_start=int(s), entry_stop=int(s + chunk))
        runs = ak.to_numpy(a["run_number"])
        pb = np.array([SPEC[int(r)].p_beam if int(r) in SPEC else np.nan for r in runs])
        pb_tr = ak.broadcast_arrays(pb, a["mom0"])[0]
        sel = (a["nhtrack"] >= 15) & (abs(a["helix_dz"]) < 0.05)
        if not all_tracks:
            sel = sel & (abs(a["mom0"] - pb_tr) < np.maximum(0.12, 0.2 * pb_tr))
        q = ak.broadcast_arrays(a["charge"][sel], a["hitlayer"][sel])[0]
        ly = ak.to_numpy(ak.flatten(a["hitlayer"][sel], axis=None)).astype(int)
        rx = ak.to_numpy(ak.flatten(a["residual_x"][sel], axis=None))
        rz = ak.to_numpy(ak.flatten(a["residual_z"][sel], axis=None))
        hx = ak.to_numpy(ak.flatten(a["hitpos_x"][sel], axis=None))
        hz = ak.to_numpy(ak.flatten(a["hitpos_z"][sel], axis=None))
        rc = np.rint(ak.to_numpy(ak.flatten(a["track_cluster_row_center"][sel], axis=None))).astype(int)
        qq = ak.to_numpy(ak.flatten(q, axis=None))
        phi = np.arctan2(hx, hz - Z_TARGET)
        rt = rx * np.cos(phi) - rz * np.sin(phi)
        ok = np.isfinite(rt) & (ly >= 0) & (ly < N_LAYER)
        ly, rt, qq, rc = ly[ok], rt[ok], qq[ok], rc[ok]
        np.add.at(lay_sum, ly, qq * rt)
        np.add.at(lay_n, ly, 1)
        key = ly * 1000 + rc
        uk, inv = np.unique(key, return_inverse=True)
        sums = np.bincount(inv, weights=rt)
        cnts = np.bincount(inv)
        for k, sm, c in zip(uk, sums, cnts):
            pad_sum[k] = pad_sum.get(k, 0.0) + sm
            pad_n[k] = pad_n.get(k, 0) + int(c)
        rt_all.append(rt[np.abs(rt) < 5])
    pad_mean = np.array([pad_sum[k] / pad_n[k] for k in pad_sum if pad_n[k] >= 50])
    rt_all = np.concatenate(rt_all)
    return lay_sum / np.maximum(lay_n, 1), pad_mean, rt_all


def fit_eps_b(med: Dict[int, float]):
    rows = [(SPEC[r].p_beam, SPEC[r].expect_charge, m) for r, m in med.items()
            if SPEC[r].p_beam >= MIN_P_FIT and np.isfinite(m)]
    if len(rows) < 3:
        return np.nan, np.nan
    p = np.array([x[0] for x in rows])
    q = np.array([x[1] for x in rows], dtype=float)
    pm = np.array([x[2] for x in rows])
    y = 1.0 / pm - 1.0 / p
    A = np.vstack([1.0 / p, q]).T
    (eps, b), *_ = np.linalg.lstsq(A, y, rcond=None)
    return -eps, b  # eps>0 means measured |p| larger than beam


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="LABEL=path")
    ap.add_argument("--tag", default="")
    ap.add_argument("--img", action="store_true",
                    help="write to results/img/run00000/tpcpos_iter (curated final only); default: scratch")
    ap.add_argument("--htof", action="append", default=[], help="LABEL=path (multi-pion)")
    ap.add_argument("--res-chunks", type=int, default=20)
    ap.add_argument("--res-chunk-size", type=int, default=5000)
    args = ap.parse_args()

    labels, paths = [], []
    for s in args.inputs:
        lab, p = s.split("=", 1)
        labels.append(lab)
        paths.append(Path(p))

    res = {}
    for lab, p in zip(labels, paths):
        print(f"[eval] {lab}: {p}", flush=True)
        mom = load_momentum(p)
        med = {r: float(np.median(v)) if len(v) > 50 else np.nan for r, v in mom.items()}
        lay, padm, rt = load_residual(p, args.res_chunks, args.res_chunk_size)
        eps, b = fit_eps_b(med)
        res[lab] = dict(mom=mom, med=med, lay=lay, padm=padm, rt=rt, eps=eps, b=b)

    out_dir = (RESULTS_IMG / "run00000" / "tpcpos_iter") if args.img else SCRATCH_OUT
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = "tpcpos_iter_eval" + (f"_{args.tag}" if args.tag else "")
    SCRATCH_OUT.mkdir(parents=True, exist_ok=True)
    pdf_path, txt_path = out_dir / f"{stem}.pdf", SCRATCH_OUT / f"{stem}.txt"  # numbers stay in tmp

    lines = []
    lines.append("TPCPOS iteration QA (beam-through, inclusive helix)")
    for lab, p in zip(labels, paths):
        lines.append(f"  {lab}: {p}")
    lines.append("")
    lines.append(f"fit (p_beam>={MIN_P_FIT}): 1/p_meas-1/p_beam = -eps/p_beam + q*b")
    lines.append(f"{'label':>12} {'eps':>9} {'b[1/GeV]':>9} {'|<rt>|pad med':>14} "
                 f"{'rt rms':>7} {'max|q*rt| layer':>15}")
    for lab in labels:
        d = res[lab]
        lines.append(f"{lab:>12} {d['eps']:+9.4f} {d['b']:+9.4f} {np.median(np.abs(d['padm'])):14.4f} "
                     f"{np.std(d['rt']):7.4f} {np.max(np.abs(d['lay'])):15.4f}")
    lines.append("")
    lines.append("per-run width (q84-q16)/2 of mom0/p_beam, median over runs: "
                 + "  ".join(f"{lab}={np.median([np.diff(np.percentile(v / SPEC[r].p_beam, [16, 84]))[0] / 2 for r, v in res[lab]['mom'].items() if len(v) > 50]):.4f}" for lab in labels))
    hres = {}
    for s_ in args.htof:
        hl, hp = s_.split("=", 1)
        print(f"[eval] htof {hl}: {hp}", flush=True)
        hres[hl] = load_residual(Path(hp), args.res_chunks, args.res_chunk_size, all_tracks=True)
    if hres:
        lines.append("")
        lines.append("multi-pion (HTOF6) inclusive residual, all tracks")
        lines.append(f"{'label':>12} {'|<rt>|pad med':>14} {'rt rms':>7} {'max|q*rt| layer':>15}")
        for hl, (lay, padm, rt) in hres.items():
            lines.append(f"{hl:>12} {np.median(np.abs(padm)):14.4f} {np.std(rt):7.4f} {np.max(np.abs(lay)):15.4f}")
    lines.append("")
    lines.append("charge-pair median(mom0)/p_beam  [(+) - (-)] per species/|p|")
    pairs: Dict[tuple, Dict[int, List[int]]] = {}
    for r in RUNS:
        pairs.setdefault((r.species, r.p_beam), {}).setdefault(r.expect_charge, []).append(r.run)
    head = f"{'species':>7} {'p':>6} " + " ".join(f"{lab:>22}" for lab in labels)
    lines.append(head)
    for (sp, p), byq in sorted(pairs.items(), key=lambda x: (x[0][0], -x[0][1])):
        if +1 not in byq or -1 not in byq:
            continue
        cells = []
        for lab in labels:
            med = res[lab]["med"]
            rp = np.nanmean([med.get(r, np.nan) for r in byq[+1]]) / p
            rm = np.nanmean([med.get(r, np.nan) for r in byq[-1]]) / p
            cells.append(f"{rp:.4f}/{rm:.4f} d={rp - rm:+.4f}")
        lines.append(f"{sp:>7} {p:6.3f} " + " ".join(f"{c:>22}" for c in cells))
    lines.append("")
    lines.append("per-run median(mom0)/p_beam")
    lines.append(f"{'run':>6} {'sp':>3} {'q':>3} {'p':>6} " + " ".join(f"{lab:>12}" for lab in labels))
    for r in RUNS:
        cells = [res[lab]["med"].get(r.run, np.nan) / r.p_beam for lab in labels]
        if all(np.isnan(cells)):
            continue
        lines.append(f"{r.run:6d} {r.species:>3} {r.expect_charge:+3d} {r.p_beam:6.3f} "
                     + " ".join(f"{c:12.4f}" for c in cells))
    txt = "\n".join(lines)
    txt_path.write_text(txt + "\n")
    print(txt)

    with PdfPages(pdf_path) as pdf:
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        for lab in labels:
            med = res[lab]["med"]
            for q, mk in ((+1, "o"), (-1, "s")):
                xs = [SPEC[r].p_beam for r in med if SPEC[r].expect_charge == q]
                ys = [med[r] / SPEC[r].p_beam for r in med if SPEC[r].expect_charge == q]
                o = np.argsort(xs)
                ax[0].plot(np.array(xs)[o], np.array(ys)[o], mk + "-", ms=4,
                           label=f"{lab} q={q:+d}")
        ax[0].axhline(1, c="k", lw=0.8)
        ax[0].set_xlabel("p_beam [GeV/c]")
        ax[0].set_ylabel("median mom0 / p_beam")
        ax[0].legend(fontsize=7)
        ax[0].set_title("beam-through momentum by run polarity")
        for lab in labels:
            med = res[lab]["med"]
            xs, ds = [], []
            for (sp, p), byq in pairs.items():
                if +1 in byq and -1 in byq:
                    rp = np.nanmean([med.get(r, np.nan) for r in byq[+1]]) / p
                    rm = np.nanmean([med.get(r, np.nan) for r in byq[-1]]) / p
                    xs.append(p)
                    ds.append(rp - rm)
            o = np.argsort(xs)
            ax[1].plot(np.array(xs)[o], np.array(ds)[o], "o-", ms=4,
                       label=f"{lab}: eps={res[lab]['eps']:+.4f} b={res[lab]['b']:+.4f}")
        ax[1].axhline(0, c="k", lw=0.8)
        ax[1].set_xlabel("p_beam [GeV/c]")
        ax[1].set_ylabel("ratio(q=+) - ratio(q=-)")
        ax[1].legend(fontsize=7)
        ax[1].set_title("charge-pair difference (position distortion)")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(1, 3, figsize=(17, 5))
        for lab in labels:
            d = res[lab]
            ax[0].plot(np.arange(N_LAYER), d["lay"], "o-", ms=3, label=lab)
            ax[1].hist(d["padm"], bins=np.linspace(-1, 1, 81), histtype="step",
                       label=f"{lab} med|.|={np.median(np.abs(d['padm'])):.3f}")
            ax[2].hist(d["rt"], bins=np.linspace(-3, 3, 121), histtype="step",
                       label=f"{lab} rms={np.std(d['rt']):.3f}")
        ax[0].axhline(0, c="k", lw=0.8)
        ax[0].set_xlabel("layer")
        ax[0].set_ylabel("mean q*r_t [mm]")
        ax[1].set_xlabel("pad <r_t> [mm] (n>=50)")
        ax[2].set_xlabel("r_t [mm]")
        ax[2].set_yscale("log")
        for a in ax:
            a.legend(fontsize=7)
        fig.suptitle("inclusive helix residual (sampled entries)")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    print(pdf_path)
    print(txt_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
