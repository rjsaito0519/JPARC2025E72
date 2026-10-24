#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""lambda_K = p_k / p_miss of K-p elastic tip-cut events for several TPCPOS maps.

Input per label: a directory with kpsc.root (DstKpScattering) and helix.root
(DstTPCHelixTracking) made from the same skim (entries aligned).
Tip cut: mmass>0, effective_ntTpc==2, diff_angle<0.1, close_dist<5, i_kcand>=0.
K track reach: max hitlayer of track i_kcand (>=31 = outermost layer).

Usage: tpcpos_lamk_eval.py LABEL=dir [LABEL=dir ...] [--tag NAME]
Output: /group/had/sks/Users/sryuta/tmp/tpcpos_iter/tpcpos_lamk_eval[_TAG].pdf / .txt
        (with --img the PDF goes to RESULTS_IMG/run00000/tpcpos_iter/; the .txt always stays in tmp)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import uproot

RESULTS_IMG = Path("/home/had/sryuta/JPARC2025E72/results/img")
SCRATCH_OUT = Path("/group/had/sks/Users/sryuta/tmp/tpcpos_iter")
KB = ["mmass", "diff_angle", "effective_ntTpc", "close_dist", "i_kcand", "p_k", "p_miss",
      "theta_k_lab", "p_p", "p_p_calc"]
TH_EDGES = np.linspace(0.1, 1.0, 10)
M_K = 0.493677


def load(d: Path):
    k = uproot.open(d / "kpsc.root")["kpsc"].arrays(KB, library="np")
    tip = ((k["mmass"] > 0) & (k["effective_ntTpc"] == 2) & (k["diff_angle"] < 0.1)
           & (k["close_dist"] < 5) & (k["i_kcand"] >= 0))
    lay = uproot.open(d / "helix.root")["tpc"]["hitlayer"].array(library="np")
    reach = np.full(len(tip), -1)
    for e in np.where(tip)[0]:
        i = int(k["i_kcand"][e])
        if e < len(lay) and i < len(lay[e]) and len(lay[e][i]) > 0:
            reach[e] = int(np.max(lay[e][i]))
    s = tip & (reach >= 0)
    return {"lam": k["p_k"][s] / k["p_miss"][s], "th": k["theta_k_lab"][s],
            "out": reach[s] >= 31, "mm": k["mmass"][s], "lp": k["p_p"][s] / k["p_p_calc"][s]}


def prof(x, y, edges):
    i = np.digitize(x, edges) - 1
    m = np.array([np.median(y[i == j]) if (i == j).sum() > 15 else np.nan for j in range(len(edges) - 1)])
    return 0.5 * (edges[1:] + edges[:-1]), m


def boot_err(y, n=200):
    if len(y) < 5:
        return np.nan
    rng = np.random.default_rng(1)
    return float(np.std([np.median(rng.choice(y, len(y))) for _ in range(n)]))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--tag", default="")
    ap.add_argument("--img", action="store_true",
                    help="write to results/img/run00000/tpcpos_iter (curated final only); default: scratch")
    args = ap.parse_args()
    labs = [s.split("=", 1)[0] for s in args.inputs]
    D = {s.split("=", 1)[0]: load(Path(s.split("=", 1)[1])) for s in args.inputs}

    out = (RESULTS_IMG / "run00000" / "tpcpos_iter") if args.img else SCRATCH_OUT
    out.mkdir(parents=True, exist_ok=True)
    stem = "tpcpos_lamk_eval" + (f"_{args.tag}" if args.tag else "")
    lines = [f"lambda_K = p_k/p_miss, mmass>0 skim, tip cut (median +- bootstrap); tag={args.tag}",
             f"{'label':>12} {'N':>5} {'lam_K all':>15} {'reach31':>15} {'not31':>15} "
             f"{'diff':>8} {'f31':>5} {'lam_p':>7} {'mmass med':>9} {'mmass IQR':>9}"]
    for lab in labs:
        d = D[lab]
        a, o, n = d["lam"], d["lam"][d["out"]], d["lam"][~d["out"]]
        q75, q25 = np.percentile(d["mm"], [75, 25])
        lines.append(f"{lab:>12} {len(a):5d} {np.median(a):.4f}+-{boot_err(a):.4f} "
                     f"{np.median(o):.4f}+-{boot_err(o):.4f} {np.median(n):.4f}+-{boot_err(n):.4f} "
                     f"{np.median(o) - np.median(n):+8.4f} {d['out'].mean():5.2f} {np.median(d['lp']):7.4f} "
                     f"{np.median(d['mm']):9.4f} {q75 - q25:9.4f}")
    lines.append("")
    lines.append("median lambda_K vs theta_K [rad]")
    x = 0.5 * (TH_EDGES[1:] + TH_EDGES[:-1])
    lines.append(f"{'theta':>6} " + " ".join(f"{lab:>12}" for lab in labs))
    profs = {lab: prof(D[lab]["th"], D[lab]["lam"], TH_EDGES)[1] for lab in labs}
    for j, xv in enumerate(x):
        lines.append(f"{xv:6.3f} " + " ".join(f"{profs[lab][j]:12.4f}" for lab in labs))
    SCRATCH_OUT.mkdir(parents=True, exist_ok=True)
    (SCRATCH_OUT / f"{stem}.txt").write_text("\n".join(lines) + "\n")  # numbers stay in tmp
    print("\n".join(lines))

    with PdfPages(out / f"{stem}.pdf") as pdf:
        fig, ax = plt.subplots(1, 3, figsize=(17, 5))
        for lab in labs:
            d = D[lab]
            ax[0].plot(x, profs[lab], "o-", ms=3, label=lab)
            for sel, ls, tag in ((d["out"], "-", "reach31"), (~d["out"], "--", "not31")):
                ax[1].plot(x, prof(d["th"][sel], d["lam"][sel], TH_EDGES)[1], "o" + ls, ms=3,
                           label=f"{lab} {tag}")
            ax[2].hist(d["mm"], bins=np.linspace(0.3, 0.7, 81), histtype="step", label=lab)
        for a_ in ax[:2]:
            a_.axhline(1, c="k", lw=0.8)
            a_.set_xlabel("theta_K lab [rad]")
            a_.set_ylabel("median lambda_K")
        ax[2].axvline(M_K, c="k", lw=0.8)
        ax[2].set_xlabel("missing mass [GeV/c^2]")
        for a_ in ax:
            a_.legend(fontsize=7)
        fig.suptitle(f"K-p elastic, mmass>0 skim, tip cut ({args.tag})")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    print(out / f"{stem}.pdf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
