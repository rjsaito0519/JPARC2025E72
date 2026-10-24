#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lambda -> p pi- peak from DstTPCHelixTracking (branch lambda_mass) for several TPCPOS maps.

Candidates with lambda_close_dist < 5 mm; Gaussian + linear fit in [1.09, 1.145] GeV/c^2.
Inputs are taken in the given order (e.g. the iteration order Map_0, s1, ..., s6).

Pages: peak overlay; peak shift and width vs input; peak shift vs p-pi- opening angle phi.
phi is the angle between the p and pi- momentum vectors, so it does not depend on the
momentum magnitudes; a mass that is flat in phi means consistent p and pi- momentum scales.
The phi bins are equal-population bins of the first input. chi2 is that of a common
shift over the phi bins (ndf = nbins - 1).
Usage: tpcpos_lam_eval.py LABEL=helix.root [...] [--tag NAME] [--phi-bins N]
Output: /group/had/sks/Users/sryuta/tmp/tpcpos_iter/tpcpos_lam_eval[_TAG].pdf / .txt
        (with --img the PDF goes to RESULTS_IMG/run00000/tpcpos_iter/; the .txt always stays in tmp)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import awkward as ak
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import uproot
from scipy.optimize import curve_fit

RESULTS_IMG = Path("/home/had/sryuta/JPARC2025E72/results/img")
SCRATCH_OUT = Path("/group/had/sks/Users/sryuta/tmp/tpcpos_iter")
M_L = 1.115683
LO, HI, NB = 1.09, 1.145, 55


def model(x, a, m, s, c0, c1):
    return a * np.exp(-0.5 * ((x - m) / s) ** 2) + c0 + c1 * (x - M_L)


def fit(m):
    h, e = np.histogram(m, bins=NB, range=(LO, HI))
    x = 0.5 * (e[1:] + e[:-1])
    p0 = [h.max(), M_L, 0.003, np.median(h[:5]), 0.0]
    p, cov = curve_fit(model, x, h, p0=p0, sigma=np.sqrt(np.maximum(h, 1)), maxfev=20000)
    return p, np.sqrt(np.diag(cov)), x, h


MP, MPI = 0.938272, 0.139570
BRANCHES = ["mom0", "charge", "lambda_mass", "lambda_close_dist", "lambda_mom_x", "lambda_mom_y",
            "lambda_mom_z", "lambda_track_id1", "lambda_track_id2"]


def candidates(path):
    """Lambda mass and p-pi- opening angle [rad] of the candidates with close_dist < 5 mm."""
    a = uproot.open(path)["tpc"].arrays(BRANCHES)
    sel = a["lambda_close_dist"] < 5
    i1, i2 = a["lambda_track_id1"][sel], a["lambda_track_id2"][sel]
    f = lambda x: ak.to_numpy(ak.flatten(x))
    m = f(a["lambda_mass"][sel])
    p1, p2 = f(a["mom0"][i1]), f(a["mom0"][i2])
    q1 = f(a["charge"][i1])
    big_p = f(np.sqrt(a["lambda_mom_x"][sel] ** 2 + a["lambda_mom_y"][sel] ** 2
                      + a["lambda_mom_z"][sel] ** 2))
    pp, pm = np.where(q1 > 0, p1, p2), np.where(q1 > 0, p2, p1)
    cphi = np.clip((big_p ** 2 - pp ** 2 - pm ** 2) / (2 * pp * pm), -1, 1)
    w = (m > LO) & (m < HI)
    return m[w], np.arccos(cphi[w])


def bin_peak(m):
    try:
        p, e, _, _ = fit(m)
        return p[1] - M_L, e[1]
    except (RuntimeError, ValueError):
        return np.nan, np.nan


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--tag", default="")
    ap.add_argument("--img", action="store_true",
                    help="write to results/img/run00000/tpcpos_iter (curated final only); default: scratch")
    ap.add_argument("--phi-bins", type=int, default=6)
    args = ap.parse_args()
    out = (RESULTS_IMG / "run00000" / "tpcpos_iter") if args.img else SCRATCH_OUT
    out.mkdir(parents=True, exist_ok=True)
    stem = "tpcpos_lam_eval" + (f"_{args.tag}" if args.tag else "")
    lines = ["Lambda peak, run03774 Lambda skim, close_dist<5 mm, gaus+lin fit",
             f"{'label':>12} {'N(win)':>7} {'shift [MeV]':>14} {'sigma [MeV]':>14}"]
    res, err, phib = {}, {}, {}
    edges = None
    for s in args.inputs:
        lab, path = s.split("=", 1)
        m, phi = candidates(path)
        p, e, x, h = fit(m)
        res[lab] = (p, x, h)
        err[lab] = e
        lines.append(f"{lab:>12} {len(m):7d} {1e3 * (p[1] - M_L):+7.2f}+-{1e3 * e[1]:.2f} "
                     f"{1e3 * abs(p[2]):7.2f}+-{1e3 * e[2]:.2f}")
        if edges is None:
            edges = np.quantile(phi, np.linspace(0, 1, args.phi_bins + 1))
            edges[0], edges[-1] = 0.0, np.pi
        k = np.digitize(phi, edges[1:-1])
        phib[lab] = np.array([bin_peak(m[k == i]) for i in range(args.phi_bins)])
    lines.append("")
    lines.append("peak shift [MeV] vs p-pi- opening angle phi; bin edges [rad]: "
                 + " ".join(f"{v:.3f}" for v in edges))
    lines.append(f"{'label':>12} " + " ".join(f"{'bin' + str(i):>7}" for i in range(args.phi_bins))
                 + f" {'spread':>7} {'chi2':>8}")
    chi2 = {}
    for lab, v in phib.items():
        d, e = 1e3 * v[:, 0], 1e3 * v[:, 1]
        ok = np.isfinite(d) & (e > 0)
        mean = np.sum(d[ok] / e[ok] ** 2) / np.sum(1 / e[ok] ** 2)
        chi2[lab] = float(np.sum(((d[ok] - mean) / e[ok]) ** 2))
        lines.append(f"{lab:>12} " + " ".join(f"{x:+7.2f}" for x in d)
                     + f" {np.nanmax(d) - np.nanmin(d):7.2f} {chi2[lab]:8.1f}")
    SCRATCH_OUT.mkdir(parents=True, exist_ok=True)
    (SCRATCH_OUT / f"{stem}.txt").write_text("\n".join(lines) + "\n")  # numbers stay in tmp
    print("\n".join(lines))
    with PdfPages(out / f"{stem}.pdf") as pdf:
        fig, ax = plt.subplots(figsize=(8, 5))
        xx = np.linspace(LO, HI, 300)
        for lab, (p, x, h) in res.items():
            l, = ax.step(x, h, where="mid", label=f"{lab}: {1e3 * (p[1] - M_L):+.2f} MeV, s={1e3 * abs(p[2]):.2f}")
            ax.plot(xx, model(xx, *p), c=l.get_color(), lw=0.8)
        ax.axvline(M_L, c="k", lw=0.8)
        ax.set_xlabel("m(p pi-) [GeV/c^2]")
        ax.legend(fontsize=8)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        labs = list(res)
        xi = np.arange(len(labs))
        fig, axs = plt.subplots(1, 3, figsize=(13, 4))
        axs[0].errorbar(xi, [1e3 * (res[l][0][1] - M_L) for l in labs],
                        yerr=[1e3 * err[l][1] for l in labs], fmt="o")
        axs[0].axhline(0, c="k", lw=0.8)
        axs[0].set_ylabel("peak - PDG [MeV]")
        axs[1].errorbar(xi, [1e3 * abs(res[l][0][2]) for l in labs],
                        yerr=[1e3 * err[l][2] for l in labs], fmt="o")
        axs[1].set_ylabel("sigma [MeV]")
        axs[2].plot(xi, [chi2[l] for l in labs], "o")
        axs[2].set_yscale("log")
        axs[2].set_ylabel(f"chi2 of common shift over phi (ndf {args.phi_bins - 1})")
        for ax in axs:
            ax.set_xticks(xi)
            ax.set_xticklabels(labs, rotation=45, ha="right")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 5))
        xc = 0.5 * (edges[1:] + edges[:-1])
        for j, lab in enumerate(labs):
            v = phib[lab]
            ax.errorbar(xc + (j - len(labs) / 2) * 0.004, 1e3 * v[:, 0], yerr=1e3 * v[:, 1],
                        fmt="o-", ms=3, label=f"{lab} (chi2 {chi2[lab]:.0f})")
        for b in edges[1:-1]:
            ax.axvline(b, c="0.85", lw=0.6)
        ax.axhline(0, c="k", lw=0.8)
        ax.set_xlabel("p-pi- opening angle phi [rad] (bin centre)")
        ax.set_ylabel("peak - PDG [MeV]")
        ax.legend(fontsize=7)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    print(out / f"{stem}.pdf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
