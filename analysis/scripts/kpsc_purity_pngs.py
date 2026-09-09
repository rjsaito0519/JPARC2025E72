#!/usr/bin/env python3
"""Kp scattering purity review: standalone PNG histograms via uproot only.

Mirrors a subset of the cut definitions used in
myanalysis/analysis/src/kpsc_purity_compare_pdf.cpp, but reads ROOT files
directly with uproot (no PyROOT) and writes individual high-resolution PNGs
instead of a multi-page PDF.

Outputs (into --outdir; default results/img/run{NNNNN}/ as in other scripts,
run number read from the tree's own "run_number" branch, not the filename):
  1. pid_dedx_vs_qp_nocut.png   dE/dx vs q*p, all TPC tracks with is_beam and
                                 is_accidental tracks removed (same beam/
                                 accidental definition as DstKpScattering's
                                 IsProtonCandidate/IsQMinusCandidate), no
                                 event-level cut otherwise
                                 (tree "kpsc" from DstKpScattering)
  2. pid_dedx_vs_qp_tipcut.png  dE/dx vs q*p, tip-cut selected events
                                 (K candidate + proton tracks only)
  3. mmass_kp_elastic.png       Kp elastic scattering missing mass, Kp
                                 good-event (tip) selection but NO cut on
                                 mmass itself
  4. lambda_mass.png            Lambda invariant mass, close_dist<5mm &&
                                 effective_ntTpc==2
                                 (tree "tpc" from DstTPCHelixTracking)

Usage:
  python3 kpsc_purity_pngs.py \\
      --kpsc run03778_KpScattering.root [more files...] \\
      --tpchelix run03778_TPCHelix.root [more files...] \\
      [--outdir DIR]

--kpsc and --tpchelix are independent: give either or both depending on
which plots you need. Multiple files per option are concatenated (e.g. to
combine several runs).
"""

import argparse
import sys
from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import uproot
from matplotlib.colors import LogNorm
from tqdm import tqdm

project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))

from lib import config

for k, v in config.PLOT_SETTINGS.items():
    plt.rcParams[k] = v
# PLOT_SETTINGS default font.size (22) is tuned for smaller figures; bump it
# up a bit here since these are single-panel figures.
plt.rcParams["font.size"] = 26

# --- cut / binning constants (edit here if the review criteria change) ---
# "tip" cut: same definition as kCutTip in kpsc_purity_compare_pdf.cpp.
TIP_ANGLE_MAX = 0.1        # diff_angle < 0.1 rad
TIP_CLOSE_DIST_MAX = 5.0   # close_dist < 5 mm
PID_QP_RANGE = (-1.0, 1.0)     # q*p [GeV/c]
PID_DEDX_RANGE = (0.0, 500.0)  # dE/dx [a.u.]
PID_NBINS = 150
MMASS_RANGE = (0.0, 0.8)   # GeV
MMASS_NBINS = 120
LAMBDA_CLOSE_DIST_MAX = 5.0  # mm
LAMBDA_NT_EFF = 2
LAMBDA_MASS_RANGE = (1.05, 1.25)  # GeV/c^2, display window (x-axis limits).
                                   # Confirmed on real data: no entries below
                                   # the p+pi- mass threshold (~1.078 GeV),
                                   # so the sharp left edge is real, not a
                                   # binning/margin artifact.
LAMBDA_MASS_NBINS = 200           # bins across LAMBDA_MASS_RANGE (~1 MeV/bin)
LAMBDA_MASS_CALC_MARGIN = 0.05    # GeV; extra margin histogrammed beyond the
                                   # display window on each side (same bin
                                   # width), so a tail beyond the display
                                   # window shows as the line running off the
                                   # frame instead of stopping exactly at the
                                   # window edge
# ---------------------------------------------------------------------

KPSC_BRANCHES = ["run_number", "mmass", "effective_ntTpc", "diff_angle", "close_dist",
                 "i_kcand", "i_p", "dEdx", "charge", "mom0", "is_beam", "is_accidental"]
TPCHELIX_BRANCHES = ["run_number", "lambda_mass", "lambda_close_dist", "effective_ntTpc"]


def default_outdir(run_numbers):
    """results/img/run{NNNNN}/, matching the convention used by the other
    scripts in analysis/scripts/. Run number comes from the tree's own
    "run_number" branch, not the filename (a filename is not guaranteed to
    encode the run number). When multiple runs are merged, everything is
    written into the directory of the lowest (first-taken) run number rather
    than creating a separate run-range directory per invocation."""
    runs = sorted(run_numbers)
    if not runs:
        return config.OUTPUT_DIR / "img" / "kpsc_purity_pngs"
    return config.OUTPUT_DIR / "img" / f"run{runs[0]:05d}"


def load_tree(files, treename, branches, library="ak", desc=None):
    """Read `branches` of `treename` from each file with a tqdm progress bar,
    then concatenate across files."""
    parts = [
        uproot.open(f)[treename].arrays(filter_name=branches, library=library)
        for f in tqdm(files, desc=desc or f"Reading {treename}", unit="file")
    ]
    if len(parts) == 1:
        return parts[0]
    if library == "ak":
        return ak.concatenate(parts)
    return {key: np.concatenate([p[key] for p in parts]) for key in parts[0]}


def pick_by_index(local_index, idx, *jagged_arrays):
    """Select, per event, the single jagged-array element at position `idx`.

    Equivalent to ROOT's `branch[idx]` in a TTree::Draw expression. Events
    where `idx` does not match any local index (e.g. idx < 0, "not found")
    contribute nothing, matching the intent of a missing candidate.
    """
    mask = local_index == idx
    return tuple(ak.flatten(arr[mask]) for arr in jagged_arrays)


def save_fig(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format="png", bbox_inches="tight", dpi=300, transparent=True)
    plt.close(fig)
    print(f"Wrote {path}")


def plot_pid_2d(qp, dedx, out_path):
    fig, ax = plt.subplots(figsize=(9, 7))
    h, xedges, yedges = np.histogram2d(
        ak.to_numpy(qp), ak.to_numpy(dedx),
        bins=PID_NBINS, range=[PID_QP_RANGE, PID_DEDX_RANGE],
    )
    mesh = ax.pcolormesh(
        xedges, yedges, h.T, shading="flat", cmap="viridis",
        norm=LogNorm(vmin=1, vmax=max(h.max(), 1)),
    )
    fig.colorbar(mesh, ax=ax, pad=0.01)
    ax.set_xlabel(r"$q \times p$ [GeV/$c$]")
    ax.set_ylabel("dE/dx [a.u.]")
    save_fig(fig, out_path)


def plot_1d(values, xlabel, bins, xrange, out_path, color="k", open_ends=False,
            xlim=None):
    """open_ends=True skips the vertical lines plt.hist(histtype="step") draws
    down to zero at the first/last bin, so the plot range edges don't read as
    a (nonexistent) cut boundary; the step just runs off the frame instead.

    `bins`/`xrange` is the range actually histogrammed; `xlim` (default:
    `xrange`) is the displayed window. Passing an `xrange` wider than `xlim`
    lets a tail outside the display window run off the frame for real,
    instead of the histogram simply stopping at the window edge.
    """
    fig, ax = plt.subplots(figsize=(8, 7))
    counts, edges = np.histogram(
        ak.to_numpy(values) if isinstance(values, ak.Array) else values,
        bins=bins, range=xrange,
    )
    if open_ends:
        ax.stairs(counts, edges, baseline=None, color=color, lw=1.5)
    else:
        ax.stairs(counts, edges, color=color, lw=1.5)
    ax.set_xlim(xlim if xlim is not None else xrange)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Counts")
    ax.yaxis.set_major_formatter(mticker.EngFormatter(sep=""))
    save_fig(fig, out_path)


def compute_tip_cut(arrs):
    """Kp "good event" selection (kCutTip in kpsc_purity_compare_pdf.cpp):
    Nt==2 + tight angle + close_dist, deliberately with NO cut on mmass
    itself (no mass window) so mmass distributions stay unbiased."""
    return (
        (arrs["mmass"] > 0)
        & (arrs["effective_ntTpc"] == 2)
        & (arrs["diff_angle"] < TIP_ANGLE_MAX)
        & (arrs["close_dist"] < TIP_CLOSE_DIST_MAX)
    )


def plot_pid_nocut(arrs, outdir):
    """All TPC tracks, no event-level cut, but with is_beam/is_accidental
    tracks removed (same beam/accidental condition DstKpScattering's
    IsProtonCandidate/IsQMinusCandidate use), so the "clean" PID plot isn't
    dominated by beam and accidental tracks."""
    clean = (arrs["is_beam"] == 0) & (arrs["is_accidental"] == 0)
    dedx_all = ak.flatten(arrs["dEdx"][clean])
    qp_all = ak.flatten((arrs["charge"] * arrs["mom0"])[clean])
    plot_pid_2d(qp_all, dedx_all, outdir / "pid_dedx_vs_qp_nocut.png")


def plot_pid_tipcut(arrs, outdir):
    tip_cut = compute_tip_cut(arrs)
    dEdx_c = arrs["dEdx"][tip_cut]
    qp_c = (arrs["charge"] * arrs["mom0"])[tip_cut]
    i_kcand_c = arrs["i_kcand"][tip_cut]
    i_p_c = arrs["i_p"][tip_cut]
    local_idx = ak.local_index(dEdx_c)

    dedx_k, qp_k = pick_by_index(local_idx, i_kcand_c, dEdx_c, qp_c)
    dedx_p, qp_p = pick_by_index(local_idx, i_p_c, dEdx_c, qp_c)
    dedx_sel = ak.concatenate([dedx_k, dedx_p])
    qp_sel = ak.concatenate([qp_k, qp_p])

    plot_pid_2d(qp_sel, dedx_sel, outdir / "pid_dedx_vs_qp_tipcut.png")


def plot_mmass(arrs, outdir):
    # Kp good-event selection (same tip cut as the PID plot), but no cut on
    # mmass itself so the missing-mass shape is unbiased.
    tip_cut = compute_tip_cut(arrs)
    mmass_sel = ak.to_numpy(arrs["mmass"][tip_cut])
    plot_1d(mmass_sel, "missing mass [GeV/$c^2$]",
            MMASS_NBINS, MMASS_RANGE, outdir / "mmass_kp_elastic.png")


def plot_lambda(arrs, outdir):
    lm = arrs["lambda_mass"]
    lcd = arrs["lambda_close_dist"]
    ent = arrs["effective_ntTpc"]

    cand_cut = (lcd < LAMBDA_CLOSE_DIST_MAX) & (ent == LAMBDA_NT_EFF)
    lambda_mass_sel = ak.to_numpy(ak.flatten(lm[cand_cut]))

    n_below = int(np.sum(lambda_mass_sel < LAMBDA_MASS_RANGE[0]))
    n_above = int(np.sum(lambda_mass_sel > LAMBDA_MASS_RANGE[1]))
    print(f"lambda_mass: {len(lambda_mass_sel)} candidates; "
          f"{n_below} below {LAMBDA_MASS_RANGE[0]} GeV, "
          f"{n_above} above {LAMBDA_MASS_RANGE[1]} GeV (outside display window)")

    bin_width = (LAMBDA_MASS_RANGE[1] - LAMBDA_MASS_RANGE[0]) / LAMBDA_MASS_NBINS
    calc_range = (LAMBDA_MASS_RANGE[0] - LAMBDA_MASS_CALC_MARGIN,
                  LAMBDA_MASS_RANGE[1] + LAMBDA_MASS_CALC_MARGIN)
    calc_nbins = round((calc_range[1] - calc_range[0]) / bin_width)

    plot_1d(lambda_mass_sel, r"$M(p\pi^-)$ [GeV/$c^2$]",
            calc_nbins, calc_range, outdir / "lambda_mass.png",
            color="k", open_ends=True, xlim=LAMBDA_MASS_RANGE)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--kpsc", nargs="+", default=[],
                        help="DstKpScattering output ROOT file(s) (tree 'kpsc')")
    parser.add_argument("--tpchelix", nargs="+", default=[],
                        help="DstTPCHelixTracking output ROOT file(s) (tree 'tpc')")
    parser.add_argument("--outdir", type=Path, default=None,
                        help="Output directory for PNGs "
                             "(default: <OUTPUT_DIR>/img/run{NNNNN}, run number "
                             "read from the tree; matches the other scripts here)")
    args = parser.parse_args()

    if not args.kpsc and not args.tpchelix:
        parser.error("give at least one of --kpsc or --tpchelix")

    run_numbers = set()
    kpsc_arrs = None
    tpc_arrs = None

    if args.kpsc:
        kpsc_arrs = load_tree(args.kpsc, "kpsc", KPSC_BRANCHES, library="ak",
                               desc="Reading kpsc files")
        run_numbers.update(int(r) for r in np.unique(ak.to_numpy(kpsc_arrs["run_number"])))
    else:
        print("No --kpsc file given: skipping PID / missing-mass plots.")

    if args.tpchelix:
        tpc_arrs = load_tree(args.tpchelix, "tpc", TPCHELIX_BRANCHES, library="ak",
                              desc="Reading tpc (TPCHelix) files")
        run_numbers.update(int(r) for r in np.unique(ak.to_numpy(tpc_arrs["run_number"])))
    else:
        print("No --tpchelix file given: skipping Lambda mass plot.")

    outdir = args.outdir or default_outdir(run_numbers)
    outdir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {outdir}")

    jobs = []
    if kpsc_arrs is not None:
        jobs.append(("PID dE/dx vs q×p (no cut)", lambda: plot_pid_nocut(kpsc_arrs, outdir)))
        jobs.append(("PID dE/dx vs q×p (tip cut)", lambda: plot_pid_tipcut(kpsc_arrs, outdir)))
        jobs.append(("Kp elastic missing mass", lambda: plot_mmass(kpsc_arrs, outdir)))
    if tpc_arrs is not None:
        jobs.append(("Lambda mass", lambda: plot_lambda(tpc_arrs, outdir)))

    for _name, job in tqdm(jobs, desc="Plots", unit="plot"):
        job()


if __name__ == "__main__":
    main()
