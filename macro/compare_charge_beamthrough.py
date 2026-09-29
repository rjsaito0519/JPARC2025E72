#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare TPCHelixTracking beam-through performance vs charge (run polarity).

Uses existing DstTPCHelixTracking ROOT (tree ``tpc``). Selection is kinematic
(not ``is_beam``): |dz|, mom0 near beam |p|, nhtrack.

Outputs:
  - per-run PDF/summary under RESULTS_IMG/runXXXXX/
  - one multi-page comparison PDF + summary under RESULTS_IMG/run00000/
    (``charge_bt_compare.pdf``; pages = summary → overlay → mom0 (all)
    → residual-vs-layer → GEM-section residual overlay → GEM-section table, per ± pair)

GEM section (1-4) is computed offline from ``hitpos_x`` / ``hitpos_z`` with the
same rule as ``tpc::GetSection`` (TPCPadHelper.hh). No DST branch is required.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import awkward as ak
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import uproot

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DECODE_DIR = Path("/group/had/sks/Users/sryuta/JPARC2025E72/root")
RESULTS_IMG = Path("/home/had/sryuta/JPARC2025E72/results/img")
TMP_DIR = Path("/group/had/sks/Users/sryuta/tmp")

BRANCHES = [
    "ntTpc",
    "nhtrack",
    "charge",
    "chisqr",
    "mom0",
    "helix_dz",
    "residual",
    "residual_x",
    "residual_y",
    "residual_z",
    "hitlayer",
    "hitpos_x",
    "hitpos_z",
]

# GEM sections used in overlays / summary (exclude boundary = -1)
GEM_SECTIONS = (1, 2, 3, 4)


@dataclass
class RunSpec:
    run: int
    expect_charge: int  # +1 / -1 from run polarity
    p_beam: float
    species: str
    label: str = ""

    def __post_init__(self) -> None:
        if not self.label:
            sign = "+" if self.expect_charge > 0 else "-"
            self.label = f"{self.species}{sign} {self.p_beam:.3f} GeV/c run{self.run:05d}"


@dataclass
class RunStats:
    run: int
    expect_charge: int
    p_beam: float
    species: str
    n_events: int = 0
    n_tracks_all: int = 0
    n_tracks_sel: int = 0
    n_charge_match: int = 0
    n_charge_mismatch: int = 0
    mean_chisqr: float = float("nan")
    sigma_chisqr: float = float("nan")
    mean_nhtrack: float = float("nan")
    sigma_nhtrack: float = float("nan")
    mean_mom0: float = float("nan")
    sigma_mom0: float = float("nan")
    mean_res: float = float("nan")
    sigma_res: float = float("nan")
    mean_res_x: float = float("nan")
    sigma_res_x: float = float("nan")
    mean_res_y: float = float("nan")
    sigma_res_y: float = float("nan")
    mean_res_z: float = float("nan")
    sigma_res_z: float = float("nan")
    # filled arrays for overlays (subsampled if huge)
    chisqr: np.ndarray = field(default_factory=lambda: np.array([]))
    nhtrack: np.ndarray = field(default_factory=lambda: np.array([]))
    mom0: np.ndarray = field(default_factory=lambda: np.array([]))
    residual: np.ndarray = field(default_factory=lambda: np.array([]))
    residual_x: np.ndarray = field(default_factory=lambda: np.array([]))
    residual_y: np.ndarray = field(default_factory=lambda: np.array([]))
    residual_z: np.ndarray = field(default_factory=lambda: np.array([]))
    hitlayer: np.ndarray = field(default_factory=lambda: np.array([]))
    section: np.ndarray = field(default_factory=lambda: np.array([]))
    charge: np.ndarray = field(default_factory=lambda: np.array([]))
    # per GEM section: {sec: {n, mean_res, sigma_res, mean_res_x, sigma_res_x}}
    by_section: Dict[int, dict] = field(default_factory=dict)


# Phase-1 π± (prefer available files)
PHASE1: List[Tuple[RunSpec, RunSpec]] = [
    (
        RunSpec(2489, -1, 1.000, "pi"),
        RunSpec(2516, +1, 1.000, "pi"),
    ),
    (
        RunSpec(2502, -1, 0.814, "pi"),
        RunSpec(2520, +1, 0.814, "pi"),
    ),
    (
        RunSpec(2508, -1, 0.645, "pi"),
        RunSpec(2524, +1, 0.645, "pi"),
    ),
    (
        RunSpec(2512, -1, 0.400, "pi"),  # 2511 missing
        RunSpec(2529, +1, 0.400, "pi"),
    ),
    (
        RunSpec(2514, -1, 0.300, "pi"),
        RunSpec(2537, +1, 0.300, "pi"),  # 2542 missing; 2537 short
    ),
]

# Phase-2 p̄ / p
PHASE2: List[Tuple[RunSpec, RunSpec]] = [
    (
        RunSpec(2491, -1, 1.000, "p"),
        RunSpec(2518, +1, 1.000, "p"),
    ),
    (
        RunSpec(2494, -1, 0.814, "p"),
        RunSpec(2522, +1, 0.814, "p"),
    ),
    (
        RunSpec(2504, -1, 0.645, "p"),
        RunSpec(2527, +1, 0.645, "p"),
    ),
]


def helix_path(run: int) -> Path:
    return DECODE_DIR / f"run{run:05d}" / f"run{run:05d}_TPCHelix.root"


def img_dir(run: int) -> Path:
    d = RESULTS_IMG / f"run{run:05d}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def mean_sigma(a: np.ndarray) -> Tuple[float, float]:
    if a.size == 0:
        return float("nan"), float("nan")
    return float(np.mean(a)), float(np.std(a, ddof=0))


def mom_window(p_beam: float) -> float:
    """Half-width of mom0 acceptance around |p_beam|."""
    return max(0.12, 0.20 * p_beam)


def gem_section(x: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Python port of ``tpc::GetSection(x, z)`` (TPCPadHelper.hh).

    Returns section 1-4, or -1 on the |x|==|z| boundary / non-finite input.
    """
    x = np.asarray(x, dtype=np.float64)
    z = np.asarray(z, dtype=np.float64)
    out = np.full(x.shape, -1, dtype=np.int32)
    finite = np.isfinite(x) & np.isfinite(z)
    ax = np.abs(x)
    az = np.abs(z)
    on_boundary = finite & (ax == az)
    usable = finite & ~on_boundary
    # |x| < |z|
    mask = usable & (ax < az) & (z < 0)
    out[mask] = 1
    mask = usable & (ax < az) & (z > 0)
    out[mask] = 3
    # |x| > |z|
    mask = usable & (ax > az) & (x < 0)
    out[mask] = 2
    mask = usable & (ax > az) & (x > 0)
    out[mask] = 4
    return out


def section_stats(residual: np.ndarray, residual_x: np.ndarray, section: np.ndarray) -> Dict[int, dict]:
    out: Dict[int, dict] = {}
    for sec in GEM_SECTIONS:
        mask = section == sec
        r = residual[mask]
        rx = residual_x[mask]
        m_r, s_r = mean_sigma(r)
        m_x, s_x = mean_sigma(rx)
        out[sec] = {
            "n": int(r.size),
            "mean_residual": m_r,
            "sigma_residual": s_r,
            "mean_residual_x": m_x,
            "sigma_residual_x": s_x,
        }
    return out


def select_and_fill(
    spec: RunSpec,
    max_abs_dz: float = 0.05,
    min_nhit: int = 15,
    max_store: int = 200_000,
) -> RunStats:
    path = helix_path(spec.run)
    if not path.is_file():
        raise FileNotFoundError(path)

    stats = RunStats(
        run=spec.run,
        expect_charge=spec.expect_charge,
        p_beam=spec.p_beam,
        species=spec.species,
    )

    with uproot.open(path) as f:
        tree = f["tpc"]
        stats.n_events = int(tree.num_entries)
        arrays = tree.arrays(BRANCHES, library="ak")

    nh = arrays["nhtrack"]
    charge = arrays["charge"]
    chisqr = arrays["chisqr"]
    mom0 = arrays["mom0"]
    dz = arrays["helix_dz"]
    res = arrays["residual"]
    res_x = arrays["residual_x"]
    res_y = arrays["residual_y"]
    res_z = arrays["residual_z"]
    layer = arrays["hitlayer"]
    hit_x = arrays["hitpos_x"]
    hit_z = arrays["hitpos_z"]

    stats.n_tracks_all = int(ak.sum(ak.num(nh)))

    dpm = mom_window(spec.p_beam)
    # Kinematic beam-like window only. Do NOT cut on is_beam / is_accidental.
    sel = (nh >= min_nhit) & (abs(dz) < max_abs_dz) & (abs(mom0 - spec.p_beam) < dpm)

    nh_s = ak.flatten(nh[sel], axis=None)
    ch_s = ak.flatten(charge[sel], axis=None)
    chi_s = ak.flatten(chisqr[sel], axis=None)
    mom_s = ak.flatten(mom0[sel], axis=None)

    # hit-level: keep hits belonging to selected tracks
    res_s = ak.flatten(res[sel], axis=None)
    resx_s = ak.flatten(res_x[sel], axis=None)
    resy_s = ak.flatten(res_y[sel], axis=None)
    resz_s = ak.flatten(res_z[sel], axis=None)
    lay_s = ak.flatten(layer[sel], axis=None)
    hx_s = ak.flatten(hit_x[sel], axis=None)
    hz_s = ak.flatten(hit_z[sel], axis=None)

    nh_np = np.asarray(ak.to_numpy(nh_s), dtype=np.float64)
    ch_np = np.asarray(ak.to_numpy(ch_s), dtype=np.int32)
    chi_np = np.asarray(ak.to_numpy(chi_s), dtype=np.float64)
    mom_np = np.asarray(ak.to_numpy(mom_s), dtype=np.float64)
    res_np = np.asarray(ak.to_numpy(res_s), dtype=np.float64)
    resx_np = np.asarray(ak.to_numpy(resx_s), dtype=np.float64)
    resy_np = np.asarray(ak.to_numpy(resy_s), dtype=np.float64)
    resz_np = np.asarray(ak.to_numpy(resz_s), dtype=np.float64)
    lay_np = np.asarray(ak.to_numpy(lay_s), dtype=np.float64)
    hx_np = np.asarray(ak.to_numpy(hx_s), dtype=np.float64)
    hz_np = np.asarray(ak.to_numpy(hz_s), dtype=np.float64)

    # drop non-finite residuals
    finite = np.isfinite(res_np)
    res_np = res_np[finite]
    resx_np = resx_np[finite]
    resy_np = resy_np[finite]
    resz_np = resz_np[finite]
    lay_np = lay_np[finite]
    hx_np = hx_np[finite]
    hz_np = hz_np[finite]
    sec_np = gem_section(hx_np, hz_np)

    stats.n_tracks_sel = int(nh_np.size)
    stats.n_charge_match = int(np.sum(ch_np == spec.expect_charge))
    stats.n_charge_mismatch = int(np.sum(ch_np != spec.expect_charge))

    stats.mean_chisqr, stats.sigma_chisqr = mean_sigma(chi_np)
    stats.mean_nhtrack, stats.sigma_nhtrack = mean_sigma(nh_np)
    stats.mean_mom0, stats.sigma_mom0 = mean_sigma(mom_np)
    stats.mean_res, stats.sigma_res = mean_sigma(res_np)
    stats.mean_res_x, stats.sigma_res_x = mean_sigma(resx_np)
    stats.mean_res_y, stats.sigma_res_y = mean_sigma(resy_np)
    stats.mean_res_z, stats.sigma_res_z = mean_sigma(resz_np)
    stats.by_section = section_stats(res_np, resx_np, sec_np)

    # Keep hit-level arrays aligned under the same subsample indices.
    if res_np.size <= max_store:
        hit_idx = None
    else:
        hit_idx = np.linspace(0, res_np.size - 1, max_store).astype(np.int64)

    def store_hit(a: np.ndarray) -> np.ndarray:
        if hit_idx is None:
            return a
        return a[hit_idx]

    # Keep track-level arrays aligned under the same subsample indices.
    if mom_np.size <= max_store:
        trk_idx = None
    else:
        trk_idx = np.linspace(0, mom_np.size - 1, max_store).astype(np.int64)

    def store_trk(a: np.ndarray) -> np.ndarray:
        if trk_idx is None:
            return a
        return a[trk_idx]

    stats.chisqr = store_trk(chi_np)
    stats.nhtrack = store_trk(nh_np)
    stats.mom0 = store_trk(mom_np)
    stats.residual = store_hit(res_np)
    stats.residual_x = store_hit(resx_np)
    stats.residual_y = store_hit(resy_np)
    stats.residual_z = store_hit(resz_np)
    stats.hitlayer = store_hit(lay_np)
    stats.section = store_hit(sec_np)
    stats.charge = store_trk(ch_np)
    return stats


def stats_to_dict(s: RunStats) -> dict:
    match_rate = (
        s.n_charge_match / s.n_tracks_sel if s.n_tracks_sel > 0 else float("nan")
    )
    d = {
        "run": s.run,
        "species": s.species,
        "expect_charge": s.expect_charge,
        "p_beam": s.p_beam,
        "n_events": s.n_events,
        "n_tracks_all": s.n_tracks_all,
        "n_tracks_sel": s.n_tracks_sel,
        "n_charge_match": s.n_charge_match,
        "n_charge_mismatch": s.n_charge_mismatch,
        "charge_match_rate": match_rate,
        "mean_chisqr": s.mean_chisqr,
        "sigma_chisqr": s.sigma_chisqr,
        "mean_nhtrack": s.mean_nhtrack,
        "sigma_nhtrack": s.sigma_nhtrack,
        "mean_mom0": s.mean_mom0,
        "sigma_mom0": s.sigma_mom0,
        "mean_residual": s.mean_res,
        "sigma_residual": s.sigma_res,
        "mean_residual_x": s.mean_res_x,
        "sigma_residual_x": s.sigma_res_x,
        "mean_residual_y": s.mean_res_y,
        "sigma_residual_y": s.sigma_res_y,
        "mean_residual_z": s.mean_res_z,
        "sigma_residual_z": s.sigma_res_z,
        "by_section": {str(k): v for k, v in s.by_section.items()},
    }
    for sec in GEM_SECTIONS:
        st = s.by_section.get(sec, {})
        d[f"sec{sec}_n"] = st.get("n", 0)
        d[f"sec{sec}_mean_residual"] = st.get("mean_residual", float("nan"))
        d[f"sec{sec}_sigma_residual"] = st.get("sigma_residual", float("nan"))
        d[f"sec{sec}_mean_residual_x"] = st.get("mean_residual_x", float("nan"))
        d[f"sec{sec}_sigma_residual_x"] = st.get("sigma_residual_x", float("nan"))
    return d


def write_summary_txt(path: Path, title: str, rows: List[dict]) -> None:
    lines = [title, "=" * len(title), ""]
    keys = [
        "run",
        "species",
        "expect_charge",
        "p_beam",
        "n_tracks_sel",
        "charge_match_rate",
        "mean_chisqr",
        "sigma_chisqr",
        "mean_nhtrack",
        "mean_mom0",
        "sigma_mom0",
        "mean_residual",
        "sigma_residual",
        "mean_residual_x",
        "sigma_residual_x",
        "mean_residual_y",
        "sigma_residual_y",
        "mean_residual_z",
        "sigma_residual_z",
    ]
    header = "  ".join(f"{k:>16}" for k in keys)
    lines.append(header)
    for r in rows:
        lines.append(
            "  ".join(
                f"{r.get(k, float('nan')):>16.6g}"
                if isinstance(r.get(k), float)
                else f"{r.get(k):>16}"
                for k in keys
            )
        )
    path.write_text("\n".join(lines) + "\n")


def plot_run_pdf(spec: RunSpec, s: RunStats, out_pdf: Path) -> None:
    with PdfPages(out_pdf) as pdf:
        fig, axes = plt.subplots(2, 3, figsize=(14, 8))
        fig.suptitle(spec.label + " (kinematic sel, no is_beam)")

        def hist(ax, data, bins, xlabel, title):
            if data.size == 0:
                ax.set_title(title + " (empty)")
                return
            ax.hist(data, bins=bins, histtype="step", color="C0", lw=1.5)
            ax.set_xlabel(xlabel)
            ax.set_title(title)

        hist(axes[0, 0], s.chisqr, np.linspace(0, 10, 81), "chisqr", "chisqr")
        hist(axes[0, 1], s.nhtrack, np.arange(0, 41) - 0.5, "nhtrack", "nhtrack")
        hist(
            axes[0, 2],
            s.mom0,
            np.linspace(
                spec.p_beam - mom_window(spec.p_beam),
                spec.p_beam + mom_window(spec.p_beam),
                61,
            ),
            "mom0 [GeV/c]",
            "mom0",
        )
        hist(axes[1, 0], s.residual, np.linspace(-5, 5, 101), "residual", "residual")
        hist(axes[1, 1], s.residual_x, np.linspace(-5, 5, 101), "residual_x", "residual_x")
        hist(axes[1, 2], s.residual_y, np.linspace(-5, 5, 101), "residual_y", "residual_y")

        match = s.n_charge_match / s.n_tracks_sel if s.n_tracks_sel else float("nan")
        fig.text(
            0.01,
            0.01,
            f"Nsel={s.n_tracks_sel}  charge_match={match:.3f}  "
            f"res mean/sigma={s.mean_res:.3g}/{s.sigma_res:.3g}",
            fontsize=9,
        )
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        pdf.savefig(fig)
        plt.close(fig)

        # single-run GEM section residual page
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle(spec.label + " residual by GEM section")
        bins = np.linspace(-5, 5, 101)
        for i, sec in enumerate(GEM_SECTIONS):
            ax = axes[i // 2, i % 2]
            mask = s.section == sec
            vals = s.residual[mask]
            if vals.size:
                ax.hist(vals, bins=bins, histtype="step", color="C0", lw=1.5, density=True)
            st = s.by_section.get(sec, {})
            ax.set_title(
                f"sec{sec} N={st.get('n', 0)} "
                f"μ={st.get('mean_residual', float('nan')):.3g} "
                f"σ={st.get('sigma_residual', float('nan')):.3g}"
            )
            ax.set_xlabel("residual")
        fig.tight_layout(rect=[0, 0.02, 1, 0.95])
        pdf.savefig(fig)
        plt.close(fig)


def plot_pair_page(
    pdf: PdfPages,
    minus: RunStats,
    plus: RunStats,
    p_beam: float,
    species: str,
    phase_name: str,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle(
        f"{phase_name} {species} ± beam-through |p|={p_beam:.3f} GeV/c (kinematic sel)"
    )

    def overlay(ax, a, b, bins, xlabel, title):
        if a.size:
            ax.hist(
                a,
                bins=bins,
                histtype="step",
                color="C0",
                lw=1.5,
                label=f"− run{minus.run:05d}",
                density=True,
            )
        if b.size:
            ax.hist(
                b,
                bins=bins,
                histtype="step",
                color="C3",
                lw=1.5,
                label=f"+ run{plus.run:05d}",
                density=True,
            )
        ax.set_xlabel(xlabel)
        ax.set_title(title)
        ax.legend(fontsize=8)

    overlay(axes[0, 0], minus.chisqr, plus.chisqr, np.linspace(0, 10, 81), "chisqr", "chisqr")
    overlay(axes[0, 1], minus.nhtrack, plus.nhtrack, np.arange(0, 41) - 0.5, "nhtrack", "nhtrack")
    w = mom_window(p_beam)
    overlay(
        axes[0, 2],
        minus.mom0,
        plus.mom0,
        np.linspace(p_beam - w, p_beam + w, 61),
        "mom0 [GeV/c]",
        "mom0",
    )
    overlay(axes[1, 0], minus.residual, plus.residual, np.linspace(-5, 5, 101), "residual", "residual")
    overlay(
        axes[1, 1],
        minus.residual_x,
        plus.residual_x,
        np.linspace(-5, 5, 101),
        "residual_x",
        "residual_x",
    )
    overlay(
        axes[1, 2],
        minus.residual_y,
        plus.residual_y,
        np.linspace(-5, 5, 101),
        "residual_y",
        "residual_y",
    )

    def rate(s: RunStats) -> float:
        return s.n_charge_match / s.n_tracks_sel if s.n_tracks_sel else float("nan")

    fig.text(
        0.01,
        0.01,
        f"− N={minus.n_tracks_sel} match={rate(minus):.3f} resμ/σ={minus.mean_res:.3g}/{minus.sigma_res:.3g}  |  "
        f"+ N={plus.n_tracks_sel} match={rate(plus):.3f} resμ/σ={plus.mean_res:.3g}/{plus.sigma_res:.3g}",
        fontsize=8,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    pdf.savefig(fig)
    plt.close(fig)


def plot_layer_residual_page(
    pdf: PdfPages,
    minus: RunStats,
    plus: RunStats,
    p_beam: float,
    species: str,
    phase_name: str,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    fig.suptitle(
        f"{phase_name} {species} ± residual vs layer |p|={p_beam:.3f} GeV/c"
    )

    def mean_profile(ax, s: RunStats, color: str, label: str):
        if s.hitlayer.size == 0:
            return
        layers = np.unique(s.hitlayer.astype(int))
        means, errs, xs = [], [], []
        for ly in layers:
            mask = s.hitlayer.astype(int) == ly
            vals = s.residual[mask]
            if vals.size < 20:
                continue
            m, sg = mean_sigma(vals)
            means.append(m)
            errs.append(sg / math.sqrt(vals.size))
            xs.append(ly)
        if xs:
            ax.errorbar(xs, means, yerr=errs, fmt="o", color=color, label=label, ms=3)

    def sigma_profile(ax, s: RunStats, color: str, label: str):
        if s.hitlayer.size == 0:
            return
        layers = np.unique(s.hitlayer.astype(int))
        xs, sigs = [], []
        for ly in layers:
            mask = s.hitlayer.astype(int) == ly
            vals = s.residual[mask]
            if vals.size < 20:
                continue
            xs.append(ly)
            sigs.append(mean_sigma(vals)[1])
        if xs:
            ax.plot(xs, sigs, "o-", color=color, label=label, ms=3)

    mean_profile(axes[0], minus, "C0", f"− run{minus.run:05d}")
    mean_profile(axes[0], plus, "C3", f"+ run{plus.run:05d}")
    axes[0].axhline(0, color="k", lw=0.5)
    axes[0].set_xlabel("layer")
    axes[0].set_ylabel("mean residual")
    axes[0].legend(fontsize=8)
    axes[0].set_title("residual mean vs layer")

    sigma_profile(axes[1], minus, "C0", f"− run{minus.run:05d}")
    sigma_profile(axes[1], plus, "C3", f"+ run{plus.run:05d}")
    axes[1].set_xlabel("layer")
    axes[1].set_ylabel("residual sigma")
    axes[1].legend(fontsize=8)
    axes[1].set_title("residual sigma vs layer")

    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    pdf.savefig(fig)
    plt.close(fig)


def plot_summary_page(pdf: PdfPages, rows: List[dict]) -> None:
    fig = plt.figure(figsize=(14, 8))
    fig.suptitle("charge beam-through ± comparison summary")
    ax = fig.add_subplot(111)
    ax.axis("off")
    if not rows:
        ax.text(0.05, 0.9, "no pairs", fontsize=12, family="monospace")
        pdf.savefig(fig)
        plt.close(fig)
        return

    # compact table (key columns only)
    cols = [
        "phase",
        "species",
        "p_beam",
        "run_minus",
        "run_plus",
        "N_minus",
        "N_plus",
        "match_minus",
        "match_plus",
        "res_mean_minus",
        "res_mean_plus",
        "res_sigma_minus",
        "res_sigma_plus",
        "delta_res_mean",
    ]
    header = " ".join(
        f"{c:>12}" if c not in ("phase", "species") else f"{c:>8}" for c in cols
    )
    lines = [header, "-" * len(header)]
    for r in rows:
        parts = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                parts.append(f"{v:>12.4g}" if c not in ("phase", "species") else f"{v:>8.4g}")
            else:
                width = 8 if c in ("phase", "species") else 12
                parts.append(f"{v:>{width}}")
        lines.append(" ".join(parts))
    lines.append("")
    lines.append("Selection: |helix_dz|<0.05, nhtrack>=15, |mom0-p_beam|<max(0.12,0.2*p)")
    lines.append(
        "No is_beam/is_accidental cut. Pages: summary → overlay → mom0 → layer → GEM."
    )
    ax.text(
        0.02,
        0.98,
        "\n".join(lines),
        va="top",
        ha="left",
        fontsize=8,
        family="monospace",
        transform=ax.transAxes,
    )
    pdf.savefig(fig)
    plt.close(fig)


def plot_mom0_page(
    pdf: PdfPages,
    minus: RunStats,
    plus: RunStats,
    p_beam: float,
    species: str,
    phase_name: str,
) -> None:
    """mom0 in kinematic window for all selected tracks (no flag split)."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(
        f"{phase_name} {species} mom0 (all tracks) |p|={p_beam:.3f} GeV/c"
    )
    w = mom_window(p_beam)
    bins = np.linspace(p_beam - w, p_beam + w, 61)

    def draw_all(ax, s: RunStats, color: str, title: str):
        if s.mom0.size == 0:
            ax.set_title(title + " (empty)")
            return
        ax.hist(
            s.mom0,
            bins=bins,
            histtype="step",
            color=color,
            lw=1.5,
            label=f"all N={s.mom0.size}",
        )
        ax.axvline(p_beam, color="k", ls="--", lw=0.8, label=f"p_beam={p_beam:.3g}")
        ax.set_xlabel("mom0 [GeV/c]")
        ax.set_title(title)
        ax.legend(fontsize=8)

    draw_all(
        axes[0],
        minus,
        "C0",
        f"− run{minus.run:05d}  μ={minus.mean_mom0:.3g} σ={minus.sigma_mom0:.3g}",
    )
    draw_all(
        axes[1],
        plus,
        "C3",
        f"+ run{plus.run:05d}  μ={plus.mean_mom0:.3g} σ={plus.sigma_mom0:.3g}",
    )
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    pdf.savefig(fig)
    plt.close(fig)


def plot_gem_section_overlay_page(
    pdf: PdfPages,
    minus: RunStats,
    plus: RunStats,
    p_beam: float,
    species: str,
    phase_name: str,
) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    fig.suptitle(
        f"{phase_name} {species} ± residual by GEM section |p|={p_beam:.3f} GeV/c"
    )
    bins = np.linspace(-5, 5, 101)

    def draw(ax, minus_vals, plus_vals, title):
        if minus_vals.size:
            ax.hist(
                minus_vals,
                bins=bins,
                histtype="step",
                color="C0",
                lw=1.5,
                label=f"− run{minus.run:05d}",
                density=True,
            )
        if plus_vals.size:
            ax.hist(
                plus_vals,
                bins=bins,
                histtype="step",
                color="C3",
                lw=1.5,
                label=f"+ run{plus.run:05d}",
                density=True,
            )
        ax.set_title(title)
        ax.legend(fontsize=7)

    for i, sec in enumerate(GEM_SECTIONS):
        m_mask = minus.section == sec
        p_mask = plus.section == sec
        draw(
            axes[0, i],
            minus.residual[m_mask],
            plus.residual[p_mask],
            f"residual sec{sec}",
        )
        axes[0, i].set_xlabel("residual")
        draw(
            axes[1, i],
            minus.residual_x[m_mask],
            plus.residual_x[p_mask],
            f"residual_x sec{sec}",
        )
        axes[1, i].set_xlabel("residual_x")

    fig.tight_layout(rect=[0, 0.02, 1, 0.94])
    pdf.savefig(fig)
    plt.close(fig)


def plot_gem_section_table_page(
    pdf: PdfPages,
    minus: RunStats,
    plus: RunStats,
    p_beam: float,
    species: str,
    phase_name: str,
) -> None:
    fig = plt.figure(figsize=(14, 8))
    fig.suptitle(
        f"{phase_name} {species} ± GEM-section residual summary |p|={p_beam:.3f} GeV/c"
    )
    ax = fig.add_subplot(111)
    ax.axis("off")
    lines = [
        f"{'sec':>4}  {'N-':>8}  {'N+':>8}  "
        f"{'resμ-':>10}  {'resμ+':>10}  {'Δresμ':>10}  "
        f"{'resσ-':>10}  {'resσ+':>10}  "
        f"{'resxμ-':>10}  {'resxμ+':>10}  {'Δresxμ':>10}",
        "-" * 110,
    ]
    for sec in GEM_SECTIONS:
        sm = minus.by_section.get(sec, {})
        sp = plus.by_section.get(sec, {})
        m_res = sm.get("mean_residual", float("nan"))
        p_res = sp.get("mean_residual", float("nan"))
        m_sx = sm.get("sigma_residual", float("nan"))
        p_sx = sp.get("sigma_residual", float("nan"))
        m_rx = sm.get("mean_residual_x", float("nan"))
        p_rx = sp.get("mean_residual_x", float("nan"))
        d_res = (
            p_res - m_res
            if math.isfinite(p_res) and math.isfinite(m_res)
            else float("nan")
        )
        d_rx = (
            p_rx - m_rx
            if math.isfinite(p_rx) and math.isfinite(m_rx)
            else float("nan")
        )
        lines.append(
            f"{sec:>4}  {sm.get('n', 0):>8}  {sp.get('n', 0):>8}  "
            f"{m_res:>10.4g}  {p_res:>10.4g}  {d_res:>10.4g}  "
            f"{m_sx:>10.4g}  {p_sx:>10.4g}  "
            f"{m_rx:>10.4g}  {p_rx:>10.4g}  {d_rx:>10.4g}"
        )
    lines.append("")
    lines.append("Section from hitpos_x/z via tpc::GetSection rule (boundary |x|==|z| excluded).")
    ax.text(
        0.02,
        0.98,
        "\n".join(lines),
        va="top",
        ha="left",
        fontsize=9,
        family="monospace",
        transform=ax.transAxes,
    )
    pdf.savefig(fig)
    plt.close(fig)


def process_pair(
    minus_spec: RunSpec,
    plus_spec: RunSpec,
    cache: Dict[int, RunStats],
    summary_rows: List[dict],
) -> Tuple[RunStats, RunStats]:
    for spec in (minus_spec, plus_spec):
        if spec.run not in cache:
            print(f"[info] processing run{spec.run:05d} ({spec.label})")
            st = select_and_fill(spec)
            cache[spec.run] = st
            d = stats_to_dict(st)
            summary_rows.append(d)
            out_dir = img_dir(spec.run)
            plot_run_pdf(spec, st, out_dir / f"run{spec.run:05d}_charge_bt_kin.pdf")
            write_summary_txt(
                out_dir / f"run{spec.run:05d}_charge_bt_kin_summary.txt",
                f"run{spec.run:05d} charge beam-through kinematic summary",
                [d],
            )
            (out_dir / f"run{spec.run:05d}_charge_bt_kin_summary.json").write_text(
                json.dumps(d, indent=2) + "\n"
            )
        else:
            print(f"[info] reuse cache run{spec.run:05d}")
    return cache[minus_spec.run], cache[plus_spec.run]


def run_phase(
    pairs: List[Tuple[RunSpec, RunSpec]],
    phase_name: str,
    cache: Dict[int, RunStats],
    all_rows: List[dict],
) -> Tuple[List[dict], List[Tuple[str, RunStats, RunStats, float, str]]]:
    """Return compare_rows and list of (phase, minus, plus, p_beam, species) for PDF pages."""
    compare_rows: List[dict] = []
    page_specs: List[Tuple[str, RunStats, RunStats, float, str]] = []
    for minus_spec, plus_spec in pairs:
        path_m = helix_path(minus_spec.run)
        path_p = helix_path(plus_spec.run)
        if not path_m.is_file() or not path_p.is_file():
            print(
                f"[warn] skip pair {minus_spec.run}/{plus_spec.run}: "
                f"minus={path_m.is_file()} plus={path_p.is_file()}"
            )
            continue
        m, p = process_pair(minus_spec, plus_spec, cache, all_rows)
        page_specs.append((phase_name, m, p, minus_spec.p_beam, minus_spec.species))
        dm = stats_to_dict(m)
        dp = stats_to_dict(p)
        row = {
            "phase": phase_name,
            "species": minus_spec.species,
            "p_beam": minus_spec.p_beam,
            "run_minus": m.run,
            "run_plus": p.run,
            "N_minus": m.n_tracks_sel,
            "N_plus": p.n_tracks_sel,
            "match_minus": dm["charge_match_rate"],
            "match_plus": dp["charge_match_rate"],
            "chisqr_mean_minus": m.mean_chisqr,
            "chisqr_mean_plus": p.mean_chisqr,
            "res_mean_minus": m.mean_res,
            "res_mean_plus": p.mean_res,
            "res_sigma_minus": m.sigma_res,
            "res_sigma_plus": p.sigma_res,
            "res_x_mean_minus": m.mean_res_x,
            "res_x_mean_plus": p.mean_res_x,
            "mom0_mean_minus": m.mean_mom0,
            "mom0_mean_plus": p.mean_mom0,
            "delta_res_mean": (
                p.mean_res - m.mean_res
                if math.isfinite(p.mean_res) and math.isfinite(m.mean_res)
                else float("nan")
            ),
            "delta_res_sigma": (
                p.sigma_res - m.sigma_res
                if math.isfinite(p.sigma_res) and math.isfinite(m.sigma_res)
                else float("nan")
            ),
            "by_section_minus": dm["by_section"],
            "by_section_plus": dp["by_section"],
        }
        for sec in GEM_SECTIONS:
            sm = m.by_section.get(sec, {})
            sp = p.by_section.get(sec, {})
            row[f"sec{sec}_n_minus"] = sm.get("n", 0)
            row[f"sec{sec}_n_plus"] = sp.get("n", 0)
            row[f"sec{sec}_res_mean_minus"] = sm.get("mean_residual", float("nan"))
            row[f"sec{sec}_res_mean_plus"] = sp.get("mean_residual", float("nan"))
            row[f"sec{sec}_res_sigma_minus"] = sm.get("sigma_residual", float("nan"))
            row[f"sec{sec}_res_sigma_plus"] = sp.get("sigma_residual", float("nan"))
            row[f"sec{sec}_resx_mean_minus"] = sm.get("mean_residual_x", float("nan"))
            row[f"sec{sec}_resx_mean_plus"] = sp.get("mean_residual_x", float("nan"))
        compare_rows.append(row)
    return compare_rows, page_specs


def write_combined_compare_pdf(
    out_pdf: Path,
    compare_rows: List[dict],
    page_specs: List[Tuple[str, RunStats, RunStats, float, str]],
) -> None:
    with PdfPages(out_pdf) as pdf:
        plot_summary_page(pdf, compare_rows)
        for phase_name, m, p, p_beam, species in page_specs:
            plot_pair_page(pdf, m, p, p_beam, species, phase_name)
            plot_mom0_page(pdf, m, p, p_beam, species, phase_name)
            plot_layer_residual_page(pdf, m, p, p_beam, species, phase_name)
            plot_gem_section_overlay_page(pdf, m, p, p_beam, species, phase_name)
            plot_gem_section_table_page(pdf, m, p, p_beam, species, phase_name)


def cleanup_fragmented_compare_pdfs(out_dir: Path) -> None:
    """Remove old per-pair PDFs left by earlier versions of this script."""
    patterns = (
        "charge_bt_compare_phase*.pdf",
        "charge_bt_res_vs_layer_*.pdf",
    )
    for pat in patterns:
        for path in out_dir.glob(pat):
            path.unlink()
            print(f"[info] removed old fragment {path.name}")


def write_compare_summary(path: Path, rows: List[dict]) -> None:
    if not rows:
        path.write_text("no pairs\n")
        return
    # flat numeric/str columns only (skip nested by_section dicts in the text table)
    keys = [k for k in rows[0].keys() if not isinstance(rows[0][k], dict)]
    lines = ["charge beam-through ± comparison", "=" * 40, ""]
    lines.append("  ".join(f"{k:>18}" for k in keys))
    for r in rows:
        lines.append(
            "  ".join(
                f"{r[k]:>18.6g}" if isinstance(r[k], float) else f"{r[k]:>18}"
                for k in keys
            )
        )
    lines.append("")
    lines.append("notes:")
    lines.append("- Selection: |helix_dz|<0.05, nhtrack>=15, |mom0-p_beam|<max(0.12,0.2*p)")
    lines.append("- is_beam / is_accidental are NOT used as selection cuts")
    lines.append("- delta_res_mean = plus - minus; sign-flip of residual mean suggests charge-asymmetric bias")
    lines.append("- Plots: single multi-page PDF charge_bt_compare.pdf in run00000")
    lines.append("- GEM section from hitpos_x/z (tpc::GetSection); see secN_* columns / by_section_* in JSON")
    path.write_text("\n".join(lines) + "\n")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--phase",
        choices=("1", "2", "all"),
        default="all",
        help="which phase pairs to run (default: all)",
    )
    args = parser.parse_args(argv)

    RESULTS_IMG.mkdir(parents=True, exist_ok=True)
    TMP_DIR.mkdir(parents=True, exist_ok=True)

    cache: Dict[int, RunStats] = {}
    all_rows: List[dict] = []
    compare_rows: List[dict] = []
    page_specs: List[Tuple[str, RunStats, RunStats, float, str]] = []

    if args.phase in ("1", "all"):
        rows, pages = run_phase(PHASE1, "phase1", cache, all_rows)
        compare_rows.extend(rows)
        page_specs.extend(pages)
    if args.phase in ("2", "all"):
        rows, pages = run_phase(PHASE2, "phase2", cache, all_rows)
        compare_rows.extend(rows)
        page_specs.extend(pages)

    out00000 = img_dir(0)
    cleanup_fragmented_compare_pdfs(out00000)
    combined_pdf = out00000 / "charge_bt_compare.pdf"
    write_combined_compare_pdf(combined_pdf, compare_rows, page_specs)

    write_compare_summary(out00000 / "charge_bt_compare_summary.txt", compare_rows)
    (out00000 / "charge_bt_compare_summary.json").write_text(
        json.dumps(compare_rows, indent=2) + "\n"
    )
    (out00000 / "charge_bt_all_runs_summary.json").write_text(
        json.dumps(all_rows, indent=2) + "\n"
    )

    n_pages = 1 + 5 * len(page_specs)  # summary + 5 pages per pair
    print(f"[done] combined PDF -> {combined_pdf}")
    print(f"[done] compare summary -> {out00000 / 'charge_bt_compare_summary.txt'}")
    print(f"[done] processed {len(cache)} runs, {len(compare_rows)} pairs, {n_pages} pages")
    return 0


if __name__ == "__main__":
    sys.exit(main())
