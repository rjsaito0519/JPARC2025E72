#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare TPCHelixTracking beam-through with / without TPCPOS correction.

ON  = root/runNNNNN/runNNNNN_TPCHelix.root          (KinFit_th50_1)
OFF = scratch_root/.../runNNNNN_TPCHelix_tpcpos0.root (Map_0)

Same kinematic selection as compare_charge_beamthrough.py
(|dz|<0.05, nhtrack>=15, |mom0-p_beam|<max(0.12,0.2*p); no is_beam cut).

Output under RESULTS_IMG/run00000/bt_compare/:
  tpcpos_bt_compare.pdf / _summary.txt / _summary.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from compare_charge_beamthrough import (
    GEM_SECTIONS,
    RESULTS_IMG,
    RunSpec,
    RunStats,
    img_dir,
    mean_sigma,
    mom_window,
    select_and_fill,
    stats_to_dict,
)

DECODE_ON = Path("/group/had/sks/Users/sryuta/JPARC2025E72/root")
DECODE_OFF = Path("/group/had/sks/Users/sryuta/JPARC2025E72/scratch_root")

# All 22 runs from scratch/dst_tpchelix_02489-02537_22runs_tpcpos0.yml
RUNS: List[RunSpec] = [
    RunSpec(2489, -1, 1.000, "pi"),
    RunSpec(2491, -1, 1.000, "p"),
    RunSpec(2494, -1, 0.814, "p"),
    RunSpec(2496, -1, 0.814, "p"),
    RunSpec(2500, -1, 0.814, "K"),
    RunSpec(2502, -1, 0.814, "pi"),
    RunSpec(2504, -1, 0.645, "p"),
    RunSpec(2506, -1, 0.645, "K"),
    RunSpec(2508, -1, 0.645, "pi"),
    RunSpec(2509, -1, 0.645, "pi"),
    RunSpec(2512, -1, 0.400, "pi"),
    RunSpec(2514, -1, 0.300, "pi"),
    RunSpec(2516, +1, 1.000, "pi"),
    RunSpec(2518, +1, 1.000, "p"),
    RunSpec(2520, +1, 0.814, "pi"),
    RunSpec(2522, +1, 0.814, "p"),
    RunSpec(2524, +1, 0.645, "pi"),
    RunSpec(2527, +1, 0.645, "p"),
    RunSpec(2529, +1, 0.400, "pi"),
    RunSpec(2531, +1, 0.400, "p"),
    RunSpec(2535, +1, 0.400, "pi"),
    RunSpec(2537, +1, 0.300, "pi"),
]


def path_on(run: int) -> Path:
    return DECODE_ON / f"run{run:05d}" / f"run{run:05d}_TPCHelix.root"


def path_off(run: int) -> Path:
    return DECODE_OFF / f"run{run:05d}" / f"run{run:05d}_TPCHelix_tpcpos0.root"


def plot_summary_page(pdf: PdfPages, rows: List[dict]) -> None:
    fig, ax = plt.subplots(figsize=(14, max(4.0, 0.35 * len(rows) + 1.5)))
    ax.axis("off")
    ax.set_title("TPCPOS ON (KinFit) vs OFF (Map_0) — kinematic beam-through", fontsize=12)

    cols = [
        "run",
        "species",
        "q",
        "p_beam",
        "N_on",
        "N_off",
        "χ²μ_on",
        "χ²μ_off",
        "resμ_on",
        "resμ_off",
        "resσ_on",
        "resσ_off",
        "Δresμ",
        "Δresσ",
    ]
    cell = [cols]
    for r in rows:
        cell.append(
            [
                f"{r['run']:05d}",
                r["species"],
                f"{r['expect_charge']:+d}",
                f"{r['p_beam']:.3f}",
                f"{r['N_on']}",
                f"{r['N_off']}",
                f"{r['chisqr_mean_on']:.3g}",
                f"{r['chisqr_mean_off']:.3g}",
                f"{r['res_mean_on']:.3g}",
                f"{r['res_mean_off']:.3g}",
                f"{r['res_sigma_on']:.3g}",
                f"{r['res_sigma_off']:.3g}",
                f"{r['delta_res_mean']:.3g}",
                f"{r['delta_res_sigma']:.3g}",
            ]
        )
    table = ax.table(cellText=cell, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1.0, 1.15)
    fig.text(
        0.02,
        0.01,
        "Δ = OFF − ON.  Selection: |dz|<0.05, nh>=15, |mom0−p|<max(0.12,0.2p).  No is_beam cut.",
        fontsize=8,
    )
    pdf.savefig(fig)
    plt.close(fig)


def _overlay(ax, a, b, bins, xlabel, title, run: int) -> None:
    if a.size:
        ax.hist(a, bins=bins, histtype="step", color="C0", lw=1.5, label="ON KinFit", density=True)
    if b.size:
        ax.hist(b, bins=bins, histtype="step", color="C3", lw=1.5, label="OFF Map_0", density=True)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.legend(fontsize=8)


def plot_pair_page(pdf: PdfPages, on: RunStats, off: RunStats, spec: RunSpec) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle(f"TPCPOS ON vs OFF  {spec.label}")
    _overlay(axes[0, 0], on.chisqr, off.chisqr, np.linspace(0, 10, 81), "chisqr", "chisqr", spec.run)
    _overlay(axes[0, 1], on.nhtrack, off.nhtrack, np.arange(0, 41) - 0.5, "nhtrack", "nhtrack", spec.run)
    w = mom_window(spec.p_beam)
    _overlay(
        axes[0, 2],
        on.mom0,
        off.mom0,
        np.linspace(spec.p_beam - w, spec.p_beam + w, 61),
        "mom0 [GeV/c]",
        "mom0",
        spec.run,
    )
    _overlay(axes[1, 0], on.residual, off.residual, np.linspace(-5, 5, 101), "residual", "residual", spec.run)
    _overlay(
        axes[1, 1], on.residual_x, off.residual_x, np.linspace(-5, 5, 101), "residual_x", "residual_x", spec.run
    )
    _overlay(
        axes[1, 2], on.residual_y, off.residual_y, np.linspace(-5, 5, 101), "residual_y", "residual_y", spec.run
    )
    fig.text(
        0.01,
        0.01,
        f"ON N={on.n_tracks_sel} resμ/σ={on.mean_res:.3g}/{on.sigma_res:.3g}  |  "
        f"OFF N={off.n_tracks_sel} resμ/σ={off.mean_res:.3g}/{off.sigma_res:.3g}",
        fontsize=8,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    pdf.savefig(fig)
    plt.close(fig)


def plot_mom0_page(pdf: PdfPages, on: RunStats, off: RunStats, spec: RunSpec) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    fig.suptitle(f"TPCPOS mom0 (all sel tracks)  {spec.label}")
    w = mom_window(spec.p_beam)
    bins = np.linspace(spec.p_beam - w, spec.p_beam + w, 61)
    for ax, s, tag, color in (
        (axes[0], on, "ON KinFit", "C0"),
        (axes[1], off, "OFF Map_0", "C3"),
    ):
        if s.mom0.size:
            ax.hist(s.mom0, bins=bins, histtype="step", color=color, lw=1.5)
        ax.axvline(spec.p_beam, color="k", ls="--", lw=0.8, label=f"p_beam={spec.p_beam:.3g}")
        ax.set_xlabel("mom0 [GeV/c]")
        ax.set_title(f"{tag}  μ={s.mean_mom0:.4g} σ={s.sigma_mom0:.3g} N={s.n_tracks_sel}")
        ax.legend(fontsize=8)
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    pdf.savefig(fig)
    plt.close(fig)


def plot_layer_residual_page(pdf: PdfPages, on: RunStats, off: RunStats, spec: RunSpec) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    fig.suptitle(f"TPCPOS residual vs layer  {spec.label}")

    def mean_profile(ax, s: RunStats, color: str, label: str) -> None:
        if s.hitlayer.size == 0:
            return
        xs, means, errs = [], [], []
        for ly in np.unique(s.hitlayer.astype(int)):
            vals = s.residual[s.hitlayer.astype(int) == ly]
            if vals.size < 20:
                continue
            m, sg = mean_sigma(vals)
            xs.append(ly)
            means.append(m)
            errs.append(sg / math.sqrt(vals.size))
        if xs:
            ax.errorbar(xs, means, yerr=errs, fmt="o", color=color, label=label, ms=3)

    def sigma_profile(ax, s: RunStats, color: str, label: str) -> None:
        if s.hitlayer.size == 0:
            return
        xs, sigs = [], []
        for ly in np.unique(s.hitlayer.astype(int)):
            vals = s.residual[s.hitlayer.astype(int) == ly]
            if vals.size < 20:
                continue
            xs.append(ly)
            sigs.append(mean_sigma(vals)[1])
        if xs:
            ax.plot(xs, sigs, "o-", color=color, label=label, ms=3)

    mean_profile(axes[0], on, "C0", "ON KinFit")
    mean_profile(axes[0], off, "C3", "OFF Map_0")
    axes[0].axhline(0, color="k", lw=0.5)
    axes[0].set_xlabel("layer")
    axes[0].set_ylabel("mean residual")
    axes[0].legend(fontsize=8)

    sigma_profile(axes[1], on, "C0", "ON KinFit")
    sigma_profile(axes[1], off, "C3", "OFF Map_0")
    axes[1].set_xlabel("layer")
    axes[1].set_ylabel("σ(residual)")
    axes[1].legend(fontsize=8)
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    pdf.savefig(fig)
    plt.close(fig)


def plot_gem_section_overlay_page(pdf: PdfPages, on: RunStats, off: RunStats, spec: RunSpec) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(f"TPCPOS residual by GEM section  {spec.label}")
    bins = np.linspace(-5, 5, 101)
    for i, sec in enumerate(GEM_SECTIONS):
        ax = axes[i // 2, i % 2]
        for s, tag, color in ((on, "ON", "C0"), (off, "OFF", "C3")):
            vals = s.residual[s.section == sec]
            if vals.size:
                ax.hist(vals, bins=bins, histtype="step", color=color, lw=1.5, density=True, label=tag)
        st_on = on.by_section.get(sec, {})
        st_off = off.by_section.get(sec, {})
        ax.set_title(
            f"sec{sec}  ON N={st_on.get('n', 0)} μ={st_on.get('mean_residual', float('nan')):.3g}  |  "
            f"OFF N={st_off.get('n', 0)} μ={st_off.get('mean_residual', float('nan')):.3g}"
        )
        ax.set_xlabel("residual")
        ax.legend(fontsize=8)
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    pdf.savefig(fig)
    plt.close(fig)


def plot_gem_section_table_page(pdf: PdfPages, on: RunStats, off: RunStats, spec: RunSpec) -> None:
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.axis("off")
    ax.set_title(f"TPCPOS GEM-section summary  {spec.label}", fontsize=11)
    header = ["sec", "N_on", "N_off", "resμ_on", "resμ_off", "resσ_on", "resσ_off", "resxμ_on", "resxμ_off"]
    cell = [header]
    for sec in GEM_SECTIONS:
        a = on.by_section.get(sec, {})
        b = off.by_section.get(sec, {})
        cell.append(
            [
                str(sec),
                str(a.get("n", 0)),
                str(b.get("n", 0)),
                f"{a.get('mean_residual', float('nan')):.4g}",
                f"{b.get('mean_residual', float('nan')):.4g}",
                f"{a.get('sigma_residual', float('nan')):.4g}",
                f"{b.get('sigma_residual', float('nan')):.4g}",
                f"{a.get('mean_residual_x', float('nan')):.4g}",
                f"{b.get('mean_residual_x', float('nan')):.4g}",
            ]
        )
    table = ax.table(cellText=cell, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.1, 1.4)
    pdf.savefig(fig)
    plt.close(fig)


def process_run(spec: RunSpec) -> Tuple[RunStats, RunStats, dict]:
    print(f"[info] processing run{spec.run:05d} ({spec.label})")
    on = select_and_fill(spec, path=path_on(spec.run))
    off = select_and_fill(spec, path=path_off(spec.run))
    d_on = stats_to_dict(on)
    d_off = stats_to_dict(off)
    row = {
        "run": spec.run,
        "species": spec.species,
        "expect_charge": spec.expect_charge,
        "p_beam": spec.p_beam,
        "label": spec.label,
        "N_on": on.n_tracks_sel,
        "N_off": off.n_tracks_sel,
        "chisqr_mean_on": on.mean_chisqr,
        "chisqr_mean_off": off.mean_chisqr,
        "res_mean_on": on.mean_res,
        "res_mean_off": off.mean_res,
        "res_sigma_on": on.sigma_res,
        "res_sigma_off": off.sigma_res,
        "res_x_mean_on": on.mean_res_x,
        "res_x_mean_off": off.mean_res_x,
        "mom0_mean_on": on.mean_mom0,
        "mom0_mean_off": off.mean_mom0,
        "delta_res_mean": (
            off.mean_res - on.mean_res
            if math.isfinite(off.mean_res) and math.isfinite(on.mean_res)
            else float("nan")
        ),
        "delta_res_sigma": (
            off.sigma_res - on.sigma_res
            if math.isfinite(off.sigma_res) and math.isfinite(on.sigma_res)
            else float("nan")
        ),
        "by_section_on": d_on["by_section"],
        "by_section_off": d_off["by_section"],
        "on": d_on,
        "off": d_off,
    }
    for sec in GEM_SECTIONS:
        so = on.by_section.get(sec, {})
        sf = off.by_section.get(sec, {})
        row[f"sec{sec}_n_on"] = so.get("n", 0)
        row[f"sec{sec}_n_off"] = sf.get("n", 0)
        row[f"sec{sec}_res_mean_on"] = so.get("mean_residual", float("nan"))
        row[f"sec{sec}_res_mean_off"] = sf.get("mean_residual", float("nan"))
        row[f"sec{sec}_res_sigma_on"] = so.get("sigma_residual", float("nan"))
        row[f"sec{sec}_res_sigma_off"] = sf.get("sigma_residual", float("nan"))
    return on, off, row


def write_summary(path: Path, rows: List[dict]) -> None:
    lines = [
        "TPCPOS ON (KinFit) vs OFF (Map_0) beam-through comparison",
        "=" * 60,
        "",
        f"{'run':>6} {'sp':>3} {'q':>3} {'p':>6} {'N_on':>8} {'N_off':>8} "
        f"{'resμ_on':>10} {'resμ_off':>10} {'resσ_on':>10} {'resσ_off':>10} "
        f"{'Δresμ':>10} {'Δresσ':>10}",
    ]
    for r in rows:
        lines.append(
            f"{r['run']:6d} {r['species']:>3} {r['expect_charge']:+3d} {r['p_beam']:6.3f} "
            f"{r['N_on']:8d} {r['N_off']:8d} "
            f"{r['res_mean_on']:10.4g} {r['res_mean_off']:10.4g} "
            f"{r['res_sigma_on']:10.4g} {r['res_sigma_off']:10.4g} "
            f"{r['delta_res_mean']:10.4g} {r['delta_res_sigma']:10.4g}"
        )
    lines += [
        "",
        "notes:",
        "- ON  = root/.../runNNNNN_TPCHelix.root (TPCPOS KinFit_th50_1)",
        "- OFF = scratch_root/.../runNNNNN_TPCHelix_tpcpos0.root (TPCPOS Map_0)",
        "- Selection: |helix_dz|<0.05, nhtrack>=15, |mom0-p_beam|<max(0.12,0.2*p)",
        "- is_beam / is_accidental are NOT used as selection cuts",
        "- Δ = OFF − ON",
        "- Plots: tpcpos_bt_compare.pdf in run00000/bt_compare",
    ]
    path.write_text("\n".join(lines) + "\n")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runs",
        nargs="*",
        type=int,
        default=None,
        help="subset of run numbers (default: all 22)",
    )
    args = parser.parse_args(argv)

    specs = RUNS
    if args.runs:
        want = set(args.runs)
        specs = [s for s in RUNS if s.run in want]
        missing = want - {s.run for s in specs}
        if missing:
            print(f"[warn] unknown runs ignored: {sorted(missing)}")

    RESULTS_IMG.mkdir(parents=True, exist_ok=True)
    out00000 = img_dir(0) / "bt_compare"
    out00000.mkdir(parents=True, exist_ok=True)

    rows: List[dict] = []
    pairs: List[Tuple[RunSpec, RunStats, RunStats]] = []
    for spec in specs:
        pon, poff = path_on(spec.run), path_off(spec.run)
        if not pon.is_file() or not poff.is_file():
            print(f"[warn] skip run{spec.run:05d}: on={pon.is_file()} off={poff.is_file()}")
            continue
        on, off, row = process_run(spec)
        rows.append(row)
        pairs.append((spec, on, off))

    out_pdf = out00000 / "tpcpos_bt_compare.pdf"
    with PdfPages(out_pdf) as pdf:
        plot_summary_page(pdf, rows)
        for spec, on, off in pairs:
            plot_pair_page(pdf, on, off, spec)
            plot_mom0_page(pdf, on, off, spec)
            plot_layer_residual_page(pdf, on, off, spec)
            plot_gem_section_overlay_page(pdf, on, off, spec)
            plot_gem_section_table_page(pdf, on, off, spec)

    write_summary(out00000 / "tpcpos_bt_compare_summary.txt", rows)
    # JSON: drop nested full stats if huge — keep row dicts as built
    (out00000 / "tpcpos_bt_compare_summary.json").write_text(json.dumps(rows, indent=2) + "\n")

    n_pages = 1 + 5 * len(pairs)
    print(f"[done] combined PDF -> {out_pdf}")
    print(f"[done] compare summary -> {out00000 / 'tpcpos_bt_compare_summary.txt'}")
    print(f"[done] processed {len(pairs)} runs, {n_pages} pages")
    return 0


if __name__ == "__main__":
    sys.exit(main())
