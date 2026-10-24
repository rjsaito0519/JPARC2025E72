#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TPCPOS Phase 2: inclusive helix residual QA from beam-through TPCHelix (Map_0).

formula_id: helix_inclusive_tangent_phi_v1
  r_t = residual_x * cos(phi) - residual_z * sin(phi)
  phi = atan2(hitpos_x, hitpos_z - Z_TARGET)
  Spec: docs/tpc/TPCPOS_RESIDUAL_HELIX.ja.md

Input: standard ``root/runNNNNN/runNNNNN_TPCHelix.root`` only (no BcOut).
Selection matches Phase 1 / compare_charge_beamthrough (kinematic; no is_beam).

Outputs under RESULTS_IMG/run00000/tpcpos_bt/residual/:
  tpcpos_bt_residual.pdf, _summary.txt, _summary.json
Optional HTOF calib runs on diagnostic pages only.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import awkward as ak
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import uproot

DECODE_DIR = Path("/group/had/sks/Users/sryuta/JPARC2025E72/root")
RESULTS_IMG = Path("/home/had/sryuta/JPARC2025E72/results/img")

FORMULA_ID = "helix_inclusive_tangent_phi_v1"
Z_TARGET = -143.0
N_LAYER = 32
N_PHI = 12
Y_EDGES = np.array([-300.0, -200.0, -100.0, 0.0, 100.0, 200.0], dtype=np.float64)
N_Y = len(Y_EDGES) - 1

BRANCHES = [
    "event_number",
    "nhtrack",
    "charge",
    "mom0",
    "helix_dz",
    "residual",
    "residual_x",
    "residual_y",
    "residual_z",
    "hitlayer",
    "hitpos_x",
    "hitpos_y",
    "hitpos_z",
    "calpos_x",
    "calpos_y",
    "calpos_z",
]


@dataclass
class RunSpec:
    run: int
    expect_charge: int
    p_beam: float
    species: str


def helix_path(run: int) -> Path:
    return DECODE_DIR / f"run{run:05d}" / f"run{run:05d}_TPCHelix.root"


def mom_window(p_beam: float) -> float:
    return max(0.12, 0.20 * p_beam)


BT_RUNS: List[RunSpec] = [
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

HTOF_RUNS: List[RunSpec] = [
    RunSpec(2601, -1, 1.000, "pi"),
    RunSpec(2602, -1, 1.000, "pi"),
    RunSpec(2603, -1, 1.000, "pi"),
    RunSpec(2604, -1, 1.000, "pi"),
    RunSpec(2606, -1, 1.000, "pi"),
    RunSpec(2607, -1, 1.000, "pi"),
]

OUT_DIR = RESULTS_IMG / "run00000" / "tpcpos_bt" / "residual"
SIDE_DIR = Path("/group/had/sks/Users/sryuta/tmp/tpcpos_bt/residual")  # txt/json summaries (not in results/img)


@dataclass
class BinAgg:
    n_hit: int = 0
    sum_qrt: float = 0.0
    sum_qrt2: float = 0.0
    sum_rt: float = 0.0
    sum_rt2: float = 0.0
    sum_qrx: float = 0.0
    sum_ry: float = 0.0
    sum_check_dx: float = 0.0  # |residual_x - (hit_x - cal_x)|

    def add_hits(
        self,
        n: int,
        qrt: np.ndarray,
        rt: np.ndarray,
        qrx: np.ndarray,
        ry: np.ndarray,
        check_dx: np.ndarray,
    ) -> None:
        self.n_hit += int(n)
        if n == 0:
            return
        self.sum_qrt += float(np.sum(qrt))
        self.sum_qrt2 += float(np.sum(qrt * qrt))
        self.sum_rt += float(np.sum(rt))
        self.sum_rt2 += float(np.sum(rt * rt))
        self.sum_qrx += float(np.sum(qrx))
        self.sum_ry += float(np.sum(ry))
        self.sum_check_dx += float(np.sum(check_dx))

    def mean_std(self, sum1: float, sum2: float) -> Tuple[float, float]:
        if self.n_hit == 0:
            return float("nan"), float("nan")
        m = sum1 / self.n_hit
        v = max(0.0, sum2 / self.n_hit - m * m)
        return m, math.sqrt(v)


def bin_indices(
    layer: np.ndarray,
    phi_deg: np.ndarray,
    y: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    ly = np.clip(layer.astype(np.int32), 0, N_LAYER - 1)
    ip = np.floor(phi_deg / (360.0 / N_PHI)).astype(np.int32)
    ip = np.clip(ip, 0, N_PHI - 1)
    iy = np.digitize(y, Y_EDGES) - 1
    iy = np.clip(iy, 0, N_Y - 1)
    return ly, ip, iy


def linear_bin(ly: np.ndarray, ip: np.ndarray, iy: np.ndarray) -> np.ndarray:
    return ly * (N_PHI * N_Y) + ip * N_Y + iy


def tangent_residual(
    hx: np.ndarray,
    hz: np.ndarray,
    rx: np.ndarray,
    rz: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (phi_deg, r_t) for formula_id helix_inclusive_tangent_phi_v1."""
    phi_rad = np.arctan2(hx, hz - Z_TARGET)
    phi_deg = np.degrees(phi_rad) % 360.0
    r_t = rx * np.cos(phi_rad) - rz * np.sin(phi_rad)
    return phi_deg, r_t


def process_run(
    spec: RunSpec,
    max_events: Optional[int],
    agg: List[BinAgg],
) -> Dict[str, object]:
    path = helix_path(spec.run)
    if not path.is_file():
        raise FileNotFoundError(path)

    with uproot.open(path) as f:
        tree = f["tpc"]
        nent = int(tree.num_entries)
        if max_events is not None:
            nent = min(nent, max_events)
        arrays = tree.arrays(BRANCHES, entry_start=0, entry_stop=nent, library="ak")

    nh = arrays["nhtrack"]
    charge = arrays["charge"]
    mom0 = arrays["mom0"]
    dz = arrays["helix_dz"]
    layer = arrays["hitlayer"]
    hx = arrays["hitpos_x"]
    hy = arrays["hitpos_y"]
    hz = arrays["hitpos_z"]
    cx = arrays["calpos_x"]
    rx = arrays["residual_x"]
    ry = arrays["residual_y"]
    rz = arrays["residual_z"]

    dpm = mom_window(spec.p_beam)
    sel = (nh >= 15) & (abs(dz) < 0.05) & (abs(mom0 - spec.p_beam) < dpm)

    n_tracks_sel = int(ak.sum(ak.num(nh[sel])))
    n_events = nent

    layer_sel = layer[sel]
    hx_sel = hx[sel]
    hy_sel = hy[sel]
    hz_sel = hz[sel]
    cx_sel = cx[sel]
    rx_sel = rx[sel]
    ry_sel = ry[sel]
    rz_sel = rz[sel]
    ch_per_track = charge[sel]

    ch_t = ak.flatten(ak.broadcast_arrays(ch_per_track, layer_sel)[0], axis=None)

    layer_np = ak.to_numpy(ak.flatten(layer_sel, axis=None))
    hx_np = ak.to_numpy(ak.flatten(hx_sel, axis=None))
    hy_np = ak.to_numpy(ak.flatten(hy_sel, axis=None))
    hz_np = ak.to_numpy(ak.flatten(hz_sel, axis=None))
    cx_np = ak.to_numpy(ak.flatten(cx_sel, axis=None))
    rx_np = ak.to_numpy(ak.flatten(rx_sel, axis=None))
    ry_np = ak.to_numpy(ak.flatten(ry_sel, axis=None))
    rz_np = ak.to_numpy(ak.flatten(rz_sel, axis=None))
    ch_np = ak.to_numpy(ch_t)

    phi_deg, r_t = tangent_residual(hx_np, hz_np, rx_np, rz_np)
    qrt = ch_np * r_t
    qrx = ch_np * rx_np
    check_dx = np.abs(rx_np - (hx_np - cx_np))

    ly, ip, iy = bin_indices(layer_np, phi_deg, hy_np)
    lb = linear_bin(ly, ip, iy)

    n_hit = layer_np.size
    for b in range(N_LAYER * N_PHI * N_Y):
        m = lb == b
        if not np.any(m):
            continue
        agg[b].add_hits(
            int(np.sum(m)),
            qrt[m],
            r_t[m],
            qrx[m],
            ry_np[m],
            check_dx[m],
        )

    return {
        "run": spec.run,
        "species": spec.species,
        "p_beam": spec.p_beam,
        "expect_charge": spec.expect_charge,
        "n_events": n_events,
        "n_tracks_sel": n_tracks_sel,
        "n_hits_sel": int(n_hit),
        "path": str(path),
        "mean_abs_resx_minus_hit_cal": float(np.mean(check_dx)) if n_hit else float("nan"),
    }


def build_grids(
    agg: List[BinAgg],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    hits = np.zeros((N_LAYER, N_PHI, N_Y), dtype=np.int64)
    qrt_mean = np.full((N_LAYER, N_PHI, N_Y), np.nan)
    rt_mean = np.full((N_LAYER, N_PHI, N_Y), np.nan)
    qrx_mean = np.full((N_LAYER, N_PHI, N_Y), np.nan)
    for ly in range(N_LAYER):
        for ip in range(N_PHI):
            for iy in range(N_Y):
                b = ly * N_PHI * N_Y + ip * N_Y + iy
                hits[ly, ip, iy] = agg[b].n_hit
                if agg[b].n_hit == 0:
                    continue
                qrt_mean[ly, ip, iy] = agg[b].sum_qrt / agg[b].n_hit
                rt_mean[ly, ip, iy] = agg[b].sum_rt / agg[b].n_hit
                qrx_mean[ly, ip, iy] = agg[b].sum_qrx / agg[b].n_hit
    hits_lp = hits.sum(axis=2)
    qrt_layer = np.full(N_LAYER, np.nan)
    rt_layer = np.full(N_LAYER, np.nan)
    for ly in range(N_LAYER):
        n = 0
        sq = 0.0
        sr = 0.0
        for ip in range(N_PHI):
            for iy in range(N_Y):
                b = ly * N_PHI * N_Y + ip * N_Y + iy
                n += agg[b].n_hit
                sq += agg[b].sum_qrt
                sr += agg[b].sum_rt
        if n > 0:
            qrt_layer[ly] = sq / n
            rt_layer[ly] = sr / n
    return hits_lp, qrt_mean, rt_mean, qrx_mean, qrt_layer, rt_layer


def aggregate_runs(
    specs: List[RunSpec],
    max_events: Optional[int],
) -> Tuple[List[Dict[str, object]], List[BinAgg], Tuple]:
    n_bins = N_LAYER * N_PHI * N_Y
    agg = [BinAgg() for _ in range(n_bins)]
    rows: List[Dict[str, object]] = []
    for spec in specs:
        rows.append(process_run(spec, max_events, agg))
    grids = build_grids(agg)
    return rows, agg, grids


def write_summary(
    path: Path,
    run_rows: List[Dict[str, object]],
    agg: List[BinAgg],
    qrt_layer: np.ndarray,
    label: str,
) -> None:
    total_hits = sum(a.n_hit for a in agg)
    check_sum = sum(a.sum_check_dx for a in agg)
    mean_check = check_sum / total_hits if total_hits else float("nan")
    lines = [
        f"TPCPOS Phase 2 residual QA ({label})",
        f"formula_id={FORMULA_ID}",
        "ref=inclusive helix (calpos); measured=hitpos after TPCPOS; residual=hit-cal",
        f"total_hits_sel={total_hits}",
        f"mean_|residual_x-(hitpos_x-calpos_x)|={mean_check:.3e} (expect ~0)",
        "",
        "Mean q*r_t vs layer:",
    ]
    for ly in range(N_LAYER):
        lines.append(f"  layer {ly:02d}: {qrt_layer[ly]:+.4f}")
    lines.append("")
    lines.append("Per-run:")
    for r in run_rows:
        lines.append(
            f"  run{r['run']:05d} {r['species']} p={r['p_beam']:.3f} "
            f"events={r['n_events']} tracks_sel={r['n_tracks_sel']} hits_sel={r['n_hits_sel']}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def heatmap(ax: plt.Axes, grid: np.ndarray, title: str, vmin=None, vmax=None, cmap="RdBu_r") -> None:
    im = ax.imshow(
        grid,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_xlabel("phi bin")
    ax.set_ylabel("layer")
    ax.set_title(title)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def main() -> int:
    parser = argparse.ArgumentParser(description="TPCPOS helix residual QA (Phase 2)")
    parser.add_argument("--max-events", type=int, default=None, help="limit events per run (test)")
    parser.add_argument("--skip-htof", action="store_true", help="skip HTOF diagnostic runs")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    bt_rows, bt_agg, grids = aggregate_runs(BT_RUNS, args.max_events)
    hits_lp, qrt_mean, rt_mean, qrx_mean, qrt_layer, rt_layer = grids

    meta = {
        "phase": 2,
        "formula_id": FORMULA_ID,
        "selection": "kinematic no is_beam",
        "reference": "inclusive helix calpos",
        "measured": "hitpos after TPCPOS",
        "residual_sign": "hit - cal (measured - reference); map update negative",
        "n_phi": N_PHI,
        "y_edges_mm": Y_EDGES.tolist(),
        "z_target_mm": Z_TARGET,
        "beam_through_runs": [r.run for r in BT_RUNS],
        "per_run": bt_rows,
        "qrt_mean_vs_layer": [None if math.isnan(x) else x for x in qrt_layer.tolist()],
    }

    pdf_path = OUT_DIR / "tpcpos_bt_residual.pdf"
    with PdfPages(pdf_path) as pdf:
        fig, ax = plt.subplots(figsize=(10, 8))
        qrt_2d = np.nanmean(qrt_mean, axis=2)
        heatmap(
            ax,
            qrt_2d,
            f"Mean q*r_t ({FORMULA_ID}); layer x phi, avg y",
            vmin=-0.5,
            vmax=0.5,
        )
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(10, 8))
        rt_2d = np.nanmean(rt_mean, axis=2)
        heatmap(
            ax,
            rt_2d,
            "Mean r_t (sign restored; map update uses -r_t); layer x phi",
            vmin=-0.5,
            vmax=0.5,
        )
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(10, 6))
        layers = np.arange(N_LAYER)
        ax.plot(layers, qrt_layer, "o-", label="mean q*r_t")
        ax.plot(layers, rt_layer, "s--", label="mean r_t")
        ax.axhline(0.0, color="k", lw=0.6)
        ax.set_xlabel("layer")
        ax.set_ylabel("mm")
        ax.set_title("Layer profile (beam-through, inclusive helix residual)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(10, 8))
        qrx_2d = np.nanmean(qrx_mean, axis=2)
        heatmap(
            ax,
            qrx_2d,
            "Phase-1-style mean q*residual_x (diag; not Phase-2 formal)",
            vmin=-0.5,
            vmax=0.5,
        )
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(10, 8))
        heatmap(ax, hits_lp.astype(float), "Hits used (layer x phi)", cmap="viridis")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        if not args.skip_htof:
            ht_rows, _, h_grids = aggregate_runs(HTOF_RUNS, args.max_events)
            h_qrt = np.nanmean(h_grids[1], axis=2)
            meta["htof_runs"] = [r.run for r in HTOF_RUNS]
            meta["htof_per_run"] = ht_rows
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            heatmap(axes[0], h_qrt, "HTOF calib (diag): mean q*r_t", vmin=-0.5, vmax=0.5)
            heatmap(axes[1], h_grids[0].astype(float), "HTOF calib (diag): hits", cmap="viridis")
            fig.suptitle("Diagnostic only — not used for map fit")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    SIDE_DIR.mkdir(parents=True, exist_ok=True)
    write_summary(
        SIDE_DIR / "tpcpos_bt_residual_summary.txt",
        bt_rows,
        bt_agg,
        qrt_layer,
        "beam-through",
    )
    (SIDE_DIR / "tpcpos_bt_residual_summary.json").write_text(
        json.dumps(meta, indent=2) + "\n",
        encoding="utf-8",
    )

    print(f"Wrote {pdf_path}")
    print(f"Wrote {SIDE_DIR / 'tpcpos_bt_residual_summary.txt'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
