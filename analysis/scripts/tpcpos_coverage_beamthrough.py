#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TPCPOS Phase 1: layer x phi x y coverage from beam-through TPCHelix (Map_0).

Input: standard ``root/runNNNNN/runNNNNN_TPCHelix.root`` (no _tpcpos0 suffix).
Selection matches compare_charge_beamthrough (kinematic; no is_beam).

Outputs under RESULTS_IMG/run00000/tpcpos_bt/coverage/:
  tpcpos_bt_coverage.pdf, _summary.txt, _summary.json
Optional HTOF calib runs (2601-2604, 2606-2607) on diagnostic pages only.
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

Z_TARGET = -143.0
N_LAYER = 32
N_PHI = 12
# Map_0 y grid centers: -250,-150,-50,50,150 with dy=100 -> bin edges between centers
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

OUT_DIR = RESULTS_IMG / "run00000" / "tpcpos_bt" / "coverage"
SIDE_DIR = Path("/group/had/sks/Users/sryuta/tmp/tpcpos_bt/coverage")  # txt/json summaries (not in results/img)


@dataclass
class BinAgg:
    n_hit: int = 0
    n_track: int = 0
    sum_qres_x: float = 0.0
    sum_qres_x2: float = 0.0
    sum_res_y: float = 0.0
    sum_res_y2: float = 0.0

    def add_hits(
        self,
        n: int,
        qres_x: np.ndarray,
        res_y: np.ndarray,
    ) -> None:
        self.n_hit += int(n)
        if n == 0:
            return
        self.sum_qres_x += float(np.sum(qres_x))
        self.sum_qres_x2 += float(np.sum(qres_x * qres_x))
        self.sum_res_y += float(np.sum(res_y))
        self.sum_res_y2 += float(np.sum(res_y * res_y))

    def mean_std_qres_x(self) -> Tuple[float, float]:
        if self.n_hit == 0:
            return float("nan"), float("nan")
        m = self.sum_qres_x / self.n_hit
        v = max(0.0, self.sum_qres_x2 / self.n_hit - m * m)
        return m, math.sqrt(v)

    def mean_std_res_y(self) -> Tuple[float, float]:
        if self.n_hit == 0:
            return float("nan"), float("nan")
        m = self.sum_res_y / self.n_hit
        v = max(0.0, self.sum_res_y2 / self.n_hit - m * m)
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


def process_run(
    spec: RunSpec,
    max_events: Optional[int],
    track_keys_per_bin: List[set],
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
    rx = arrays["residual_x"]
    ry = arrays["residual_y"]
    ev = arrays["event_number"]

    dpm = mom_window(spec.p_beam)
    sel = (nh >= 15) & (abs(dz) < 0.05) & (abs(mom0 - spec.p_beam) < dpm)

    n_tracks_sel = int(ak.sum(ak.num(nh[sel])))
    n_events = nent

    ev_per_track = ak.broadcast_arrays(ev, nh)[0][sel]
    tr_per_track = ak.local_index(nh)[sel]

    layer_sel = layer[sel]
    hx_sel = hx[sel]
    hy_sel = hy[sel]
    hz_sel = hz[sel]
    rx_sel = rx[sel]
    ry_sel = ry[sel]
    ch_per_track = charge[sel]

    ev_t = ak.flatten(ak.broadcast_arrays(ev_per_track, layer_sel)[0], axis=None)
    tr_t = ak.flatten(ak.broadcast_arrays(tr_per_track, layer_sel)[0], axis=None)
    ch_t = ak.flatten(ak.broadcast_arrays(ch_per_track, layer_sel)[0], axis=None)

    layer_f = ak.flatten(layer_sel, axis=None)
    hx_f = ak.flatten(hx_sel, axis=None)
    hy_f = ak.flatten(hy_sel, axis=None)
    hz_f = ak.flatten(hz_sel, axis=None)
    rx_f = ak.flatten(rx_sel, axis=None)
    ry_f = ak.flatten(ry_sel, axis=None)

    layer_np = ak.to_numpy(layer_f)
    hx_np = ak.to_numpy(hx_f)
    hy_np = ak.to_numpy(hy_f)
    hz_np = ak.to_numpy(hz_f)
    rx_np = ak.to_numpy(rx_f)
    ry_np = ak.to_numpy(ry_f)
    ch_np = ak.to_numpy(ch_t)
    ev_np = ak.to_numpy(ev_t)
    tr_np = ak.to_numpy(tr_t)

    phi_deg = np.degrees(np.arctan2(hx_np, hz_np - Z_TARGET)) % 360.0
    ly, ip, iy = bin_indices(layer_np, phi_deg, hy_np)
    lb = linear_bin(ly, ip, iy)
    qrx = ch_np * rx_np

    n_hit = layer_np.size
    for b in range(N_LAYER * N_PHI * N_Y):
        m = lb == b
        if not np.any(m):
            continue
        agg[b].add_hits(int(np.sum(m)), qrx[m], ry_np[m])
        keys = track_keys_per_bin[b]
        for e, t in zip(ev_np[m], tr_np[m], strict=False):
            keys.add((spec.run, int(e), int(t)))

    return {
        "run": spec.run,
        "species": spec.species,
        "p_beam": spec.p_beam,
        "expect_charge": spec.expect_charge,
        "n_events": n_events,
        "n_tracks_sel": n_tracks_sel,
        "n_hits_sel": int(n_hit),
        "path": str(path),
    }


def finalize_tracks(track_keys_per_bin: List[set], agg: List[BinAgg]) -> None:
    for b, keys in enumerate(track_keys_per_bin):
        agg[b].n_track = len(keys)


def heatmap_layer_phi(hit_grid: np.ndarray, title: str, ax: plt.Axes) -> None:
    im = ax.imshow(
        hit_grid,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap="viridis",
    )
    ax.set_xlabel("phi bin")
    ax.set_ylabel("layer")
    ax.set_title(title)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="hits")


def heatmap_layer_y(hit_grid: np.ndarray, title: str, ax: plt.Axes) -> None:
    im = ax.imshow(
        hit_grid,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap="viridis",
    )
    ax.set_xlabel("y bin")
    ax.set_ylabel("layer")
    ylabels = [f"{Y_EDGES[i]:.0f},{Y_EDGES[i+1]:.0f}" for i in range(N_Y)]
    ax.set_xticks(range(N_Y))
    ax.set_xticklabels(ylabels, rotation=45, ha="right", fontsize=7)
    ax.set_title(title)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="hits")


def build_grids(agg: List[BinAgg]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    hits = np.zeros((N_LAYER, N_PHI, N_Y), dtype=np.int64)
    tracks = np.zeros((N_LAYER, N_PHI, N_Y), dtype=np.int64)
    qrx_mean = np.full((N_LAYER, N_PHI, N_Y), np.nan)
    ry_mean = np.full((N_LAYER, N_PHI, N_Y), np.nan)
    for ly in range(N_LAYER):
        for ip in range(N_PHI):
            for iy in range(N_Y):
                b = ly * N_PHI * N_Y + ip * N_Y + iy
                hits[ly, ip, iy] = agg[b].n_hit
                tracks[ly, ip, iy] = agg[b].n_track
                m, _ = agg[b].mean_std_qres_x()
                qrx_mean[ly, ip, iy] = m
                my, _ = agg[b].mean_std_res_y()
                ry_mean[ly, ip, iy] = my
    hits_lp = hits.sum(axis=2)
    hits_ly = hits.sum(axis=1)
    return hits_lp, hits_ly, qrx_mean, tracks


def write_summary(
    path: Path,
    run_rows: List[Dict[str, object]],
    hits_lp: np.ndarray,
    hits_ly: np.ndarray,
    tracks: np.ndarray,
    label: str,
) -> None:
    total_hits = int(hits_lp.sum())
    total_tracks = int(tracks.sum())
    thin = int(np.sum(hits_lp < 100))
    lines = [
        f"TPCPOS Phase 1 coverage ({label})",
        f"total_hits_sel={total_hits}",
        f"total_track-bin entries (sum over bins)={total_tracks}",
        f"layer-phi cells with hits<100: {thin} / {N_LAYER * N_PHI}",
        "",
        "Per-run:",
    ]
    for r in run_rows:
        lines.append(
            f"  run{r['run']:05d} {r['species']} p={r['p_beam']:.3f} "
            f"events={r['n_events']} tracks_sel={r['n_tracks_sel']} hits_sel={r['n_hits_sel']}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def aggregate_runs(
    specs: List[RunSpec],
    max_events: Optional[int],
    _label: str,
) -> Tuple[List[Dict[str, object]], List[BinAgg], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n_bins = N_LAYER * N_PHI * N_Y
    agg = [BinAgg() for _ in range(n_bins)]
    track_keys = [set() for _ in range(n_bins)]
    rows: List[Dict[str, object]] = []
    for spec in specs:
        rows.append(process_run(spec, max_events, track_keys, agg))
    finalize_tracks(track_keys, agg)
    hits_lp, hits_ly, qrx_mean, tracks3 = build_grids(agg)
    return rows, agg, hits_lp, hits_ly, qrx_mean, tracks3


def main() -> int:
    parser = argparse.ArgumentParser(description="TPCPOS beam-through coverage (Phase 1)")
    parser.add_argument("--max-events", type=int, default=None, help="limit events per run (test)")
    parser.add_argument("--skip-htof", action="store_true", help="skip HTOF diagnostic runs")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    bt_rows, _, hits_lp, hits_ly, qrx_mean, tracks3 = aggregate_runs(
        BT_RUNS, args.max_events, "beam-through"
    )

    meta = {
        "phase": 1,
        "selection": "kinematic no is_beam",
        "n_phi": N_PHI,
        "y_edges_mm": Y_EDGES.tolist(),
        "z_target_mm": Z_TARGET,
        "beam_through_runs": [r.run for r in BT_RUNS],
        "per_run": bt_rows,
    }

    pdf_path = OUT_DIR / "tpcpos_bt_coverage.pdf"
    with PdfPages(pdf_path) as pdf:
        fig, ax = plt.subplots(figsize=(10, 8))
        heatmap_layer_phi(hits_lp, "Beam-through: hits (layer x phi)", ax)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(10, 8))
        heatmap_layer_y(hits_ly, "Beam-through: hits (layer x y bin)", ax)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(10, 8))
        qrx_2d = np.nanmean(qrx_mean, axis=2)
        im = ax.imshow(qrx_2d, aspect="auto", origin="lower", cmap="RdBu_r", vmin=-0.5, vmax=0.5)
        ax.set_title(
            "Mean charge*residual_x (hit-helix diag only; not Phase-2 BcOut residual)"
        )
        ax.set_xlabel("phi bin")
        ax.set_ylabel("layer")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        if not args.skip_htof:
            ht_rows, _, h_lp, h_ly, _, _ = aggregate_runs(
                HTOF_RUNS, args.max_events, "HTOF-calib diagnostic"
            )
            meta["htof_runs"] = [r.run for r in HTOF_RUNS]
            meta["htof_per_run"] = ht_rows
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            heatmap_layer_phi(h_lp, "HTOF calib (diag): hits layer x phi", axes[0])
            heatmap_layer_y(h_ly, "HTOF calib (diag): hits layer x y", axes[1])
            fig.suptitle("Diagnostic only — not used for map fit")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    SIDE_DIR.mkdir(parents=True, exist_ok=True)
    write_summary(
        SIDE_DIR / "tpcpos_bt_coverage_summary.txt",
        bt_rows,
        hits_lp,
        hits_ly,
        tracks3,
        "beam-through",
    )
    (SIDE_DIR / "tpcpos_bt_coverage_summary.json").write_text(
        json.dumps(meta, indent=2) + "\n",
        encoding="utf-8",
    )

    print(f"Wrote {pdf_path}")
    print(f"Wrote {SIDE_DIR / 'tpcpos_bt_coverage_summary.txt'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
