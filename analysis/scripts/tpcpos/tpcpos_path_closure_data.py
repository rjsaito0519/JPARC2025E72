#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TPCPOS Phase 3: data-path closure — hitpos difference vs Map_0.

Compares small-sample TPCHelix ROOT produced with:
  Map_0 / const_vx1 (+1,0) / ydep_vx (vx=iy).

Outputs under RESULTS_IMG/run00000/tpcpos_bt/path_closure/.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import awkward as ak
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import uproot

RESULTS_IMG = Path("/home/had/sryuta/JPARC2025E72/results/img")
OUT_DIR = RESULTS_IMG / "run00000" / "tpcpos_bt" / "path_closure"
SIDE_DIR = Path("/group/had/sks/Users/sryuta/tmp/tpcpos_bt/path_closure")  # txt/json summaries (not in results/img)
DEFAULT_RECON = Path("/group/had/sks/Users/sryuta/tmp/debug/tpcpos/recon")
STD_MAP0 = Path(
    "/group/had/sks/Users/sryuta/JPARC2025E72/root/run02489/run02489_TPCHelix.root"
)

BRANCHES = [
    "event_number",
    "nhtrack",
    "hitlayer",
    "hitpos_x",
    "hitpos_y",
    "hitpos_z",
]

Y0 = -250.0
DY = 100.0


def flatten_hits(path: Path, max_events: int | None = None) -> Dict[str, np.ndarray]:
    with uproot.open(path) as f:
        tree = f["tpc"]
        nent = int(tree.num_entries)
        if max_events is not None:
            nent = min(nent, max_events)
        arrays = tree.arrays(BRANCHES, entry_start=0, entry_stop=nent, library="ak")

    nh = arrays["nhtrack"]
    ev = arrays["event_number"]
    # keep tracks with any hits
    layer = arrays["hitlayer"]
    hx = arrays["hitpos_x"]
    hy = arrays["hitpos_y"]
    hz = arrays["hitpos_z"]

    ev_per = ak.broadcast_arrays(ev, nh)[0]
    tr_per = ak.local_index(nh)
    ev_h = ak.flatten(ak.broadcast_arrays(ev_per, layer)[0], axis=None)
    tr_h = ak.flatten(ak.broadcast_arrays(tr_per, layer)[0], axis=None)
    return {
        "event": ak.to_numpy(ev_h).astype(np.int64),
        "track": ak.to_numpy(tr_h).astype(np.int64),
        "layer": ak.to_numpy(ak.flatten(layer, axis=None)),
        "x": ak.to_numpy(ak.flatten(hx, axis=None)),
        "y": ak.to_numpy(ak.flatten(hy, axis=None)),
        "z": ak.to_numpy(ak.flatten(hz, axis=None)),
        "n_events": nent,
    }


def match_hits(
    a: Dict[str, np.ndarray],
    b: Dict[str, np.ndarray],
) -> Tuple[np.ndarray, ...]:
    """Match hits by (event, track, hit-index within track).

    Assumes identical track/hit ordering for same events (same DST binary).
    Falls back to sorting by (event, track) and zip within equal counts.
    """
    # Build per-(event,track) hit lists via stable order already present.
    key_a = a["event"] * 1_000_000 + a["track"]
    key_b = b["event"] * 1_000_000 + b["track"]
    # Index maps: for each key, list of row indices
    from collections import defaultdict

    ia: Dict[int, List[int]] = defaultdict(list)
    ib: Dict[int, List[int]] = defaultdict(list)
    for i, k in enumerate(key_a):
        ia[int(k)].append(i)
    for i, k in enumerate(key_b):
        ib[int(k)].append(i)

    idx_a: List[int] = []
    idx_b: List[int] = []
    n_key_mismatch = 0
    for k, lista in ia.items():
        listb = ib.get(k)
        if listb is None:
            n_key_mismatch += 1
            continue
        n = min(len(lista), len(listb))
        if len(lista) != len(listb):
            n_key_mismatch += 1
        idx_a.extend(lista[:n])
        idx_b.extend(listb[:n])
    return (
        np.asarray(idx_a, dtype=np.int64),
        np.asarray(idx_b, dtype=np.int64),
        n_key_mismatch,
    )


def median_mad(x: np.ndarray) -> Tuple[float, float]:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan"), float("nan")
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    return med, mad


def main() -> int:
    parser = argparse.ArgumentParser(description="TPCPOS data path closure QA")
    parser.add_argument("--recon-dir", type=Path, default=DEFAULT_RECON)
    parser.add_argument("--std-map0", type=Path, default=STD_MAP0)
    parser.add_argument("--max-events", type=int, default=None)
    args = parser.parse_args()

    paths = {
        "map0": args.recon_dir / "run02489_TPCHelix_map0.root",
        "const": args.recon_dir / "run02489_TPCHelix_const_vx1.root",
        "ydep": args.recon_dir / "run02489_TPCHelix_ydep.root",
    }
    for p in paths.values():
        if not p.is_file():
            raise FileNotFoundError(p)

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    data = {k: flatten_hits(p, args.max_events) for k, p in paths.items()}
    map0 = data["map0"]
    const = data["const"]
    ydep = data["ydep"]

    ia_c, ib_c, mism_c = match_hits(map0, const)
    ia_y, ib_y, mism_y = match_hits(map0, ydep)

    dx_c = const["x"][ib_c] - map0["x"][ia_c]
    dy_c = const["y"][ib_c] - map0["y"][ia_c]
    dz_c = const["z"][ib_c] - map0["z"][ia_c]
    y_c = map0["y"][ia_c]

    dx_y = ydep["x"][ib_y] - map0["x"][ia_y]
    dy_y = ydep["y"][ib_y] - map0["y"][ia_y]
    dz_y = ydep["z"][ib_y] - map0["z"][ia_y]
    y_y = map0["y"][ia_y]

    med_dx_c, mad_dx_c = median_mad(dx_c)
    med_dz_c, mad_dz_c = median_mad(dz_c)
    med_dy_c, _ = median_mad(dy_c)

    # Expected ydep: at lattice, vx=iy with iy=round((y-y0)/dy) for exact nodes;
    # between nodes linear interp. For #22: y < y0 -> iy1=0 -> vx=0.
    below = y_y < Y0
    med_dx_below, _ = median_mad(dx_y[below]) if np.any(below) else (float("nan"), float("nan"))
    # mid band around y=-50 (iy~2)
    mid = (y_y > -100.0) & (y_y < 0.0)
    med_dx_mid, _ = median_mad(dx_y[mid]) if np.any(mid) else (float("nan"), float("nan"))

    # Optional: compare tmp Map_0 vs standard helix (first N events)
    std_note = "skipped (file missing)"
    std_med_dx = float("nan")
    if args.std_map0.is_file():
        n_cmp = min(map0["n_events"], 100)
        std = flatten_hits(args.std_map0, n_cmp)
        # restrict map0 to first n_cmp events already
        ia_s, ib_s, mism_s = match_hits(map0, std)
        dxs = std["x"][ib_s] - map0["x"][ia_s]
        std_med_dx, _ = median_mad(dxs)
        std_note = (
            f"tmp Map_0 vs std helix: matched={len(ia_s)} key_mismatches~={mism_s} "
            f"median(dx)={std_med_dx:.4g} mm (expect ~0 if same Map_0 path)"
        )

    meta = {
        "phase": 3,
        "run": 2489,
        "n_events_tmp": map0["n_events"],
        "const_minus_map0": {
            "n_matched_hits": int(len(ia_c)),
            "key_mismatches": int(mism_c),
            "median_dx_mm": med_dx_c,
            "mad_dx_mm": mad_dx_c,
            "median_dy_mm": med_dy_c,
            "median_dz_mm": med_dz_c,
            "mad_dz_mm": mad_dz_c,
            "expect_dx_mm": 1.0,
            "expect_dz_mm": 0.0,
        },
        "ydep_minus_map0": {
            "n_matched_hits": int(len(ia_y)),
            "key_mismatches": int(mism_y),
            "median_dx_y_below_y0_mm": med_dx_below,
            "n_hits_y_below_y0": int(np.sum(below)),
            "median_dx_y_mid_m100_0_mm": med_dx_mid,
            "n_hits_y_mid": int(np.sum(mid)),
            "issue22_note": (
                "GetCorrectionVector: if y<y0 then iy1=iy2=0 (clamp to bottom node). "
                "For ydep map bottom node has vx=0, so Delta x ~0 below y0."
            ),
        },
        "std_map0_compare": std_note,
        "paths": {k: str(v) for k, v in paths.items()},
    }

    lines = [
        "TPCPOS Phase 3 data-path closure (run02489, max_loop~100)",
        f"n_events={map0['n_events']}",
        "",
        "const_vx1 - Map_0:",
        f"  matched_hits={len(ia_c)}  key_mismatches~={mism_c}",
        f"  median dx={med_dx_c:.4f} mm (expect ~+1)  mad={mad_dx_c:.4f}",
        f"  median dy={med_dy_c:.4f} mm (expect ~0)",
        f"  median dz={med_dz_c:.4f} mm (expect ~0)  mad={mad_dz_c:.4f}",
        "",
        "ydep_vx - Map_0:",
        f"  matched_hits={len(ia_y)}  key_mismatches~={mism_y}",
        f"  y < y0={Y0}: n={int(np.sum(below))} median dx={med_dx_below:.4f} mm "
        f"(#22 clamp -> expect ~0 for this map)",
        f"  y in (-100,0): n={int(np.sum(mid))} median dx={med_dx_mid:.4f} mm "
        f"(expect ~2 at lattice iy=2)",
        "",
        "Issue #22 observation:",
        "  Code clamps iy1<0 to bottom bin (no extrapolation).",
        "  With ydep map (vx=0 at y=-250), hits with y<y0 show Delta x ~ 0,",
        "  i.e. same correction as bottom node — not linear extrapolation below y0.",
        "",
        std_note,
        "",
        "G4 policy: truth only; no TPCPOS / no intentional degradation (Phase 3).",
    ]
    SIDE_DIR.mkdir(parents=True, exist_ok=True)
    (SIDE_DIR / "tpcpos_path_closure_summary.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    (SIDE_DIR / "tpcpos_path_closure_summary.json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8"
    )

    pdf_path = OUT_DIR / "tpcpos_path_closure.pdf"
    with PdfPages(pdf_path) as pdf:
        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        axes[0].hist(dx_c, bins=80, range=(-0.5, 2.5), histtype="step", color="C0")
        axes[0].axvline(1.0, color="k", ls="--", lw=0.8)
        axes[0].set_title(f"const-Map0 dx  med={med_dx_c:.3f}")
        axes[0].set_xlabel("dx [mm]")
        axes[1].hist(dz_c, bins=80, range=(-1, 1), histtype="step", color="C1")
        axes[1].axvline(0.0, color="k", ls="--", lw=0.8)
        axes[1].set_title(f"const-Map0 dz  med={med_dz_c:.3f}")
        axes[1].set_xlabel("dz [mm]")
        axes[2].hist(dy_c, bins=80, range=(-1, 1), histtype="step", color="C2")
        axes[2].set_title(f"const-Map0 dy  med={med_dy_c:.3f}")
        axes[2].set_xlabel("dy [mm]")
        fig.suptitle("Constant map (vx=+1): expect dx=+1, dz=dy=0")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 6))
        # subsample for scatter
        rng = np.random.default_rng(0)
        if len(y_y) > 20000:
            sel = rng.choice(len(y_y), 20000, replace=False)
        else:
            sel = np.arange(len(y_y))
        ax.scatter(y_y[sel], dx_y[sel], s=2, alpha=0.3, c="C0")
        ax.axvline(Y0, color="r", ls="--", label=f"y0={Y0}")
        ax.set_xlabel("hitpos_y Map_0 [mm]")
        ax.set_ylabel("dx = ydep - Map_0 [mm]")
        ax.set_title("y-dependent map: dx vs y (#22 clamp below y0)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(dx_y[below], bins=60, histtype="step", label=f"y<y0 n={np.sum(below)}")
        ax.hist(dx_y[mid], bins=60, histtype="step", label=f"mid n={np.sum(mid)}")
        ax.set_xlabel("dx [mm]")
        ax.set_title("ydep-Map0 dx by y region")
        ax.legend()
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    print(f"Wrote {pdf_path}")
    print(f"Wrote {SIDE_DIR / 'tpcpos_path_closure_summary.txt'}")
    print(f"const median dx={med_dx_c:.4f} dz={med_dz_c:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
