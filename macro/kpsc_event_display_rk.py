#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
K-p scattering event display with the beam trajectory (real data: RK beam VPs; G4: beam used for the
kinematics, the separate beam_prm particle, and truth directions).

Builds on check_tpc_helix_track_3d.py (same directory): helix parametrisation, theta range, TPC frame
and the display convention (matplotlib X, Y, Z) = (local z, local x, local y) are imported from there,
so the helix drawing is identical to that macro.

Inputs
  --helix : DstTPCHelixTracking(-Geant4) output (tree "tpc")
  --kpsc  : DstKpScattering(-Geant4) output (tree "kpsc") of the same events
  --g4    : G4 files (DstKpScatteringGeant4 + DstTPCHelixTrackingGeant4 with truth branches)

What is drawn (3D + top view z-x + side view z-y)
  - all TPC helices (grey), the selected proton track i_p (red) and kaon candidate i_kcand (blue),
    is_beam tracks dashed; hits of the selected tracks as points
  - reconstructed vertex (black star), LH2 cylinder (R = 40 mm, |y| <= 50 mm, centre z = -143 mm)
  - data: RK beam trajectory from DstKpScattering (xvp, yvp, zvp), orange
  - G4  : beam used for the kinematics (beam_used_*) through the vertex (orange),
          beam_prm_* straight line from its generation point (magenta dashed; separate particle),
          truth proton / kaon directions from the truth vertex (dotted)
  - text: p_meas, p_calc (elastic two-body from |p_beam| and the lab angle), dk = 1/p_meas - 1/p_calc

Usage
  # one event (data)
  python3 kpsc_event_display_rk.py --helix /group/.../run03774_TPCHelix.root \\
      --kpsc /group/.../root/run03774/run03774_KpScattering_lh2off.root --run 3774 --event 12345 --out ev.png
  # one event by kpsc entry (G4: event_number is not unique)
  python3 kpsc_event_display_rk.py --g4 --helix ... --kpsc ... --entry 25 --out ev_g4.png
  # N elastic-enriched LH2-core events into one PDF (data)
  python3 kpsc_event_display_rk.py --helix ... --kpsc ... --n 12 --out data_events.pdf
  # G4 (honest chain)
  python3 kpsc_event_display_rk.py --g4 --helix .../helix_B1T_truth_n100k.root \\
      --kpsc .../kpsc_B1T_hon_off_s0_n100k.root --n 12 --out g4_events.pdf

Selection for --n (same as the LH2 eloss studies):
  mmass>0, effective_ntTpc==2, close_dist<5, PID (|nsigma_k_sel|<2, !(|nsigma_pi_sel|<3 && nsigma_k_sel<-2),
  nsigma_p_sel>-1.5), LH2 core vertex (rho_xz<40, |y|<50 mm), diff_angle<0.1, ||dphi|-pi|<0.15.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from typing import Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_tpc_helix_track_3d as h3d  # noqa: E402  (helix geometry and TPC frame)

M_K = 0.493677
M_P = 0.938272
Z_TGT = h3d.TPC_Z_TARGET
LH2_R = 40.0
LH2_H = 50.0


# ---------------------------------------------------------------------------- kinematics
def p_two_body(pb: float, cos_t: float, m3: float, mmiss: float, pmeas: float) -> float:
    """|p_3| of the detected particle (mass m3) at cos(theta_lab) for K- p elastic (target at rest).
    4(S^2-c^2) x^2 - 4 A c x + 4 S^2 m3^2 - A^2 = 0; the physical root closest to pmeas."""
    S = math.sqrt(pb * pb + M_K * M_K) + M_P
    c = pb * cos_t
    A = S * S + m3 * m3 - pb * pb - mmiss * mmiss
    qa, qb, qc = 4.0 * (S * S - c * c), -4.0 * A * c, 4.0 * S * S * m3 * m3 - A * A
    disc = qb * qb - 4.0 * qa * qc
    if not (disc >= 0.0) or abs(qa) < 1e-18:
        return float("nan")
    best = float("nan")
    for x in ((-qb + math.sqrt(disc)) / (2 * qa), (-qb - math.sqrt(disc)) / (2 * qa)):
        if x > 1e-6 and A + 2.0 * c * x > 0.0:
            if not math.isfinite(best) or abs(x - pmeas) < abs(best - pmeas):
                best = x
    return best


# ---------------------------------------------------------------------------- I/O
def _open_tree(path: str, name: str):
    import uproot

    f = uproot.open(path)
    return f[name]


def _has(tree, b: str) -> bool:
    try:
        return b in tree.keys()
    except Exception:
        return False


KPSC_SCALARS = ["event_number", "mmass", "diff_angle", "close_dist", "effective_ntTpc", "nsigma_p_sel", "nsigma_k_sel",
                "nsigma_pi_sel", "vtx_x", "vtx_y", "vtx_z", "delta_phi_scat", "i_p", "i_kcand", "p_p", "p_k",
                "theta_p_lab", "theta_k_lab", "d5_momentum", "px_p", "py_p", "pz_p"]
KPSC_SCALARS_DATA = ["run_number"]
KPSC_SCALARS_G4 = ["beam_used_px", "beam_used_py", "beam_used_pz", "p_p_true", "p_k_true",
                   "px_p_true", "py_p_true", "pz_p_true", "vtx_x_true", "vtx_y_true", "vtx_z_true",
                   "beam_vtx_px", "beam_vtx_py", "beam_vtx_pz"]


def load_kpsc_scalars(path: str, g4: bool) -> Dict[str, np.ndarray]:
    t = _open_tree(path, "kpsc")
    names = KPSC_SCALARS + (KPSC_SCALARS_G4 if g4 else KPSC_SCALARS_DATA)
    names = [n for n in names if _has(t, n)]
    arr = t.arrays(names, library="np")
    if "run_number" not in arr:
        arr["run_number"] = np.zeros(len(arr["event_number"]), dtype=np.int64)
    return arr


def elastic_mask(a: Dict[str, np.ndarray]) -> np.ndarray:
    rho = np.hypot(a["vtx_x"], a["vtx_z"] - Z_TGT)
    m = (a["mmass"] > 0) & (a["effective_ntTpc"] == 2) & (a["close_dist"] < 5)
    m &= (a["nsigma_k_sel"] > -2) & (a["nsigma_k_sel"] < 2)
    m &= ~((np.abs(a["nsigma_pi_sel"]) < 3) & (a["nsigma_k_sel"] < -2)) & (a["nsigma_p_sel"] > -1.5)
    m &= (rho < LH2_R) & (np.abs(a["vtx_y"]) < LH2_H)
    m &= (a["diff_angle"] < 0.1) & (np.abs(np.abs(a["delta_phi_scat"]) - math.pi) < 0.15)
    return m


class HelixFile:
    """Random access to the helix tree.

    The kpsc tree is written entry by entry from the helix tree, so the same entry index is used when the
    event numbers agree at that entry (G4 event_number is not unique in the combined samples). Otherwise the
    (run, event) index is used, which requires unique event numbers (real data)."""

    BR = ["ntTpc", "nhtrack", "is_beam", "charge", "mom0", "pid", "helix_cx", "helix_cy", "helix_z0", "helix_r", "helix_dz",
          "helix_t", "hitpos_x", "hitpos_y", "hitpos_z"]
    BR_OPT = ["helix_theta_min", "helix_theta_max"]
    BR_G4 = ["beam_prm_px", "beam_prm_py", "beam_prm_pz", "beam_prm_x", "beam_prm_y", "beam_prm_z",
             "prm_k_px", "prm_k_py", "prm_k_pz"]

    def __init__(self, path: str, g4: bool):
        self.t = _open_tree(path, "tpc")
        self.g4 = g4
        ev = self.t["event_number"].array(library="np")
        run = self.t["run_number"].array(library="np") if (not g4 and _has(self.t, "run_number")) else np.zeros(len(ev), dtype=np.int64)
        self.ev = ev
        self.run = run
        self.index = {(int(r), int(e)): i for i, (r, e) in enumerate(zip(run, ev))}
        self.branches = [b for b in self.BR + self.BR_OPT + (self.BR_G4 if g4 else []) if _has(self.t, b)]

    def get(self, run: int, event: int, kpsc_entry: Optional[int] = None) -> Optional[dict]:
        i = None
        if kpsc_entry is not None and kpsc_entry < len(self.ev) and int(self.ev[kpsc_entry]) == int(event) \
                and (self.g4 or int(self.run[kpsc_entry]) == int(run)):
            i = kpsc_entry
        if i is None:
            if self.g4:
                return None  # event_number is not unique in G4; only entry-aligned lookup is safe
            i = self.index.get((int(run), int(event)))
        if i is None:
            return None
        a = self.t.arrays(self.branches, entry_start=i, entry_stop=i + 1, library="ak")
        import awkward as ak

        out = {}
        for b in self.branches:
            out[b] = ak.to_list(a[b][0])
        out["_entry"] = i
        return out


def kpsc_vp(path: str, entry: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    t = _open_tree(path, "kpsc")
    if not _has(t, "xvp"):
        return np.array([]), np.array([]), np.array([])
    a = t.arrays(["xvp", "yvp", "zvp"], entry_start=entry, entry_stop=entry + 1, library="np")
    x, y, z = np.asarray(a["xvp"][0], float), np.asarray(a["yvp"][0], float), np.asarray(a["zvp"][0], float)
    o = np.argsort(z)
    return x[o], y[o], z[o]


# ---------------------------------------------------------------------------- drawing helpers
def helix_points(hx: dict, i: int, n: int = 200) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    rng = None
    if "helix_theta_min" in hx and i < len(hx["helix_theta_min"]):
        a, b = float(hx["helix_theta_min"][i]), float(hx["helix_theta_max"][i])
        if math.isfinite(a) and math.isfinite(b) and abs(b - a) > 1e-12:
            lo, hi = min(a, b), max(a, b)
            m = (hi - lo) * 0.02 + 5e-4
            rng = (lo - m, hi + m)
    if rng is None and i < len(hx["helix_t"]):
        rng = h3d.helix_theta_draw_range(hx["helix_t"][i])
    r = float(hx["helix_r"][i])
    if rng is None or not math.isfinite(r):
        return None
    th = np.linspace(rng[0], rng[1], max(n, h3d.helix_polyline_sample_count(rng[0], rng[1])))
    return h3d.helix_xyz(float(hx["helix_cx"][i]), float(hx["helix_cy"][i]), float(hx["helix_z0"][i]), r, float(hx["helix_dz"][i]), th)


def line_through(p0: Tuple[float, float, float], d: Tuple[float, float, float], z0: float, z1: float):
    """Straight line through p0 with direction d, sampled from z0 to z1 (local coordinates)."""
    if abs(d[2]) < 1e-12:
        return None
    zs = np.linspace(z0, z1, 50)
    t = (zs - p0[2]) / d[2]
    return p0[0] + t * d[0], p0[1] + t * d[1], zs


def ray(p0, d, length: float = 250.0):
    n = math.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
    if not n > 0:
        return None
    s = np.linspace(0, length, 30)
    return p0[0] + s * d[0] / n, p0[1] + s * d[1] / n, p0[2] + s * d[2] / n


class Canvas:
    """3D + top view (z-x) + side view (z-y) + text panel."""

    def __init__(self, title: str):
        self.fig = plt.figure(figsize=(16, 10))
        self.ax3 = self.fig.add_subplot(2, 2, 1, projection="3d")
        self.axt = self.fig.add_subplot(2, 2, 2)
        self.axs = self.fig.add_subplot(2, 2, 4)
        self.axx = self.fig.add_subplot(2, 2, 3)
        self.axx.axis("off")
        self.fig.suptitle(title, fontsize=13)
        h3d._draw_tpc_frame(self.ax3)
        self.ax3.set_xlim(*h3d.DISPLAY_XLIM); self.ax3.set_ylim(*h3d.DISPLAY_YLIM); self.ax3.set_zlim(*h3d.DISPLAY_ZLIM)
        self.ax3.set_xlabel("z [mm]"); self.ax3.set_ylabel("x [mm]"); self.ax3.set_zlabel("y [mm]")
        self.ax3.view_init(elev=h3d.HELIX_VIEW_ELEV, azim=h3d.HELIX_VIEW_AZIM)
        for ax, lab in ((self.axt, "x [mm]"), (self.axs, "y [mm]")):
            ax.set_xlim(-320, 320); ax.set_ylim(-320, 320); ax.set_aspect("equal")
            ax.set_xlabel("z [mm]"); ax.set_ylabel(lab); ax.grid(alpha=0.3)
        self.axt.set_title("top view (z-x)", fontsize=10)
        self.axs.set_title("side view (z-y)", fontsize=10)
        th = np.linspace(0, 2 * math.pi, 100)
        self.axt.plot(Z_TGT + LH2_R * np.cos(th), LH2_R * np.sin(th), color="0.3", lw=1)
        self.axs.plot([Z_TGT - LH2_R, Z_TGT + LH2_R, Z_TGT + LH2_R, Z_TGT - LH2_R, Z_TGT - LH2_R],
                      [-LH2_H, -LH2_H, LH2_H, LH2_H, -LH2_H], color="0.3", lw=1)
        # TPC outer boundary (octagon inscribed radius ~ 250 mm around the target is drawn in 3D; 2D: circle)
        self.axt.plot(Z_TGT + 250 * np.cos(th), 250 * np.sin(th), color="0.8", lw=0.8, ls=":")

    def curve(self, x, y, z, **kw):
        X, Y, Z = h3d.tpc_local_to_display_vec(x, y, z)
        self.ax3.plot(X, Y, Z, **kw)
        kw2 = {k: v for k, v in kw.items() if k != "label"}
        self.axt.plot(z, x, **kw)
        self.axs.plot(z, y, **kw2)

    def points(self, x, y, z, **kw):
        X, Y, Z = h3d.tpc_local_to_display_vec(x, y, z)
        self.ax3.scatter(X, Y, Z, **kw)
        self.axt.scatter(z, x, **kw)
        self.axs.scatter(z, y, **kw)

    def text(self, lines: List[str]):
        self.axx.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", family="monospace", fontsize=10.5, transform=self.axx.transAxes)

    def finish(self):
        self.axt.legend(loc="upper left", fontsize=8, framealpha=0.8)
        self.fig.tight_layout(rect=(0, 0, 1, 0.96))
        return self.fig


# ---------------------------------------------------------------------------- one event
def draw_event(k: Dict[str, np.ndarray], ik: int, hx: dict, g4: bool, kpsc_path: str, label: str):
    run, ev = int(k["run_number"][ik]), int(k["event_number"][ik])
    ip, ikc = int(k["i_p"][ik]), int(k["i_kcand"][ik])
    vtx = (float(k["vtx_x"][ik]), float(k["vtx_y"][ik]), float(k["vtx_z"][ik]))
    c = Canvas(f"{label}   run {run}  event {ev}   (kpsc entry {ik}, helix entry {hx['_entry']})")
    if hx["_entry"] != ik and len(hx["helix_r"]) <= max(int(k["i_p"][ik]), int(k["i_kcand"][ik])):
        print(f"warning: helix entry {hx['_entry']} has fewer tracks than the kpsc track indices")
    # helices
    nt = len(hx["helix_r"])
    for i in range(nt):
        pts = helix_points(hx, i)
        if pts is None:
            continue
        if i == ip:
            kw = dict(color="tab:red", lw=2.4, label=f"proton track (i_p={ip})")
        elif i == ikc:
            kw = dict(color="tab:blue", lw=2.4, label=f"kaon candidate (i_kcand={ikc})")
        elif int(hx["is_beam"][i]) == 1:
            kw = dict(color="0.35", lw=1.2, ls="--", label="is_beam track")
        else:
            kw = dict(color="0.6", lw=1.0)
        c.curve(*pts, **kw)
        if i in (ip, ikc):
            c.points(np.asarray(hx["hitpos_x"][i]), np.asarray(hx["hitpos_y"][i]), np.asarray(hx["hitpos_z"][i]),
                     s=10, color="tab:red" if i == ip else "tab:blue", alpha=0.8)
    # vertex
    c.points(np.array([vtx[0]]), np.array([vtx[1]]), np.array([vtx[2]]), s=90, marker="*", color="k", label="vertex (reco)")
    # beam
    if not g4:
        x, y, z = kpsc_vp(kpsc_path, ik)
        if len(z):
            c.curve(x, y, z, color="tab:orange", lw=2.0, label="RK beam (xvp, yvp, zvp)")
    else:
        bu = (float(k["beam_used_px"][ik]), float(k["beam_used_py"][ik]), float(k["beam_used_pz"][ik]))
        ln = line_through(vtx, bu, -320.0, vtx[2])
        if ln is not None:
            c.curve(*ln, color="tab:orange", lw=2.0, label="beam used (through vtx)")
        if "beam_prm_px" in hx and math.isfinite(float(hx["beam_prm_px"])):
            p0 = (float(hx["beam_prm_x"]), float(hx["beam_prm_y"]), float(hx["beam_prm_z"]))
            d = (float(hx["beam_prm_px"]), float(hx["beam_prm_py"]), float(hx["beam_prm_pz"]))
            ln2 = line_through(p0, d, -320.0, 320.0)
            if ln2 is not None:
                c.curve(*ln2, color="tab:purple", lw=1.4, ls="--", label="beam_prm (straight, separate particle)")
        vt = (float(k["vtx_x_true"][ik]), float(k["vtx_y_true"][ik]), float(k["vtx_z_true"][ik]))
        rp = ray(vt, (float(k["px_p_true"][ik]), float(k["py_p_true"][ik]), float(k["pz_p_true"][ik])))
        if rp is not None:
            c.curve(*rp, color="tab:red", lw=1.2, ls=":", label="truth proton direction")
        if "prm_k_px" in hx and math.isfinite(float(hx["prm_k_px"])):
            rk = ray(vt, (float(hx["prm_k_px"]), float(hx["prm_k_py"]), float(hx["prm_k_pz"])))
            if rk is not None:
                c.curve(*rk, color="tab:blue", lw=1.2, ls=":", label="truth kaon direction")
    # numbers
    pb = float(k["d5_momentum"][ik])
    pp, pk = float(k["p_p"][ik]), float(k["p_k"][ik])
    pcp = p_two_body(pb, math.cos(float(k["theta_p_lab"][ik])), M_P, M_K, pp)
    pck = p_two_body(pb, math.cos(float(k["theta_k_lab"][ik])), M_K, M_P, pk)
    def dk(pm, pc):
        return 1.0 / pm - 1.0 / pc if (pm > 0 and pc > 0) else float("nan")
    lines = [f"{label}", "",
             f"|p_beam| (d5_momentum) = {pb:.4f} GeV/c",
             f"mmass - m_K = {float(k['mmass'][ik]) - M_K:+.4f} GeV",
             f"diff_angle  = {float(k['diff_angle'][ik]):.4f} rad",
             f"close_dist  = {float(k['close_dist'][ik]):.2f} mm",
             f"vertex (x,y,z) = ({vtx[0]:.1f}, {vtx[1]:.1f}, {vtx[2]:.1f}) mm", "",
             f"proton: theta_lab = {float(k['theta_p_lab'][ik]):.3f} rad",
             f"   p_meas = {pp:.4f}  p_calc = {pcp:.4f} GeV/c   dk = {dk(pp, pcp):+.3f} (GeV/c)^-1",
             f"kaon  : theta_lab = {float(k['theta_k_lab'][ik]):.3f} rad",
             f"   p_meas = {pk:.4f}  p_calc = {pck:.4f} GeV/c   dk = {dk(pk, pck):+.3f} (GeV/c)^-1"]
    if ip < len(hx["nhtrack"]) and ikc < len(hx["nhtrack"]):
        lines += ["", f"nhits: proton {hx['nhtrack'][ip]}, kaon {hx['nhtrack'][ikc]};  helix_r: proton {hx['helix_r'][ip]:.0f}, kaon {hx['helix_r'][ikc]:.0f} mm"]
    if g4:
        lines += ["", f"truth: p_p = {float(k['p_p_true'][ik]):.4f}, p_K = {float(k['p_k_true'][ik]):.4f} GeV/c",
                  "orange: beam used for the kinematics (incoming K- at the vertex, truth),",
                  "purple dashed: beam_prm (beam entry of the combined sample, z=-1150 mm)",
                  "dotted: truth proton / kaon directions from the truth vertex"]
    else:
        lines += ["", "orange: RK beam trajectory (HSTrack VPs, from D5 at z=-1300.9 mm)"]
    c.text(lines)
    return c.finish()


# ---------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--helix", required=True)
    ap.add_argument("--kpsc", required=True)
    ap.add_argument("--g4", action="store_true")
    ap.add_argument("--run", type=int, default=None)
    ap.add_argument("--event", type=int, default=None)
    ap.add_argument("--entry", type=int, default=None, help="kpsc tree entry (use for G4, where event_number is not unique)")
    ap.add_argument("--n", type=int, default=0, help="number of elastic-enriched events to draw into one PDF")
    ap.add_argument("--skip", type=int, default=0, help="skip the first N selected events")
    ap.add_argument("--label", default=None)
    ap.add_argument("--out", required=True, help=".png (one event) or .pdf")
    a = ap.parse_args(argv)
    label = a.label or ("G4 (honest chain)" if a.g4 else "real data")
    k = load_kpsc_scalars(a.kpsc, a.g4)
    hx = HelixFile(a.helix, a.g4)
    if a.n > 0:
        idx = np.nonzero(elastic_mask(k))[0][a.skip:]
        print(f"selected {len(idx) + a.skip} elastic-enriched events; drawing {min(a.n, len(idx))}")
        with PdfPages(a.out) as pdf:
            nd = 0
            for ik in idx:
                if nd >= a.n:
                    break
                h = hx.get(int(k["run_number"][ik]), int(k["event_number"][ik]), int(ik))
                if h is None:
                    continue
                fig = draw_event(k, int(ik), h, a.g4, a.kpsc, label)
                pdf.savefig(fig); plt.close(fig); nd += 1
        print(f"wrote {a.out} ({nd} events)")
        return
    if a.entry is not None:
        if not (0 <= a.entry < len(k["event_number"])):
            sys.exit("--entry out of range")
        ik = int(a.entry)
    else:
        if a.event is None:
            ap.error("give --event (and --run for data), --entry, or --n")
        sel = np.nonzero((k["event_number"] == a.event) & ((k["run_number"] == a.run) if (a.run is not None and not a.g4) else True))[0]
        if len(sel) == 0:
            sys.exit("event not found in kpsc")
        if len(sel) > 1:
            print(f"note: event_number {a.event} appears {len(sel)} times; drawing the first (kpsc entry {int(sel[0])}); use --entry to choose")
        ik = int(sel[0])
    h = hx.get(int(k["run_number"][ik]), int(k["event_number"][ik]), ik)
    if h is None:
        sys.exit("event not found in helix file")
    fig = draw_event(k, ik, h, a.g4, a.kpsc, label)
    fig.savefig(a.out, dpi=150)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
