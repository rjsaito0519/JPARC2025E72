#!/usr/bin/env python3
"""Approximate chamber (dx, dy) from DCGEO Ofs [mm].

Per detector (BLC1a/b, BLC2a/b):
  mu_U = mean Ofs of 4 U/UP planes
  mu_V = mean Ofs of 4 V/VP planes
  TA_U = mean TA of those U planes  (from DCGEO)
  TA_V = mean TA of those V planes  (from DCGEO)

Solve (RA=0):
  mu_U = cos(TA_U)*dx + sin(TA_U)*dy
  mu_V = cos(TA_V)*dx + sin(TA_V)*dy

Reference for Delta Ofs:
  zero   (default) — use Ofs as-is (vs Ofs=0)
  survey — Ofs - (cos(TA)*X + sin(TA)*Y)  [layer TA from DCGEO]
  FILE   — Ofs - Ofs_ref from another DCGEO

Local plane-by-plane scatter is ignored (U/V averages only).

Example:
  python3 invert_dcgeo_ofs_xy.py param/DCGEO/e72/DCGeomParam_run03773_K
  python3 invert_dcgeo_ofs_xy.py DCGeomParam_run03773_K --ref survey
  python3 invert_dcgeo_ofs_xy.py DCGeomParam_run03773_K --ref param/DCGEO/DCGeomParam_e72_example
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path

DETECTORS = ("BLC1a", "BLC1b", "BLC2a", "BLC2b")


def parse_blc_rows(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open() as src:
        for line in src:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 13:
                continue
            name = parts[1]
            if not name.startswith("BLC") or "-" not in name:
                continue
            try:
                det, plane = name.split("-", 1)
                rows.append(
                    {
                        "name": name,
                        "det": det,
                        "plane": plane,
                        "X": float(parts[2]),
                        "Y": float(parts[3]),
                        "TA": float(parts[5]),
                        "Ofs": float(parts[12]),
                    }
                )
            except ValueError:
                continue
    return rows


def is_u(plane: str) -> bool:
    return plane.startswith("U")


def is_v(plane: str) -> bool:
    return plane.startswith("V")


def delta_ofs(row: dict, ref_ofs: dict[str, float] | None, ref_mode: str) -> float:
    if ref_mode == "zero":
        return row["Ofs"]
    if ref_mode == "survey":
        ta = math.radians(row["TA"])
        return row["Ofs"] - (math.cos(ta) * row["X"] + math.sin(ta) * row["Y"])
    assert ref_ofs is not None
    if row["name"] not in ref_ofs:
        raise KeyError(f"ref missing layer {row['name']}")
    return row["Ofs"] - ref_ofs[row["name"]]


def approx_dx_dy(mu_u: float, mu_v: float, ta_u_deg: float, ta_v_deg: float):
    """Solve mu = cos(TA)*dx + sin(TA)*dy for U and V means."""
    cu = math.cos(math.radians(ta_u_deg))
    su = math.sin(math.radians(ta_u_deg))
    cv = math.cos(math.radians(ta_v_deg))
    sv = math.sin(math.radians(ta_v_deg))
    det = cu * sv - su * cv
    if abs(det) < 1e-12:
        raise ValueError(f"singular TA_U={ta_u_deg}, TA_V={ta_v_deg}")
    dx = (mu_u * sv - su * mu_v) / det
    dy = (cu * mu_v - mu_u * cv) / det
    return dx, dy


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Approximate BLC (dx,dy) from averaged U/V Ofs (TA from DCGEO)."
    )
    parser.add_argument("dcgeo", type=Path, help="DCGEO parameter file")
    parser.add_argument(
        "--ref",
        default="zero",
        help="zero | survey | path to reference DCGEO (default: zero)",
    )
    args = parser.parse_args()

    if not args.dcgeo.is_file():
        print(f"[Error] not found: {args.dcgeo}", file=sys.stderr)
        return 1

    ref_mode = args.ref
    ref_ofs: dict[str, float] | None = None
    if ref_mode not in ("zero", "survey"):
        ref_path = Path(ref_mode)
        if not ref_path.is_file():
            print(f"[Error] ref DCGEO not found: {ref_path}", file=sys.stderr)
            return 1
        ref_ofs = {r["name"]: r["Ofs"] for r in parse_blc_rows(ref_path)}
        ref_mode = "file"

    rows = parse_blc_rows(args.dcgeo)
    by_det: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_det[r["det"]].append(r)

    print(f"DCGEO: {args.dcgeo}")
    print(f"ref:   {args.ref}")
    print("approx: mu_U/V = cos(TA)*dx + sin(TA)*dy  (TA_U/V = mean from DCGEO)")
    print(
        f"{'det':6} {'TA_U':>7} {'TA_V':>7} {'mu_U':>9} {'mu_V':>9} "
        f"{'dx':>9} {'dy':>9}  nU nV"
    )
    print("-" * 72)

    for det in DETECTORS:
        planes = by_det.get(det, [])
        if not planes:
            print(f"{det:6}  (no layers)")
            continue
        u_rows = [r for r in planes if is_u(r["plane"])]
        v_rows = [r for r in planes if is_v(r["plane"])]
        try:
            us = [delta_ofs(r, ref_ofs, ref_mode) for r in u_rows]
            vs = [delta_ofs(r, ref_ofs, ref_mode) for r in v_rows]
        except KeyError as exc:
            print(f"[Error] {exc}", file=sys.stderr)
            return 1
        if len(us) != 4 or len(vs) != 4:
            print(
                f"[Warn] {det}: expected 4 U and 4 V, got nU={len(us)} nV={len(vs)}",
                file=sys.stderr,
            )
        if not us or not vs:
            print(f"{det:6}  (missing U or V)")
            continue
        mu_u = sum(us) / len(us)
        mu_v = sum(vs) / len(vs)
        ta_u = sum(r["TA"] for r in u_rows) / len(u_rows)
        ta_v = sum(r["TA"] for r in v_rows) / len(v_rows)
        try:
            dx, dy = approx_dx_dy(mu_u, mu_v, ta_u, ta_v)
        except ValueError as exc:
            print(f"[Error] {det}: {exc}", file=sys.stderr)
            return 1
        print(
            f"{det:6} {ta_u:7.2f} {ta_v:7.2f} {mu_u:9.4f} {mu_v:9.4f} "
            f"{dx:9.4f} {dy:9.4f}  {len(us):2d} {len(vs):2d}"
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
