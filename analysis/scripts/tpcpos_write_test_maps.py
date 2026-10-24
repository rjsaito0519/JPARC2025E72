#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Write TPCPOS test maps for Phase 3 path closure (tmp only; not param/).

Reads Map_0 lattice and writes:
  - const_vx1: (vx,vz)=(+1,0) everywhere
  - ydep_vx: vx = iy (mm) at lattice index iy=0..n_y-1, vz=0
"""

from __future__ import annotations

import argparse
from pathlib import Path

MAP0 = Path("/home/had/sryuta/analyzer/JPARC2025E72/param/TPCPOS/TPCPositionCorrectionMap_0")
OUT_DIR = Path("/group/had/sks/Users/sryuta/tmp/debug/tpcpos/maps")


def rewrite(path_in: Path, path_out: Path, mode: str) -> None:
    lines = path_in.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise RuntimeError(f"empty map: {path_in}")
    header = lines[0].split()
    n_y = int(header[0])
    y0 = float(header[1])
    dy = float(header[2])
    out = [lines[0]]
    n_body = 0
    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        layer, row = int(parts[0]), int(parts[1])
        y = float(parts[2])
        iy = int(round((y - y0) / dy))
        if mode == "const_vx1":
            vx, vz = 1.0, 0.0
        elif mode == "ydep_vx":
            vx, vz = float(iy), 0.0
        else:
            raise ValueError(mode)
        out.append(f"{layer}   {row}   {y:g}   {vx:g}   {vz:g}")
        n_body += 1
    path_out.write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"Wrote {path_out}  n_y={n_y} body_lines={n_body}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--map0", type=Path, default=MAP0)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rewrite(
        args.map0,
        args.out_dir / "TPCPositionCorrectionMap_test_const_vx1",
        "const_vx1",
    )
    rewrite(
        args.map0,
        args.out_dir / "TPCPositionCorrectionMap_test_ydep_vx",
        "ydep_vx",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
