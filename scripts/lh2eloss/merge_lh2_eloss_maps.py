#!/usr/bin/env python3
"""Combine single-particle LH2 eloss maps into one Lh2ElossMan map file.

Inputs are 5-column maps written by kpsc_g4_lh2_eloss_map
(particle_id exit_id p L Loss, plus Target* / VesselRadius header keys).
The header keys must agree between all inputs; each particle_id may come
from one input only. The output keeps every input's comment block, so the
build provenance of each species stays in the file.

Usage:
  merge_lh2_eloss_maps.py -o Lh2ElossMap_e72 Lh2ElossMap_proton Lh2ElossMap_pion Lh2ElossMap_kaon
"""
import argparse
import sys

KEYS = ("TargetRadius", "TargetHalfHeight", "TargetCenterX", "TargetCenterY",
        "TargetCenterZ", "VesselRadius")
NAMES = {0: "proton", 1: "pion", 2: "kaon"}


def read_map(path):
    header, comments, rows = {}, [], []
    for n, line in enumerate(open(path), 1):
        body = line.split("#")[0].split()
        if not body:
            if line.lstrip().startswith("#"):
                comments.append(line.rstrip("\n"))
            continue
        if body[0] in KEYS:
            header[body[0]] = float(body[1])
            continue
        if len(body) != 5 or not body[0].isdigit():
            sys.exit(f"{path}:{n}: not a 5-column grid line: {line.strip()}")
        pid, exit_id = int(body[0]), int(body[1])
        if pid not in NAMES or exit_id not in (0, 1):
            sys.exit(f"{path}:{n}: unknown particle_id/exit_id: {line.strip()}")
        rows.append((pid, exit_id, body[2], body[3], body[4]))
    missing = [k for k in KEYS if k not in header]
    if missing:
        sys.exit(f"{path}: missing header keys {missing}")
    if not rows:
        sys.exit(f"{path}: no grid lines")
    return header, comments, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    ref, owner, blocks = None, {}, []
    for path in args.inputs:
        header, comments, rows = read_map(path)
        if ref is None:
            ref = header
        elif any(abs(header[k] - ref[k]) > 1e-9 for k in KEYS):
            sys.exit(f"{path}: header {header} differs from {args.inputs[0]} {ref}")
        for pid in sorted({r[0] for r in rows}):
            if pid in owner:
                sys.exit(f"particle {NAMES[pid]} in both {owner[pid]} and {path}")
            owner[pid] = path
        blocks.append((path, comments, rows))
    if 0 not in owner:
        sys.exit("no proton grid (particle_id 0); Lh2ElossMan requires it")

    with open(args.output, "w") as f:
        f.write("#\n#  LH2 eloss map, merged by merge_lh2_eloss_maps.py from:\n")
        for path, _, rows in blocks:
            names = ",".join(NAMES[p] for p in sorted({r[0] for r in rows}))
            f.write(f"#    {names}: {path}\n")
        f.write("#  columns: particle_id exit_id p[MeV/c] L[mm] Loss[MeV/c]\n")
        f.write("#  particle_id: 0=proton 1=pion 2=kaon   exit_id: 0=side 1=cap\n#\n")
        for k in KEYS:
            f.write(f"{k:<18} {ref[k]}    # mm\n")
        for path, comments, rows in blocks:
            f.write("\n#______________________________________________________________________________\n")
            f.write(f"#  from {path}\n")
            for c in comments:
                if c.strip("#_ ").strip():
                    f.write(c + "\n")
            for r in rows:
                f.write("  ".join(str(x) for x in r) + "\n")
    print(f"wrote {args.output}: " + ", ".join(f"{NAMES[p]} <- {owner[p]}" for p in sorted(owner)))


if __name__ == "__main__":
    main()
