#!/usr/bin/env python3
"""Project BLC chamber survey displacements onto DCGEO Ofs values.

Survey input is a rigid translation (dx, dy [, dz]) of the chamber center in
the FF / experiment coordinate system [mm].  For each BLC layer the script
uses the existing TA, RA1, RA2 to build the local-s axis (kS2s in
DCGeomRecord::CalcVectors) and adds

    delta_ofs = scale * (dsdx*dx + dsdy*dy + dsdz*dz)

to the current Ofs column only.  X, Y, Z, L, TA, and RA are never modified.

Recommended workflow
--------------------
1. apply_dcgeo_offset.py  -- first-pass survey correction (this script)
2. tracking / D5 validation
3. update_param.py residual  -- per-layer fine tuning (Ofs column, start_col=12)

Notes
-----
* Prefer per-detector offsets via --offset DET:DX,DY (dz optional, default 0).
* TA is read per layer from DCGEO (nominal U/V tilt may be ~±45 deg but need not
  be exact, e.g. 44 deg after survey); projection uses each layer's actual value.
* RA enters only through the s-axis projection; it is not written back to DCGEO.
* With RA1=RA2=0 (typical BLC), dsdz=0 so dz does not change Ofs; adjust Z via
  other geometry if needed.
* Try scale=+1 first; if residuals worsen, retry with --scale -1.

See also: docs/dcgeo_tracking_geometry.md (Ofs / TA / RA in tracking).
"""

import argparse
import math
import os
import sys
from dataclasses import dataclass

BLC_DETECTORS = ("BLC1a", "BLC1b", "BLC2a", "BLC2b")

_DSDZ_EPS = 1e-12


@dataclass
class DCGeoRecord:
    line_index: int
    fields: list[str]
    raw_line: str
    separator: str
    line_ending: str

    @property
    def name(self):
        return self.fields[1]

    @property
    def detector(self):
        return self.name.split("-", 1)[0]

    @property
    def ta(self):
        return float(self.fields[5])

    @property
    def ra1(self):
        return float(self.fields[6])

    @property
    def ra2(self):
        return float(self.fields[7])

    @property
    def ofs(self):
        return float(self.fields[12])


@dataclass
class OffsetResult:
    record: DCGeoRecord
    dx: float
    dy: float
    dz: float
    dsdx: float
    dsdy: float
    dsdz: float
    delta_ofs: float

    @property
    def new_ofs(self):
        return self.record.ofs + self.delta_ofs


def is_blc_layer(name):
    if "-" not in name:
        return False
    detector = name.split("-", 1)[0]
    return detector in BLC_DETECTORS


def split_dcgeo_line(line):
    body = line.rstrip("\r\n")
    line_ending = line[len(body):]
    if not line_ending:
        line_ending = "\n"

    if "\t" in body:
        fields = body.split("\t")
        separator = "\t"
    else:
        fields = body.split()
        separator = " "
    return fields, separator, line_ending


def calc_s_axis(ta_deg, ra1_deg, ra2_deg):
    """Return the local-s axis in global coordinates.

    This is the kS2s definition used by DCGeomRecord::CalcVectors().
    """

    ta = math.radians(ta_deg)
    ra1 = math.radians(ra1_deg)
    ra2 = math.radians(ra2_deg)

    ct0 = math.cos(ta)
    st0 = math.sin(ta)
    ct1 = math.cos(ra1)
    st1 = math.sin(ra1)
    ct2 = math.cos(ra2)
    st2 = math.sin(ra2)

    dsdx = ct0 * ct2 + st0 * st1 * st2
    dsdy = st0 * ct1
    dsdz = -ct0 * st2 + st0 * st1 * ct2
    return dsdx, dsdy, dsdz


def parse_dcgeo(path):
    with open(path, "r") as src:
        lines = src.readlines()

    records = []
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        fields, separator, line_ending = split_dcgeo_line(line)
        if len(fields) < 13:
            continue
        if not is_blc_layer(fields[1]):
            continue
        try:
            float(fields[5])
            float(fields[6])
            float(fields[7])
            float(fields[12])
        except ValueError:
            continue
        records.append(DCGeoRecord(i, fields, line, separator, line_ending))

    return lines, records


def parse_detector_offset(text):
    try:
        detector, values = text.split(":", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "offset must be DET:DX,DY[,DZ], e.g. BLC1a:0.5,-0.3"
        ) from exc

    detector = detector.strip()
    if detector not in BLC_DETECTORS:
        raise argparse.ArgumentTypeError(
            f"detector must be one of {', '.join(BLC_DETECTORS)}"
        )

    parts = [p.strip() for p in values.split(",")]
    if len(parts) not in (2, 3):
        raise argparse.ArgumentTypeError("offset values must be DX,DY or DX,DY,DZ")
    try:
        dx = float(parts[0])
        dy = float(parts[1])
        dz = float(parts[2]) if len(parts) == 3 else 0.0
    except ValueError as exc:
        raise argparse.ArgumentTypeError("offset values must be numbers") from exc
    return detector, (dx, dy, dz)


def select_detectors(args):
    selected = set(args.detector or [])
    if args.bcin:
        selected.update(("BLC1a", "BLC1b"))
    if args.bcout:
        selected.update(("BLC2a", "BLC2b"))
    for detector, _ in args.offset:
        selected.add(detector)
    if not selected:
        selected.update(BLC_DETECTORS)
    return tuple(det for det in BLC_DETECTORS if det in selected)


def build_offsets(args):
    offsets = {det: (args.dx, args.dy, args.dz) for det in BLC_DETECTORS}
    for detector, offset in args.offset:
        offsets[detector] = offset
    return offsets


def compute_results(records, detectors, offsets, scale):
    selected = set(detectors)
    results = []
    for record in records:
        if record.detector not in selected:
            continue
        dx, dy, dz = offsets[record.detector]
        dsdx, dsdy, dsdz = calc_s_axis(record.ta, record.ra1, record.ra2)
        delta_ofs = scale * (dsdx * dx + dsdy * dy + dsdz * dz)
        results.append(OffsetResult(record, dx, dy, dz, dsdx, dsdy, dsdz, delta_ofs))
    return results


def format_float(value):
    return f"{value:.8f}".rstrip("0").rstrip(".")


def format_ofs(value, reference_str=None):
    if reference_str and "." in reference_str:
        ref_decimals = len(reference_str.split(".", 1)[1])
        decimals = max(ref_decimals, 12)
        formatted = f"{value:.{decimals}f}"
    else:
        formatted = f"{value:.12g}"
    if "." in formatted:
        formatted = formatted.rstrip("0").rstrip(".")
    return formatted


def rebuild_line(record, new_ofs_str):
    fields = record.fields[:]
    fields[12] = new_ofs_str
    return record.separator.join(fields) + record.line_ending


def print_results(results):
    if not results:
        return

    header = (
        "Name", "TA", "RA1", "RA2", "dx", "dy", "dz",
        "dsdx", "dsdy", "dsdz", "old_ofs", "delta", "new_ofs",
    )
    rows = []
    for res in results:
        rec = res.record
        rows.append(
            (
                rec.name,
                format_float(rec.ta),
                format_float(rec.ra1),
                format_float(rec.ra2),
                format_float(res.dx),
                format_float(res.dy),
                format_float(res.dz),
                format_float(res.dsdx),
                format_float(res.dsdy),
                format_float(res.dsdz),
                format_float(rec.ofs),
                format_float(res.delta_ofs),
                format_float(res.new_ofs),
            )
        )

    widths = [max(len(str(row[i])) for row in (header, *rows)) for i in range(len(header))]
    fmt = "  ".join(f"{{:<{width}}}" for width in widths)
    print(fmt.format(*header))
    print(fmt.format(*("-" * width for width in widths)))
    for row in rows:
        print(fmt.format(*row))


def print_warnings(results, offsets, detectors, write_path=None):
    if not results:
        print("No BLC records matched the selection.", file=sys.stderr)
        return False

    max_dz_input = max(abs(offsets[det][2]) for det in detectors)
    if max_dz_input > _DSDZ_EPS:
        max_dsdz = max(abs(res.dsdz) for res in results)
        if max_dsdz < _DSDZ_EPS:
            print(
                "Warning: dz is non-zero but |dsdz| is ~0 for all selected layers; "
                "with RA1=RA2=0, dz does not affect Ofs.",
                file=sys.stderr,
            )

    if write_path:
        by_det = {}
        for res in results:
            by_det.setdefault(res.record.detector, 0.0)
            by_det[res.record.detector] += res.delta_ofs
        parts = [f"{det}: sum(delta_ofs)={by_det[det]:.6f}" for det in sorted(by_det)]
        print(
            f"Writing {len(results)} layer(s) to {write_path}: " + ", ".join(parts),
            file=sys.stderr,
        )

    return True


def write_dcgeo(input_lines, output_path, results):
    by_line = {res.record.line_index: res for res in results}
    output_lines = list(input_lines)
    for line_index, res in by_line.items():
        rec = res.record
        new_ofs_str = format_ofs(res.new_ofs, rec.fields[12])
        output_lines[line_index] = rebuild_line(rec, new_ofs_str)

    with open(output_path, "w") as dst:
        dst.writelines(output_lines)


def _assert_close(actual, expected, tol=1e-12):
    if abs(actual - expected) >= tol:
        raise AssertionError(f"expected {expected}, got {actual}")


def _expected_s_axis_ra0(ta_deg):
    """Reference s-axis for RA1=RA2=0 (any TA, not only ±45 deg)."""
    ta = math.radians(ta_deg)
    return math.cos(ta), math.sin(ta), 0.0


def _run_self_test():
    # calc_s_axis uses cos/sin only; no special-case constants for nominal tilts.
    for ta_deg in (-45, 135, 45, 44):
        exp_dsdx, exp_dsdy, exp_dsdz = _expected_s_axis_ra0(ta_deg)
        dsdx, dsdy, dsdz = calc_s_axis(ta_deg, 0, 0)
        _assert_close(dsdx, exp_dsdx)
        _assert_close(dsdy, exp_dsdy)
        _assert_close(dsdz, exp_dsdz)

    dsdx, dsdy, _ = calc_s_axis(44, 0, 0)
    _assert_close(dsdx * 1.0 + dsdy * 0.0, math.cos(math.radians(44)))

    print("self-test passed")


def make_arg_parser():
    epilog = """\
examples:
  %(prog)s param/DCGEO/e72/DCGeomParam_run03772_Pi \\
      --offset BLC1a:0.5,-0.3 --bcin
  %(prog)s param/DCGEO/e72/DCGeomParam_run03772_Pi \\
      --offset BLC1a:0.5,-0.3 --bcin --write DCGeomParam_updated
  %(prog)s --self-test

After survey correction, fine-tune per-layer Ofs with:
  update_param.py <run> residual <suffix>
"""
    parser = argparse.ArgumentParser(
        description=(
            "Project BLC chamber-center survey offsets (dx, dy) onto DCGEO Ofs. "
            "Only the Ofs column is updated; TA/RA are used for projection only."
        ),
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "dcgeo",
        nargs="?",
        help="input DCGeomParam file (not required with --self-test)",
    )
    parser.add_argument(
        "--dx",
        type=float,
        default=0.0,
        help="default global X offset [mm] for all BLC detectors (prefer --offset)",
    )
    parser.add_argument(
        "--dy",
        type=float,
        default=0.0,
        help="default global Y offset [mm] for all BLC detectors (prefer --offset)",
    )
    parser.add_argument(
        "--dz",
        type=float,
        default=0.0,
        help="default global Z offset [mm]; ineffective when RA1=RA2=0 (prefer --offset)",
    )
    parser.add_argument(
        "--offset",
        action="append",
        type=parse_detector_offset,
        default=[],
        metavar="DET:DX,DY[,DZ]",
        help="per-detector survey offset at chamber center, e.g. BLC1a:0.5,-0.3",
    )
    parser.add_argument(
        "--detector",
        action="append",
        choices=BLC_DETECTORS,
        help="limit target detector; can be given multiple times",
    )
    parser.add_argument("--bcin", action="store_true", help="target BLC1a and BLC1b")
    parser.add_argument("--bcout", action="store_true", help="target BLC2a and BLC2b")
    parser.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help="multiply the projected offset before adding to Ofs; use -1 to invert sign",
    )
    parser.add_argument("--write", metavar="PATH", help="write updated DCGeomParam to PATH")
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="run internal consistency checks for calc_s_axis and exit",
    )
    return parser


def main():
    parser = make_arg_parser()
    args = parser.parse_args()

    if args.self_test:
        _run_self_test()
        return

    if not args.dcgeo:
        parser.error("dcgeo is required unless --self-test is given")

    if not os.path.exists(args.dcgeo):
        parser.error(f"input file does not exist: {args.dcgeo}")

    lines, records = parse_dcgeo(args.dcgeo)
    detectors = select_detectors(args)
    offsets = build_offsets(args)
    results = compute_results(records, detectors, offsets, args.scale)

    print_results(results)
    ok = print_warnings(results, offsets, detectors, write_path=args.write)
    if not ok:
        sys.exit(1)

    if args.write:
        write_dcgeo(lines, args.write, results)
        print(f"\nWrote updated DCGEO: {args.write}")


if __name__ == "__main__":
    main()
