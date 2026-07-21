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
1. apply_dcgeo_offset.py DCGEO --write   -- survey: X,Y columns -> Ofs (overwrites file)
2. tracking / D5 validation
3. apply_dcgeo_offset.py DCGEO --skip-xy --bcin --resiy 0.2 --write
   -- extra correction from CLI only (ignores X,Y columns)
4. update_param.py residual  -- per-layer fine tuning (Ofs column, start_col=12)

Notes
-----
* By default each layer's X,Y columns are read as survey displacement [mm] from
  nominal (global FF coordinates) and projected onto Ofs via TA/RA.
* Use --skip-xy for follow-up runs: only --resix/--resiy/--residual-offset apply.
* Re-running without --skip-xy on the same file double-counts into Ofs; keep a .bk.
* --offset DET:DX,DY overrides file X,Y for that detector/layer when not --skip-xy.
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
BCIN_DETECTORS = ("BLC1a", "BLC1b")
BCOUT_DETECTORS = ("BLC2a", "BLC2b")

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

    @property
    def survey_x(self):
        return float(self.fields[2])

    @property
    def survey_y(self):
        return float(self.fields[3])


@dataclass
class OffsetResult:
    record: DCGeoRecord
    dx_survey: float
    dy_survey: float
    dz_survey: float
    dx_residual: float
    dy_residual: float
    dz_residual: float
    dsdx: float
    dsdy: float
    dsdz: float
    delta_ofs_survey: float
    delta_ofs_residual: float
    delta_ofs: float

    @property
    def dx(self):
        return self.dx_survey + self.dx_residual

    @property
    def dy(self):
        return self.dy_survey + self.dy_residual

    @property
    def dz(self):
        return self.dz_survey + self.dz_residual

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
            float(fields[2])
            float(fields[3])
            float(fields[5])
            float(fields[6])
            float(fields[7])
            float(fields[12])
        except ValueError:
            continue
        records.append(DCGeoRecord(i, fields, line, separator, line_ending))

    return lines, records


def _parse_xyz_values(values, label):
    parts = [p.strip() for p in values.split(",")]
    if len(parts) not in (2, 3):
        raise argparse.ArgumentTypeError(f"{label} values must be DX,DY or DX,DY,DZ")
    try:
        dx = float(parts[0])
        dy = float(parts[1])
        dz = float(parts[2]) if len(parts) == 3 else 0.0
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"{label} values must be numbers") from exc
    return dx, dy, dz


def parse_displacement_spec(text, label):
    try:
        key, values = text.split(":", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"{label} must be KEY:DX,DY[,DZ], e.g. BLC1a:0.5,-0.3 or BLC1a-U1:0.1,0"
        ) from exc

    key = key.strip()
    xyz = _parse_xyz_values(values, label)
    if key in BLC_DETECTORS:
        return ("detector", key, xyz)
    if is_blc_layer(key):
        return ("layer", key, xyz)
    raise argparse.ArgumentTypeError(
        f"{label} key must be a BLC detector ({', '.join(BLC_DETECTORS)}) "
        f"or layer name, e.g. BLC1a-U1"
    )


def parse_detector_offset(text):
    return parse_displacement_spec(text, "offset")


def parse_residual_offset(text):
    return parse_displacement_spec(text, "residual-offset")


class _SelectBcInAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        namespace.bcin = True
        namespace._residual_group = "bcin"


class _SelectBcOutAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        namespace.bcout = True
        namespace._residual_group = "bcout"


class _ScopedResidualAction(argparse.Action):
    _AXIS = {"--resix": 0, "--resiy": 1, "--resiz": 2}

    def __call__(self, parser, namespace, values, option_string=None):
        group = getattr(namespace, "_residual_group", None)
        if group is None:
            parser.error(f"{option_string} must follow --bcin or --bcout")
        scoped = getattr(namespace, "scoped_residual", None)
        if scoped is None:
            scoped = {"bcin": [None, None, None], "bcout": [None, None, None]}
            namespace.scoped_residual = scoped
        scoped[group][self._AXIS[option_string]] = values


def _default_scoped_residual():
    return {"bcin": [None, None, None], "bcout": [None, None, None]}


def _scoped_residual_set(scoped):
    return any(v is not None for vals in scoped.values() for v in vals)


def select_detectors(args):
    scoped = getattr(args, "scoped_residual", _default_scoped_residual())
    selected = set(args.detector or [])
    if args.bcin or any(v is not None for v in scoped["bcin"]):
        selected.update(BCIN_DETECTORS)
    if args.bcout or any(v is not None for v in scoped["bcout"]):
        selected.update(BCOUT_DETECTORS)
    for kind, key, _ in args.offset + args.residual_offset:
        if kind == "detector":
            selected.add(key)
        else:
            selected.add(key.split("-", 1)[0])
    if not selected:
        selected.update(BLC_DETECTORS)
    return tuple(det for det in BLC_DETECTORS if det in selected)


def build_override_maps(items):
    by_detector = {}
    by_layer = {}
    for kind, key, xyz in items:
        if kind == "detector":
            by_detector[key] = xyz
        else:
            by_layer[key] = xyz
    return by_detector, by_layer


def resolve_survey_displacement(record, skip_xy, override_det, override_layer, global_xyz):
    if record.name in override_layer:
        return override_layer[record.name]
    if record.detector in override_det:
        return override_det[record.detector]
    if skip_xy:
        return global_xyz
    return record.survey_x, record.survey_y, 0.0


def build_survey_maps(args):
    return build_override_maps(args.offset)


def apply_group_overrides(by_detector, detectors, dx, dy, dz):
    if dx is None and dy is None and dz is None:
        return
    for det in detectors:
        old = by_detector.get(det, (0.0, 0.0, 0.0))
        by_detector[det] = (
            dx if dx is not None else old[0],
            dy if dy is not None else old[1],
            dz if dz is not None else old[2],
        )


def build_residual_maps(args):
    scoped = getattr(args, "scoped_residual", _default_scoped_residual())
    by_detector, by_layer = build_override_maps(args.residual_offset)
    apply_group_overrides(by_detector, BCIN_DETECTORS, *scoped["bcin"])
    apply_group_overrides(by_detector, BCOUT_DETECTORS, *scoped["bcout"])
    return by_detector, by_layer


def resolve_residual_displacement(record, residual_det, residual_layer):
    if record.name in residual_layer:
        return residual_layer[record.name]
    if record.detector in residual_det:
        return residual_det[record.detector]
    return 0.0, 0.0, 0.0


def resolve_displacement(record, by_detector, by_layer, default_xyz):
    if record.name in by_layer:
        return by_layer[record.name]
    return by_detector.get(record.detector, default_xyz)


def project_displacement(dx, dy, dz, dsdx, dsdy, dsdz, scale):
    return scale * (dsdx * dx + dsdy * dy + dsdz * dz)


def compute_results(records, detectors, skip_xy, survey_maps, residual_maps, global_xyz, scale):
    selected = set(detectors)
    survey_by_det, survey_by_layer = survey_maps
    residual_by_det, residual_by_layer = residual_maps
    results = []
    for record in records:
        if record.detector not in selected:
            continue
        dx_s, dy_s, dz_s = resolve_survey_displacement(
            record, skip_xy, survey_by_det, survey_by_layer, global_xyz
        )
        dx_r, dy_r, dz_r = resolve_residual_displacement(
            record, residual_by_det, residual_by_layer
        )
        dsdx, dsdy, dsdz = calc_s_axis(record.ta, record.ra1, record.ra2)
        delta_survey = project_displacement(dx_s, dy_s, dz_s, dsdx, dsdy, dsdz, scale)
        delta_residual = project_displacement(dx_r, dy_r, dz_r, dsdx, dsdy, dsdz, scale)
        results.append(
            OffsetResult(
                record,
                dx_s,
                dy_s,
                dz_s,
                dx_r,
                dy_r,
                dz_r,
                dsdx,
                dsdy,
                dsdz,
                delta_survey,
                delta_residual,
                delta_survey + delta_residual,
            )
        )
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
        "Name", "TA", "RA1", "RA2",
        "dx_srv", "dy_srv", "dx_res", "dy_res", "dx", "dy", "dz",
        "d_srv", "d_res", "delta", "old_ofs", "new_ofs",
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
                format_float(res.dx_survey),
                format_float(res.dy_survey),
                format_float(res.dx_residual),
                format_float(res.dy_residual),
                format_float(res.dx),
                format_float(res.dy),
                format_float(res.dz),
                format_float(res.delta_ofs_survey),
                format_float(res.delta_ofs_residual),
                format_float(res.delta_ofs),
                format_float(rec.ofs),
                format_float(res.new_ofs),
            )
        )

    widths = [max(len(str(row[i])) for row in (header, *rows)) for i in range(len(header))]
    fmt = "  ".join(f"{{:<{width}}}" for width in widths)
    print(fmt.format(*header))
    print(fmt.format(*("-" * width for width in widths)))
    for row in rows:
        print(fmt.format(*row))


def print_warnings(results, detectors, write_path=None):
    if not results:
        print("No BLC records matched the selection.", file=sys.stderr)
        return False

    max_dz_input = max(abs(res.dz) for res in results)
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
    elif any(abs(res.delta_ofs) > _DSDZ_EPS for res in results):
        print(
            "Dry-run only (use --write to update Ofs). "
            "Re-applying file X,Y without --skip-xy will double-count Ofs.",
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

    rec = DCGeoRecord(
        0,
        ["0", "BLC1a-U1"] + ["0"] * 11,
        "",
        "\t",
        "\n",
    )

    dsdx, dsdy, _ = calc_s_axis(44, 0, 0)
    _assert_close(dsdx * 1.0 + dsdy * 0.0, math.cos(math.radians(44)))

    rec.fields[2] = "0.5"
    rec.fields[3] = "0.3"
    rec.fields[5] = "-45"
    rec.fields[6] = "0"
    rec.fields[7] = "0"
    rec.fields[12] = "0"
    out = compute_results(
        [rec],
        ("BLC1a",),
        False,
        ({}, {}),
        ({"BLC1a": (0.0, 0.5, 0.0)}, {"BLC1a-U1": (0.2, 0.0, 0.0)}),
        (0.0, 0.0, 0.0),
        1.0,
    )[0]
    _assert_close(out.dx_survey, 0.5)
    _assert_close(out.dy_survey, 0.3)
    _assert_close(out.dx_residual, 0.2)
    _assert_close(out.dy_residual, 0.0)
    _assert_close(out.delta_ofs, out.delta_ofs_survey + out.delta_ofs_residual)

    out_skip = compute_results(
        [rec],
        ("BLC1a",),
        True,
        ({}, {}),
        ({"BLC1a": (0.0, 0.2, 0.0)}, {}),
        (0.0, 0.0, 0.0),
        1.0,
    )[0]
    _assert_close(out_skip.dx_survey, 0.0)
    _assert_close(out_skip.dy_survey, 0.0)
    _assert_close(out_skip.dy_residual, 0.2)

    survey_maps = (
        {"BLC1a": (1.0, 0.0, 0.0), "BLC1b": (0.0, 0.0, 0.0),
         "BLC2a": (0.0, 0.0, 0.0), "BLC2b": (0.0, 0.0, 0.0)},
        {},
    )
    residual_maps = (
        {"BLC1a": (0.0, 0.5, 0.0), "BLC1b": (0.0, 0.0, 0.0),
         "BLC2a": (0.0, 0.0, 0.0), "BLC2b": (0.0, 0.0, 0.0)},
        {"BLC1a-U1": (0.2, 0.0, 0.0)},
    )
    out2 = compute_results(
        [rec],
        ("BLC1a",),
        True,
        survey_maps,
        residual_maps,
        (0.0, 0.0, 0.0),
        1.0,
    )[0]
    _assert_close(out2.dx_survey, 1.0)
    _assert_close(out2.dy_survey, 0.0)
    _assert_close(out2.dx_residual, 0.2)
    _assert_close(out2.dy_residual, 0.0)
    _assert_close(out2.delta_ofs, out2.delta_ofs_survey + out2.delta_ofs_residual)

    bcin_maps = build_residual_maps(
        argparse.Namespace(
            residual_offset=[],
            scoped_residual={
                "bcin": [None, 0.3, None],
                "bcout": [0.4, None, None],
            },
        )
    )
    _assert_close(bcin_maps[0]["BLC1a"][1], 0.3)
    _assert_close(bcin_maps[0]["BLC2a"][0], 0.4)
    _assert_close(bcin_maps[0]["BLC2a"][1], 0.0)

    print("self-test passed")


def make_arg_parser():
    epilog = """\
examples:
  %(prog)s param/DCGEO/e72/DCGeomParam_run03772_Pi --write
  %(prog)s param/DCGEO/e72/DCGeomParam_run03772_Pi \\
      --skip-xy --bcin --resiy 0.2 --write
  %(prog)s param/DCGEO/e72/DCGeomParam_run03772_Pi \\
      --write DCGeomParam_copy
  %(prog)s param/DCGEO/e72/DCGeomParam_run03772_Pi \\
      --skip-xy --bcin --resiy 0.2 --bcout --resix 0.1
  %(prog)s --self-test

After survey correction, fine-tune per-layer Ofs with:
  update_param.py <run> residual <suffix>
"""
    parser = argparse.ArgumentParser(
        description=(
            "Project BLC survey displacements from DCGEO X,Y columns (default) "
            "and/or CLI residual offsets onto Ofs.  Only Ofs is written; TA/RA "
            "are used for projection."
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
        help="with --skip-xy only: uniform X offset [mm] when no --offset",
    )
    parser.add_argument(
        "--dy",
        type=float,
        default=0.0,
        help="with --skip-xy only: uniform Y offset [mm] when no --offset",
    )
    parser.add_argument(
        "--dz",
        type=float,
        default=0.0,
        help="with --skip-xy only: uniform Z offset [mm] when no --offset",
    )
    parser.add_argument(
        "--offset",
        action="append",
        type=parse_detector_offset,
        default=[],
        metavar="DET:DX,DY[,DZ]",
        help="override file X,Y for detector/layer (--skip-xy: sole survey source)",
    )
    parser.add_argument(
        "--skip-xy",
        action="store_true",
        help="ignore X,Y columns; apply only CLI residuals (--resix/--resiy/etc.)",
    )
    parser.add_argument(
        "--residual-offset",
        action="append",
        type=parse_residual_offset,
        default=[],
        metavar="KEY:DX,DY[,DZ]",
        help=(
            "optional per-detector/layer analysis residual; "
            "prefer --bcin/--bcout with --resix/--resiy"
        ),
    )
    parser.add_argument(
        "--resix",
        action=_ScopedResidualAction,
        type=float,
        metavar="MM",
        help="uniform X residual [mm] for preceding --bcin or --bcout",
    )
    parser.add_argument(
        "--resiy",
        action=_ScopedResidualAction,
        type=float,
        metavar="MM",
        help="uniform Y residual [mm] for preceding --bcin or --bcout",
    )
    parser.add_argument(
        "--resiz",
        action=_ScopedResidualAction,
        type=float,
        metavar="MM",
        help="uniform Z residual [mm] for preceding --bcin or --bcout",
    )
    parser.add_argument(
        "--detector",
        action="append",
        choices=BLC_DETECTORS,
        help="limit target detector; can be given multiple times",
    )
    parser.add_argument(
        "--bcin",
        action=_SelectBcInAction,
        nargs=0,
        help="target BLC1a and BLC1b",
    )
    parser.add_argument(
        "--bcout",
        action=_SelectBcOutAction,
        nargs=0,
        help="target BLC2a and BLC2b",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help="multiply the projected offset before adding to Ofs; use -1 to invert sign",
    )
    parser.add_argument(
        "--write",
        nargs="?",
        const="",
        default=None,
        metavar="PATH",
        help="write updated Ofs; omit PATH to overwrite input DCGEO",
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="run internal consistency checks for calc_s_axis and exit",
    )
    return parser


def resolve_write_path(dcgeo, write_arg):
    if write_arg is None:
        return None
    return dcgeo if write_arg == "" else write_arg


def main():
    parser = make_arg_parser()
    parser.set_defaults(bcin=False, bcout=False, _residual_group=None)
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
    survey_maps = build_survey_maps(args)
    residual_maps = build_residual_maps(args)
    global_xyz = (args.dx, args.dy, args.dz)
    results = compute_results(
        records,
        detectors,
        args.skip_xy,
        survey_maps,
        residual_maps,
        global_xyz,
        args.scale,
    )

    print_results(results)
    write_path = resolve_write_path(args.dcgeo, args.write)
    ok = print_warnings(results, detectors, write_path=write_path)
    if not ok:
        sys.exit(1)

    if write_path:
        write_dcgeo(lines, write_path, results)
        if write_path == args.dcgeo:
            print(f"\nWrote updated DCGEO (in place): {write_path}")
        else:
            print(f"\nWrote updated DCGEO: {write_path}")


if __name__ == "__main__":
    main()
