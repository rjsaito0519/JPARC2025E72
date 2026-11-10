"""Beam run condition lookup for calibration sessions.

Reads myanalysis calibration/config/beam_run_condition.yml (copy of analyzer
BeamRunCondition ranges). Optional UserParam peek for BeamParticle / BeamMom /
BeamCharge. Does not write analyzer tables or UserParam.
"""
from __future__ import annotations

from pathlib import Path

import yaml

from . import config

PARTICLES = ("pion", "kaon", "proton")
PARTICLE_FROM_CODE = {0: "pion", 1: "kaon", 2: "proton"}
SOURCE_NONE = "none"
SOURCE_TABLE = "run_table"
SOURCE_USER_PARAM = "user_param"
SOURCE_USER_INPUT = "user_input"

DEFAULT_TABLE = (
    Path(__file__).resolve().parents[2] / "calibration" / "config" / "beam_run_condition.yml"
)


#_____________________________________________________________________________
def default_table_path():
    return DEFAULT_TABLE


#_____________________________________________________________________________
def load_table(path=None):
    path = Path(path) if path is not None else default_table_path()
    if not path.is_file():
        raise FileNotFoundError(f"beam run table not found: {path}")
    doc = yaml.safe_load(path.read_text())
    if not isinstance(doc, dict) or not isinstance(doc.get("ranges"), list):
        raise ValueError(f"beam run table must have a ranges list: {path}")
    ranges = []
    for i, row in enumerate(doc["ranges"]):
        if not isinstance(row, dict):
            raise ValueError(f"ranges[{i}] must be a mapping")
        for key in ("first", "last", "particle", "momentum_gev_c", "charge"):
            if key not in row:
                raise ValueError(f"ranges[{i}] missing {key}")
        particle = str(row["particle"]).lower()
        if particle not in PARTICLES:
            raise ValueError(f"ranges[{i}] bad particle {particle!r}")
        first = int(row["first"])
        last = int(row["last"])
        if first > last:
            raise ValueError(f"ranges[{i}] first>last")
        mom = float(row["momentum_gev_c"])
        charge = int(row["charge"])
        if charge not in (-1, 1):
            raise ValueError(f"ranges[{i}] bad charge {charge}")
        ranges.append({
            "first": first,
            "last": last,
            "particle": particle,
            "momentum_gev_c": mom,
            "charge": charge,
        })
    return path, ranges


#_____________________________________________________________________________
def find_in_ranges(ranges, run):
    """Return the covering range row or None (no neighbor guess)."""
    for row in ranges:
        if row["first"] <= run <= row["last"]:
            return dict(row)
    return None


#_____________________________________________________________________________
def user_param_paths(run):
    """Candidate UserParam paths (first existing wins for peek)."""
    sub = config.SUB_DIR
    return [
        config.PARAM_DIR / "USER" / sub / f"UserParam_run{run:05d}",
        config.PARAM_DIR / "USER" / f"UserParam_run{run:05d}",
    ]


#_____________________________________________________________________________
def _parse_user_param_file(path):
    """Return dict of Beam* keys found (raw strings/ints/floats)."""
    out = {}
    for raw in path.read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.replace("\t", " ").split()
        if len(parts) < 2:
            continue
        key, value = parts[0], parts[1]
        if key in ("BeamParticle", "BeamMom", "BeamCharge"):
            out[key] = value
    return out


#_____________________________________________________________________________
def peek_user_param(run):
    """Return overrides from UserParam if present (does not write)."""
    found = {}
    used_path = None
    for path in user_param_paths(run):
        if not path.is_file():
            continue
        raw = _parse_user_param_file(path)
        if not raw:
            continue
        used_path = path
        if "BeamParticle" in raw:
            code = int(float(raw["BeamParticle"]))
            if code not in PARTICLE_FROM_CODE:
                raise ValueError(f"{path}: bad BeamParticle {raw['BeamParticle']!r}")
            found["particle"] = PARTICLE_FROM_CODE[code]
        if "BeamMom" in raw:
            mom = float(raw["BeamMom"])
            if not (mom > 0):
                raise ValueError(f"{path}: bad BeamMom {raw['BeamMom']!r}")
            found["momentum_gev_c"] = mom
        if "BeamCharge" in raw:
            charge = int(float(raw["BeamCharge"]))
            if charge not in (-1, 1):
                raise ValueError(f"{path}: bad BeamCharge {raw['BeamCharge']!r}")
            found["charge"] = charge
        break
    return found, used_path


#_____________________________________________________________________________
def resolve(run, table_path=None):
    """Resolve beam condition. Priority: UserParam > table. No neighbor guess.

    Returns dict with particle/momentum/charge and *_source fields.
    """
    _, ranges = load_table(table_path)
    table_hit = find_in_ranges(ranges, run)
    user_hit, user_path = peek_user_param(run)

    particle = None
    particle_source = SOURCE_NONE
    momentum = None
    momentum_source = SOURCE_NONE
    charge = None
    charge_source = SOURCE_NONE

    if table_hit:
        particle = table_hit["particle"]
        particle_source = SOURCE_TABLE
        momentum = table_hit["momentum_gev_c"]
        momentum_source = SOURCE_TABLE
        charge = table_hit["charge"]
        charge_source = SOURCE_TABLE

    if "particle" in user_hit:
        particle = user_hit["particle"]
        particle_source = SOURCE_USER_PARAM
    if "momentum_gev_c" in user_hit:
        momentum = user_hit["momentum_gev_c"]
        momentum_source = SOURCE_USER_PARAM
    if "charge" in user_hit:
        charge = user_hit["charge"]
        charge_source = SOURCE_USER_PARAM

    return {
        "run": run,
        "particle": particle,
        "particle_source": particle_source,
        "momentum_gev_c": momentum,
        "momentum_source": momentum_source,
        "charge": charge,
        "charge_source": charge_source,
        "user_param_path": str(user_path) if user_path else None,
        "table_hit": table_hit is not None,
    }


#_____________________________________________________________________________
def upsert_run(run, particle, momentum_gev_c, charge, table_path=None):
    """Insert or update a single-run row in the myanalysis table.

    - Same values already covering the run: no write, returns 'unchanged'
    - Existing first==last==run row: overwrite
    - Run inside a multi-run range with different values: raise ValueError
    - Otherwise append first=last=run
    """
    particle = str(particle).lower()
    if particle not in PARTICLES:
        raise ValueError(f"bad particle {particle!r}")
    momentum_gev_c = float(momentum_gev_c)
    if not (momentum_gev_c > 0):
        raise ValueError(f"bad momentum {momentum_gev_c}")
    charge = int(charge)
    if charge not in (-1, 1):
        raise ValueError(f"bad charge {charge}")

    path, ranges = load_table(table_path)
    wanted = {
        "first": run,
        "last": run,
        "particle": particle,
        "momentum_gev_c": momentum_gev_c,
        "charge": charge,
    }

    for i, row in enumerate(ranges):
        if not (row["first"] <= run <= row["last"]):
            continue
        same = (
            row["particle"] == particle
            and abs(row["momentum_gev_c"] - momentum_gev_c) < 1e-9
            and row["charge"] == charge
        )
        if same:
            return path, "unchanged"
        if row["first"] == row["last"] == run:
            ranges[i] = wanted
            _write_ranges(path, ranges)
            return path, "updated"
        raise ValueError(
            f"run {run} is inside range {row['first']}-{row['last']} with different "
            f"values; edit {path} by hand"
        )

    ranges.append(wanted)
    ranges.sort(key=lambda r: (r["first"], r["last"]))
    _write_ranges(path, ranges)
    return path, "inserted"


#_____________________________________________________________________________
def _write_ranges(path, ranges):
    lines = [
        "# Copy of analyzer src/BeamRunCondition.cc kBeamRunRanges (same content).",
        "# When analyzer table is regenerated, refresh this file to match.",
        "# calib_session --write-table may append single-run rows here only (not the .cc).",
        "ranges:",
    ]
    for row in ranges:
        lines.append(
            "  - {"
            f"first: {row['first']}, last: {row['last']}, "
            f"particle: {row['particle']}, "
            f"momentum_gev_c: {row['momentum_gev_c']}, "
            f"charge: {row['charge']}"
            "}"
        )
    path.write_text("\n".join(lines) + "\n")
