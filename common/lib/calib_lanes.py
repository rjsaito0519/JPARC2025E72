"""Calib dispatch lane split (hodo / bcout / bcin / d5)."""
from __future__ import annotations

import os
from typing import Any, Optional

LANE_HODO = "hodo"
LANE_BCOUT = "bcout"
LANE_BCIN = "bcin"
LANE_D5 = "d5"
LANE_ORDER = (LANE_HODO, LANE_BCOUT, LANE_BCIN, LANE_D5)
LANE_ENV = "E72_CALIB_LANE"

HODO_NAMES = frozenset({
    "hodo", "hodo_hdprm", "hodo_t0", "hodo_hdphc", "hodo_hdphc_adopt",
})
DC_NAMES = frozenset({
    "dc_chamber", "dc_t0", "dc_drift", "dc_ch_resi", "dc_reso", "dc_resi",
})
D5_NAMES = frozenset({"d5_align", "d5_pass"})


#_____________________________________________________________________________
def active_lane() -> Optional[str]:
    """Lane id from env (lane worker), or None for coordinator / legacy."""
    text = str(os.environ.get(LANE_ENV) or "").strip()
    return text or None


#_____________________________________________________________________________
def step_name_opts(raw) -> tuple[str, dict[str, Any]]:
    if isinstance(raw, str):
        return raw, {}
    if isinstance(raw, dict) and len(raw) == 1:
        name, opts = next(iter(raw.items()))
        return str(name), dict(opts or {})
    return "?", {}


#_____________________________________________________________________________
def lane_of_step(raw) -> Optional[str]:
    """Map one session/dispatch step to a lane id."""
    name, opts = step_name_opts(raw)
    if name in HODO_NAMES:
        return LANE_HODO
    if name in D5_NAMES:
        return LANE_D5
    if name in DC_NAMES:
        if opts.get("bcin"):
            return LANE_BCIN
        return LANE_BCOUT
    return None


#_____________________________________________________________________________
def split_steps_by_lane(steps: list) -> dict[str, list]:
    """Preserve order within each lane; omit empty lanes."""
    out: dict[str, list] = {lane: [] for lane in LANE_ORDER}
    for raw in steps:
        lane = lane_of_step(raw)
        if lane is None:
            continue
        out[lane].append(raw)
    return {lane: rows for lane, rows in out.items() if rows}


#_____________________________________________________________________________
def d5_ready(lane_status: dict[str, str], lane_steps: dict[str, list]) -> bool:
    """True when d5 may start: every present DC lane is done."""
    for lane in (LANE_BCOUT, LANE_BCIN):
        if lane not in lane_steps:
            continue
        if lane_status.get(lane) != "done":
            return False
    return True


#_____________________________________________________________________________
def parallel_ready_lanes(lane_steps: dict[str, list]) -> list[str]:
    """Lanes that can start immediately (not d5)."""
    return [lane for lane in (LANE_HODO, LANE_BCOUT, LANE_BCIN) if lane in lane_steps]
