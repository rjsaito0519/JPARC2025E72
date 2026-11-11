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
# Only in-process work conflicts with a new start of the same lane.
RUNNING_STATUSES = frozenset({"running"})
HELD_STATUSES = frozenset({"held"})
# Legacy alias: treat as "blocks --force restart", not attach/resume.
BUSY_STATUSES = frozenset({"running", "held"})
# Coordinator / attach / D5 wait: worker finished its steps (human may still approve).
COMPLETE_FOR_DEPS_STATUSES = frozenset({
    "done",
    "awaiting_final_approval",
    "approved",
    "finished",
})
# Attach only starts lanes with no prior status in this set.
ATTACH_BLOCKED_STATUSES = frozenset({
    "running", "held", "stopped", "done", "failed",
    "awaiting_final_approval", "approved", "finished",
})
# Coordinator local map: no more child process for this lane.
COORD_TERMINAL_STATUSES = frozenset({
    "done", "failed", "stopped", "held",
})

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
def lane_complete_for_deps(status: Optional[str]) -> bool:
    """True when a lane's work is finished for D5 / attach skip / prior_lanes.

    Includes awaiting_final_approval (worker done, human approval pending).
    Does not treat held/stopped/failed as complete.
    """
    return str(status or "") in COMPLETE_FOR_DEPS_STATUSES


#_____________________________________________________________________________
def normalize_status_for_deps(status: Optional[str]) -> str:
    """Map disk status to coordinator-facing value used in deps checks."""
    text = str(status or "")
    if lane_complete_for_deps(text):
        return "done"
    return text


#_____________________________________________________________________________
def d5_ready(lane_status: dict[str, str], lane_steps: dict[str, list]) -> bool:
    """True when d5 may start: every present DC lane finished its work."""
    for lane in (LANE_BCOUT, LANE_BCIN):
        if lane not in lane_steps:
            continue
        if not lane_complete_for_deps(lane_status.get(lane)):
            return False
    return True


#_____________________________________________________________________________
def parallel_ready_lanes(lane_steps: dict[str, list]) -> list[str]:
    """Lanes that can start immediately (not d5)."""
    return [lane for lane in (LANE_HODO, LANE_BCOUT, LANE_BCIN) if lane in lane_steps]


#_____________________________________________________________________________
def requested_lanes(steps: list) -> list[str]:
    """Lane ids implied by session/CLI steps, in LANE_ORDER."""
    split = split_steps_by_lane(steps)
    return [lane for lane in LANE_ORDER if lane in split]


#_____________________________________________________________________________
def is_busy_status(status: Optional[str]) -> bool:
    return str(status or "") in BUSY_STATUSES


#_____________________________________________________________________________
def is_running_status(status: Optional[str]) -> bool:
    return str(status or "") in RUNNING_STATUSES


#_____________________________________________________________________________
def is_held_status(status: Optional[str]) -> bool:
    return str(status or "") in HELD_STATUSES


#_____________________________________________________________________________
def busy_lanes(status_map: dict[str, str]) -> set[str]:
    """Lanes that are running or held (blocks full --force restart)."""
    return {
        str(lane)
        for lane, status in (status_map or {}).items()
        if is_busy_status(status)
    }


#_____________________________________________________________________________
def running_lanes(status_map: dict[str, str]) -> set[str]:
    """Lanes currently executing (cannot start/attach the same lane)."""
    return {
        str(lane)
        for lane, status in (status_map or {}).items()
        if is_running_status(status)
    }


#_____________________________________________________________________________
def held_lanes(status_map: dict[str, str]) -> set[str]:
    """Lanes waiting on decide (resume OK while other lanes run)."""
    return {
        str(lane)
        for lane, status in (status_map or {}).items()
        if is_held_status(status)
    }


#_____________________________________________________________________________
def stopped_lanes(status_map: dict[str, str]) -> set[str]:
    return {
        str(lane)
        for lane, status in (status_map or {}).items()
        if str(status or "") == "stopped"
    }


#_____________________________________________________________________________
def free_requested_lanes(
    want: list[str] | set[str], status_map: dict[str, str],
) -> list[str]:
    """Wanted lanes with no prior work (attach candidates).

    Stopped / held / complete / failed lanes need continue/resume, not attach.
    """
    want_set = {str(x) for x in want}
    out = []
    for lane in LANE_ORDER:
        if lane not in want_set:
            continue
        status = str((status_map or {}).get(lane) or "")
        if status in ATTACH_BLOCKED_STATUSES:
            continue
        out.append(lane)
    return out


#_____________________________________________________________________________
def overlap_lanes(want: list[str] | set[str], busy: set[str]) -> list[str]:
    """Busy lanes that the request also wants, in LANE_ORDER."""
    want_set = {str(x) for x in want}
    return [lane for lane in LANE_ORDER if lane in want_set and lane in busy]
