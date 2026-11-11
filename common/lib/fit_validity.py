"""Finite / fit-status checks before writing calibration parameters."""
from __future__ import annotations

import math
from typing import Any, Iterable, Optional


# ROOT Fit() composite status: 4000 is an IMPROVE-only failure.
# MIGRAD / MINOS / HESSE succeeded; finite parameter and statistics checks remain required.
PHC_FIT_STATUS_OK = 0
PHC_FIT_STATUSES_USABLE = frozenset({0, 4000})


#_____________________________________________________________________________
def is_finite_number(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


#_____________________________________________________________________________
def all_finite(values: Iterable[Any]) -> bool:
    return all(is_finite_number(v) for v in values)


#_____________________________________________________________________________
def phc_fit_status_ok(status: Optional[Any]) -> bool:
    """Allow successful primary fits, including an IMPROVE-only failure."""
    if status is None:
        return True
    try:
        return int(status) in PHC_FIT_STATUSES_USABLE
    except (TypeError, ValueError):
        return False
