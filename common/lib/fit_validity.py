"""Finite / fit-status checks before writing calibration parameters."""
from __future__ import annotations

import math
from typing import Any, Iterable, Optional


# ROOT TMinuit Fit() return: 0 = OK (see common/src/hodo.cpp write_phc_tree).
PHC_FIT_STATUS_OK = 0


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
    """True when status is missing (legacy tree) or equals Fit OK."""
    if status is None:
        return True
    try:
        return int(status) == PHC_FIT_STATUS_OK
    except (TypeError, ValueError):
        return False
