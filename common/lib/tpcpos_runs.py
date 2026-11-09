# -*- coding: utf-8 -*-
"""Beam-through TPCPOS run table (run, charge, beam momentum, species).

Shared by analysis/scripts/tpcpos/*.py. Formerly defined in
macro/compare_tpcpos_beamthrough.py (PDF compare scripts now under archive/).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List


@dataclass
class RunSpec:
    run: int
    expect_charge: int  # +1 / -1 from run polarity
    p_beam: float
    species: str
    label: str = ""

    def __post_init__(self) -> None:
        if not self.label:
            sign = "+" if self.expect_charge > 0 else "-"
            self.label = f"{self.species}{sign} {self.p_beam:.3f} GeV/c run{self.run:05d}"


# All 22 runs from scratch/dst_tpchelix_02489-02537_22runs_tpcpos0.yml
RUNS: List[RunSpec] = [
    RunSpec(2489, -1, 1.000, "pi"),
    RunSpec(2491, -1, 1.000, "p"),
    RunSpec(2494, -1, 0.814, "p"),
    RunSpec(2496, -1, 0.814, "p"),
    RunSpec(2500, -1, 0.814, "K"),
    RunSpec(2502, -1, 0.814, "pi"),
    RunSpec(2504, -1, 0.645, "p"),
    RunSpec(2506, -1, 0.645, "K"),
    RunSpec(2508, -1, 0.645, "pi"),
    RunSpec(2509, -1, 0.645, "pi"),
    RunSpec(2512, -1, 0.400, "pi"),
    RunSpec(2514, -1, 0.300, "pi"),
    RunSpec(2516, +1, 1.000, "pi"),
    RunSpec(2518, +1, 1.000, "p"),
    RunSpec(2520, +1, 0.814, "pi"),
    RunSpec(2522, +1, 0.814, "p"),
    RunSpec(2524, +1, 0.645, "pi"),
    RunSpec(2527, +1, 0.645, "p"),
    RunSpec(2529, +1, 0.400, "pi"),
    RunSpec(2531, +1, 0.400, "p"),
    RunSpec(2535, +1, 0.400, "pi"),
    RunSpec(2537, +1, 0.300, "pi"),
]
