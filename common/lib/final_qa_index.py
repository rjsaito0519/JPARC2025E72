"""Index existing calib QA PDFs for final approval (#25). Does not rename files."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Optional

from . import config

_FAMILY_PATTERNS = {
    "hodo": (
        "hodo_*.pdf",
        "*_HDPRM_*.pdf",
        "*_PHC_*.pdf",
        "*_T0_Offset_*.pdf",
    ),
    "dc": (
        "BLC*_*.pdf",
    ),
    "d5": (
        "d5_pos_review_run*_*.pdf",
        "d5_plane_resi_qa_run*_*.pdf",
    ),
}


#_____________________________________________________________________________
def img_dir(run: int) -> Path:
    return Path(config.OUTPUT_DIR) / "img" / f"run{int(run):05d}"


#_____________________________________________________________________________
def _particle_tag(particle: Optional[str]) -> Optional[str]:
    if particle == "kaon":
        return "K"
    if particle == "pion":
        return "Pi"
    return None


#_____________________________________________________________________________
def list_qa_pdfs(
    run: int,
    *,
    particle: Optional[str] = None,
    families: Optional[Iterable[str]] = None,
) -> list[dict[str, Any]]:
    """List known QA PDFs under results/img/runXXXXX (no rename).

    Optional filters: calib particle (kaon/pion) and families (hodo/dc/d5).
    """
    root = img_dir(run)
    if not root.is_dir():
        return []
    fams = list(families) if families is not None else ["hodo", "dc", "d5"]
    tag = _particle_tag(particle)
    run_i = int(run)
    seen: set[Path] = set()
    rows: list[dict[str, Any]] = []
    for family in fams:
        patterns = _FAMILY_PATTERNS.get(family) or ()
        for pat in patterns:
            # Prefix run-scoped names where the pattern is relative.
            if pat.startswith("d5_"):
                glob_pat = pat.replace("run*", f"run{run_i:05d}", 1)
            else:
                glob_pat = f"run{run_i:05d}_{pat}"
            for path in sorted(root.glob(glob_pat)):
                if path in seen or not path.is_file():
                    continue
                if tag is not None and f"_{tag}.pdf" not in path.name:
                    continue
                seen.add(path)
                rows.append({
                    "name": path.name,
                    "path": str(path),
                    "mtime": path.stat().st_mtime,
                    "size": path.stat().st_size,
                    "family": family,
                })
    rows.sort(key=lambda r: r["name"])
    return rows


#_____________________________________________________________________________
def index_summary(
    run: int,
    *,
    particle: Optional[str] = None,
    families: Optional[Iterable[str]] = None,
) -> dict[str, Any]:
    pdfs = list_qa_pdfs(run, particle=particle, families=families)
    return {
        "run": int(run),
        "img_dir": str(img_dir(run)),
        "n_pdf": len(pdfs),
        "pdfs": pdfs,
        "particle": particle,
        "families": list(families) if families is not None else None,
    }
