"""Index existing calib QA PDFs for final approval (#25). Does not rename files."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from . import config


#_____________________________________________________________________________
def img_dir(run: int) -> Path:
    return Path(config.OUTPUT_DIR) / "img" / f"run{int(run):05d}"


#_____________________________________________________________________________
def list_qa_pdfs(run: int) -> list[dict[str, Any]]:
    """List known QA PDFs under results/img/runXXXXX (no rename)."""
    root = img_dir(run)
    if not root.is_dir():
        return []
    patterns = (
        f"run{int(run):05d}_hodo_*.pdf",
        f"run{int(run):05d}_*_HDPRM_*.pdf",
        f"run{int(run):05d}_*_PHC_*.pdf",
        f"run{int(run):05d}_*_T0_Offset_*.pdf",
        f"run{int(run):05d}_BLC*_*.pdf",
        f"d5_pos_review_run{int(run):05d}_*.pdf",
        f"d5_plane_resi_qa_run{int(run):05d}_*.pdf",
    )
    seen: set[Path] = set()
    rows: list[dict[str, Any]] = []
    for pat in patterns:
        for path in sorted(root.glob(pat)):
            if path in seen or not path.is_file():
                continue
            seen.add(path)
            rows.append({
                "name": path.name,
                "path": str(path),
                "mtime": path.stat().st_mtime,
                "size": path.stat().st_size,
            })
    rows.sort(key=lambda r: r["name"])
    return rows


#_____________________________________________________________________________
def index_summary(run: int) -> dict[str, Any]:
    pdfs = list_qa_pdfs(run)
    return {
        "run": int(run),
        "img_dir": str(img_dir(run)),
        "n_pdf": len(pdfs),
        "pdfs": pdfs,
    }
