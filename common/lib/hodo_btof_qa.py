"""Draw CBtof0 QA pages via bin/CBtof0_QA (ROOT) for hdphc books."""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Optional

from . import hodo_btof_monitor, hodo_qa_pdf, runner


#_____________________________________________________________________________
def btof0_pdf_path(run: int, *, kaon: bool = False) -> Path:
    suffix = hodo_qa_pdf.particle_suffix(kaon)
    return hodo_qa_pdf.img_dir(run) / f"run{int(run):05d}_hodo_btof0_{suffix}.pdf"


#_____________________________________________________________________________
def draw_btof0_pdf(
    run: int,
    *,
    kaon: bool = False,
    root_path: Path | None = None,
    out_path: Path | None = None,
) -> tuple[bool, str, Optional[Path]]:
    """Run CBtof0_QA. Returns (ok, msg, path). out_path is ignored (bin path)."""
    del out_path  # binary writes the canonical img path from paths.h / run
    root = Path(root_path) if root_path is not None else hodo_btof_monitor.hodo_root(run)
    if not root.is_file():
        return False, f"Hodo ROOT missing: {root}", None
    binary = runner.BIN_DIR / "CBtof0_QA"
    if not binary.is_file():
        return False, "CBtof0_QA binary missing (run build.sh)", None

    suffix = hodo_qa_pdf.particle_suffix(kaon)
    cmd = [str(binary), str(root), suffix, str(int(run))]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=120, check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, f"CBtof0_QA failed: {type(exc).__name__}: {exc}", None

    out = btof0_pdf_path(run, kaon=kaon)
    if result.returncode != 0 or not out.is_file():
        err = (result.stderr or result.stdout or "").strip()
        if len(err) > 400:
            err = err[-400:]
        return False, f"CBtof0_QA failed (rc={result.returncode}): {err}", None
    note = (result.stdout or "").strip().splitlines()
    msg = note[-1] if note else f"wrote {out}"
    return True, msg, out


#_____________________________________________________________________________
def ensure_btof0_pdf(
    run: int,
    *,
    kaon: bool = False,
    root_path: Path | None = None,
    force: bool = False,
) -> tuple[bool, str, Optional[Path]]:
    """Return existing btof0 PDF or draw a new one via ROOT."""
    out = btof0_pdf_path(run, kaon=kaon)
    root = (
        Path(root_path) if root_path is not None
        else hodo_btof_monitor.hodo_root(run)
    )
    if out.is_file() and not force:
        if root.is_file() and out.stat().st_mtime >= root.stat().st_mtime:
            return True, f"reuse {out.name}", out
    return draw_btof0_pdf(
        run, kaon=kaon, root_path=root_path, out_path=out,
    )
