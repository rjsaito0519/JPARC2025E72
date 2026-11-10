"""Merge per-detector Hodo QA PDFs into one book per stage (page concat only).

Uses pdfunite / pdfinfo (poppler). Does not re-draw pages. Source PDFs are kept.
"""
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

from . import config

MERGE_MODES = ("hdprm", "t0", "hdphc")
PDFUNITE = shutil.which("pdfunite")
PDFINFO = shutil.which("pdfinfo")


#_____________________________________________________________________________
def img_dir(run):
    return Path(config.OUTPUT_DIR) / "img" / f"run{run:05d}"


#_____________________________________________________________________________
def particle_suffix(kaon):
    return "K" if kaon else "Pi"


#_____________________________________________________________________________
def output_name(run, mode, kaon, ftof=False):
    s = particle_suffix(kaon)
    mid = f"{mode}_ftof" if ftof and mode in ("hdprm", "hdphc") else mode
    return f"run{run:05d}_hodo_{mid}_{s}.pdf"


#_____________________________________________________________________________
def source_paths(run, mode, kaon, ftof=False):
    """Ordered list of expected source PDF paths (no glob)."""
    if mode not in MERGE_MODES:
        raise ValueError(f"unsupported mode {mode!r}")
    s = particle_suffix(kaon)
    base = img_dir(run)
    if mode == "t0":
        return [base / f"run{run:05d}_T0_Offset_{s}.pdf"]
    dets = config.hodo_detectors(mode, ftof)
    tag = "HDPRM" if mode == "hdprm" else "PHC"
    return [base / f"run{run:05d}_{det}_{tag}_{s}.pdf" for det in dets]


#_____________________________________________________________________________
def pdf_page_count(path):
    """Return page count via pdfinfo, or raise RuntimeError."""
    if not PDFINFO:
        raise RuntimeError("pdfinfo not found (poppler-utils)")
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(str(path))
    result = subprocess.run(
        [PDFINFO, str(path)], capture_output=True, text=True, timeout=30,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"pdfinfo failed for {path}: {result.stderr.strip() or result.stdout.strip()}"
        )
    match = re.search(r"^Pages:\s*(\d+)\s*$", result.stdout, re.M)
    if not match:
        raise RuntimeError(f"pdfinfo: no Pages line for {path}")
    return int(match.group(1))


#_____________________________________________________________________________
def merge(run, mode, kaon=False, ftof=False, dry_run=False, session=None):
    """Merge stage PDFs. Returns (ok: bool, message: str, out_path: Path|None).

    When session is set (auto path from run_hodo), require sources to differ from
    the begin() snapshot before uniting. Manual CLI leaves session=None.
    """
    if mode not in MERGE_MODES:
        return False, f"unsupported mode {mode!r}", None
    if ftof and mode == "t0":
        ftof = False

    sources = source_paths(run, mode, kaon, ftof)
    out = img_dir(run) / output_name(run, mode, kaon, ftof)

    if session is not None:
        from . import qa_session
        ok, message, _hashes = qa_session.verify_sources(session)
        if not ok:
            return False, message, out
    else:
        missing = [p for p in sources if not p.is_file()]
        if missing:
            names = ", ".join(p.name for p in missing)
            return False, f"missing PDF(s): {names}", out

    page_counts = []
    for path in sources:
        try:
            page_counts.append(pdf_page_count(path))
        except (RuntimeError, FileNotFoundError) as e:
            return False, f"unreadable PDF: {e}", out
    total_pages = sum(page_counts)

    if dry_run:
        listing = "; ".join(f"{p.name}({n})" for p, n in zip(sources, page_counts))
        return True, f"dry-run: {listing} -> {out.name} ({total_pages} pages)", out

    if not PDFUNITE:
        return False, "pdfunite not found (poppler-utils)", out

    out.parent.mkdir(parents=True, exist_ok=True)
    if len(sources) == 1:
        shutil.copy2(sources[0], out)
    else:
        cmd = [PDFUNITE, *[str(p) for p in sources], str(out)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            return False, (
                f"pdfunite failed: {result.stderr.strip() or result.stdout.strip()}"
            ), out

    try:
        out_pages = pdf_page_count(out)
    except RuntimeError as e:
        return False, f"merged PDF unreadable: {e}", out
    if out_pages != total_pages:
        return False, (
            f"page count mismatch: sources={total_pages} merged={out_pages} ({out})"
        ), out

    listing = "+".join(str(n) for n in page_counts)
    return True, f"wrote {out} ({listing}={out_pages} pages)", out
