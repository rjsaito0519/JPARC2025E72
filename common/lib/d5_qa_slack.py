"""Post existing D5 QA PDFs to the run's Slack thread (no merge).

Filenames match d5_pos_review / d5_plane_resi_qa:
  d5_pos_review_runXXXXX_{S}.pdf
  d5_plane_resi_qa_runXXXXX_{S}.pdf
"""
from __future__ import annotations

from pathlib import Path

from . import config
from .local_notify_config import load_slack_config
from . import slack_notify

KINDS = ("pos_review", "plane_resi_qa")
_PREFIX = {
    "pos_review": "d5_pos_review",
    "plane_resi_qa": "d5_plane_resi_qa",
}


#_____________________________________________________________________________
def particle_suffix(kaon: bool) -> str:
    return "K" if kaon else "Pi"


#_____________________________________________________________________________
def pdf_name(run: int, kind: str, kaon: bool) -> str:
    if kind not in _PREFIX:
        raise ValueError(f"unsupported kind {kind!r}")
    return f"{_PREFIX[kind]}_run{run:05d}_{particle_suffix(kaon)}.pdf"


#_____________________________________________________________________________
def source_path(run: int, kind: str, *, kaon: bool = False) -> Path:
    if kind not in KINDS:
        raise ValueError(f"unsupported kind {kind!r}")
    base = Path(config.OUTPUT_DIR) / "img" / f"run{run:05d}"
    return base / pdf_name(run, kind, kaon)


#_____________________________________________________________________________
def post_d5_qa(
    run: int,
    kind: str,
    *,
    kaon: bool = False,
    config_path: Path | None = None,
    dry_run: bool = False,
) -> tuple[bool, str]:
    """Post one D5 QA PDF for (run, kind). Missing PDF -> failure (unless dry-run)."""
    if kind not in KINDS:
        return False, f"unsupported kind {kind!r}"
    path = source_path(run, kind, kaon=kaon)
    if not dry_run and not path.is_file():
        return False, f"missing PDF: {path.name}"

    suffix = particle_suffix(kaon)
    mode_key = f"d5_{kind}"
    caption = f"run{run:05d} d5 {kind} ({suffix})"
    return slack_notify.post_pdf(
        run,
        mode_key,
        path,
        caption,
        particle_suffix=suffix,
        config_path=config_path,
        dry_run=dry_run,
    )


#_____________________________________________________________________________
def try_post_after_d5(
    run: int,
    kind: str,
    *,
    kaon: bool = False,
) -> tuple[bool, str]:
    """Best-effort post from run_d5. Skip if Slack config missing/disabled."""
    _, status = load_slack_config()
    if status in ("missing", "disabled"):
        return True, f"slack skip: {status}"
    return post_d5_qa(run, kind, kaon=kaon, dry_run=False)
