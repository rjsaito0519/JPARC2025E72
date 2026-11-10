"""Post existing DC (BLC) QA PDFs to the run's Slack thread (no merge).

Filenames match BLC_TDC / BLC_DRFT / BLC_residual / BLC_pull:
  runXXXXX_BLC{1|2}_{TDC|DRIFT|residual|pull}_{S}.pdf
BcIn -> BLC1, BcOut -> BLC2. reso --all posts both.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

from . import config
from .local_notify_config import load_slack_config
from . import qa_session, slack_notify

MODES = ("t0", "drift", "resi", "reso")
_TAG = {
    "t0": "TDC",
    "drift": "DRIFT",
    "resi": "residual",
    "reso": "pull",
}


#_____________________________________________________________________________
def particle_suffix(kaon: bool) -> str:
    return "K" if kaon else "Pi"


#_____________________________________________________________________________
def blc_ids(bcin: bool = False, all_chambers: bool = False) -> list[str]:
    if all_chambers:
        return ["1", "2"]
    return ["1" if bcin else "2"]


#_____________________________________________________________________________
def pdf_name(run: int, blc: str, mode: str, kaon: bool) -> str:
    if mode not in _TAG:
        raise ValueError(f"unsupported mode {mode!r}")
    return f"run{run:05d}_BLC{blc}_{_TAG[mode]}_{particle_suffix(kaon)}.pdf"


#_____________________________________________________________________________
def source_paths(
    run: int,
    mode: str,
    *,
    kaon: bool = False,
    bcin: bool = False,
    all_chambers: bool = False,
) -> list[Path]:
    if mode not in MODES:
        raise ValueError(f"unsupported mode {mode!r}")
    if all_chambers and mode != "reso":
        all_chambers = False
    base = Path(config.OUTPUT_DIR) / "img" / f"run{run:05d}"
    return [
        base / pdf_name(run, blc, mode, kaon)
        for blc in blc_ids(bcin=bcin, all_chambers=all_chambers)
    ]


#_____________________________________________________________________________
def post_dc_qa(
    run: int,
    mode: str,
    *,
    kaon: bool = False,
    bcin: bool = False,
    all_chambers: bool = False,
    config_path: Path | None = None,
    dry_run: bool = False,
    manifest: dict[str, Any] | None = None,
    pdf_paths: list[Path] | None = None,
) -> tuple[bool, str]:
    """Post chamber QA PDF(s) for (run, mode). Missing PDF -> failure message."""
    if mode not in MODES:
        return False, f"unsupported mode {mode!r}"
    paths = list(pdf_paths) if pdf_paths is not None else source_paths(
        run, mode, kaon=kaon, bcin=bcin, all_chambers=all_chambers,
    )
    missing = [p.name for p in paths if not p.is_file()]
    if missing and not dry_run:
        return False, f"missing PDF(s): {', '.join(missing)}"

    suffix = particle_suffix(kaon)
    chamber = "all" if all_chambers and mode == "reso" else ("bcin" if bcin else "bcout")
    mode_key = f"dc_{mode}_{chamber}"
    notes = []
    for path in paths:
        caption = f"run{run:05d} dc {mode} {path.name}"
        ok, msg = slack_notify.post_pdf(
            run,
            mode_key,
            path,
            caption,
            particle_suffix=suffix,
            config_path=config_path,
            dry_run=dry_run,
            manifest=manifest,
        )
        if not ok:
            return False, msg
        notes.append(msg)
    return True, "; ".join(notes)


#_____________________________________________________________________________
def try_post_after_dcprm(
    run: int,
    mode: str,
    *,
    kaon: bool = False,
    bcin: bool = False,
    all_chambers: bool = False,
    session: dict[str, Any] | None = None,
) -> tuple[bool, str]:
    """Best-effort post from run_dcprm. Skip if Slack config missing/disabled."""
    _, status = load_slack_config()
    if status in ("missing", "disabled"):
        return True, f"slack skip: {status}"
    manifest = None
    if session is not None:
        ok, msg, hashes = qa_session.verify_sources(session)
        if not ok:
            return False, msg
        paths = [Path(p) for p in session["expected"]]
        manifest_path = qa_session.write_manifest(
            session, source_hashes=hashes, output_pdfs=paths,
        )
        manifest = qa_session.load_manifest(manifest_path)
    return post_dc_qa(
        run, mode, kaon=kaon, bcin=bcin, all_chambers=all_chambers, dry_run=False,
        manifest=manifest,
    )
