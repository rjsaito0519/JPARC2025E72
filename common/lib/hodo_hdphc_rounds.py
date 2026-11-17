"""Snapshot / journal / adopt helpers for multi-pass hdphc."""
from __future__ import annotations

import json
import hashlib
import math
import os
import tempfile
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from . import config, qa_session
from .phc_comparison import compare_candidates

DEFAULT_MAX_ROUNDS = 10
DEFAULT_PATIENCE = 3


#_____________________________________________________________________________
def rounds_dir(run: int) -> Path:
    return Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}" / "hdphc_rounds"


#_____________________________________________________________________________
def journal_path(run: int) -> Path:
    return rounds_dir(run) / "journal.jsonl"


#_____________________________________________________________________________
def summary_path(run: int) -> Path:
    return rounds_dir(run) / "summary.txt"


#_____________________________________________________________________________
def round_dir(run: int, round_i: int) -> Path:
    return rounds_dir(run) / f"round_{int(round_i):02d}"


#_____________________________________________________________________________
def hdphc_param_path(run: int, suffix: str) -> Path:
    return (
        Path(config.PARAM_DIR)
        / "HDPHC"
        / config.SUB_DIR
        / f"HodoPHCParam_run{int(run):05d}_{suffix}"
    )


#_____________________________________________________________________________
def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


#_____________________________________________________________________________
def clear_rounds(run: int, *, reset_range_confirmed: bool = False) -> None:
    """Reset hdphc round journal/snapshots for a new multi-pass series.

    Keeps session PHC range_confirmed by default so Go → round1 refit does not
    re-open the range ask. Pass reset_range_confirmed=True only for an explicit
    fresh start (restart / clear_dispatch already clears via clear_phc_ask_state).
    """
    path = rounds_dir(run)
    if path.is_dir():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)
    if reset_range_confirmed:
        from . import phc_range_gate

        phc_range_gate.clear_range_confirmed(int(run))


#_____________________________________________________________________________
def snapshot_round(
    run: int,
    round_i: int,
    suffix: str,
    score: dict[str, Any],
    *,
    input_path: Path | None = None,
    qa_paths=None,
) -> dict[str, Any]:
    """Copy current HDPHC param (and optional fit QA PDFs) into round_XX/."""
    dest = round_dir(run, round_i)
    dest.mkdir(parents=True, exist_ok=True)
    src = hdphc_param_path(run, suffix)
    if not src.is_file():
        raise FileNotFoundError(f"HDPHC param missing: {src}")
    shutil.copy2(src, dest / src.name)
    param_sha = qa_session.file_sha256(src)
    tok = None
    if input_path is not None and Path(input_path).is_file():
        tok = qa_session.input_token(input_path)
    qa_files: list[dict[str, Any]] = []
    if qa_paths:
        qa_dir = dest / "qa"
        qa_dir.mkdir(parents=True, exist_ok=True)
        for path in qa_paths:
            path = Path(path)
            if not path.is_file():
                qa_files.append({
                    "name": path.name,
                    "sha": None,
                    "missing": True,
                })
                continue
            sealed = qa_dir / path.name
            shutil.copy2(path, sealed)
            qa_files.append({
                "name": path.name,
                "sha": qa_session.file_sha256(sealed),
                "path": str(sealed),
                "missing": False,
            })
    entry = {
        "round": int(round_i),
        "width": score.get("width"),
        "n_entries": score.get("n_entries"),
        "ok": bool(score.get("ok")),
        "input_token": tok,
        "param_sha": param_sha,
        "param_file": src.name,
        "qa_kind": "fit",
        "qa_files": qa_files,
        "at": _utc_now(),
        "per_det": score.get("per_det") or {},
        "scoring_context": score.get("scoring_context"),
    }
    if score.get("timings"):
        entry["timings"] = score["timings"]
    with journal_path(run).open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry, sort_keys=True) + "\n")
    return entry


#_____________________________________________________________________________
def load_journal(run: int) -> list[dict[str, Any]]:
    path = journal_path(run)
    if not path.is_file():
        return []
    rows = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        entry = json.loads(line)
        round_i = int(entry["round"])
        # Refit of a round invalidates every later snapshot's dependency chain.
        rows = {index: row for index, row in rows.items() if index < round_i}
        rows[round_i] = entry
    return [rows[index] for index in sorted(rows)]


#_____________________________________________________________________________
def journal_signature(journal):
    return hashlib.sha256(json.dumps(journal, sort_keys=True).encode()).hexdigest()


#_____________________________________________________________________________
def request_finish(run, *, hold_id, user_id=None):
    """Pin an operator's decision to finish this exact, successfully scored series."""
    from .calib_state_io import write_json
    journal = load_journal(run)
    if not journal or [int(row["round"]) for row in journal] != list(range(1, len(journal) + 1)):
        raise ValueError("PHC の連続した評価記録がありません。未完了ラウンドを確認してください。")
    if any(not row.get("ok") or not width_ok(row.get("width")) for row in journal):
        raise ValueError("PHC の幅評価が未完了なので、正常終了できません。")
    chosen = pick_round(journal, "auto")
    selected = next(row for row in journal if int(row["round"]) == chosen)
    snapshot = round_dir(run, chosen) / selected["param_file"]
    if not snapshot.is_file() or qa_session.file_sha256(snapshot) != selected.get("param_sha"):
        raise ValueError("最良ラウンドのパラメータ保存を確認できません。")
    record = {"at": _utc_now(), "hold_id": hold_id, "user_id": user_id,
              "completed_round": int(journal[-1]["round"]), "best_round": chosen,
              "journal_signature": journal_signature(journal)}
    write_json(rounds_dir(run) / "finish_request.json", record)
    return record


#_____________________________________________________________________________
def finish_requested(run, journal):
    """Refitting or resetting the series invalidates an earlier finish request."""
    try:
        record = json.loads((rounds_dir(run) / "finish_request.json").read_text())
    except (OSError, ValueError):
        return False
    return bool(journal) and record.get("journal_signature") == journal_signature(journal)


#_____________________________________________________________________________
def adaptive_progress(journal, *, patience=DEFAULT_PATIENCE, max_rounds=DEFAULT_MAX_ROUNDS):
    """Track best width; keep updating live params until patience or the cap."""
    if patience < 1 or max_rounds < 1:
        raise ValueError("PHC patience and round cap must be positive")
    best_width = None
    best_round = None
    without_improvement = 0
    comparisons = compare_candidates(journal)
    for entry in comparisons:
        width = entry.get("width")
        if not entry["comparable"] or not width_ok(width):
            continue
        width = float(width)
        if best_width is None or width < best_width:
            best_width = width
            best_round = int(entry["round"])
            without_improvement = 0
        else:
            without_improvement += 1
    completed = max((int(entry["round"]) for entry in journal), default=0)
    reason = None
    if completed >= max_rounds:
        reason = "round_limit"
    elif best_round is not None and without_improvement >= patience:
        reason = "no_improvement"
    return {
        "finished": reason is not None, "reason": reason,
        "completed": completed, "best_round": best_round, "best_width": best_width,
        "without_improvement": without_improvement,
        "comparisons": comparisons,
        "review_required": any(not row["comparable"] for row in comparisons),
    }


#_____________________________________________________________________________
def pick_round(journal: list[dict[str, Any]], adopt: str) -> int:
    """Return 1-based round index. adopt is 'auto' or an integer string."""
    if not journal:
        raise ValueError("empty hdphc journal")
    if adopt != "auto":
        k = int(adopt)
        rounds = {int(e["round"]) for e in journal}
        if k not in rounds:
            raise ValueError(f"adopt round {k} not in journal {sorted(rounds)}")
        return k
    candidates = [row for row in compare_candidates(journal) if row["comparable"]]
    if not candidates:
        raise ValueError("PHC auto adopt: no comparable candidate; review summary or restart/re-evaluate the series")
    return min(candidates, key=lambda row: row["width"])["round"]


#_____________________________________________________________________________
def width_ok(value) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


#_____________________________________________________________________________
def restore_round(run: int, round_i: int, suffix: str) -> Path:
    """Copy round snapshot back to the live HDPHC param path."""
    src = round_dir(run, round_i) / f"HodoPHCParam_run{int(run):05d}_{suffix}"
    if not src.is_file():
        raise FileNotFoundError(f"round snapshot missing: {src}")
    dest = hdphc_param_path(run, suffix)
    dest.parent.mkdir(parents=True, exist_ok=True)
    selected = next((entry for entry in load_journal(run)
                     if int(entry["round"]) == int(round_i)), None)
    if not selected or qa_session.file_sha256(src) != selected.get("param_sha"):
        raise ValueError(f"round {round_i} parameter snapshot does not match its journal")
    from .param import exclusive_param_lock
    with exclusive_param_lock(dest):
        descriptor, temporary = tempfile.mkstemp(prefix=dest.name + ".", dir=dest.parent)
        os.close(descriptor)
        try:
            shutil.copyfile(src, temporary)
            shutil.copymode(src, temporary)
            os.replace(temporary, dest)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)
    from .calib_state_io import write_json
    write_json(rounds_dir(run) / "adoption_validation.json",
               {"round": int(round_i), "matches": None, "status": "pending_decode", "restored_at": _utc_now()})
    return dest


#_____________________________________________________________________________
def score_differences(expected, actual):
    """Compare the same saved PHC candidate after a fresh decode; no quality cuts."""
    differences = []
    if not expected.get("ok") or not actual.get("ok"):
        differences.append("score is not valid")
    def compare(label, left, right):
        finite_values = width_ok(left) and width_ok(right)
        same_value = finite_values and math.isclose(float(left), float(right), rel_tol=1e-8, abs_tol=1e-10)
        if not same_value:
            differences.append(f"{label}: recorded={left}, decoded={right}")
    if expected.get("scoring_context") is not None and expected.get("scoring_context") != actual.get("scoring_context"):
        differences.append("scoring context changed")
    for key in ("width", "n_entries"):
        compare(key, expected.get(key), actual.get(key))
    for detector, values in (expected.get("per_det") or {}).items():
        current = (actual.get("per_det") or {}).get(detector) or {}
        empty_detector = (values.get("n_seg") == current.get("n_seg") == 0
                          and values.get("n_entries") == current.get("n_entries") == 0
                          and not width_ok(values.get("width")) and not width_ok(current.get("width")))
        for key in ("width", "n_entries", "n_seg"):
            if key != "width" or not empty_detector:
                compare(f"{detector}.{key}", values.get(key), current.get(key))
        if "segments" in values:
            recorded_segments = {segment["histogram"]: segment for segment in values["segments"]}
            current_segments = {segment["histogram"]: segment for segment in current.get("segments", [])}
            if recorded_segments.keys() != current_segments.keys():
                differences.append(f"{detector}: histogram set changed")
            for name, segment in recorded_segments.items():
                decoded_segment = current_segments.get(name) or {}
                for key in ("n_entries", "n_total"):
                    compare(f"{name}.{key}", segment.get(key), decoded_segment.get(key))
                for key in ("n_tof_flow", "n_de_flow", "n_both_flow"):
                    if key in segment:
                        compare(f"{name}.{key}", segment[key], decoded_segment.get(key))
                if "axes_sha" in segment and segment["axes_sha"] != decoded_segment.get("axes_sha"):
                    differences.append(f"{name}: histogram axes changed")
                empty_segment = (segment.get("n_entries") == decoded_segment.get("n_entries") == 0
                                 and segment.get("width") is None and decoded_segment.get("width") is None)
                if not empty_segment:
                    compare(f"{name}.width", segment.get("width"), decoded_segment.get("width"))
    return differences


#_____________________________________________________________________________
def write_summary(run: int, adopted: Optional[int] = None) -> Path:
    rows = load_journal(run)
    lines = [
        f"hdphc rounds run{int(run):05d}",
        f"{'round':>5}  {'width':>10}  {'n_entries':>12}  {'ok':>3}  param_sha",
    ]
    comparisons = {row["round"]: row for row in compare_candidates(rows)}
    for entry in rows:
        mark = " *" if adopted is not None and int(entry["round"]) == int(adopted) else ""
        width = entry.get("width")
        wtxt = f"{float(width):.4f}" if width_ok(width) else "nan"
        lines.append(
            f"{int(entry['round']):5d}  {wtxt:>10}  {float(entry.get('n_entries') or 0):12.0f}  "
            f"{'Y' if entry.get('ok') else 'N':>3}  {entry.get('param_sha')}{mark}"
        )
        comparison = comparisons[int(entry["round"])]
        label = "comparable" if comparison["comparable"] else "REVIEW"
        lines.append(f"  fixed population: {label}, width={comparison['width']}, "
                     f"reference=round {comparison['reference_round']}, segments={comparison['required_segments']}")
        lines.append(f"  missing={comparison['missing_segments']}, "
                     f"statistics decreased={comparison['statistics_decreased']}, "
                     f"flow increased={comparison['flow_increased']}, "
                     f"reference minimum entries={comparison['reference_min_entries']}")
        lines.extend("    " + reason for reason in comparison["reasons"])
        for segment in comparison["segments"]:
            lines.append(f"    {segment['histogram']}: entries={segment['n_entries']} "
                         f"(reference {segment['reference_entries']}), total={segment['n_total']} "
                         f"(reference {segment['reference_total']}), "
                         f"flow={segment.get('flow_fraction')} "
                         f"(reference {segment.get('reference_flow_fraction')}), "
                         f"TOF-only flow={segment.get('n_tof_flow')} "
                         f"(reference {segment.get('reference_n_tof_flow')}), "
                         f"dE-only flow={segment.get('n_de_flow')} "
                         f"(reference {segment.get('reference_n_de_flow')})")
    comparison_path = rounds_dir(run) / "comparison.json"
    comparison_path.write_text(json.dumps(list(comparisons.values()), indent=2) + "\n", encoding="utf-8")
    if adopted is not None:
        lines.append(f"adopted: round {int(adopted)}")
    path = summary_path(run)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
