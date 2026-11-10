"""Final QA approval subject / journal for calib_dispatch (#25).

Process complete → awaiting_final_approval (not auto-approved).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional

from . import config, qa_session

STATUS_AWAITING = "awaiting_final_approval"
STATUS_APPROVED = "approved"
# Legacy process-complete status; treat as awaiting.
STATUS_FINISHED_LEGACY = "finished"
COMPLETE_STATUSES = frozenset({
    STATUS_AWAITING, STATUS_APPROVED, STATUS_FINISHED_LEGACY,
})


#_____________________________________________________________________________
def approval_journal_path(run: int) -> Path:
    return Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}" / "approval_journal.jsonl"


#_____________________________________________________________________________
def append_approval_journal(run: int, event: dict[str, Any]) -> None:
    path = approval_journal_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    row = dict(event)
    if "at" not in row:
        from datetime import datetime, timezone
        row["at"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(row, sort_keys=True) + "\n")


#_____________________________________________________________________________
def effective_status(st: Optional[dict[str, Any]]) -> Optional[str]:
    if not st:
        return None
    status = st.get("status")
    if status == STATUS_FINISHED_LEGACY:
        return STATUS_AWAITING
    return status if isinstance(status, str) else None


#_____________________________________________________________________________
def is_process_complete(st: Optional[dict[str, Any]]) -> bool:
    """True when the step chain has run to the end (awaiting or approved)."""
    if not st:
        return False
    if st.get("finished") is True:
        return True
    return effective_status(st) in (STATUS_AWAITING, STATUS_APPROVED, STATUS_FINISHED_LEGACY)


#_____________________________________________________________________________
def latest_step_rows(steps_hist: list) -> list[dict[str, Any]]:
    by_index: dict[int, dict[str, Any]] = {}
    for entry in steps_hist or []:
        if not isinstance(entry, dict):
            continue
        idx = entry.get("index")
        if not isinstance(idx, int):
            continue
        by_index[idx] = {
            "index": idx,
            "step": entry.get("step"),
            "rc": entry.get("rc"),
        }
    return [by_index[i] for i in sorted(by_index)]


#_____________________________________________________________________________
def _param_sha(run: int, category: str, name_prefix: str, suffix: str) -> Optional[str]:
    path = (
        Path(config.PARAM_DIR) / category / config.SUB_DIR
        / f"{name_prefix}_run{int(run):05d}_{suffix}"
    )
    if not path.is_file():
        return None
    return qa_session.file_sha256(path)


#_____________________________________________________________________________
def build_subject(run: int, st: dict[str, Any]) -> dict[str, Any]:
    """Snapshot of the version under final review."""
    suffix = "K" if st.get("calib_particle") == "kaon" or st.get("kaon") is True else "Pi"
    # Prefer session calib_particle on state.
    particle = st.get("calib_particle")
    if particle == "kaon":
        suffix = "K"
    elif particle == "pion":
        suffix = "Pi"

    summary = (
        Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}"
        / "hdphc_rounds" / "summary.txt"
    )
    hdphc_summary_sha = None
    if summary.is_file():
        hdphc_summary_sha = qa_session.file_sha256(summary)

    from . import channel_flags
    open_flags = channel_flags.open_needs_review(run)
    subject = {
        "plan": st.get("plan"),
        "steps": latest_step_rows(st.get("steps") or []),
        "hold": st.get("hold") if isinstance(st.get("hold"), dict) else None,
        "unverified": list(st.get("unverified") or []),
        "open_channel_flags": [
            {
                "mode": r.get("mode"),
                "det": r.get("det"),
                "channel": r.get("channel"),
            }
            for r in open_flags
        ],
        "hdphc_summary_sha": hdphc_summary_sha,
        "hdphc_param_sha": _param_sha(run, "HDPHC", "HodoPHCParam", suffix),
        "hdprm_param_sha": _param_sha(run, "HDPRM", "HodoParam", suffix),
        "suffix": suffix,
    }
    return subject


#_____________________________________________________________________________
def subject_sha(subject: dict[str, Any]) -> str:
    blob = json.dumps(subject, sort_keys=True, separators=(",", ":"))
    return qa_session.sha256_hex(blob, n=16)


#_____________________________________________________________________________
def approval_summary(st: Optional[dict[str, Any]]) -> str:
    if not st:
        return "approval=none"
    status = effective_status(st)
    bits = [f"status={status}"]
    approval = st.get("approval") if isinstance(st.get("approval"), dict) else None
    if approval:
        sha = approval.get("subject_sha") or ""
        bits.append(f"approval.sha={sha[:8]}… via={approval.get('via')!r}")
    elif status == STATUS_AWAITING:
        bits.append("approval=pending")
    else:
        bits.append("approval=none")
    return " ".join(bits)


#_____________________________________________________________________________
def invalidate_approval(run: int, reason: str) -> bool:
    """If process-complete and approved, drop to awaiting. Return True if changed."""
    from datetime import datetime, timezone

    path = (
        Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}"
        / "calib_dispatch_state.json"
    )
    if not path.is_file():
        return False
    try:
        st = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    if not isinstance(st, dict) or not is_process_complete(st):
        return False
    if effective_status(st) != STATUS_APPROVED and not st.get("approval"):
        # Already awaiting (or legacy finished without approval).
        if effective_status(st) == STATUS_FINISHED_LEGACY:
            st["status"] = STATUS_AWAITING
            st["updated_at"] = datetime.now(timezone.utc).strftime(
                "%Y-%m-%dT%H:%M:%SZ"
            )
            path.write_text(
                json.dumps(st, indent=2, sort_keys=True) + "\n", encoding="utf-8",
            )
            return True
        return False
    st["status"] = STATUS_AWAITING
    st["finished"] = True
    st.pop("approval", None)
    st["updated_at"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    path.write_text(json.dumps(st, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    append_approval_journal(run, {
        "event": "invalidate",
        "reason": reason,
    })
    return True


#_____________________________________________________________________________
def awaiting_action_blocks(run: int, sha: str) -> list[dict[str, Any]]:
    """Block Kit Approve/Revoke buttons for the awaiting note."""
    short = (sha or "")[:8]
    run_s = str(int(run))
    return [
        {
            "type": "section",
            "text": {
                "type": "mrkdwn",
                "text": (
                    f"*run{int(run):05d}* awaiting final QA approval "
                    f"(subject_sha=`{short}…`)"
                ),
            },
        },
        {
            "type": "actions",
            "block_id": f"calib_awaiting_{run_s}",
            "elements": [
                {
                    "type": "button",
                    "action_id": "calib_approve_run",
                    "text": {"type": "plain_text", "text": "Approve"},
                    "style": "primary",
                    "value": run_s,
                },
                {
                    "type": "button",
                    "action_id": "calib_revoke_run",
                    "text": {"type": "plain_text", "text": "Revoke"},
                    "value": run_s,
                },
            ],
        },
    ]


#_____________________________________________________________________________
def notify_awaiting_slack(run: int, sha: str) -> tuple[bool, str]:
    """Short thread note that final approval is needed (with Block Kit buttons)."""
    from . import slack_notify

    from . import final_qa_index

    short = (sha or "")[:8]
    n_pdf = final_qa_index.index_summary(int(run)).get("n_pdf", 0)
    text = (
        f"run{int(run):05d}: awaiting final QA approval "
        f"(subject_sha={short}…, QA PDFs={n_pdf}). "
        f"/calib approve {int(run)}  "
        f"(Web: final_approval_web port 8766)"
    )
    return slack_notify.try_notify_failure(
        int(run),
        f"awaiting_final_approval:{short}",
        text,
        blocks=awaiting_action_blocks(int(run), sha),
    )


#_____________________________________________________________________________
def dispatch_state_path(run: int) -> Path:
    return (
        Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}"
        / "calib_dispatch_state.json"
    )


#_____________________________________________________________________________
def write_blocked(run: int, mode: str | None = None) -> tuple[bool, str]:
    """True if param writes must be refused due to hold or open channel flags.

    Blocks when status=held or hold_requested. Optional mode also blocks when
    channel_flags has open needs_review for that mode (or any mode if mode is
    None). Missing/unreadable state does not block. awaiting_final_approval /
    approved do not block by themselves.
    """
    path = dispatch_state_path(run)
    if path.is_file():
        try:
            st = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            st = None
        if isinstance(st, dict):
            if st.get("hold_requested") is True:
                reason = st.get("hold_reason") or "hold_requested"
                return True, f"hold_requested ({reason})"
            if st.get("status") == "held":
                hold = st.get("hold") if isinstance(st.get("hold"), dict) else {}
                reason = hold.get("reason") or "held"
                return True, f"status=held ({reason})"
    from . import channel_flags
    return channel_flags.write_blocked_by_flags(run, mode=mode)
