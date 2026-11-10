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
def particle_suffix(st: dict[str, Any]) -> tuple[str, str]:
    """Return (particle_label, param_suffix) from plan/state."""
    plan = st.get("plan") if isinstance(st.get("plan"), dict) else {}
    particle = plan.get("calib_particle") or st.get("calib_particle")
    if particle == "kaon" or st.get("kaon") is True or plan.get("kaon") is True:
        return "kaon", "K"
    if particle == "pion":
        return "pion", "Pi"
    return "pion", "Pi"


#_____________________________________________________________________________
def plan_step_names(st: dict[str, Any]) -> set[str]:
    plan = st.get("plan") if isinstance(st.get("plan"), dict) else {}
    names: set[str] = set()
    for entry in plan.get("steps") or []:
        if isinstance(entry, dict) and entry.get("name"):
            names.add(str(entry["name"]))
    return names


#_____________________________________________________________________________
def plan_families(st: dict[str, Any]) -> set[str]:
    names = plan_step_names(st)
    families: set[str] = set()
    if "hodo" in names:
        families.add("hodo")
    if "dc_resi" in names or "dc_chamber" in names:
        families.add("dc")
    if "d5_align" in names:
        families.add("d5")
    return families


#_____________________________________________________________________________
def _param_entry(
    run: int,
    category: str,
    name_prefix: str,
    suffix: str,
    *,
    ext: str = "",
    key: str | None = None,
) -> dict[str, Any]:
    path = (
        Path(config.PARAM_DIR) / category / config.SUB_DIR
        / f"{name_prefix}_run{int(run):05d}_{suffix}{ext}"
    )
    missing = not path.is_file()
    return {
        "key": key or f"{category}/{path.name}",
        "name": path.name,
        "path": str(path),
        "sha": None if missing else qa_session.file_sha256(path),
        "missing": missing,
    }


#_____________________________________________________________________________
def _decode_entry(path: Path, key: str) -> dict[str, Any]:
    path = Path(path)
    missing = not path.is_file()
    return {
        "key": key,
        "name": path.name,
        "path": str(path),
        "token": None if missing else qa_session.input_token(path),
        "missing": missing,
    }


#_____________________________________________________________________________
def _dc_chambers_from_plan(st: dict[str, Any]) -> set[str]:
    """Chamber tags needed by dc_* steps (BcIn / BcOut)."""
    plan = st.get("plan") if isinstance(st.get("plan"), dict) else {}
    chambers: set[str] = set()
    for entry in plan.get("steps") or []:
        if not isinstance(entry, dict):
            continue
        if entry.get("name") not in ("dc_resi", "dc_chamber"):
            continue
        opts = entry.get("opts") if isinstance(entry.get("opts"), dict) else {}
        chambers.add("BcIn" if opts.get("bcin") else "BcOut")
    return chambers


#_____________________________________________________________________________
def _latest_qa_manifests(
    run: int,
    particle: str,
    families: set[str],
) -> tuple[list[dict[str, Any]], list[str]]:
    """Latest qa_manifest per (family, mode) for this particle; missing families listed."""
    root = qa_session.manifest_dir(int(run))
    by_key: dict[tuple[str, str], dict[str, Any]] = {}
    if root.is_dir():
        for path in sorted(root.glob("*.json")):
            try:
                data = qa_session.load_manifest(path)
            except (OSError, ValueError, TypeError):
                continue
            meta = data.get("meta") if isinstance(data.get("meta"), dict) else {}
            family = str(meta.get("family") or "")
            if family not in families:
                continue
            meta_particle = meta.get("particle")
            if meta_particle is not None and str(meta_particle) != particle:
                continue
            mode = str(meta.get("mode") or data.get("label") or "unknown")
            key = (family, mode)
            prev = by_key.get(key)
            finished = str(data.get("finished_at") or "")
            if prev is None or finished >= str(prev.get("finished_at") or ""):
                by_key[key] = {
                    "family": family,
                    "mode": mode,
                    "exec_id": data.get("exec_id"),
                    "manifest": path.name,
                    "finished_at": data.get("finished_at"),
                    "input_token": data.get("input_token"),
                    "output_shas": [
                        {"name": o.get("name"), "sha": o.get("sha")}
                        for o in (data.get("outputs") or [])
                        if isinstance(o, dict)
                    ],
                    "param_shas": [
                        {
                            "name": p.get("name"),
                            "sha": p.get("sha"),
                            "role": p.get("role"),
                        }
                        for p in (data.get("params") or [])
                        if isinstance(p, dict)
                    ],
                }
    rows = [by_key[k] for k in sorted(by_key)]
    present = {row["family"] for row in rows}
    missing = [f"qa_manifest:{fam}" for fam in sorted(families) if fam not in present]
    return rows, missing


#_____________________________________________________________________________
def build_subject(run: int, st: dict[str, Any]) -> dict[str, Any]:
    """Snapshot of the version under final review (plan-gated params/decode/QA)."""
    particle, suffix = particle_suffix(st)
    families = plan_families(st)
    run_i = int(run)
    run_tag = f"run{run_i:05d}"

    params: list[dict[str, Any]] = []
    decodes: list[dict[str, Any]] = []
    missing: list[str] = []

    if "hodo" in families:
        for entry in (
            _param_entry(run_i, "HDPRM", "HodoParam", suffix, key="HDPRM"),
            _param_entry(run_i, "HDPHC", "HodoPHCParam", suffix, key="HDPHC"),
        ):
            params.append(entry)
            if entry["missing"]:
                missing.append(f"param:{entry['key']}")
        hodo_root = Path(config.DATA_DIR) / f"{run_tag}_Hodo.root"
        dec = _decode_entry(hodo_root, "Hodo")
        decodes.append(dec)
        if dec["missing"]:
            missing.append("decode:Hodo")

    if "dc" in families:
        for entry in (
            _param_entry(run_i, "DCTDC", "DCTdcParam", suffix, key="DCTDC"),
            _param_entry(
                run_i, "DCDRFT", "DCDriftParam", suffix, ext=".root", key="DCDRFT",
            ),
            _param_entry(run_i, "DCGEO", "DCGeomParam", suffix, key="DCGEO"),
        ):
            params.append(entry)
            if entry["missing"]:
                missing.append(f"param:{entry['key']}")
        chambers = _dc_chambers_from_plan(st) or {"BcOut"}
        for chamber in sorted(chambers):
            path = Path(config.DATA_DIR) / f"{run_tag}_{chamber}.root"
            dec = _decode_entry(path, chamber)
            decodes.append(dec)
            if dec["missing"]:
                missing.append(f"decode:{chamber}")

    if "d5" in families:
        geo = _param_entry(run_i, "DCGEO", "DCGeomParam", suffix, key="DCGEO")
        if not any(p.get("key") == "DCGEO" for p in params):
            params.append(geo)
            if geo["missing"]:
                missing.append("param:DCGEO")
        d5_normal = (
            Path(config.DECODE_DIR) / run_tag / f"{run_tag}_D5_{suffix}.root"
        )
        d5_ref = (
            Path(config.SCRATCH_DIR) / run_tag
            / f"{run_tag}_D5_{suffix}_bc1reso.root"
        )
        for key, path in (("D5", d5_normal), ("D5_ref", d5_ref)):
            dec = _decode_entry(path, key)
            decodes.append(dec)
            if dec["missing"]:
                missing.append(f"decode:{key}")

    summary = (
        Path(config.OUTPUT_DIR) / "calib" / run_tag / "hdphc_rounds" / "summary.txt"
    )
    hdphc_summary_sha = (
        qa_session.file_sha256(summary) if summary.is_file() else None
    )
    if "hodo" in families and hdphc_summary_sha is None:
        missing.append("hdphc_summary")

    qa_manifests, qa_missing = _latest_qa_manifests(run_i, particle, families)
    missing.extend(qa_missing)

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
        "families": sorted(families),
        "particle": particle,
        "suffix": suffix,
        "params": params,
        "decodes": decodes,
        "qa_manifests": qa_manifests,
        "missing": missing,
        "hdphc_summary_sha": hdphc_summary_sha,
        # Flat aliases kept for older readers / journal diffs.
        "hdphc_param_sha": _param_sha(run_i, "HDPHC", "HodoPHCParam", suffix),
        "hdprm_param_sha": _param_sha(run_i, "HDPRM", "HodoParam", suffix),
    }
    return subject


#_____________________________________________________________________________
def subject_missing(subject: dict[str, Any]) -> list[str]:
    rows = subject.get("missing") or []
    return [str(x) for x in rows]


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
                    # Pin the viewed subject version: run:subject_sha
                    "value": f"{run_s}:{sha}",
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
def notify_awaiting_slack(
    run: int,
    sha: str,
    st: dict[str, Any] | None = None,
) -> tuple[bool, str]:
    """Short thread note that final approval is needed (with Block Kit buttons)."""
    from . import final_qa_index, slack_notify

    short = (sha or "")[:8]
    particle = None
    families = None
    if isinstance(st, dict):
        particle, _suffix = particle_suffix(st)
        families = sorted(plan_families(st))
    n_pdf = final_qa_index.index_summary(
        int(run), particle=particle, families=families,
    ).get("n_pdf", 0)
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
