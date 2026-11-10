"""Thin wrappers so Slack (or other UIs) can drive calib_dispatch safely.

No arbitrary shell: only status / stop / resume / approve / revoke against
a run's session YAML.
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from . import config

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DISPATCH_SCRIPT = PROJECT_ROOT / "calibration" / "scripts" / "calib_dispatch.py"
OPS = ("status", "stop", "resume", "approve", "revoke")
DEDUP_SECONDS = 5


#_____________________________________________________________________________
def session_yml(run: int) -> Path:
    return Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}" / "calib_session.yml"


#_____________________________________________________________________________
def ops_state_path(run: int) -> Path:
    return Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}" / "slack_bot_ops.json"


#_____________________________________________________________________________
def resume_log_path(run: int) -> Path:
    return Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}" / "slack_bot_resume.log"


#_____________________________________________________________________________
def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


#_____________________________________________________________________________
def parse_calib_text(text: str) -> tuple[Optional[str], Optional[int], str]:
    """Parse 'status 3774' / 'approve 12'. Returns (op, run, error)."""
    parts = (text or "").strip().split()
    if len(parts) != 2:
        return None, None, "usage: /calib status|stop|resume|approve|revoke <run>"
    op = parts[0].lower()
    if op not in OPS:
        return None, None, (
            f"unsupported op {parts[0]!r}; "
            f"use status|stop|resume|approve|revoke"
        )
    try:
        run = int(parts[1])
    except ValueError:
        return None, None, f"run must be an integer, got {parts[1]!r}"
    if run < 1:
        return None, None, "run must be >= 1"
    return op, run, ""


#_____________________________________________________________________________
def user_channel_allowed(
    bot_cfg: dict[str, Any], user_id: str, channel_id: str,
) -> bool:
    users = bot_cfg.get("allowed_user_ids") or []
    channels = bot_cfg.get("allowed_channel_ids") or []
    return user_id in users and channel_id in channels


#_____________________________________________________________________________
def check_dedup(
    run: int, op: str, user_id: str, *, window_s: int = DEDUP_SECONDS,
) -> tuple[bool, str]:
    """Return (is_duplicate, message)."""
    path = ops_state_path(run)
    if not path.is_file():
        return False, ""
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False, ""
    if not isinstance(data, dict):
        return False, ""
    if data.get("op") != op:
        return False, ""
    try:
        last = float(data.get("epoch") or 0)
    except (TypeError, ValueError):
        return False, ""
    if time.time() - last > window_s:
        return False, ""
    return True, f"dedup: {op} run{run:05d} within {window_s}s (ignored)"


#_____________________________________________________________________________
def record_op(run: int, op: str, user_id: str) -> None:
    path = ops_state_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "op": op,
        "user_id": user_id,
        "epoch": time.time(),
        "at": _utc_now(),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


#_____________________________________________________________________________
def run_dispatch(
    command: str,
    run: int,
    *,
    background_resume: bool = True,
    timeout_s: float = 30.0,
    user_id: Optional[str] = None,
) -> tuple[bool, str]:
    """Run calib_dispatch status|stop|resume|approve|revoke for session YAML."""
    if command not in OPS:
        return False, f"unsupported command {command!r}"
    yml = session_yml(run)
    if not yml.is_file():
        return False, (
            f"no session YAML for run{run:05d}: {yml.name} "
            f"(run calib_session.py init first)"
        )
    if not DISPATCH_SCRIPT.is_file():
        return False, f"missing dispatcher: {DISPATCH_SCRIPT}"
    cmd = [sys.executable, str(DISPATCH_SCRIPT), command, str(yml)]
    if command == "approve":
        cmd += ["--via", "slack"]
        if user_id:
            cmd += ["--user-id", str(user_id)]

    if command == "resume" and background_resume:
        log_path = resume_log_path(run)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            log_file = log_path.open("a", encoding="utf-8")
        except OSError as exc:
            return False, f"cannot open resume log: {type(exc).__name__}"
        try:
            log_file.write(f"\n--- resume {_utc_now()} ---\n")
            log_file.flush()
            try:
                proc = subprocess.Popen(
                    cmd,
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
            except OSError as exc:
                return False, f"resume spawn failed: {type(exc).__name__}"
        finally:
            # Child keeps its own fd; close the parent handle.
            log_file.close()
        return True, (
            f"resume started run{run:05d} pid={proc.pid}; "
            f"follow {log_path.name} / calib_dispatch_state.json"
        )

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout_s,
        )
    except subprocess.TimeoutExpired:
        return False, f"{command} timed out after {timeout_s:.0f}s"
    except OSError as exc:
        return False, f"{command} failed: {type(exc).__name__}"
    text = ((result.stdout or "") + (result.stderr or "")).strip()
    if command == "status":
        summary = _status_summary_line(run)
        if summary:
            text = (summary + "\n" + text) if text else summary
    if len(text) > 1800:
        text = "…\n" + text[-1800:]
    if not text:
        text = f"{command} exit={result.returncode}"
    return result.returncode == 0, text


#_____________________________________________________________________________
def _status_summary_line(run: int) -> str:
    """Prefix Slack status with hold + approval summary from dispatch state."""
    from . import calib_approval

    path = Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}" / "calib_dispatch_state.json"
    if not path.is_file():
        return ""
    try:
        st = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return ""
    if not isinstance(st, dict):
        return ""
    hold = st.get("hold") if isinstance(st.get("hold"), dict) else None
    bits = [
        f"status={calib_approval.effective_status(st) or st.get('status')}",
        f"hold_requested={bool(st.get('hold_requested'))}",
    ]
    if hold:
        bits.append(
            f"hold.source={hold.get('source')} "
            f"decision={hold.get('decision')!r} "
            f"reason={hold.get('reason')!r}"
        )
    else:
        bits.append("hold=none")
    bits.append(calib_approval.approval_summary(st))
    return " ".join(bits)
