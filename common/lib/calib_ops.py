"""Thin wrappers so Slack (or other UIs) can drive calib_dispatch safely.

No arbitrary shell: defined ops only against a run's session YAML (or
extra_hdphc_iter for hdphc-again).
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import yaml

from . import config

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DISPATCH_SCRIPT = PROJECT_ROOT / "calibration" / "scripts" / "calib_dispatch.py"
EXTRA_HDPHC_SCRIPT = PROJECT_ROOT / "calibration" / "scripts" / "extra_hdphc_iter.py"
OPS = (
    "status", "stop", "resume", "approve", "revoke",
    "hold", "decide", "hdphc-again", "flag", "clear",
)
DECIDE_ACTIONS = ("refit", "continue")
DEDUP_SECONDS = 5
USAGE = (
    "usage: /calib status|stop|resume|approve|revoke|hdphc-again <run> | "
    "hold <run> <reason...> | decide <run> refit|continue | "
    "flag <run> <mode> <det> <channel> <reason...> | "
    "clear <run> <mode> <det> <channel>"
)


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
def hdphc_again_log_path(run: int) -> Path:
    return (
        Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}"
        / "slack_bot_hdphc_again.log"
    )


#_____________________________________________________________________________
def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


#_____________________________________________________________________________
def parse_calib_text(
    text: str,
) -> tuple[Optional[str], Optional[int], str, dict[str, Any]]:
    """Parse slash text. Returns (op, run, error, extras).

    extras may include reason (hold) or action (decide).
    """
    parts = (text or "").strip().split()
    if not parts:
        return None, None, USAGE, {}
    op = parts[0].lower()
    if op not in OPS:
        return None, None, f"unsupported op {parts[0]!r}; {USAGE}", {}

    if op == "hold":
        if len(parts) < 3:
            return None, None, "usage: /calib hold <run> <reason...>", {}
        try:
            run = int(parts[1])
        except ValueError:
            return None, None, f"run must be an integer, got {parts[1]!r}", {}
        if run < 1:
            return None, None, "run must be >= 1", {}
        reason = " ".join(parts[2:]).strip()
        if not reason:
            return None, None, "hold requires a non-empty reason", {}
        return op, run, "", {"reason": reason}

    if op == "decide":
        if len(parts) != 3:
            return None, None, "usage: /calib decide <run> refit|continue", {}
        try:
            run = int(parts[1])
        except ValueError:
            return None, None, f"run must be an integer, got {parts[1]!r}", {}
        if run < 1:
            return None, None, "run must be >= 1", {}
        action = parts[2].lower()
        if action not in DECIDE_ACTIONS:
            return None, None, f"decide action must be refit|continue, got {parts[2]!r}", {}
        return op, run, "", {"action": action}

    if op == "flag":
        if len(parts) < 6:
            return None, None, (
                "usage: /calib flag <run> <mode> <det> <channel> <reason...>"
            ), {}
        try:
            run = int(parts[1])
        except ValueError:
            return None, None, f"run must be an integer, got {parts[1]!r}", {}
        if run < 1:
            return None, None, "run must be >= 1", {}
        reason = " ".join(parts[5:]).strip()
        if not reason:
            return None, None, "flag requires a non-empty reason", {}
        return op, run, "", {
            "mode": parts[2],
            "det": parts[3],
            "channel": parts[4],
            "reason": reason,
        }

    if op == "clear":
        if len(parts) != 5:
            return None, None, (
                "usage: /calib clear <run> <mode> <det> <channel>"
            ), {}
        try:
            run = int(parts[1])
        except ValueError:
            return None, None, f"run must be an integer, got {parts[1]!r}", {}
        if run < 1:
            return None, None, "run must be >= 1", {}
        return op, run, "", {
            "mode": parts[2],
            "det": parts[3],
            "channel": parts[4],
        }

    if len(parts) != 2:
        return None, None, USAGE, {}
    try:
        run = int(parts[1])
    except ValueError:
        return None, None, f"run must be an integer, got {parts[1]!r}", {}
    if run < 1:
        return None, None, "run must be >= 1", {}
    return op, run, "", {}


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
def session_wants_kaon(run: int) -> bool:
    path = session_yml(run)
    if not path.is_file():
        return False
    try:
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError):
        return False
    if not isinstance(doc, dict):
        return False
    if doc.get("kaon") is True:
        return True
    particle = str(doc.get("calib_particle") or "").lower()
    return particle == "kaon"


#_____________________________________________________________________________
def _spawn_background(cmd: list[str], log_path: Path, label: str, run: int) -> tuple[bool, str]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        log_file = log_path.open("a", encoding="utf-8")
    except OSError as exc:
        return False, f"cannot open {label} log: {type(exc).__name__}"
    try:
        log_file.write(f"\n--- {label} {_utc_now()} ---\n")
        log_file.flush()
        try:
            proc = subprocess.Popen(
                cmd,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
        except OSError as exc:
            return False, f"{label} spawn failed: {type(exc).__name__}"
    finally:
        log_file.close()
    return True, (
        f"{label} started run{run:05d} pid={proc.pid}; "
        f"follow {log_path.name}"
    )


#_____________________________________________________________________________
def run_dispatch(
    command: str,
    run: int,
    *,
    background_resume: bool = True,
    timeout_s: float = 30.0,
    user_id: Optional[str] = None,
    reason: Optional[str] = None,
    action: Optional[str] = None,
    mode: Optional[str] = None,
    det: Optional[str] = None,
    channel: Optional[str] = None,
) -> tuple[bool, str]:
    """Run a defined calib op for the run's session YAML (or extra hdphc)."""
    if command not in OPS:
        return False, f"unsupported command {command!r}"
    yml = session_yml(run)
    if not yml.is_file():
        return False, (
            f"no session YAML for run{run:05d}: {yml.name} "
            f"(run calib_session.py init first)"
        )

    if command == "hdphc-again":
        if not EXTRA_HDPHC_SCRIPT.is_file():
            return False, f"missing script: {EXTRA_HDPHC_SCRIPT}"
        cmd = [sys.executable, str(EXTRA_HDPHC_SCRIPT), str(run)]
        if session_wants_kaon(run):
            cmd.append("--kaon")
        if background_resume:
            return _spawn_background(
                cmd, hdphc_again_log_path(run), "hdphc-again", run,
            )
        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=timeout_s,
            )
        except subprocess.TimeoutExpired:
            return False, f"hdphc-again timed out after {timeout_s:.0f}s"
        except OSError as exc:
            return False, f"hdphc-again failed: {type(exc).__name__}"
        text = ((result.stdout or "") + (result.stderr or "")).strip()
        if len(text) > 1800:
            text = "…\n" + text[-1800:]
        if not text:
            text = f"hdphc-again exit={result.returncode}"
        return result.returncode == 0, text

    if not DISPATCH_SCRIPT.is_file():
        return False, f"missing dispatcher: {DISPATCH_SCRIPT}"
    dispatch_cmd = {
        "flag": "flag-channel",
        "clear": "clear-channel",
    }.get(command, command)
    cmd = [sys.executable, str(DISPATCH_SCRIPT), dispatch_cmd, str(yml)]
    if command == "approve":
        cmd += ["--via", "slack"]
        if user_id:
            cmd += ["--user-id", str(user_id)]
    if command == "hold":
        if not reason or not str(reason).strip():
            return False, "hold requires a reason"
        cmd += ["--reason", str(reason).strip()]
    if command == "decide":
        if action not in DECIDE_ACTIONS:
            return False, f"decide requires action refit|continue, got {action!r}"
        cmd += ["--action", str(action)]
    if command in ("flag", "clear"):
        if not mode or not det or not channel:
            return False, f"{command} requires mode det channel"
        cmd += ["--mode", str(mode), "--det", str(det), "--channel", str(channel)]
        if command == "flag":
            if not reason or not str(reason).strip():
                return False, "flag requires a reason"
            cmd += ["--reason", str(reason).strip()]
        elif reason:
            cmd += ["--reason", str(reason).strip()]

    if command == "resume" and background_resume:
        return _spawn_background(cmd, resume_log_path(run), "resume", run)

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
