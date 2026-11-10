"""Cluster-wide single Socket Mode calib Slack bot (shared heartbeat + lock).

cw01–08 share OUTPUT_DIR. Only one bot process may run; ensure_slack_bot
starts it on the calling host via tmux when the heartbeat is missing/stale.
"""
from __future__ import annotations

import atexit
import fcntl
import json
import os
import socket
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from . import config, local_notify_config

TMUX_SESSION = "e72_calib_slack_bot"
HEARTBEAT_TTL_S = 90.0
HEARTBEAT_INTERVAL_S = 30.0
ENSURE_WAIT_S = 12.0
PROJECT_ROOT = Path(__file__).resolve().parents[2]
BOT_SCRIPT = PROJECT_ROOT / "calibration" / "scripts" / "calib_slack_bot.py"

_heartbeat_stop: Optional[threading.Event] = None


#_____________________________________________________________________________
def slack_bot_dir() -> Path:
    return Path(config.OUTPUT_DIR) / "calib" / "slack_bot"


#_____________________________________________________________________________
def heartbeat_path() -> Path:
    return slack_bot_dir() / "heartbeat.json"


#_____________________________________________________________________________
def ensure_lock_path() -> Path:
    return slack_bot_dir() / "ensure.lock"


#_____________________________________________________________________________
def bot_log_path() -> Path:
    return slack_bot_dir() / "bot.log"


#_____________________________________________________________________________
def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


#_____________________________________________________________________________
def _parse_utc(stamp: str) -> Optional[datetime]:
    text = (stamp or "").strip()
    if not text:
        return None
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    try:
        dt = datetime.fromisoformat(text)
    except ValueError:
        return None
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt.astimezone(timezone.utc)


#_____________________________________________________________________________
def _hostname() -> str:
    return socket.gethostname().split(".")[0]


#_____________________________________________________________________________
def read_heartbeat() -> Optional[dict[str, Any]]:
    path = heartbeat_path()
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return data if isinstance(data, dict) else None


#_____________________________________________________________________________
def write_heartbeat(
    *,
    host: Optional[str] = None,
    pid: Optional[int] = None,
    tmux_session: str = TMUX_SESSION,
    started_at: Optional[str] = None,
) -> dict[str, Any]:
    slack_bot_dir().mkdir(parents=True, exist_ok=True)
    prev = read_heartbeat() or {}
    now = _utc_now()
    payload = {
        "host": host or _hostname(),
        "pid": int(pid if pid is not None else os.getpid()),
        "tmux_session": str(tmux_session or TMUX_SESSION),
        "started_at": started_at or prev.get("started_at") or now,
        "updated_at": now,
    }
    path = heartbeat_path()
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tmp.replace(path)
    return payload


#_____________________________________________________________________________
def clear_heartbeat_if_owner() -> None:
    data = read_heartbeat()
    if not data:
        return
    if str(data.get("host") or "") != _hostname():
        return
    try:
        if int(data.get("pid") or -1) != os.getpid():
            return
    except (TypeError, ValueError):
        return
    path = heartbeat_path()
    try:
        path.unlink(missing_ok=True)
    except OSError:
        pass


#_____________________________________________________________________________
def _pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False
    return True


#_____________________________________________________________________________
def _tmux_available() -> bool:
    try:
        result = subprocess.run(
            ["tmux", "-V"],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0


#_____________________________________________________________________________
def _tmux_has_session(session: str = TMUX_SESSION) -> bool:
    try:
        result = subprocess.run(
            ["tmux", "has-session", "-t", session],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0


#_____________________________________________________________________________
def _tmux_kill_session(session: str = TMUX_SESSION) -> None:
    try:
        subprocess.run(
            ["tmux", "kill-session", "-t", session],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        pass


#_____________________________________________________________________________
def heartbeat_is_fresh(
    data: Optional[dict[str, Any]],
    *,
    now: Optional[datetime] = None,
    ttl_s: float = HEARTBEAT_TTL_S,
) -> bool:
    """True if shared heartbeat looks like a live cluster bot."""
    if not isinstance(data, dict):
        return False
    updated = _parse_utc(str(data.get("updated_at") or ""))
    if updated is None:
        return False
    now = now or datetime.now(timezone.utc)
    age = (now - updated).total_seconds()
    if age < 0:
        age = 0.0
    if age > float(ttl_s):
        return False
    host = str(data.get("host") or "")
    try:
        pid = int(data.get("pid") or -1)
    except (TypeError, ValueError):
        return False
    if host == _hostname():
        if not _pid_alive(pid):
            return False
    return True


#_____________________________________________________________________________
def _with_ensure_lock(fn):
    slack_bot_dir().mkdir(parents=True, exist_ok=True)
    with ensure_lock_path().open("a+", encoding="utf-8") as lock_file:
        fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        try:
            return fn()
        finally:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)


#_____________________________________________________________________________
def _spawn_bot_tmux() -> tuple[bool, str]:
    log_path = bot_log_path()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    if _tmux_has_session(TMUX_SESSION):
        _tmux_kill_session(TMUX_SESSION)
    # Login-less pane: cd repo and run bot; append to shared log.
    inner = (
        f"cd {PROJECT_ROOT} && "
        f"{sys.executable} {BOT_SCRIPT} "
        f">> {log_path} 2>&1"
    )
    try:
        result = subprocess.run(
            ["tmux", "new-session", "-d", "-s", TMUX_SESSION, inner],
            capture_output=True,
            text=True,
            timeout=15,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, f"tmux spawn failed: {type(exc).__name__}"
    if result.returncode != 0:
        err = (result.stderr or result.stdout or "").strip()
        return False, f"tmux new-session failed: {err or result.returncode}"
    return True, f"tmux session {TMUX_SESSION} started; log={log_path}"


#_____________________________________________________________________________
def _spawn_bot_popen() -> tuple[bool, str]:
    log_path = bot_log_path()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        log_file = log_path.open("a", encoding="utf-8")
    except OSError as exc:
        return False, f"cannot open bot log: {type(exc).__name__}"
    try:
        log_file.write(f"\n--- ensure popen {_utc_now()} ---\n")
        log_file.flush()
        try:
            proc = subprocess.Popen(
                [sys.executable, str(BOT_SCRIPT)],
                cwd=str(PROJECT_ROOT),
                stdout=log_file,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
        except OSError as exc:
            return False, f"popen spawn failed: {type(exc).__name__}"
    finally:
        log_file.close()
    return True, f"bot started pid={proc.pid}; log={log_path}"


#_____________________________________________________________________________
def _wait_for_local_heartbeat(*, timeout_s: float = ENSURE_WAIT_S) -> bool:
    deadline = time.time() + float(timeout_s)
    host = _hostname()
    while time.time() < deadline:
        data = read_heartbeat()
        if heartbeat_is_fresh(data) and str((data or {}).get("host") or "") == host:
            return True
        time.sleep(0.4)
    return False


#_____________________________________________________________________________
def bot_config_gap_message(status: str) -> str:
    """Human-readable gap list for slack_local.yml Socket Mode setup."""
    data, yaml_status = local_notify_config.load_yaml_config(
        local_notify_config.SLACK_LOCAL,
    )
    gaps: list[str] = []
    if yaml_status != "ok" or data is None:
        return (
            f"slack bot config not ready ({status or yaml_status}); "
            f"fix calibration/config/slack_local.yml "
            f"(enabled:true, bot_enabled:true, token, app_token, "
            f"allowed_user_ids). See wiki Calibration-Automation "
            f"Socket Mode section."
        )
    if data.get("bot_enabled") is not True:
        gaps.append("bot_enabled: true")
    if not str(data.get("token") or "").strip():
        gaps.append("token (xoxb-…)")
    if not str(data.get("app_token") or "").strip():
        gaps.append("app_token (xapp-…, Socket Mode connections:write)")
    if not local_notify_config._id_list(data.get("allowed_user_ids")):
        gaps.append("allowed_user_ids: [U…]")
    if not str(data.get("channel_id") or "").strip() and not local_notify_config._id_list(
        data.get("allowed_channel_ids"),
    ):
        gaps.append("channel_id or allowed_channel_ids")
    if gaps:
        return (
            f"slack bot config incomplete ({status}): need "
            + ", ".join(gaps)
            + " in calibration/config/slack_local.yml"
        )
    return (
        f"slack bot config not ready ({status}); "
        f"check calibration/config/slack_local.yml"
    )


#_____________________________________________________________________________
def ensure_slack_bot(*, required: bool = False) -> tuple[bool, str]:
    """Ensure exactly one cluster Slack bot; start on this host if needed.

    Returns (ok, message).
    When required=False, missing bot config → ok=True skip (PDF notify only).
    When required=True, missing config or spawn failure → ok=False.
    """
    _cfg, status = local_notify_config.load_slack_bot_config()
    if status != "ok":
        msg = bot_config_gap_message(status)
        if required:
            return False, msg
        return True, f"slack bot ensure skipped ({status}); {msg}"
    if not BOT_SCRIPT.is_file():
        return False, f"missing bot script: {BOT_SCRIPT}"

    data = read_heartbeat()
    if heartbeat_is_fresh(data):
        return True, (
            f"slack bot already running on {data.get('host')} "
            f"pid={data.get('pid')}"
        )

    def _locked() -> tuple[bool, str]:
        data2 = read_heartbeat()
        if heartbeat_is_fresh(data2):
            return True, (
                f"slack bot already running on {data2.get('host')} "
                f"pid={data2.get('pid')}"
            )
        if _tmux_available():
            ok_spawn, msg_spawn = _spawn_bot_tmux()
        else:
            ok_spawn, msg_spawn = _spawn_bot_popen()
        if not ok_spawn:
            return False, msg_spawn
        if not _wait_for_local_heartbeat():
            return False, (
                f"bot spawn reported ok but heartbeat not fresh: {msg_spawn}"
            )
        data3 = read_heartbeat() or {}
        return True, (
            f"slack bot ensured on {data3.get('host')} "
            f"pid={data3.get('pid')}; {msg_spawn}"
        )

    return _with_ensure_lock(_locked)


#_____________________________________________________________________________
def start_heartbeat_loop(
    *,
    interval_s: float = HEARTBEAT_INTERVAL_S,
    tmux_session: str = TMUX_SESSION,
) -> None:
    """Write heartbeat now and keep refreshing until process exit."""
    global _heartbeat_stop
    write_heartbeat(tmux_session=tmux_session)
    atexit.register(clear_heartbeat_if_owner)
    if _heartbeat_stop is not None:
        return
    stop = threading.Event()
    _heartbeat_stop = stop

    def _loop() -> None:
        while not stop.wait(float(interval_s)):
            try:
                write_heartbeat(tmux_session=tmux_session)
            except OSError:
                pass

    thread = threading.Thread(
        target=_loop, name="calib-slack-bot-heartbeat", daemon=True,
    )
    thread.start()
