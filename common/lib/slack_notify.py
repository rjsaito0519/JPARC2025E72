"""Post calibration QA PDFs to a per-run Slack thread.

Reads calibration/config/slack_local.yml (gitignored). No env vars, no slack_sdk.
Uses chat.postMessage and the current file upload API (getUploadURLExternal).
"""
from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Optional

from . import config, hodo_qa_pdf
from .local_notify_config import load_slack_config

API_BASE = "https://slack.com/api"
STATE_NAME = "slack_thread.json"


#_____________________________________________________________________________
def _redact(text: str) -> str:
    """Strip token-like / webhook-like substrings from error text."""
    out = text
    for needle in ("xoxb-", "xoxp-", "xoxa-", "discord.com/api/webhooks"):
        if needle in out:
            out = "request failed (secret redacted)"
            break
    return out


#_____________________________________________________________________________
def state_path(run: int) -> Path:
    return Path(config.OUTPUT_DIR) / "calib" / f"run{run:05d}" / STATE_NAME


#_____________________________________________________________________________
def load_state(run: int) -> dict[str, Any]:
    path = state_path(run)
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


#_____________________________________________________________________________
def save_state(run: int, data: dict[str, Any]) -> None:
    path = state_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


#_____________________________________________________________________________
def _api(
    method: str,
    token: str,
    payload: dict[str, Any],
    *,
    as_json: bool = True,
) -> tuple[bool, dict[str, Any], str]:
    """Call Slack Web API. Returns (ok, body, message). Never echoes the token."""
    url = f"{API_BASE}/{method}"
    headers = {"Authorization": f"Bearer {token}"}
    if as_json:
        body = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json; charset=utf-8"
    else:
        body = urllib.parse.urlencode(payload).encode("utf-8")
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    req = urllib.request.Request(url, data=body, headers=headers, method="POST")
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            raw = resp.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        detail = ""
        try:
            detail = exc.read().decode("utf-8", errors="replace")
        except Exception:  # noqa: BLE001
            pass
        return False, {}, _redact(f"HTTP {exc.code}: {detail or exc.reason}")
    except urllib.error.URLError as exc:
        return False, {}, _redact(f"URL error: {exc.reason}")
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError:
        return False, {}, "invalid JSON from Slack"
    if not isinstance(parsed, dict):
        return False, {}, "unexpected Slack response"
    if not parsed.get("ok"):
        err = str(parsed.get("error") or "unknown_error")
        return False, parsed, f"Slack {method}: {err}"
    return True, parsed, "ok"


#_____________________________________________________________________________
def ensure_run_thread(
    run: int,
    text: str,
    *,
    config_path: Path | None = None,
    dry_run: bool = False,
) -> tuple[bool, str, Optional[str]]:
    """Ensure a parent message for the run; return (ok, message, thread_ts)."""
    cfg, status = load_slack_config(config_path)
    if status != "ok" or cfg is None:
        return False, f"slack config: {status}", None
    if dry_run:
        return True, f"dry-run: ensure thread run{run:05d}", "dry-run-ts"

    state = load_state(run)
    existing = state.get("thread_ts")
    if isinstance(existing, str) and existing:
        return True, f"reuse thread_ts={existing}", existing

    ok, body, msg = _api(
        "chat.postMessage",
        cfg["token"],
        {"channel": cfg["channel_id"], "text": text},
    )
    if not ok:
        return False, msg, None
    ts = body.get("ts")
    if not isinstance(ts, str) or not ts:
        return False, "chat.postMessage: no ts", None
    state["thread_ts"] = ts
    state["channel_id"] = cfg["channel_id"]
    state.setdefault("posts", {})
    save_state(run, state)
    return True, f"created thread_ts={ts}", ts


#_____________________________________________________________________________
def _upload_pdf(
    token: str,
    channel_id: str,
    thread_ts: str,
    pdf_path: Path,
    caption: str,
) -> tuple[bool, str]:
    pdf_path = Path(pdf_path)
    if not pdf_path.is_file():
        return False, f"PDF not found: {pdf_path.name}"
    data = pdf_path.read_bytes()
    ok, body, msg = _api(
        "files.getUploadURLExternal",
        token,
        {"filename": pdf_path.name, "length": len(data)},
        as_json=False,
    )
    if not ok:
        return False, msg
    upload_url = body.get("upload_url")
    file_id = body.get("file_id")
    if not isinstance(upload_url, str) or not isinstance(file_id, str):
        return False, "getUploadURLExternal: missing upload_url/file_id"

    put_req = urllib.request.Request(
        upload_url,
        data=data,
        method="PUT",
        headers={"Content-Type": "application/pdf"},
    )
    try:
        with urllib.request.urlopen(put_req, timeout=120) as resp:
            if resp.status and resp.status >= 400:
                return False, f"upload PUT HTTP {resp.status}"
    except urllib.error.HTTPError as exc:
        return False, _redact(f"upload PUT HTTP {exc.code}")
    except urllib.error.URLError as exc:
        return False, _redact(f"upload PUT: {exc.reason}")

    complete: dict[str, Any] = {
        "files": [{"id": file_id, "title": pdf_path.name}],
        "channel_id": channel_id,
        "thread_ts": thread_ts,
    }
    if caption:
        complete["initial_comment"] = caption
    ok, _, msg = _api("files.completeUploadExternal", token, complete)
    if not ok:
        return False, msg
    return True, f"uploaded {pdf_path.name}"


#_____________________________________________________________________________
def post_key(mode: str, particle_suffix: str, pdf_name: str) -> str:
    return f"{mode}|{particle_suffix}|{pdf_name}"


#_____________________________________________________________________________
def post_pdf(
    run: int,
    mode: str,
    pdf_path: Path,
    caption: str,
    *,
    particle_suffix: str = "Pi",
    config_path: Path | None = None,
    dry_run: bool = False,
) -> tuple[bool, str]:
    """Upload a PDF into the run's Slack thread. Records the post in state."""
    cfg, status = load_slack_config(config_path)
    if status != "ok" or cfg is None:
        return False, f"slack config: {status}"

    pdf_path = Path(pdf_path)
    key = post_key(mode, particle_suffix, pdf_path.name)
    parent_text = f"calib QA run{run:05d}"

    if dry_run:
        return True, f"dry-run: post {pdf_path.name} key={key}"

    ok, msg, thread_ts = ensure_run_thread(
        run, parent_text, config_path=config_path, dry_run=False,
    )
    if not ok or not thread_ts:
        return False, msg

    ok, msg = _upload_pdf(
        cfg["token"], cfg["channel_id"], thread_ts, pdf_path, caption,
    )
    if not ok:
        return False, msg

    state = load_state(run)
    posts = state.setdefault("posts", {})
    if not isinstance(posts, dict):
        posts = {}
        state["posts"] = posts
    posts[key] = {"pdf": pdf_path.name, "mode": mode, "particle": particle_suffix}
    state["thread_ts"] = thread_ts
    state["channel_id"] = cfg["channel_id"]
    save_state(run, state)
    return True, msg


#_____________________________________________________________________________
def post_hodo_qa(
    run: int,
    mode: str,
    *,
    kaon: bool = False,
    ftof: bool = False,
    pdf_path: Path | None = None,
    config_path: Path | None = None,
    dry_run: bool = False,
) -> tuple[bool, str]:
    """Post the merged Hodo stage PDF for (run, mode)."""
    if mode not in hodo_qa_pdf.MERGE_MODES:
        return False, f"unsupported mode {mode!r}"
    suffix = hodo_qa_pdf.particle_suffix(kaon)
    if pdf_path is None:
        pdf_path = hodo_qa_pdf.img_dir(run) / hodo_qa_pdf.output_name(
            run, mode, kaon, ftof,
        )
    mid = f"{mode}_ftof" if ftof and mode in ("hdprm", "hdphc") else mode
    caption = f"run{run:05d} hodo {mid} ({suffix})"
    return post_pdf(
        run,
        mode,
        Path(pdf_path),
        caption,
        particle_suffix=suffix,
        config_path=config_path,
        dry_run=dry_run,
    )


#_____________________________________________________________________________
def try_post_hodo_after_merge(
    run: int,
    mode: str,
    *,
    kaon: bool = False,
    ftof: bool = False,
    pdf_path: Path | None = None,
) -> tuple[bool, str]:
    """Best-effort post used from run_hodo. Skip quietly if config missing/disabled."""
    _, status = load_slack_config()
    if status in ("missing", "disabled"):
        return True, f"slack skip: {status}"
    return post_hodo_qa(
        run, mode, kaon=kaon, ftof=ftof, pdf_path=pdf_path, dry_run=False,
    )
