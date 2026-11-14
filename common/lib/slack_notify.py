"""Post calibration QA PDFs to a per-run Slack thread.

Reads calibration/config/slack_local.yml (gitignored). No env vars, no slack_sdk.
Uses chat.postMessage and the current file upload API (getUploadURLExternal).
"""
from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from . import config, hodo_qa_pdf
from .calib_state_io import locked_path, write_json
from .local_notify_config import load_slack_config

API_BASE = "https://slack.com/api"
STATE_NAME = "slack_thread.json"

# Transport / I/O failures that must not escape into calibration callers.
# HTTPError is a URLError (and OSError); keep HTTPError first for status detail.
_NOTIFY_TRANSPORT_ERRORS = (
    TimeoutError,
    ConnectionError,
    BrokenPipeError,
    urllib.error.HTTPError,
    urllib.error.URLError,
    OSError,
)


#_____________________________________________________________________________
def _redact(text: str) -> str:
    """Strip token-like / webhook-like / signed-URL substrings from error text."""
    out = str(text)
    for needle in (
        "xoxb-", "xoxp-", "xoxa-",
        "discord.com/api/webhooks",
        "files.slack.com",
        "upload/v1/",
    ):
        if needle in out:
            return "request failed (secret redacted)"
    return out


#_____________________________________________________________________________
def _exc_message(prefix: str, exc: BaseException) -> str:
    return _redact(f"{prefix}: {type(exc).__name__}")


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
def save_state(run: int, data: dict[str, Any]) -> bool:
    """Write state. Returns False on I/O failure (does not raise)."""
    path = state_path(run)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        write_json(path, data)
        return True
    except OSError:
        return False


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
    except _NOTIFY_TRANSPORT_ERRORS as exc:
        return False, {}, _exc_message(f"Slack {method}", exc)
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
def _utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


#_____________________________________________________________________________
def _archive_active_thread(state: dict[str, Any]) -> None:
    """Move the active channel thread into history; clear posts/failures."""
    old_ts = state.get("thread_ts")
    old_channel = state.get("channel_id")
    if not (isinstance(old_ts, str) and old_ts) and old_channel is None:
        state["posts"] = {}
        state["failures"] = {}
        state.pop("thread_ts", None)
        state.pop("channel_id", None)
        return
    entry = {
        "channel_id": old_channel,
        "thread_ts": old_ts if isinstance(old_ts, str) else None,
        "posts": state.get("posts") if isinstance(state.get("posts"), dict) else {},
        "failures": (
            state.get("failures") if isinstance(state.get("failures"), dict) else {}
        ),
        "closed_at": _utc_now(),
    }
    history = state.get("history")
    if not isinstance(history, list):
        history = []
    history.append(entry)
    state["history"] = history
    state["posts"] = {}
    state["failures"] = {}
    state.pop("thread_ts", None)
    state.pop("channel_id", None)


#_____________________________________________________________________________
def ensure_run_thread(
    run: int,
    text: str,
    *,
    config_path: Path | None = None,
    dry_run: bool = False,
) -> tuple[bool, str, Optional[str]]:
    """Ensure a parent message for the run; return (ok, message, thread_ts)."""
    if dry_run:
        return True, f"dry-run: ensure thread run{run:05d}", "dry-run-ts"
    try:
        with locked_path(state_path(run)):
            cfg, status = load_slack_config(config_path)
            if status != "ok" or cfg is None:
                return False, f"slack config: {status}", None

            channel_id = cfg["channel_id"]
            state = load_state(run)
            existing = state.get("thread_ts")
            if (
                isinstance(existing, str)
                and existing
                and state.get("channel_id") == channel_id
            ):
                return True, f"reuse thread_ts={existing}", existing

            if isinstance(existing, str) and existing:
                _archive_active_thread(state)

            ok, body, msg = _api(
                "chat.postMessage",
                cfg["token"],
                {"channel": channel_id, "text": text},
            )
            if not ok:
                return False, msg, None
            ts = body.get("ts")
            if not isinstance(ts, str) or not ts:
                return False, "chat.postMessage: no ts", None
            state["thread_ts"] = ts
            state["channel_id"] = channel_id
            state.setdefault("posts", {})
            state.setdefault("failures", {})
            if not save_state(run, state):
                return False, "state save failed", None
            return True, f"created thread_ts={ts}", ts
    except _NOTIFY_TRANSPORT_ERRORS as exc:
        return False, _exc_message("ensure_run_thread", exc), None


#_____________________________________________________________________________
def _upload_pdf(
    token: str,
    channel_id: str,
    thread_ts: str,
    pdf_path: Path,
    caption: str,
) -> tuple[bool, str]:
    pdf_path = Path(pdf_path)
    try:
        if not pdf_path.is_file():
            return False, f"PDF not found: {pdf_path.name}"
        data = pdf_path.read_bytes()
    except OSError as exc:
        return False, _exc_message("PDF read", exc)

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

    # Slack expects POST of the raw bytes to upload_url (not PUT).
    upload_req = urllib.request.Request(
        upload_url,
        data=data,
        method="POST",
        headers={"Content-Type": "application/octet-stream"},
    )
    try:
        with urllib.request.urlopen(upload_req, timeout=180) as resp:
            if resp.status and resp.status >= 400:
                return False, f"upload POST HTTP {resp.status}"
    except urllib.error.HTTPError as exc:
        return False, _redact(f"upload POST HTTP {exc.code}")
    except _NOTIFY_TRANSPORT_ERRORS as exc:
        return False, _exc_message("upload POST", exc)

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
    manifest: dict[str, Any] | None = None,
) -> tuple[bool, str]:
    """Upload a PDF into the run's Slack thread. Records the post in state."""
    pdf_path = Path(pdf_path)
    key = post_key(mode, particle_suffix, pdf_path.name)
    parent_text = f"calib QA run{run:05d}"
    if manifest and manifest.get("exec_id"):
        caption = f"{caption}\n結果ID: {str(manifest['exec_id']).split('_')[-1]}"

    caption = f"【QA PDF】{caption}"
    if dry_run:
        return True, f"dry-run: post {pdf_path.name} key={key}"

    try:
        with locked_path(state_path(run)):
            cfg, status = load_slack_config(config_path)
            if status != "ok" or cfg is None:
                return False, f"slack config: {status}"

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
            entry = {"pdf": pdf_path.name, "mode": mode, "particle": particle_suffix}
            if manifest:
                from . import qa_session
                entry.update(qa_session.manifest_post_fields(manifest, pdf_path))
            posts[key] = entry
            state["thread_ts"] = thread_ts
            state["channel_id"] = cfg["channel_id"]
            if not save_state(run, state):
                return False, "state save failed"
            return True, msg
    except _NOTIFY_TRANSPORT_ERRORS as exc:
        return False, _exc_message("post_pdf", exc)


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
    manifest: dict[str, Any] | None = None,
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
    meta = (manifest or {}).get("meta") or {}
    if mode == "hdphc":
        round_i = meta.get("hdphc_round")
        total = meta.get("hdphc_rounds")
        adopted = meta.get("adopted_round")
        denominator = f"/{total}" if total else ""
        if adopted is not None:
            caption += f"\n【採用後 QA】第{adopted}{denominator}ラウンドを採用。"
        elif round_i is not None:
            caption += f"\n【PHC 第{round_i}ラウンド" + (f"・最大{total}回" if total else "") + "】"
            if meta.get("fit_rc") == 0:
                caption += "fit 結果です。この後、更新後 decode と幅評価を行い、完了後に round の承認ボタンを表示します。"
            elif meta.get("fit_rc") == 4:
                caption += "有効な更新結果なし・未適用で停止。round は未完了です。"
            elif meta.get("fit_rc") == 3:
                caption += "要確認です。後続のレンジ確認または確認ボタンから判断してください。"
    return post_pdf(
        run,
        mode,
        Path(pdf_path),
        caption,
        particle_suffix=suffix,
        config_path=config_path,
        dry_run=dry_run,
        manifest=manifest,
    )


#_____________________________________________________________________________
def try_post_hodo_after_merge(
    run: int,
    mode: str,
    *,
    kaon: bool = False,
    ftof: bool = False,
    pdf_path: Path | None = None,
    manifest: dict[str, Any] | None = None,
) -> tuple[bool, str]:
    """Best-effort post used from run_hodo. Skip quietly if config missing/disabled."""
    try:
        _, status = load_slack_config()
        if status in ("missing", "disabled"):
            return True, f"slack skip: {status}"
        return post_hodo_qa(
            run, mode, kaon=kaon, ftof=ftof, pdf_path=pdf_path, dry_run=False,
            manifest=manifest,
        )
    except _NOTIFY_TRANSPORT_ERRORS as exc:
        return False, _exc_message("try_post_hodo_after_merge", exc)


#_____________________________________________________________________________
def status_color(text):
    """Consistent color bars; labels still carry meaning without color."""
    if any(word in text for word in ("エラー", "失敗", "failed", "FAIL", "ERROR")):
        return "#D64545"
    if any(word in text for word in ("確認待ち", "承認待ち", "要確認", "CHECK", "SETTLED", "held", "取り消")):
        return "#E3A008"
    if any(word in text for word in ("停止", "止める", "stopped", "開始待ち", "状態記録なし")):
        return "#7A869A"
    if any(word in text for word in ("承認済み", "自動処理完了", "【選択済み】最終承認", "CONVERGED", "収束")):
        return "#2E9D63"
    return "#2F80ED"


#_____________________________________________________________________________
def colored_notice(text, blocks=None):
    content = review_blocks(blocks) if blocks else [
        {"type": "section", "text": {"type": "mrkdwn", "text": text}},
    ]
    if content and content[0].get("type") == "header":
        heading, content = content[0], content[1:]
    else:
        lines = str((content[0].get("text") or {}).get("text") or text).splitlines()
        title = lines[0].replace("*", "")[:150] if len(lines) > 1 else "E72 Calibration"
        heading = {"type": "header", "text": {"type": "plain_text", "text": title}}
        if len(lines) > 1:
            content = [{"type": "section", "text": {"type": "mrkdwn", "text": "\n".join(lines[1:])}}, *content[1:]]
    return {"blocks": [heading], "attachments": [{"color": status_color(text), "blocks": content}]}


#_____________________________________________________________________________
def review_blocks(blocks):
    """Keep the stage prominent, explanations readable, and actions last."""
    if not blocks or blocks[0].get("type") != "section":
        return blocks
    if not any(block.get("type") == "actions" for block in blocks):
        return blocks
    source = (blocks[0].get("text") or {}).get("text") or ""
    lines = source.splitlines()
    if not lines:
        return blocks
    heading = lines[0].replace("*", "")[:150]
    body, status, references = [], [], []
    for line in lines[1:]:
        if line.startswith("対象PDF:") or any(line.startswith(prefix) for prefix in ("step ", "stale_input:", "qa_delivery_failed:", "write_blocked:", "phc_round_missing:")):
            references.append(line)
        elif line.startswith("【"):
            status.append(line)
        else:
            body.append(line)
    rendered = [{"type": "header", "text": {"type": "plain_text", "text": heading}}]
    detail = "\n\n".join(part for part in ("\n".join(status), "\n".join(body)) if part)
    if detail:
        rendered.append({"type": "section", "text": {"type": "mrkdwn", "text": detail}})
    if references:
        rendered.append({"type": "context", "elements": [{"type": "mrkdwn", "text": "\n".join(references)}]})
    rendered.extend(block for block in blocks[1:] if block.get("type") != "actions")
    rendered.extend(block for block in blocks[1:] if block.get("type") == "actions")
    return rendered


#_____________________________________________________________________________
def notify_failure(
    run: int,
    key: str,
    text: str,
    *,
    config_path: Path | None = None,
    dry_run: bool = False,
    blocks: list | None = None,
) -> tuple[bool, str, bool]:
    """Post a short clear-failure note to the run thread.

    Returns (ok, message, posted). Same key+text already recorded → dedup (posted=False).
    No quality tracking beyond that. Does not change calibration exit codes.
    Optional Block Kit `blocks` are attached to chat.postMessage.
    """
    if dry_run:
        return True, f"dry-run: failure {key}", True

    try:
        with locked_path(state_path(run)):
            cfg, status = load_slack_config(config_path)
            if status in ("missing", "disabled"):
                return True, f"slack skip: {status}", False
            if status != "ok" or cfg is None:
                return False, f"slack config: {status}", False

            parent = f"calib QA run{run:05d}"
            ok, msg, thread_ts = ensure_run_thread(
                run, parent, config_path=config_path, dry_run=False,
            )
            if not ok or not thread_ts:
                return False, msg, False

            # Dedup only against the active channel after ensure_run_thread.
            state = load_state(run)
            failures = state.setdefault("failures", {})
            if not isinstance(failures, dict):
                failures = {}
                state["failures"] = failures
            prev = failures.get(key)
            if isinstance(prev, dict) and prev.get("text") == text:
                return True, f"dedup: {key}", False

            payload = {
                "channel": cfg["channel_id"],
                "thread_ts": thread_ts,
                "text": text,
            }
            payload.update(colored_notice(text, blocks))
            ok, body, msg = _api(
                "chat.postMessage",
                cfg["token"],
                payload,
            )
            if not ok:
                return False, msg, False

            failures[key] = {
                "text": text,
                "ts": body.get("ts") if isinstance(body.get("ts"), str) else None,
                "color_layout": 2,
            }
            state["thread_ts"] = thread_ts
            state["channel_id"] = cfg["channel_id"]
            if not save_state(run, state):
                return False, "state save failed", False
            return True, f"notified failure {key}", True
    except _NOTIFY_TRANSPORT_ERRORS as exc:
        return False, _exc_message("notify_failure", exc), False


#_____________________________________________________________________________
def try_notify_failure(
    run: int, key: str, text: str, *, blocks: list | None = None,
) -> tuple[bool, str]:
    """Best-effort failure note for run_*. Skip if Slack config missing/disabled."""
    try:
        ok, msg, _posted = notify_failure(
            run, key, text, dry_run=False, blocks=blocks,
        )
        return ok, msg
    except _NOTIFY_TRANSPORT_ERRORS as exc:
        return False, _exc_message("try_notify_failure", exc)


#_____________________________________________________________________________
def update_failure_message(
    run: int,
    key: str,
    text: str,
    *,
    blocks: list | None = None,
    config_path: Path | None = None,
) -> tuple[bool, str]:
    """chat.update an earlier failure note (same key) when its ts is known."""
    try:
        with locked_path(state_path(run)):
            cfg, status = load_slack_config(config_path)
            if status in ("missing", "disabled"):
                return True, f"slack skip: {status}"
            if status != "ok" or cfg is None:
                return False, f"slack config: {status}"

            state = load_state(run)
            failures = state.get("failures")
            if not isinstance(failures, dict):
                return False, f"no failure note for key {key}"
            prev = failures.get(key)
            if not isinstance(prev, dict):
                return False, f"no failure note for key {key}"
            ts = prev.get("ts")
            if not isinstance(ts, str) or not ts:
                return False, f"no message ts for key {key}"

            payload: dict[str, Any] = {
                "channel": cfg["channel_id"],
                "ts": ts,
                "text": text,
            }
            payload.update(colored_notice(text, blocks))
            ok, _body, msg = _api("chat.update", cfg["token"], payload)
            if not ok:
                return False, msg
            failures[key] = {"text": text, "ts": ts, "color_layout": 2}
            state["failures"] = failures
            if not save_state(run, state):
                return False, "state save failed"
            return True, f"updated failure {key}"
    except _NOTIFY_TRANSPORT_ERRORS as exc:
        return False, _exc_message("update_failure_message", exc)


#_____________________________________________________________________________
def try_update_failure_message(
    run: int, key: str, text: str, *, blocks: list | None = None,
) -> tuple[bool, str]:
    try:
        return update_failure_message(
            int(run), key, text, blocks=blocks,
        )
    except _NOTIFY_TRANSPORT_ERRORS as exc:
        return False, _exc_message("try_update_failure_message", exc)
