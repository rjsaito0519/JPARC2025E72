"""Maintain a read-only calibration summary in each existing Slack parent message."""
from __future__ import annotations

import atexit
import hashlib
import json
import re
import threading
from datetime import datetime, timedelta, timezone
from pathlib import Path

from . import calib_approval, config, slack_notify
from .calib_state_io import locked_path

POLL_SECONDS = 30.0
LANES = ("hodo", "bcin", "bcout", "d5")
LABELS = {"hodo": "Hodo", "bcin": "BcIn / BLC1", "bcout": "BcOut / BLC2", "d5": "D5"}
_stop = None


#_____________________________________________________________________________
def read_object(path):
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return value if isinstance(value, dict) else {}


#_____________________________________________________________________________
def phase_of(state):
    status = calib_approval.effective_status(state)
    hold = state.get("hold") or {}
    if status == "approved":
        return "approved", "✅ 承認済み"
    if status in ("awaiting_final_approval", "done"):
        return "awaiting", "☑ 処理完了・最終承認待ち"
    if status == "held":
        if hold.get("decision") == "stop":
            return "stopped", "⏸ 停止中（ユーザー操作）"
        if hold.get("decision") in ("continue", "refit"):
            return "running", "▶ 再開要求を受け付け済み"
        if state.get("rc") == 3 or hold.get("source") == "exit_check":
            return "check", "⚠ 要確認・停止中"
        return "waiting", "🟡 確認待ち"
    if status == "running":
        if state.get("stop_requested"):
            return "running", "▶ 停止要求中（現在の処理終了待ち）"
        return "running", "▶ 実行中"
    if status in ("stopped", "failed"):
        if state.get("rc") == 1:
            return "stopped", "❌ エラーで停止"
        return "stopped", "⏸ 停止中"
    return "pending", "開始待ち" if status else "状態記録なし"


#_____________________________________________________________________________
def overall_phase(phases):
    if not phases:
        return "状態記録なし"
    if all(p == "approved" for p in phases):
        return "✅ 承認済み"
    if "check" in phases:
        return "⚠ 要確認あり" + ("・ほかのレーンは実行中" if "running" in phases else "")
    if "waiting" in phases:
        return "🟡 確認待ちあり" + ("・ほかのレーンは実行中" if "running" in phases else "")
    if "running" in phases:
        return "▶ 実行中" + ("・停止中のレーンあり" if "stopped" in phases else "")
    if "stopped" in phases:
        return "⏸ 停止中" + ("・処理完了のレーンあり" if "awaiting" in phases else "")
    if "awaiting" in phases and all(p in ("awaiting", "approved") for p in phases):
        return "☑ 処理完了・最終承認待ち"
    return "開始待ち"


#_____________________________________________________________________________
def progress_tail(directory):
    """Bound I/O even when a long calibration has a large journal."""
    try:
        with (directory / "progress.jsonl").open("rb") as stream:
            stream.seek(0, 2)
            stream.seek(max(0, stream.tell() - 65536))
            lines = stream.read().decode("utf-8", errors="replace").splitlines()
    except OSError:
        return []
    rows = []
    for line in lines:
        try:
            row = json.loads(line)
        except ValueError:
            continue
        if isinstance(row, dict):
            rows.append(row)
    return rows


#_____________________________________________________________________________
def summary_text(run):
    directory = Path(config.OUTPUT_DIR) / "calib" / f"run{run:05d}"
    parent = read_object(directory / "calib_dispatch_state.json")
    states = {lane: read_object(directory / "lanes" / lane / "calib_dispatch_state.json") for lane in LANES}
    states = {lane: state for lane, state in states.items() if state}
    # Actual lane records take precedence over a coordinator's cached status map.
    for lane, status in (parent.get("lanes") or {}).items():
        if lane in LANES and lane not in states:
            states[lane] = {"status": status}
    if not states and parent:
        states = {"calib": parent}
    global_approved = (parent.get("status") == "approved" and bool(states)
                       and all(phase_of(state)[0] in ("awaiting", "approved") for state in states.values()))
    phases = []
    lines = []
    stamps = [str(parent["updated_at"])] if global_approved and parent.get("updated_at") else []
    latest = {}
    for row in progress_tail(directory):
        if row.get("event") == "substage":
            latest[row.get("lane") or "calib"] = row
    for lane in (*LANES, "calib"):
        state = states.get(lane)
        if not state:
            continue
        phase, label = phase_of(state)
        if global_approved:
            phase, label = "approved", "✅ 承認済み"
        phases.append(phase)
        title = calib_approval.stage_title(str(state.get("step") or ""), state.get("opts"))
        if lane == "d5" and phase in ("awaiting", "approved", "stopped"):
            try:
                marker = (directory / "d5_converged.txt").read_text().strip()
                count = int(marker.removeprefix("pass="))
                title = f"第{count}周回で収束"
            except (OSError, ValueError):
                pass
        line = f"*{LABELS.get(lane, '較正')}*：{label}" + (f"｜{title}" if title else "")
        row = latest.get(lane)
        if phase == "running" and row and str(row.get("at") or "") >= str(state.get("updated_at") or ""):
            stage = str(row.get("stage") or "")
            if stage:
                prefix = stage.split("[", 1)[0]
                detail = {"ref": "基準の準備", "pos": "位置補正", "bcin_resi": "BLC1 residual補正",
                          "resi": "residual補正", "reso": "分解能補正", "hdprm": "ADC/TDC補正",
                          "hdphc": "PHC更新", "t0": "時間補正"}.get(prefix, prefix)
                iteration = re.search(r"(?:pass=|\[)(\d+)", stage)
                if iteration:
                    detail += f"・反復{iteration[1]}"
                line += f"（直近: {detail}）"
            if row.get("at"):
                stamps.append(str(row["at"]))
        if phase == "check":
            reason = str((state.get("hold") or {}).get("reason") or "")
            explanations = {"stale_input:": "入力の再decodeが必要", "qa_delivery_failed:": "QAの生成・送信を確認",
                            "write_blocked:": "更新の保留理由を確認", "phc_round_missing:": "前ラウンドの評価が未完了",
                            "phc_range_proposal:": "PHCの使用範囲を確認", "phc_no_update:": "有効なPHC更新なし",
                            "phc_score_invalid:": "PHCの幅評価が未完了"}
            for marker, explanation in explanations.items():
                if marker in reason:
                    line += f" — {explanation}"
                    break
        lines.append(line)
        if state.get("updated_at"):
            stamps.append(str(state["updated_at"]))
    heading = f"*E72 Calibration · run{run:05d}* — {overall_phase(phases)}"
    lines.insert(0, heading)
    if any(p in ("waiting", "check", "awaiting") for p in phases):
        lines.append("操作は、このスレッド内の対応する最新の確認メッセージから行ってください。")
    if stamps:
        stamp = max(stamps)
        try:
            stamp = datetime.fromisoformat(stamp.replace("Z", "+00:00")).astimezone(
                timezone(timedelta(hours=9))).strftime("%m/%d %H:%M JST")
        except ValueError:
            pass
        lines.append(f"状態の最終更新: {stamp}")
    return "\n".join(lines)


#_____________________________________________________________________________
def summary_blocks(run, text):
    """Use Slack's layout blocks instead of a dense paragraph."""
    lines = text.splitlines()
    overall = lines[0].split(" — ", 1)[-1] if lines else "状態記録なし"
    fields = {}
    context = []
    for line in lines[1:]:
        if line.startswith("*") and "*：" in line:
            name, content = line[1:].split("*：", 1)
            content, _, reason = content.partition(" — ")
            phase, _, stage = content.partition("｜")
            phase = phase.replace("処理完了・最終承認待ち", "最終承認待ち").replace("停止中（ユーザー操作）", "停止（ユーザー）")
            fields[name] = f"*{name}:* {phase}" + (f"｜{stage}" if stage else "")
            if reason:
                context.append(f"{name}: {reason}")
        elif line.startswith("状態の最終更新:"):
            context.append(line.replace("状態の最終更新:", "状態更新:"))
        elif line.startswith("操作は"):
            context.append("操作はスレッド内の最新の確認メッセージから。")
    blocks = [
        {"type": "header", "text": {"type": "plain_text", "text": f"run{run:05d} · 校正状況"}},
        {"type": "section", "text": {"type": "mrkdwn", "text": f"*全体*　{overall}"}},
        {"type": "divider"},
    ]
    if fields:
        ordered = [fields.pop(name) for name in ("Hodo", "BcIn / BLC1", "BcOut / BLC2", "D5") if name in fields]
        ordered.extend(fields.values())
        blocks.append({"type": "section", "text": {"type": "mrkdwn", "text": "\n".join(ordered)}})
    if context:
        blocks.append({"type": "context", "elements": [{"type": "mrkdwn", "text": "\n".join(context)}]})
    return blocks


#_____________________________________________________________________________
def parent_layout(run, text):
    blocks = summary_blocks(run, text)
    heading, attachments, context = [], [], []
    for block in blocks:
        body = (block.get("text") or {}).get("text", "")
        if block.get("type") == "section" and any(body.startswith(f"*{name}:") for name in (*LABELS.values(), "較正")):
            for row in body.splitlines():
                attachments.append({"color": slack_notify.status_color(row),
                                    "blocks": [{"type": "section", "text": {"type": "mrkdwn", "text": row}}]})
        elif block.get("type") == "context":
            context.append(block)
        else:
            heading.append(block)
    if context:
        if attachments:
            attachments[-1]["blocks"].extend(context)
        else:
            heading.extend(context)
    return {"blocks": heading, "attachments": attachments}


#_____________________________________________________________________________
def refresh_review_colors(run):
    """Restyle only still-active holds, including gates from older workers."""
    directory = Path(config.OUTPUT_DIR) / "calib" / f"run{run:05d}"
    for path in [directory / "calib_dispatch_state.json", *sorted(directory.glob("lanes/*/calib_dispatch_state.json"))]:
        if not path.is_file():
            continue
        with locked_path(path):
            state = read_object(path)
            hold = state.get("hold") or {}
            if state.get("status") != "held" or hold.get("decision") is not None:
                continue
            source, hold_id = hold.get("source"), hold.get("id")
            if not hold_id or source not in ("step_gate", "exit_check"):
                continue
            notes = slack_notify.load_state(run).get("failures") or {}
            keys = [key for key in notes if key.endswith(":" + str(hold_id))]
            if "phc_range_proposal:" in str(hold.get("reason") or ""):
                from . import phc_range_gate
                keys = [phc_range_gate.ask_slack_key(run)]
            for key in keys:
                note = notes.get(key) or {}
                if not note.get("ts") or note.get("color_layout") == 2:
                    continue
                if key.startswith("phc_range_ask"):
                    proposal = phc_range_gate.read_proposal(run) or {}
                    blocks = phc_range_gate.phc_range_action_blocks(run, proposal.get("proposals") or [])
                elif source == "step_gate":
                    blocks = calib_approval.step_gate_action_blocks(run, hold["step"], hold["index"], str(hold_id))
                    blocks[0]["text"]["text"] = note.get("text") or blocks[0]["text"]["text"]
                else:
                    blocks = calib_approval.exit_check_action_blocks(
                        run, hold["step"], hold["index"], str(hold_id),
                        lane=str(hold.get("lane") or ""), reason=str(hold.get("reason") or ""),
                    )
                ok, message = slack_notify.update_failure_message(run, key, note["text"], blocks=blocks)
                if not ok:
                    print(f"[Warning] Slack review color run{run:05d}: {message}", flush=True)


#_____________________________________________________________________________
def refresh_parent(run, *, config_path=None):
    """Update only the existing parent. Never create a thread or change calibration."""
    path = slack_notify.state_path(run)
    with locked_path(path.with_name("slack_parent_status.json")):
        text = summary_text(run)
        layout = parent_layout(run, text)
        signature = hashlib.sha256(json.dumps({"text": text, **layout}, ensure_ascii=False, sort_keys=True).encode()).hexdigest()
        cfg, status = slack_notify.load_slack_config(config_path)
        if status in ("missing", "disabled"):
            return True, f"slack skip: {status}"
        if status != "ok" or cfg is None:
            return False, f"slack config: {status}"
        with locked_path(path):
            state = slack_notify.load_state(run)
            ts = state.get("thread_ts")
            channel = state.get("channel_id")
            if not ts or channel != cfg["channel_id"]:
                return True, "skip: no existing thread in current channel"
            previous = state.get("header_status") or {}
            if (previous.get("signature"), previous.get("thread_ts"), previous.get("channel_id")) == (signature, ts, channel):
                return True, "unchanged"
        # Keep the PDF/notification state lock free while waiting for Slack.
        ok, _, message = slack_notify._api("chat.update", cfg["token"], {
            "channel": channel, "ts": ts, "text": text,
            **layout,
        })
        if not ok:
            return False, message
        with locked_path(path):
            current = slack_notify.load_state(run)
            if current.get("thread_ts") != ts or current.get("channel_id") != channel:
                return True, "thread changed during update; next cycle will refresh it"
            current["header_status"] = {"signature": signature, "thread_ts": ts, "channel_id": channel,
                                        "text": text, "updated_at": datetime.now(timezone.utc).isoformat()}
            if not slack_notify.save_state(run, current):
                return False, "header status save failed"
        return True, "updated"


#_____________________________________________________________________________
def start_status_loop(*, config_path=None, interval_s=POLL_SECONDS):
    """One daemon in the cluster bot also observes already-running dispatchers."""
    global _stop
    if _stop is not None:
        return _stop
    stop = threading.Event()
    _stop = stop

    def refresh():
        while not stop.is_set():
            base = Path(config.OUTPUT_DIR) / "calib"
            for path in sorted(base.glob("run*/slack_thread.json")):
                if stop.is_set():
                    return
                try:
                    run = int(path.parent.name.removeprefix("run"))
                    ok, message = refresh_parent(run, config_path=config_path)
                    refresh_review_colors(run)
                    if not ok:
                        print(f"[Warning] Slack status run{run:05d}: {message}", flush=True)
                    if message != "unchanged" and not message.startswith(("skip:", "slack skip:")):
                        stop.wait(1.3)
                except Exception as exc:
                    print(f"[Warning] Slack status: {type(exc).__name__}", flush=True)
            stop.wait(interval_s)

    threading.Thread(target=refresh, name="calib-slack-status", daemon=True).start()
    atexit.register(stop.set)
    return stop
