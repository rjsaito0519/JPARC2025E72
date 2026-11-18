"""PHC usable-segment range ask (Slack Custom then Go).

When auto unused-seg ratio is high, ask for hand ranges on the whole hdphc
detector set (PDF first, then Slack). Custom writes phc_conf and refreshes
the same Slack message; Go refits when every det has a hand range.
Hand ranges live in calibration/scripts/phc_conf.py.
"""
from __future__ import annotations

import ast
import json
import os
import tempfile
import re
from pathlib import Path
from typing import Any, Optional

from . import config
from .param import exclusive_param_lock

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PHC_CONF_PATH = PROJECT_ROOT / "calibration" / "scripts" / "phc_conf.py"
PROPOSAL_NAME = "phc_range_proposal.json"
CONFIRMED_NAME = "phc_range_confirmed.json"
# Stable Slack failure key so chat.update can refresh the same message.
SLACK_ASK_KEY = "phc_range_ask"
# Single-segment (etc.): no usable-range ask / Custom choice.
RANGE_ASK_SKIP = frozenset({"T1"})


#_____________________________________________________________________________
def needs_range_ask(det: str) -> bool:
    return str(det or "") not in RANGE_ASK_SKIP

#_____________________________________________________________________________
def calib_run_dir(run: int) -> Path:
    return Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}"


#_____________________________________________________________________________
def proposal_path(run: int) -> Path:
    return calib_run_dir(run) / PROPOSAL_NAME


#_____________________________________________________________________________
def confirmed_path(run: int) -> Path:
    return calib_run_dir(run) / CONFIRMED_NAME


#_____________________________________________________________________________
def conf_key(run: int, suffix: str, det: str) -> str:
    return f"{int(run):05d}_{suffix}_{det.lower()}"


#_____________________________________________________________________________
def propose_from_good(
    ch_list: list[int], bad_list: list[int],
) -> Optional[tuple[int, int]]:
    """Inclusive [first, last] from good channels; None if none remain."""
    bad = {int(c) for c in bad_list}
    good = sorted(int(c) for c in ch_list if int(c) not in bad)
    if not good:
        return None
    return int(good[0]), int(good[-1])


#_____________________________________________________________________________
def write_proposal(run: int, payload: dict[str, Any]) -> Path:
    path = proposal_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return path


#_____________________________________________________________________________
def read_proposal(run: int) -> Optional[dict[str, Any]]:
    path = proposal_path(run)
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return data if isinstance(data, dict) else None


#_____________________________________________________________________________
def pending_phc_ask(run: int) -> bool:
    """True when a PHC range-ask proposal is on disk for this run."""
    data = read_proposal(run)
    if not data:
        return False
    rows = data.get("proposals")
    return isinstance(rows, list) and bool(rows)


#_____________________________________________________________________________
def clear_proposal(run: int) -> bool:
    """Remove on-disk PHC ask proposal. Returns True if a file was removed."""
    path = proposal_path(run)
    if not path.is_file():
        return False
    path.unlink()
    return True


#_____________________________________________________________________________
def is_range_confirmed(run: int) -> bool:
    """True after Go in this calib session (survives resume, not restart)."""
    path = confirmed_path(run)
    if not path.is_file():
        return False
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return isinstance(data, dict) and data.get("confirmed") is True


#_____________________________________________________________________________
def set_range_confirmed(run: int) -> Path:
    """Mark PHC hand ranges confirmed for this session (after Go)."""
    from datetime import datetime, timezone

    path = confirmed_path(run)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "run": int(run),
        "confirmed": True,
        "at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return path


#_____________________________________________________________________________
def clear_range_confirmed(run: int) -> bool:
    """Clear session confirm flag so the next hdphc re-asks. Returns if removed."""
    path = confirmed_path(run)
    if not path.is_file():
        return False
    path.unlink()
    return True


#_____________________________________________________________________________
def retire_proposal_flags(run: int) -> int:
    """Invalidate automatic markers belonging to the retired range proposal."""
    from . import channel_flags
    retired = 0
    for row in channel_flags.open_needs_review(run, mode="hdphc"):
        if row.get("source") != "auto" or not str(row.get("reason") or "").startswith("phc_range_proposal:"):
            continue
        channel_flags.clear_channel(
            run, mode="hdphc", det=row["det"], channel=row["channel"],
            reason="PHC proposal retired by restart; new session requires range confirmation",
        )
        retired += 1
    return retired


#_____________________________________________________________________________
def clear_phc_ask_state(run: int) -> None:
    """Drop proposal, its automatic review markers, and confirmed on restart."""
    clear_proposal(run)
    clear_range_confirmed(run)
    retire_proposal_flags(run)


#_____________________________________________________________________________
def row_has_suggest(row: dict[str, Any]) -> bool:
    return row.get("first") is not None and row.get("last") is not None


#_____________________________________________________________________________
def _hand_range_nodes(text: str, path: Path):
    """Parse literal limits_dict without executing configuration code."""
    try:
        tree = ast.parse(text, filename=str(path))
    except SyntaxError as error:
        raise ValueError(f"Invalid PHC configuration: {path}") from error
    dictionaries = []
    for node in tree.body:
        names = node.targets if isinstance(node, ast.Assign) else (
            [node.target] if isinstance(node, ast.AnnAssign) else []
        )
        if any(isinstance(name, ast.Name) and name.id == "limits_dict" for name in names):
            dictionaries.append(node.value)
    if len(dictionaries) != 1 or not isinstance(dictionaries[0], ast.Dict):
        raise ValueError(f"PHC configuration requires one literal limits_dict: {path}")
    mapping = dictionaries[0]
    # Updates or indexed assignments would differ from a literal declaration.
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)) and node.value is mapping:
            continue
        if any(isinstance(name, ast.Name) and name.id == "limits_dict"
               for name in ast.walk(node)):
            raise ValueError(f"Dynamic limits_dict updates are not supported: {path}")
    entries = {}
    for key_node, value_node in zip(mapping.keys, mapping.values):
        try:
            key = ast.literal_eval(key_node)
            value = ast.literal_eval(value_node)
        except (ValueError, TypeError) as error:
            raise ValueError(f"PHC ranges must be literal values: {path}") from error
        if not isinstance(key, str) or key in entries:
            raise ValueError(f"Invalid or duplicate PHC range key: {path}")
        if value is not None and (
            not isinstance(value, (list, tuple)) or len(value) != 2
            or any(type(bound) is not int for bound in value)
            or value[0] > value[1]
        ):
            raise ValueError(f"Invalid PHC range for {key}: {path}")
        entries[key] = (None if value is None else tuple(value), value_node)
    return mapping, entries


#_____________________________________________________________________________
def _source_offset(data: bytes, node, *, end=False) -> int:
    lines = data.splitlines(keepends=True)
    line = node.end_lineno if end else node.lineno
    column = node.end_col_offset if end else node.col_offset
    return sum(len(value) for value in lines[:line-1]) + column


#_____________________________________________________________________________
def read_hand_range(
    run: int,
    suffix: str,
    det: str,
    *,
    path: Path | None = None,
) -> Optional[tuple[int, int]]:
    """Read the explicit range; invalid configuration is not an automatic choice."""
    conf_path = Path(path) if path is not None else PHC_CONF_PATH
    if not conf_path.is_file():
        return None
    _, entries = _hand_range_nodes(conf_path.read_text(encoding="utf-8"), conf_path)
    return entries.get(conf_key(run, suffix, det), (None, None))[0]


#_____________________________________________________________________________
def format_proposal_line(row: dict[str, Any]) -> str:
    """One det line; ranges in code spans (Slack has no colored text)."""
    det = str(row.get("det") or "?")
    cur_f, cur_l = row.get("current_first"), row.get("current_last")
    if cur_f is not None and cur_l is not None:
        return (
            f":large_green_circle: *`{det}`*  "
            f"*range* `{cur_f}`–`{cur_l}`"
        )
    if row_has_suggest(row):
        return (
            f":large_orange_circle: *`{det}`*  "
            f"_hint_ `{row.get('first')}`–`{row.get('last')}`  "
            f"— Custom で確定"
        )
    return (
        f":red_circle: *`{det}`*  "
        f"_range unset_  — Custom で指定"
    )

#_____________________________________________________________________________
def upsert_phc_conf(
    run: int,
    suffix: str,
    det: str,
    first: int,
    last: int,
    *,
    path: Path | None = None,
) -> str:
    """Persist an explicit user range without replacing unrelated settings."""
    conf_path = (Path(path) if path is not None else PHC_CONF_PATH).resolve()
    if not conf_path.is_file():
        raise FileNotFoundError(f"phc_conf not found: {conf_path}")
    if int(first) > int(last):
        raise ValueError(f"first ({first}) > last ({last})")
    key = conf_key(run, suffix, det)
    value = f"[{int(first):2d}, {int(last):2d}]".encode("utf-8")
    with exclusive_param_lock(conf_path):
        text = conf_path.read_text(encoding="utf-8")
        mapping, entries = _hand_range_nodes(text, conf_path)
        data = text.encode("utf-8")
        if key in entries:
            node = entries[key][1]
            start = _source_offset(data, node)
            end = _source_offset(data, node, end=True)
            updated = data[:start] + value + data[end:]
        else:
            close = _source_offset(data, mapping, end=True) - 1
            if mapping.values:
                end = _source_offset(data, mapping.values[-1], end=True)
                remainder = b"\n".join(line.split(b"#", 1)[0]
                                         for line in data[end:close].splitlines())
                if b"," not in remainder:
                    data = data[:end] + b"," + data[end:]
                    close += 1
            line_start = data.rfind(b"\n", 0, close) + 1
            indentation = data[line_start:close]
            entry = json.dumps(key).encode("utf-8") + b" : " + value + b",\n"
            if not indentation.strip():
                updated = data[:line_start] + indentation + b"    " + entry + data[line_start:]
            else:
                indentation = b" " * mapping.col_offset
                updated = data[:close] + b"\n" + indentation + b"    " + entry + indentation + data[close:]
        _hand_range_nodes(updated.decode("utf-8"), conf_path)
        if updated == text.encode("utf-8"):
            return key
        temporary_path = None
        try:
            with tempfile.NamedTemporaryFile(dir=conf_path.parent, prefix=conf_path.name + ".",
                                             suffix=".tmp", delete=False) as stream:
                temporary_path = Path(stream.name)
                os.fchmod(stream.fileno(), conf_path.stat().st_mode & 0o777)
                stream.write(updated)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary_path, conf_path)
        finally:
            if temporary_path is not None and temporary_path.exists():
                temporary_path.unlink()
    return key


#_____________________________________________________________________________
def phc_range_action_blocks(run: int, proposals: list[dict[str, Any]]) -> list:
    """Block Kit Custom + Go (refit when every det has a hand range)."""
    run_s = str(int(run))
    lines = [format_proposal_line(row) for row in proposals if isinstance(row, dict)]
    body = "\n".join(lines) if lines else "(no proposals)"
    result = str((read_proposal(run) or {}).get("qa_exec_id") or "").split("_")[-1]
    if result:
        body += f"\n対象PDFの結果ID: {result}"
    return [
        {
            "type": "section",
            "text": {
                "type": "mrkdwn",
                "text": (
                    f"*run{int(run):05d}* 【確認待ち】PHC 使用範囲\n"
                    "直前のPHC PDFを確認してください。このラウンドの更新はまだ適用していません。\n"
                    f"{body}\n"
                    "範囲を変更する場合は『範囲を編集』、この範囲で再計算する場合は『範囲を確定して再実行』を選んでください。"
                ),
            },
        },
        {
            "type": "actions",
            "block_id": f"calib_phc_range_{run_s}",
            "elements": [
                {
                    "type": "button",
                    "action_id": "calib_phc_custom",
                    "text": {"type": "plain_text", "text": "範囲を編集"},
                    "value": run_s,
                },
                {
                    "type": "button",
                    "action_id": "calib_phc_go",
                    "text": {"type": "plain_text", "text": "範囲を確定して再実行"},
                    "style": "primary",
                    "value": run_s,
                },
            ],
        },
    ]

#_____________________________________________________________________________
def _placeholder_range(
    proposals: list[dict[str, Any]] | None,
) -> tuple[str, str]:
    if not isinstance(proposals, list):
        return "0", "31"
    for row in proposals:
        if not isinstance(row, dict):
            continue
        if row_has_suggest(row):
            return str(row["first"]), str(row["last"])
        cur_f, cur_l = row.get("current_first"), row.get("current_last")
        if cur_f is not None and cur_l is not None:
            return str(cur_f), str(cur_l)
    return "0", "31"


#_____________________________________________________________________________
def _det_select_options(
    proposals: list[dict[str, Any]] | None,
) -> list[dict[str, Any]]:
    """Slack static_select options from the ask list (fallback: beam set)."""
    dets: list[str] = []
    if isinstance(proposals, list):
        for row in proposals:
            if isinstance(row, dict) and row.get("det"):
                det = str(row["det"])
                if det and det not in dets:
                    dets.append(det)
    if not dets:
        dets = ["BHT", "BH2", "CVC"]
    options = []
    for det in dets:
        options.append({
            "text": {"type": "plain_text", "text": det},
            "value": det,
        })
    return options


#_____________________________________________________________________________
def phc_custom_modal_view(run: int, proposals: list[dict[str, Any]] | None = None) -> dict:
    """Modal view payload for Custom PHC range (det select / first / last)."""
    run_i = int(run)
    options = _det_select_options(proposals)
    first_ph, last_ph = _placeholder_range(proposals)
    det_element: dict[str, Any] = {
        "type": "static_select",
        "action_id": "det",
        "placeholder": {"type": "plain_text", "text": "Select detector"},
        "options": options,
    }
    if options:
        det_element["initial_option"] = options[0]
    return {
        "type": "modal",
        "callback_id": "calib_phc_custom_modal",
        "private_metadata": str(run_i),
        "title": {"type": "plain_text", "text": f"PHC range {run_i:05d}"},
        "submit": {"type": "plain_text", "text": "Apply"},
        "close": {"type": "plain_text", "text": "Cancel"},
        "blocks": [
            {
                "type": "input",
                "block_id": "phc_det",
                "label": {"type": "plain_text", "text": "Detector"},
                "element": det_element,
            },
            {
                "type": "input",
                "block_id": "phc_first",
                "label": {"type": "plain_text", "text": "First segment"},
                "element": {
                    "type": "plain_text_input",
                    "action_id": "first",
                    "placeholder": {"type": "plain_text", "text": first_ph},
                },
            },
            {
                "type": "input",
                "block_id": "phc_last",
                "label": {"type": "plain_text", "text": "Last segment"},
                "element": {
                    "type": "plain_text_input",
                    "action_id": "last",
                    "placeholder": {"type": "plain_text", "text": last_ph},
                },
            },
        ],
    }


#_____________________________________________________________________________
def sync_proposal_currents(
    run: int,
    *,
    path: Path | None = None,
) -> Optional[dict[str, Any]]:
    """Refresh current_* on the proposal from phc_conf; rewrite file."""
    data = read_proposal(run)
    if not data:
        return None
    suffix = str(data.get("suffix") or "")
    rows = data.get("proposals")
    if not isinstance(rows, list) or not rows:
        return data
    for row in rows:
        if not isinstance(row, dict):
            continue
        det = str(row.get("det") or "")
        if not det or not suffix:
            continue
        hand = read_hand_range(run, suffix, det, path=path)
        if hand is None:
            row["current_first"] = None
            row["current_last"] = None
        else:
            row["current_first"], row["current_last"] = hand
    write_proposal(run, data)
    return data


#_____________________________________________________________________________
def ask_slack_key(run: int) -> str:
    proposal = read_proposal(run) or {}
    result = str(proposal.get("qa_exec_id") or "")
    return f"{SLACK_ASK_KEY}:{result}" if result else SLACK_ASK_KEY


#_____________________________________________________________________________
def notify_phc_range_slack(
    run: int, proposals: list[dict[str, Any]],
) -> tuple[bool, str]:
    from . import slack_notify

    n = len(proposals)
    dets = ",".join(str(p.get("det") or "?") for p in proposals)
    text = (
        f"run{int(run):05d}: PHC range ask for {n} detector(s) [{dets}]. "
        f"Custom then Go, or /calib phc-range {int(run)} <det> <first> <last>"
    )
    return slack_notify.try_notify_failure(
        int(run),
        ask_slack_key(int(run)),
        text,
        blocks=phc_range_action_blocks(int(run), proposals),
    )


#_____________________________________________________________________________
def update_phc_ask_slack(run: int) -> tuple[bool, str]:
    """Rewrite the existing PHC ask Slack message from on-disk proposal+conf."""
    from . import slack_notify

    data = sync_proposal_currents(run)
    if not data:
        return False, "no phc ask proposal on disk"
    rows = data.get("proposals")
    if not isinstance(rows, list) or not rows:
        return False, "empty phc ask proposal"
    n = len(rows)
    dets = ",".join(str(p.get("det") or "?") for p in rows if isinstance(p, dict))
    text = (
        f"run{int(run):05d}: PHC range ask for {n} detector(s) [{dets}]. "
        f"Custom then Go"
    )
    return slack_notify.try_update_failure_message(
        int(run),
        ask_slack_key(int(run)),
        text,
        blocks=phc_range_action_blocks(int(run), rows),
    )


#_____________________________________________________________________________
def notify_pending_phc_ask(run: int, *, qa_exec_id: str | None = None) -> tuple[bool, str]:
    """Post Slack PHC ask from on-disk proposal; no-op if none."""
    data = sync_proposal_currents(run) or read_proposal(run)
    if not data:
        return True, "no phc ask pending"
    if qa_exec_id:
        data["qa_exec_id"] = qa_exec_id
        write_proposal(run, data)
    rows = data.get("proposals")
    if not isinstance(rows, list) or not rows:
        return True, "no phc ask pending"
    return notify_phc_range_slack(int(run), rows)


#_____________________________________________________________________________
def conf_keys_on_disk(*, path: Path | None = None) -> set[str]:
    conf_path = Path(path) if path is not None else PHC_CONF_PATH
    if not conf_path.is_file():
        return set()
    _, entries = _hand_range_nodes(conf_path.read_text(encoding="utf-8"), conf_path)
    return {key for key, (value, _) in entries.items() if value is not None}


#_____________________________________________________________________________
def missing_proposal_dets(
    run: int,
    suffix: str,
    proposals: list[dict[str, Any]],
    *,
    path: Path | None = None,
) -> list[str]:
    """Detectors in the ask list that still lack a hand range in phc_conf."""
    keys = conf_keys_on_disk(path=path)
    missing = []
    for row in proposals:
        det = str(row.get("det") or "")
        if not det:
            continue
        if conf_key(run, suffix, det) not in keys:
            missing.append(det)
    return missing
