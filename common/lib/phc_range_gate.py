"""PHC usable-segment range ask (Slack Custom then Go).

When auto unused-seg ratio is high, ask for hand ranges on the whole hdphc
detector set (PDF first, then Slack). Custom writes phc_conf and refreshes
the same Slack message; Go refits when every det has a hand range.
Hand ranges live in calibration/scripts/phc_conf.py.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Optional

from . import config

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PHC_CONF_PATH = PROJECT_ROOT / "calibration" / "scripts" / "phc_conf.py"
PROPOSAL_NAME = "phc_range_proposal.json"
# Stable Slack failure key so chat.update can refresh the same message.
SLACK_ASK_KEY = "phc_range_ask"
# Single-segment (etc.): no usable-range ask / Custom choice.
RANGE_ASK_SKIP = frozenset({"T1"})


#_____________________________________________________________________________
def needs_range_ask(det: str) -> bool:
    return str(det or "") not in RANGE_ASK_SKIP

#_____________________________________________________________________________
def proposal_path(run: int) -> Path:
    return (
        Path(config.OUTPUT_DIR) / "calib" / f"run{int(run):05d}" / PROPOSAL_NAME
    )


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
def row_has_suggest(row: dict[str, Any]) -> bool:
    return row.get("first") is not None and row.get("last") is not None


#_____________________________________________________________________________
def read_hand_range(
    run: int,
    suffix: str,
    det: str,
    *,
    path: Path | None = None,
) -> Optional[tuple[int, int]]:
    """Inclusive hand range from phc_conf.py, or None if unset."""
    conf_path = Path(path) if path is not None else PHC_CONF_PATH
    if not conf_path.is_file():
        return None
    try:
        text = conf_path.read_text(encoding="utf-8")
    except OSError:
        return None
    key = conf_key(run, suffix, det)
    match = re.search(
        rf'"{re.escape(key)}"\s*:\s*\[\s*(-?\d+)\s*,\s*(-?\d+)\s*\]',
        text,
    )
    if not match:
        return None
    return int(match.group(1)), int(match.group(2))


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
    """Insert or replace limits_dict entry. Returns the conf key written."""
    conf_path = Path(path) if path is not None else PHC_CONF_PATH
    if not conf_path.is_file():
        raise FileNotFoundError(f"phc_conf not found: {conf_path}")
    if int(first) > int(last):
        raise ValueError(f"first ({first}) > last ({last})")
    key = conf_key(run, suffix, det)
    new_line = f'    "{key}" : [{int(first):2d}, {int(last):2d}],\n'
    text = conf_path.read_text(encoding="utf-8")
    pat = re.compile(
        rf'^[ \t]*"{re.escape(key)}"\s*:\s*\[[^\]]*\]\s*,?[ \t]*\n',
        re.MULTILINE,
    )
    if pat.search(text):
        text = pat.sub(new_line, text, count=1)
    else:
        # Insert before the closing brace of limits_dict.
        close = text.rfind("}")
        if close < 0:
            raise ValueError(f"no closing brace in {conf_path}")
        insert_at = close
        # Keep a blank line before } when the previous line is not blank.
        prefix = text[:insert_at].rstrip() + "\n\n"
        text = prefix + new_line + text[insert_at:]
    conf_path.write_text(text, encoding="utf-8")
    return key


#_____________________________________________________________________________
def phc_range_action_blocks(run: int, proposals: list[dict[str, Any]]) -> list:
    """Block Kit Custom + Go (refit when every det has a hand range)."""
    run_s = str(int(run))
    lines = [format_proposal_line(row) for row in proposals if isinstance(row, dict)]
    body = "\n".join(lines) if lines else "(no proposals)"
    return [
        {
            "type": "section",
            "text": {
                "type": "mrkdwn",
                "text": (
                    f"*run{int(run):05d}* PHC segment range\n"
                    f"{body}"
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
                    "text": {"type": "plain_text", "text": "Custom range"},
                    "value": run_s,
                },
                {
                    "type": "button",
                    "action_id": "calib_phc_go",
                    "text": {"type": "plain_text", "text": "Go（refit）"},
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
        SLACK_ASK_KEY,
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
        SLACK_ASK_KEY,
        text,
        blocks=phc_range_action_blocks(int(run), rows),
    )


#_____________________________________________________________________________
def notify_pending_phc_ask(run: int) -> tuple[bool, str]:
    """Post Slack PHC ask from on-disk proposal; no-op if none."""
    data = sync_proposal_currents(run) or read_proposal(run)
    if not data:
        return True, "no phc ask pending"
    rows = data.get("proposals")
    if not isinstance(rows, list) or not rows:
        return True, "no phc ask pending"
    return notify_phc_range_slack(int(run), rows)


#_____________________________________________________________________________
def conf_keys_on_disk(*, path: Path | None = None) -> set[str]:
    conf_path = Path(path) if path is not None else PHC_CONF_PATH
    if not conf_path.is_file():
        return set()
    try:
        text = conf_path.read_text(encoding="utf-8")
    except OSError:
        return set()
    return set(re.findall(r'"(\d{5}_[^"]+)"\s*:', text))


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
