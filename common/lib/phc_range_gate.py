"""PHC usable-segment range proposal (Slack Accept / Custom).

When auto unused-seg ratio is high, propose min–max of good segments instead
of a hard anomaly stop. Hand ranges live in calibration/scripts/phc_conf.py.
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
    """Block Kit Accept button; Custom is slash-command only."""
    run_s = str(int(run))
    lines = []
    for row in proposals:
        det = str(row.get("det") or "?")
        first = row.get("first")
        last = row.get("last")
        bad = row.get("bad")
        total = row.get("total")
        lines.append(
            f"• `{det}` → [{first}, {last}] "
            f"(unused {bad}/{total})"
        )
    body = "\n".join(lines) if lines else "(no proposals)"
    return [
        {
            "type": "section",
            "text": {
                "type": "mrkdwn",
                "text": (
                    f"*run{int(run):05d}* PHC usable-segment range proposal\n"
                    f"{body}\n"
                    f"Accept = write all to `phc_conf.py` and refit hdphc. "
                    f"Custom: `/calib phc-range {int(run)} <det> <first> <last>`"
                ),
            },
        },
        {
            "type": "actions",
            "block_id": f"calib_phc_range_{run_s}",
            "elements": [
                {
                    "type": "button",
                    "action_id": "calib_phc_accept",
                    "text": {"type": "plain_text", "text": "Accept 提案レンジ"},
                    "style": "primary",
                    "value": run_s,
                },
            ],
        },
    ]


#_____________________________________________________________________________
def notify_phc_range_slack(
    run: int, proposals: list[dict[str, Any]],
) -> tuple[bool, str]:
    from . import slack_notify

    n = len(proposals)
    dets = ",".join(str(p.get("det") or "?") for p in proposals)
    text = (
        f"run{int(run):05d}: PHC range proposal for {n} detector(s) [{dets}]. "
        f"Accept or /calib phc-range {int(run)} <det> <first> <last>"
    )
    return slack_notify.try_notify_failure(
        int(run),
        f"phc_range_proposal:{dets}",
        text,
        blocks=phc_range_action_blocks(int(run), proposals),
    )


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
    """Detectors in the proposal that still lack a hand range in phc_conf."""
    keys = conf_keys_on_disk(path=path)
    missing = []
    for row in proposals:
        det = str(row.get("det") or "")
        if not det:
            continue
        if conf_key(run, suffix, det) not in keys:
            missing.append(det)
    return missing
