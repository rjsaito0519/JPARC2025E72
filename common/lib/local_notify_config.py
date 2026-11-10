"""Load gitignored local YAML notify configs (Slack / Discord).

No environment variables. Missing file, enabled:false, or missing keys → skip.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_DIR = PROJECT_ROOT / "calibration" / "config"

SLACK_LOCAL = CONFIG_DIR / "slack_local.yml"
DISCORD_LOCAL = CONFIG_DIR / "discord_local.yml"
OWNER_ONLY_MODE = 0o600


#_____________________________________________________________________________
def ensure_owner_only(path: Path) -> None:
    """Restrict a secret config file to owner read/write (chmod 0600).

    Group/other bits are cleared when present. Best-effort: permission errors are ignored.
    """
    path = Path(path)
    if not path.is_file():
        return
    try:
        mode = path.stat().st_mode & 0o777
        if mode & 0o077:
            path.chmod(OWNER_ONLY_MODE)
    except OSError:
        return


#_____________________________________________________________________________
def load_yaml_config(path: Path) -> tuple[Optional[dict[str, Any]], str]:
    """Return (data, status). status is 'ok', 'missing', 'disabled', or an error reason.

    On non-ok, data is None. Does not echo secret-looking values in the reason.
    Tightens file mode to 0600 before reading when group/other bits are set.
    """
    path = Path(path)
    if not path.is_file():
        return None, "missing"
    ensure_owner_only(path)
    try:
        raw = path.read_text(encoding="utf-8")
        data = yaml.safe_load(raw)
    except Exception as exc:  # noqa: BLE001
        return None, f"unreadable: {type(exc).__name__}"
    if data is None:
        data = {}
    if not isinstance(data, dict):
        return None, "not a mapping"
    enabled = data.get("enabled", False)
    if enabled is False or enabled is None:
        return None, "disabled"
    if enabled is not True:
        return None, "enabled must be true/false"
    return data, "ok"


#_____________________________________________________________________________
def load_slack_config(path: Path | None = None) -> tuple[Optional[dict[str, str]], str]:
    """Return ({token, channel_id}, status). Requires enabled:true and non-empty keys."""
    data, status = load_yaml_config(path or SLACK_LOCAL)
    if status != "ok" or data is None:
        return None, status
    token = str(data.get("token") or "").strip()
    channel_id = str(data.get("channel_id") or "").strip()
    if not token or not channel_id:
        return None, "missing token or channel_id"
    return {"token": token, "channel_id": channel_id}, "ok"


#_____________________________________________________________________________
def load_discord_config(path: Path | None = None) -> tuple[Optional[dict[str, str]], str]:
    """Return ({webhook_url}, status). Requires enabled:true and non-empty webhook_url."""
    data, status = load_yaml_config(path or DISCORD_LOCAL)
    if status != "ok" or data is None:
        return None, status
    webhook_url = str(data.get("webhook_url") or "").strip()
    if not webhook_url:
        return None, "missing webhook_url"
    return {"webhook_url": webhook_url}, "ok"
