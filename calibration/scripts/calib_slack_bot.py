#!/usr/bin/env python3
"""Slack Socket Mode bot for calib_dispatch / extra_hdphc_iter ops.

Usage:
  pip install slack-bolt   # optional dependency
  # fill calibration/config/slack_local.yml (bot_enabled, app_token, allowlists)
  python3 calibration/scripts/calib_slack_bot.py

Slash command (register /calib in the Slack app):
  /calib status <run>
  /calib stop <run>
  /calib resume <run>
  /calib approve <run>
  /calib revoke <run>
  /calib hold <run> <reason...>
  /calib decide <run> refit|continue
  /calib hdphc-again <run>
  /calib phc-accept <run>
  /calib phc-range <run> <det> <first> <last>
  /calib flag <run> <mode> <det> <channel> <reason...>
  /calib clear <run> <mode> <det> <channel>

Block Kit:
  - Approve / Revoke on awaiting-final (calib_approve_run / calib_revoke_run)
  - Proceed / Refit / Stop on step gate
    (calib_step_proceed / calib_step_refit / calib_step_stop)
  - Accept PHC range proposal (calib_phc_accept)

Does not run start or arbitrary shell. resume, proceed, retry,
phc-accept (after conf write), and hdphc-again start resume work in
the background.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import (  # noqa: E402
    calib_ops,
    local_notify_config,
    runner,
    slack_bot_ensure,
)


#_____________________________________________________________________________
def handle_calib_request(bot_cfg, user_id, channel_id, text):
    """Core handler used by Bolt and by unit tests (no Slack I/O).

    Returns (ok, reply_text, ephemeral).
    ephemeral True → tell only the invoking user.
    """
    if not calib_ops.user_channel_allowed(bot_cfg, user_id, channel_id):
        return False, "not allowed for this user/channel", True

    op, run, err, extras = calib_ops.parse_calib_text(text)
    if err:
        return False, err, True

    dup, dup_msg = calib_ops.check_dedup(run, op, user_id)
    if dup:
        return True, dup_msg, True

    ok, message = calib_ops.run_dispatch(
        op,
        run,
        user_id=user_id,
        reason=extras.get("reason"),
        action=extras.get("action"),
        mode=extras.get("mode"),
        det=extras.get("det"),
        channel=extras.get("channel"),
        subject_sha=extras.get("subject_sha"),
        first=extras.get("first"),
        last=extras.get("last"),
        hold_id=extras.get("hold_id"),
    )
    calib_ops.record_op(run, op, user_id)
    prefix = f"run{run:05d} {op}: "
    return ok, prefix + message, False


#_____________________________________________________________________________
def build_app(bot_cfg):
    try:
        from slack_bolt import App
    except ImportError as exc:
        runner.fail(
            "slack-bolt is required for the ops bot "
            f"(pip install slack-bolt): {exc}"
        )

    app = App(token=bot_cfg["token"])

    @app.command("/calib")
    def calib_command(ack, body, respond):
        ack()
        user_id = str(body.get("user_id") or "")
        channel_id = str(body.get("channel_id") or "")
        text = str(body.get("text") or "")
        ok, message, ephemeral = handle_calib_request(
            bot_cfg, user_id, channel_id, text,
        )
        respond(text=message, response_type="ephemeral" if ephemeral else "in_channel")
        if not ok and not ephemeral:
            runner.warn(message)

    def _button_action(op):
        def _handler(ack, body, respond):
            ack()
            user = body.get("user") or {}
            if isinstance(user, dict):
                user_id = str(user.get("id") or "")
            else:
                user_id = str(user or "")
            channel = body.get("channel") or {}
            if isinstance(channel, dict):
                channel_id = str(
                    channel.get("id") or body.get("channel_id") or ""
                )
            else:
                channel_id = str(
                    body.get("channel_id") or channel or ""
                )
            actions = body.get("actions") or []
            value = ""
            if actions and isinstance(actions[0], dict):
                value = str(actions[0].get("value") or "")
            text = f"{op} {value}".strip()
            ok, message, ephemeral = handle_calib_request(
                bot_cfg, user_id, channel_id, text,
            )
            respond(
                text=message,
                response_type="ephemeral" if ephemeral else "in_channel",
            )
            if not ok and not ephemeral:
                runner.warn(message)
        return _handler

    app.action("calib_approve_run")(_button_action("approve"))
    app.action("calib_revoke_run")(_button_action("revoke"))
    app.action("calib_step_proceed")(_button_action("proceed"))
    app.action("calib_step_refit")(_button_action("retry"))
    app.action("calib_step_stop")(_button_action("gate-stop"))
    app.action("calib_phc_accept")(_button_action("phc-accept"))

    return app


#_____________________________________________________________________________
def main():
    p = runner.Parser(__doc__)
    p.add_argument(
        "--config",
        type=Path,
        default=None,
        help="slack_local.yml path (default: calibration/config/slack_local.yml)",
    )
    a = p.parse_args()

    bot_cfg, status = local_notify_config.load_slack_bot_config(a.config)
    if status != "ok" or bot_cfg is None:
        runner.fail(f"slack bot config: {status}")

    try:
        from slack_bolt.adapter.socket_mode import SocketModeHandler
    except ImportError as exc:
        runner.fail(
            "slack-bolt is required for the ops bot "
            f"(pip install slack-bolt): {exc}"
        )

    app = build_app(bot_cfg)
    slack_bot_ensure.start_heartbeat_loop()
    runner.info(
        "calib_slack_bot listening (Socket Mode); "
        f"users={len(bot_cfg['allowed_user_ids'])} "
        f"channels={len(bot_cfg['allowed_channel_ids'])}; "
        f"heartbeat={slack_bot_ensure.heartbeat_path()}"
    )
    SocketModeHandler(app, bot_cfg["app_token"]).start()


if __name__ == "__main__":
    main()
