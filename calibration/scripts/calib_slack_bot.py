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
  /calib phc-range <run> <det> <first> <last>
  /calib phc-go <run>
  /calib flag <run> <mode> <det> <channel> <reason...>
  /calib clear <run> <mode> <det> <channel>

Block Kit:
  - Approve / Revoke on awaiting-final (calib_approve_run / calib_revoke_run)
  - Proceed / Refit / Stop on step gate / EXIT_CHECK
    (calib_step_proceed / calib_step_refit / calib_step_stop)
  - Custom PHC range modal (calib_phc_custom → calib_phc_custom_modal)
  - Go PHC refit when all ranges set (calib_phc_go)

Does not run start or arbitrary shell. resume, proceed, retry,
phc-go (after conf write), and hdphc-again start resume work in
the background when no in-process waiter is present.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import (  # noqa: E402
    calib_ops,
    local_notify_config,
    phc_range_gate,
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
def _modal_values(view):
    """Extract plain_text_input / static_select values from a Slack view state."""
    state = (view or {}).get("state") or {}
    values = state.get("values") or {}
    out = {}
    for block_id, actions in values.items():
        if not isinstance(actions, dict):
            continue
        for action_id, payload in actions.items():
            if not isinstance(payload, dict):
                continue
            selected = payload.get("selected_option")
            if isinstance(selected, dict) and selected.get("value") is not None:
                text = str(selected.get("value") or "").strip()
            else:
                text = str(payload.get("value") or "").strip()
            out[action_id] = text
            out[block_id] = text
    return out


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
    app.action("calib_phc_go")(_button_action("phc-go"))

    @app.action("calib_phc_custom")
    def phc_custom_open(ack, body, client):
        ack()
        user = body.get("user") or {}
        user_id = str(user.get("id") or "") if isinstance(user, dict) else ""
        channel = body.get("channel") or {}
        if isinstance(channel, dict):
            channel_id = str(channel.get("id") or body.get("channel_id") or "")
        else:
            channel_id = str(body.get("channel_id") or "")
        runner.info(
            f"phc_custom click user={user_id} channel={channel_id}"
        )
        if not calib_ops.user_channel_allowed(bot_cfg, user_id, channel_id):
            runner.warn(
                f"phc_custom denied user={user_id} channel={channel_id}"
            )
            client.chat_postEphemeral(
                channel=channel_id or user_id,
                user=user_id,
                text="not allowed for this user/channel",
            )
            return
        actions = body.get("actions") or []
        run_s = ""
        if actions and isinstance(actions[0], dict):
            run_s = str(actions[0].get("value") or "").split(":", 1)[0]
        try:
            run = int(run_s)
        except (TypeError, ValueError):
            runner.warn(f"phc_custom bad run value: {run_s!r}")
            client.chat_postEphemeral(
                channel=channel_id or user_id,
                user=user_id,
                text=f"bad run value for PHC custom: {run_s!r}",
            )
            return
        proposal = phc_range_gate.read_proposal(run) or {}
        rows = proposal.get("proposals") if isinstance(proposal, dict) else None
        n_rows = len(rows) if isinstance(rows, list) else 0
        trigger = str(body.get("trigger_id") or "")
        view = phc_range_gate.phc_custom_modal_view(run, rows)
        runner.info(
            f"phc_custom views_open run={run:05d} rows={n_rows} "
            f"trigger={'yes' if trigger else 'no'}"
        )
        try:
            client.views_open(trigger_id=trigger, view=view)
        except Exception as exc:
            runner.warn(f"phc_custom views_open failed run={run:05d}: {exc}")
            try:
                client.chat_postEphemeral(
                    channel=channel_id or user_id,
                    user=user_id,
                    text=f"PHC Custom modal failed: {exc}",
                )
            except Exception:
                pass

    @app.view("calib_phc_custom_modal")
    def phc_custom_submit(ack, body, client):
        view = body.get("view") or {}
        vals = _modal_values(view)
        det = vals.get("det") or vals.get("phc_det") or ""
        first_s = vals.get("first") or vals.get("phc_first") or ""
        last_s = vals.get("last") or vals.get("phc_last") or ""
        runner.info(
            f"phc_custom submit det={det!r} first={first_s!r} last={last_s!r}"
        )
        try:
            run = int(str(view.get("private_metadata") or ""))
            first = int(first_s)
            last = int(last_s)
        except (TypeError, ValueError):
            ack(
                response_action="errors",
                errors={
                    "phc_first": "first/last must be integers",
                    "phc_last": "first/last must be integers",
                },
            )
            return
        if not det:
            ack(
                response_action="errors",
                errors={"phc_det": "detector is required"},
            )
            return
        if first > last:
            ack(
                response_action="errors",
                errors={"phc_last": "last must be >= first"},
            )
            return
        ack()
        user = body.get("user") or {}
        user_id = str(user.get("id") or "") if isinstance(user, dict) else ""
        ok, message = calib_ops.run_dispatch(
            "phc-range",
            run,
            user_id=user_id,
            det=det,
            first=first,
            last=last,
        )
        calib_ops.record_op(run, "phc-range", user_id)
        text = f"run{run:05d} phc-range: {message}"
        runner.info(text)
        try:
            client.chat_postMessage(channel=user_id, text=text)
        except Exception:
            runner.warn(f"phc modal DM failed; {text}")
        if not ok:
            runner.warn(text)

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
