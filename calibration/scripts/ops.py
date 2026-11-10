#!/usr/bin/env python3
"""Calibration state ops (separate from ./calib.py launcher).

Works from another terminal / cw node that shares the same paths.py:
  state lives under OUTPUT_DIR/calib/runXXXXX/.

Usage:
  ./ops.py status  <run>
  ./ops.py stop    <run>
  ./ops.py resume  <run> [--local] [--dry]
  ./ops.py hold    <run> [reason…]
  ./ops.py decide  <run> refit|continue [--note TEXT]
  ./ops.py approve <run> [--note TEXT]
  ./ops.py revoke  <run> [--reason TEXT]

Primary judgment path is Slack /calib …. stop/hold are cooperative
(honored at the next step boundary, not an immediate kill).
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import calib_dispatch as cd  # noqa: E402
import calib_session as cs  # noqa: E402
from common.lib import runner  # noqa: E402


#_____________________________________________________________________________
def resolve_cfg(run, *, local=None):
    path = cs.session_path(run)
    if not path.is_file():
        runner.fail(
            f"session not found: {path} "
            f"(start with ./calib.py {run} first)"
        )
    if local is True:
        cs.set_session_local(path, True)
    elif local is False:
        cs.set_session_local(path, False)
    _, cfg = cd.load_config(path)
    for i, raw in enumerate(cfg["steps"]):
        cd.normalize_step(raw, i)
    return path, cfg


#_____________________________________________________________________________
def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument(
        "command",
        choices=(
            "status", "stop", "resume", "hold",
            "decide", "approve", "revoke",
        ),
    )
    p.add_argument("run", type=int, help="run number")
    p.add_argument(
        "rest", nargs="*",
        help="hold: reason words; decide: refit|continue",
    )
    p.add_argument("--local", action="store_true",
                   help="resume: myrun --local")
    p.add_argument("--dry", "--dry-run", dest="dry_run", action="store_true",
                   help="resume: print iterate commands only")
    p.add_argument("--note", default=None, help="decide/approve note")
    p.add_argument("--reason", default=None, help="hold/revoke reason")
    return p.parse_args()


#_____________________________________________________________________________
def main():
    a = parse_args()
    if a.run < 1:
        runner.fail("run must be >= 1")

    if a.command == "status":
        _, cfg = resolve_cfg(a.run)
        cd.do_status(cfg)

    if a.command == "stop":
        _, cfg = resolve_cfg(a.run)
        cd.do_stop(cfg)

    if a.command == "resume":
        _, cfg = resolve_cfg(a.run, local=True if a.local else None)
        cd.do_resume(cfg, a.dry_run)

    if a.command == "hold":
        reason = a.reason
        if a.rest:
            reason = " ".join(a.rest) if reason is None else reason
        _, cfg = resolve_cfg(a.run)
        cd.do_hold(cfg, reason)

    if a.command == "decide":
        if not a.rest:
            runner.fail("decide requires refit|continue")
        action = a.rest[0]
        if action not in cd.DECIDE_ACTIONS:
            runner.fail(f"decide action must be refit|continue, got {action!r}")
        if len(a.rest) > 1 and a.note is None:
            a.note = " ".join(a.rest[1:])
        _, cfg = resolve_cfg(a.run)
        cd.do_decide(cfg, action, a.note)

    if a.command == "approve":
        _, cfg = resolve_cfg(a.run)
        cd.do_approve(cfg, a.note, "cli", None, None)

    if a.command == "revoke":
        reason = a.reason
        if a.rest and reason is None:
            reason = " ".join(a.rest)
        _, cfg = resolve_cfg(a.run)
        cd.do_revoke(cfg, reason)

    runner.fail(f"unhandled command {a.command!r}")


if __name__ == "__main__":
    main()
