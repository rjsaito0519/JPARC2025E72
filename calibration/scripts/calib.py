#!/usr/bin/env python3
"""Daily calibration launcher: ensure session and start (or resume if stopped).

Usage:
  ./calib.py <run> [--pion|--kaon] [--hodo] [--bcout] [--bcin] [--d5]
             [--local] [--dry] [--force]
  ./calib.py <run> --init [--local] [--dry]

If a prior dispatch state exists:
  - TTY: ask 続きから / 初めから
  - non-interactive: 続きから
  - --force: 初めから (not while running)

--init: prepare Hodo.root for TDC window setting (create_runlist hodo 0 +
decode if missing). Does not start a calib session.

State ops: ./ops.py status|stop|… <run>  (or Slack /calib …).
Prefer this over calib_session.py / calib_dispatch.py.
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import calib_dispatch as cd  # noqa: E402
import calib_session as cs  # noqa: E402
from common.lib import (  # noqa: E402
    calib_approval,
    calib_runlist,
    decode_run,
    hodo_btof_monitor,
    runner,
)

CHOICE_CONTINUE = "continue"
CHOICE_RESTART = "restart"


#_____________________________________________________________________________
def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int, help="run number")
    p.add_argument("--init", action="store_true",
                   help="prepare Hodo.root for TDC windows; no session/chain")
    p.add_argument("--pion", action="store_true", help="calib window pion (default)")
    p.add_argument("--kaon", action="store_true", help="calib window kaon")
    p.add_argument("--hodo", action="store_true", help="include hodo step")
    p.add_argument("--bcout", action="store_true", help="include BcOut dc_chamber")
    p.add_argument("--bcin", action="store_true", help="include BcIn dc_chamber")
    p.add_argument("--d5", action="store_true", help="include d5_align")
    p.add_argument("--local", action="store_true",
                   help="myrun --local (default: bsub)")
    p.add_argument("--dry", "--dry-run", dest="dry_run", action="store_true",
                   help="print iterate commands only")
    p.add_argument("--force", action="store_true",
                   help="restart from the beginning (rewrite session; not while running)")
    p.add_argument("--momentum", type=float, default=None,
                   help="beam momentum [GeV/c] when table/UserParam missing")
    return p.parse_args()


#_____________________________________________________________________________
def announce_chain(run, calib_particle, steps):
    """Short launch banner: calib owns create_runlist / UserParam start."""
    step_names = []
    for raw in steps:
        if isinstance(raw, dict):
            step_names.append(next(iter(raw)))
        else:
            step_names.append(str(raw))
    runner.info(
        f"calib run{run:05d} particle={calib_particle} "
        f"steps={','.join(step_names) or '(default)'}"
    )
    runner.info(
        "Each step runs create_runlist → decode → fit "
        "(UserParam / conf / runlist are created in-chain; "
        "do not hand-create runlists first)."
    )
    if "hodo" in step_names or not step_names:
        runner.info(
            "Hodo: prepare with ./calib.py <run> --init, set TDC windows, "
            "then start this chain. If windows still match the USER template, "
            "the chain stops (EXIT_CHECK); then "
            f"./ops.py decide {run} refit && ./ops.py resume {run}."
        )


#_____________________________________________________________________________
def do_init(run, *, local=False, dry_run=False):
    """Ensure Hodo.root exists for window setting; do not start calib."""
    root = hodo_btof_monitor.hodo_root(run)
    runner.info(f"calib --init run{run:05d} Hodo.root={root}")
    if root.is_file():
        runner.info(f"Hodo.root already present: {root}")
    else:
        runner.info(f"Hodo.root missing; create_runlist --hodo {run} 0 → decode")
        calib_runlist.run_create_runlist("--hodo", run, "0", dry_run=dry_run)
        runlist = calib_runlist.hodo_runlist_path(run)
        if dry_run:
            decode_run.run_decode(runlist, dry_run=True, local=local)
        else:
            if not runlist.is_file():
                runner.fail(f"hodo runlist missing after create_runlist: {runlist}")
            decode_run.run_decode(runlist, local=local)
            if not root.is_file():
                runner.fail(f"Hodo.root still missing after decode: {root}")
            runner.info(f"Hodo.root ready: {root}")
    calib_runlist.announce_tdc_window(run, after_init=True)
    sys.exit(runner.EXIT_OK)


#_____________________________________________________________________________
def format_state_summary(run, st):
    """One-block summary for the continue/restart prompt."""
    status = calib_approval.effective_status(st) or st.get("status")
    hold = st.get("hold") if isinstance(st.get("hold"), dict) else {}
    lines = [
        f"run{run:05d}: prior calib state found",
        f"  status={status} step={st.get('step')} index={st.get('index')} "
        f"rc={st.get('rc')}",
    ]
    if hold:
        lines.append(
            f"  hold.source={hold.get('source')} "
            f"decision={hold.get('decision')!r} "
            f"reason={hold.get('reason')!r}"
        )
    hist = st.get("steps") if isinstance(st.get("steps"), list) else []
    if hist:
        brief = ", ".join(
            f"{e.get('step')}(rc={e.get('rc')})"
            for e in hist[-6:] if isinstance(e, dict)
        )
        lines.append(f"  steps_hist: {brief}")
    tail = cd.read_progress_tail(run, n=5)
    if tail:
        lines.append("  progress.jsonl (tail):")
        for row in tail:
            lines.append(
                f"    {row.get('at')} {row.get('event')} "
                f"step={row.get('step')} status={row.get('status')}"
            )
    return "\n".join(lines)


#_____________________________________________________________________________
def prompt_continue_or_restart(run, st):
    """TTY: ask continue vs restart. Returns CHOICE_*."""
    print(format_state_summary(run, st), file=sys.stderr)
    print(
        "Continue from this history, or restart from the beginning?\n"
        "  [1] 続きから (continue)\n"
        "  [2] 初めから (restart)\n"
        f"(non-interactive default is continue; --force forces restart)",
        file=sys.stderr,
    )
    while True:
        try:
            raw = input(f"run{run:05d} choose [1/2]: ").strip().lower()
        except EOFError:
            runner.warn("EOF on prompt; defaulting to continue")
            return CHOICE_CONTINUE
        if raw in ("1", "c", "continue", "続き", "続きから"):
            return CHOICE_CONTINUE
        if raw in ("2", "r", "restart", "初め", "初めから"):
            return CHOICE_RESTART
        print("enter 1 (continue) or 2 (restart)", file=sys.stderr)


#_____________________________________________________________________________
def choose_action(run, st, *, force, interactive):
    """Decide continue vs restart given existing state."""
    status = calib_approval.effective_status(st) or st.get("status")
    if st.get("status") == "running" or status == "running":
        runner.fail(
            f"run{run:05d}: already running; "
            f"./ops.py stop {run} first"
        )
    if force:
        return CHOICE_RESTART
    if calib_approval.is_process_complete(st):
        if interactive and sys.stdin.isatty():
            print(format_state_summary(run, st), file=sys.stderr)
            print(
                "Process already complete "
                f"({status}). Restart from the beginning?\n"
                "  [2] 初めから (restart)\n"
                "  anything else: abort",
                file=sys.stderr,
            )
            try:
                raw = input(f"run{run:05d} choose [2 to restart]: ").strip().lower()
            except EOFError:
                raw = ""
            if raw in ("2", "r", "restart", "初め", "初めから"):
                return CHOICE_RESTART
            runner.fail(
                f"run{run:05d}: process complete ({status}); "
                f"approve/revoke via Slack/ops, or --force to restart"
            )
        runner.fail(
            f"run{run:05d}: process complete ({status}); "
            f"approve/revoke via Slack/ops, or --force to restart"
        )
    if interactive and sys.stdin.isatty():
        return prompt_continue_or_restart(run, st)
    return CHOICE_CONTINUE


#_____________________________________________________________________________
def continue_from_state(cfg, st, dry_run):
    """Resume from stopped, or decide-continue then resume from held."""
    status = st.get("status")
    if status == "stopped":
        runner.info(f"run{cfg['run']:05d}: continuing (resume from stopped)")
        cd.do_resume(cfg, dry_run)
    if status == "held":
        hold = st.get("hold") if isinstance(st.get("hold"), dict) else {}
        runner.info(
            f"run{cfg['run']:05d}: continuing from held "
            f"(source={hold.get('source')}); decide continue → resume"
        )
        cd.apply_decide(cfg["run"], "continue", note="calib.py continue")
        cd.do_resume(cfg, dry_run)
    runner.fail(
        f"run{cfg['run']:05d}: cannot continue from status={status!r}; "
        f"use ./ops.py status {cfg['run']}"
    )


#_____________________________________________________________________________
def restart_from_beginning(a, calib_particle, steps, target):
    """Clear dispatch state + progress, rewrite session, start at step 0."""
    run = a.run
    runner.info(f"run{run:05d}: restarting from the beginning")
    if not a.dry_run:
        cd.clear_dispatch_state(run)
        cd.reset_progress(run)
        cd.append_progress(run, "restart", status="cleared")
    path, _doc = cs.ensure_session(
        run,
        calib_particle=calib_particle,
        steps=steps,
        target=target,
        momentum=a.momentum,
        local=a.local,
        overwrite=True,
        prompt_if_missing=True,
    )
    if a.local:
        cs.set_session_local(path, True)
    _, cfg2 = cd.load_config(path)
    for i, raw in enumerate(cfg2["steps"]):
        cd.normalize_step(raw, i)
    cd.do_start(cfg2, a.dry_run)


#_____________________________________________________________________________
def main():
    a = parse_args()
    if a.run < 1:
        runner.fail("run must be >= 1")

    if a.init:
        do_init(a.run, local=a.local, dry_run=a.dry_run)

    calib_particle = cs.resolve_calib_particle(pion=a.pion, kaon=a.kaon)
    steps, target = cs.steps_from_flags(
        hodo=a.hodo, bcout=a.bcout, bcin=a.bcin, d5=a.d5,
    )
    announce_chain(a.run, calib_particle, steps)

    st = cd.read_state(a.run)
    interactive = (not a.dry_run)  # dry-run: no TTY ask; treat as continue/start

    if st is not None:
        choice = choose_action(
            a.run, st, force=a.force, interactive=interactive and sys.stdin.isatty(),
        )
        if choice == CHOICE_RESTART:
            restart_from_beginning(a, calib_particle, steps, target)
            return

        # continue: keep existing session unless missing
        path, _doc = cs.ensure_session(
            a.run,
            calib_particle=calib_particle,
            steps=steps,
            target=target,
            momentum=a.momentum,
            local=a.local,
            overwrite=False,
            prompt_if_missing=True,
        )
        if a.local:
            cs.set_session_local(path, True)
        _, cfg = cd.load_config(path)
        for i, raw in enumerate(cfg["steps"]):
            cd.normalize_step(raw, i)
        # Re-read in case ensure_session side effects (none expected)
        st = cd.read_state(a.run) or st
        continue_from_state(cfg, st, a.dry_run)
        return

    # No prior state → fresh start
    path, _doc = cs.ensure_session(
        a.run,
        calib_particle=calib_particle,
        steps=steps,
        target=target,
        momentum=a.momentum,
        local=a.local,
        overwrite=a.force,
        prompt_if_missing=True,
    )
    if a.local:
        cs.set_session_local(path, True)
    _, cfg = cd.load_config(path)
    for i, raw in enumerate(cfg["steps"]):
        cd.normalize_step(raw, i)
    if not a.dry_run:
        cd.reset_progress(a.run)
        cd.append_progress(a.run, "start", status="running")
    cd.do_start(cfg, a.dry_run)


if __name__ == "__main__":
    main()
