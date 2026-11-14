#!/usr/bin/env python3
"""Daily calibration launcher: ensure session and start (or resume if stopped).

Usage:
  ./calib.py <run> [--pion|--kaon] [--hodo] [--bcout] [--bcin] [--d5]
             [--local] [--dry] [--force] [--auto]
  ./calib.py <run> --init [--local] [--dry]
  ./calib.py --list calibration/caliblist/<name>.yml [--local] [--dry] [--force]

If a prior dispatch state exists:
  - TTY: ask 続きから / 初めから
  - non-interactive: 続きから
  - --force: 初めから (not while running)

--init: prepare Hodo.root for TDC window setting (create_runlist hodo 0 +
decode if missing). Does not start a calib session.

--list: run DEFAULT+RUN batch from caliblist YAML in parallel
(local: max 5 simultaneous; bsub: uncapped). Each run waits on its
own Slack gates / PHC modal as usual.

State ops: ./ops.py status|stop|… <run>  (or Slack /calib …).
Prefer this over calib_session.py / calib_dispatch.py.
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

sys.path.append(str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

CALIB_SCRIPT = Path(__file__).resolve()

import calib_dispatch as cd  # noqa: E402
import calib_session as cs  # noqa: E402
from common.lib import (  # noqa: E402
    calib_approval,
    calib_lanes,
    calib_list,
    calib_runlist,
    decode_run,
    hodo_btof_monitor,
    runner,
    slack_bot_ensure,
)

CHOICE_CONTINUE = "continue"
CHOICE_RESTART = "restart"
CHOICE_ATTACH = "attach"


#_____________________________________________________________________________
def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument(
        "run", type=int, nargs="?", default=None,
        help="run number (omit with --list)",
    )
    p.add_argument(
        "--list", type=Path, default=None,
        help="caliblist YAML (DEFAULT+RUN batch)",
    )
    p.add_argument("--list-entry", type=json.loads, default=None, help=argparse.SUPPRESS)
    p.add_argument("--init", action="store_true",
                   help="prepare Hodo.root for TDC windows; no session/chain")
    p.add_argument("--pion", action="store_true", help="calib window pion (default)")
    p.add_argument("--kaon", action="store_true", help="calib window kaon")
    p.add_argument("--hodo", action="store_true", help="include hodo step")
    p.add_argument("--bcout", action="store_true", help="include BcOut dc_chamber")
    p.add_argument("--bcin", action="store_true", help="include BcIn dc_chamber")
    p.add_argument(
        "--d5", action="store_true",
        help="include d5_align (also forces bcout+bcin)",
    )
    p.add_argument("--local", action="store_true",
                   help="myrun --local (default: bsub)")
    p.add_argument("--dry", "--dry-run", dest="dry_run", action="store_true",
                   help="print iterate commands only")
    p.add_argument("--force", action="store_true",
                   help="restart from the beginning (rewrite session; not while running)")
    p.add_argument(
        "--auto", action="store_true",
        help="run without Slack step gates (coarse hodo/dc_chamber/d5_align)",
    )
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
    if (
        any(n == "hodo" or n.startswith("hodo_") for n in step_names)
        or not step_names
    ):
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
    lanes = st.get("lanes")
    if isinstance(lanes, dict) and lanes:
        lines.append(f"  lanes={lanes}")
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
            runner.fail("no selection received; calibration was not started")
        if raw in ("1", "c", "continue", "続き", "続きから"):
            return CHOICE_CONTINUE
        if raw in ("2", "r", "restart", "初め", "初めから"):
            return CHOICE_RESTART
        print("enter 1 (continue) or 2 (restart)", file=sys.stderr)


#_____________________________________________________________________________
def restart_blockers(run, lane_status):
    """A held state without a gate waiter has no worker left to disrupt."""
    blockers = set(calib_lanes.running_lanes(lane_status))
    for lane in lane_status:
        state = cd.read_state(run, lane=lane) or {}
        if state.get("gate_waiter"):
            blockers.add(lane)
    parent = cd.read_state(run, lane="") or {}
    if parent.get("gate_waiter"):
        blockers.add("parent")
    return blockers


#_____________________________________________________________________________
def choose_action(run, st, *, force, interactive, want_lanes=None):
    """Decide continue / restart / attach given existing state."""
    status = calib_approval.effective_status(st) or st.get("status")
    lane_status = cd.collect_lane_statuses(run)
    running = calib_lanes.running_lanes(lane_status)
    held = calib_lanes.held_lanes(lane_status)
    stopped = calib_lanes.stopped_lanes(lane_status)
    busy = calib_lanes.busy_lanes(lane_status)
    want = list(want_lanes or [])
    want_running = calib_lanes.overlap_lanes(want, running) if want else []
    want_held = calib_lanes.overlap_lanes(want, held) if want else []
    want_stopped = calib_lanes.overlap_lanes(want, stopped) if want else []
    want_free = (
        calib_lanes.free_requested_lanes(want, lane_status) if want else []
    )

    if force:
        blockers = restart_blockers(run, lane_status)
        if blockers:
            runner.fail(
                f"run{run:05d}: active workers on lane(s) {sorted(blockers)}; "
                f"./ops.py stop {run} first (cannot --force while busy)"
            )
        return CHOICE_RESTART

    # True conflict: same lane already executing.
    if want_running:
        runner.fail(
            f"run{run:05d}: lane(s) {want_running} already running; "
            f"./ops.py stop {run} first"
        )

    # Unstarted lanes while other work exists → attach (no TTY prompt).
    if want_free and (running or held or stopped or lane_status):
        return CHOICE_ATTACH
    if want_free and (
        st.get("status") == "running" or status == "running"
        or st.get("mode") == "lane_coordinator"
    ):
        return CHOICE_ATTACH

    # Existing requested lanes need the same explicit TTY choice as a full run.
    if want_held or want_stopped:
        if interactive and sys.stdin.isatty():
            choice = prompt_continue_or_restart(run, st)
            if choice == CHOICE_RESTART:
                blockers = restart_blockers(run, lane_status)
                if blockers:
                    runner.fail(
                        f"run{run:05d}: active workers on lane(s) {sorted(blockers)}; "
                        "stop running work or use Slack Stop for a waiting gate first"
                    )
            return choice
        return CHOICE_CONTINUE

    if not want and (st.get("status") == "running" or status == "running"):
        runner.fail(
            f"run{run:05d}: already running; "
            f"./ops.py stop {run} first"
        )

    if calib_approval.is_process_complete(st) and not busy:
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
def _maybe_ensure_slack_bot(a, cfg):
    """Ensure cluster Slack bot unless dry-run / --auto (no interactive gates)."""
    if a.dry_run or a.auto or cfg.get("auto"):
        return
    ok_bot, bot_msg = slack_bot_ensure.ensure_slack_bot(required=True)
    if not ok_bot:
        runner.fail(f"slack bot required for interactive gates: {bot_msg}")
    runner.info(bot_msg)


#_____________________________________________________________________________
def _resume_held_lane(cfg, lane, dry_run):
    """decide continue + resume for one held lane (OK while others run)."""
    run = int(cfg["run"])
    lane_st = cd.read_state(run, lane=lane)
    if not lane_st or lane_st.get("status") != "held":
        runner.fail(f"run{run:05d}: lane [{lane}] is not held")
    hold = lane_st.get("hold") if isinstance(lane_st.get("hold"), dict) else {}
    hold_id = hold.get("id")
    runner.info(
        f"run{run:05d}: continuing held lane [{lane}] "
        f"(source={hold.get('source')}); decide continue → resume"
    )
    if dry_run:
        runner.info(f"[dry-run] would continue held lane [{lane}]")
        return
    cd.apply_decide(
        run, "continue", note="calib.py continue", hold_id=hold_id,
    )
    if lane_st.get("gate_waiter"):
        runner.info(f"lane [{lane}]: in-process waiter will continue")
        return
    # do_resume loads lanes/<lane>/session.yml for plan match.
    cd.do_resume(cfg, dry_run, lane=lane)


#_____________________________________________________________________________
def continue_from_state(cfg, st, dry_run, want_lanes=None):
    """Resume stopped/held work; lane-aware when coordinator is running."""
    run = int(cfg["run"])
    want = list(want_lanes or [])
    lane_status = cd.collect_lane_statuses(run)
    status = st.get("status")

    # Coordinator still running other lanes: resume held/stopped want, or attach.
    if (
        st.get("mode") == "lane_coordinator"
        and status == "running"
    ):
        held = calib_lanes.held_lanes(lane_status)
        stopped = calib_lanes.stopped_lanes(lane_status)
        want_held = calib_lanes.overlap_lanes(want, held) if want else []
        want_stopped = calib_lanes.overlap_lanes(want, stopped) if want else []
        want_free = (
            calib_lanes.free_requested_lanes(want, lane_status) if want else []
        )
        if want_held:
            for lane in want_held:
                _resume_held_lane(cfg, lane, dry_run)
            return
        if want_stopped:
            for lane in want_stopped:
                runner.info(f"run{run:05d}: resuming stopped lane [{lane}]")
                cd.do_resume(cfg, dry_run, lane=lane)
            return
        if want_free:
            runner.info(
                f"run{run:05d}: coordinator running; "
                f"attaching free lanes {want_free}"
            )
            cd.do_attach(cfg, dry_run, want_lanes=want_free)
            return
        runner.fail(
            f"run{run:05d}: coordinator running "
            f"(lanes={lane_status}); nothing to continue for {want or 'all'}; "
            f"use ./ops.py status {run}"
        )

    # Prefer explicit want lanes over top-level state.
    if want:
        for lane in want:
            lane_st = cd.read_state(run, lane=lane)
            if not lane_st:
                continue
            if lane_st.get("status") == "held":
                _resume_held_lane(cfg, lane, dry_run)
                return
            if lane_st.get("status") == "stopped":
                runner.info(f"run{run:05d}: resuming stopped lane [{lane}]")
                cd.do_resume(cfg, dry_run, lane=lane)
                return

    if status == "stopped":
        runner.info(f"run{run:05d}: continuing (resume from stopped)")
        cd.do_resume(cfg, dry_run)
        return
    if status == "held":
        hold = st.get("hold") if isinstance(st.get("hold"), dict) else {}
        runner.info(
            f"run{run:05d}: continuing from held "
            f"(source={hold.get('source')}); decide continue → resume"
        )
        cd.apply_decide(run, "continue", note="calib.py continue")
        cd.do_resume(cfg, dry_run)
        return
    runner.fail(
        f"run{run:05d}: cannot continue from status={status!r}; "
        f"use ./ops.py status {run}"
    )


#_____________________________________________________________________________
def restart_from_beginning(a, calib_particle, steps, target, *, max_iter=10, hdphc_iter=10):
    """Clear dispatch state + progress, rewrite session, start at step 0."""
    run = a.run
    runner.info(f"run{run:05d}: restarting from the beginning")
    if not a.dry_run:
        blockers = restart_blockers(run, cd.collect_lane_statuses(run))
        if blockers:
            runner.fail(f"run{run:05d}: active workers on lane(s) {sorted(blockers)}; cannot restart")
        cd.clear_dispatch_state(run)
        cd.reset_progress(run)
        cd.append_progress(run, "restart", status="cleared")
    path, _doc = cs.ensure_session(
        run,
        calib_particle=calib_particle,
        steps=steps,
        target=target,
        momentum=a.momentum,
        beam_particle=getattr(a, "beam_particle", None),
        charge=getattr(a, "charge", None),
        write_table=bool(getattr(a, "write_table", False)) and not a.dry_run,
        max_iter=max_iter,
        local=a.local,
        auto=a.auto,
        overwrite=True,
        prompt_if_missing=a.momentum is None,
        hdphc_iter=hdphc_iter,
    )
    if a.local:
        cs.set_session_local(path, True)
    _, cfg2 = cd.load_config(path)
    for i, raw in enumerate(cfg2["steps"]):
        cd.normalize_step(raw, i)
    _maybe_ensure_slack_bot(a, cfg2)
    cd.do_start(cfg2, a.dry_run)


#_____________________________________________________________________________
def attach_complementary_lanes(
    a, calib_particle, steps, target, *, max_iter=10, hdphc_iter=10,
):
    """Merge CLI lanes into the session and start only the free lanes."""
    run = a.run
    path = cs.session_path(run)
    if path.is_file():
        path, doc = cs.load_session_file(path)
        merged = cs.merge_session_steps(
            doc.get("steps") or [],
            steps,
            auto=bool(a.auto or doc.get("auto")),
            hdphc_iter=hdphc_iter,
            max_iter=max_iter,
        )
        doc["steps"] = merged
        doc["target"] = target or doc.get("target") or "custom"
        if a.local:
            doc["local"] = True
        cs.write_session(path, doc)
        runner.info(f"run{run:05d}: attached lanes into session {path}")
    else:
        path, doc = cs.ensure_session(
            run,
            calib_particle=calib_particle,
            steps=steps,
            target=target,
            momentum=a.momentum,
            beam_particle=getattr(a, "beam_particle", None),
            charge=getattr(a, "charge", None),
            write_table=bool(getattr(a, "write_table", False)) and not a.dry_run,
            max_iter=max_iter,
            local=a.local,
            auto=a.auto,
            overwrite=True,
            prompt_if_missing=a.momentum is None,
            hdphc_iter=hdphc_iter,
        )
    if a.local:
        cs.set_session_local(path, True)
    _, cfg2 = cd.load_config(path)
    for i, raw in enumerate(cfg2["steps"]):
        cd.normalize_step(raw, i)
    _maybe_ensure_slack_bot(a, cfg2)
    # Session may already list other lanes; attach only CLI-requested free ones.
    want = calib_lanes.requested_lanes(
        cs.expand_coarse_steps(
            steps, auto=bool(a.auto or doc.get("auto")),
            hdphc_iter=hdphc_iter, max_iter=max_iter,
        )
    )
    runner.info(
        f"run{run:05d}: attaching complementary lanes {want}"
    )
    cd.do_attach(cfg2, a.dry_run, want_lanes=want)


#_____________________________________________________________________________
def opts_to_calib_argv(run, opts, *, cli_local=False, cli_dry=False, cli_force=False) -> list[str]:
    """Pass the complete normalized entry through argv without shell interpolation."""
    entry = dict(opts)
    for name, enabled in (("local", cli_local), ("dry_run", cli_dry), ("force", cli_force)):
        if enabled:
            entry[name] = True
    return [sys.executable, str(CALIB_SCRIPT), str(int(run)),
            "--list-entry", json.dumps(entry, allow_nan=False)]


#_____________________________________________________________________________
def apply_list_entry(args):
    """Use the same single-run launcher while retaining all YAML options."""
    if not isinstance(args.list_entry, dict):
        runner.fail("caliblist entry must be a mapping")
    try:
        opts = calib_list.normalize_run_opts(args.list_entry)
        particle = cs.normalize_calib_particle(opts["calib_particle"])
        if opts.get("beam_particle") is not None:
            opts["beam_particle"] = cs.normalize_beam_particle(opts["beam_particle"])
    except (TypeError, ValueError) as error:
        runner.fail(f"invalid caliblist entry: {error}")
    for name, value in opts.items():
        setattr(args, name, value)
    args.pion, args.kaon = particle == "pion", particle == "kaon"


#_____________________________________________________________________________
def run_one_from_opts(run, opts, *, cli_local=False):
    """Start one run from caliblist-normalized opts (always fresh start)."""
    local = bool(opts.get("local")) or bool(cli_local)
    dry_run = bool(opts.get("dry_run"))
    auto = bool(opts.get("auto"))
    force = bool(opts.get("force"))
    calib_particle = cs.normalize_calib_particle(opts.get("calib_particle", "pion"))
    beam_particle = opts.get("beam_particle")
    if beam_particle is not None:
        beam_particle = cs.normalize_beam_particle(beam_particle)
    steps, target = cs.steps_from_flags(
        hodo=opts.get("hodo", False),
        bcout=opts.get("bcout", False),
        bcin=opts.get("bcin", False),
        d5=opts.get("d5", False),
    )
    announce_chain(run, calib_particle, steps)
    a = SimpleNamespace(
        run=run,
        dry_run=dry_run,
        local=local,
        auto=auto,
        force=force,
        momentum=opts.get("momentum"),
        beam_particle=beam_particle,
        charge=opts.get("charge"),
        write_table=bool(opts.get("write_table")),
    )
    restart_from_beginning(
        a, calib_particle, steps, target,
        max_iter=int(opts.get("max_iter", 10)),
        hdphc_iter=int(opts.get("hdphc_iter", 10)),
    )


#_____________________________________________________________________________
def do_list(path, *, cli_local=False, cli_dry=False, cli_force=False):
    """Run caliblist entries in parallel (local capped, bsub uncapped)."""
    entries = calib_list.load_caliblist(path)
    runner.info(
        f"caliblist {path}: {len(entries)} run(s); "
        f"parallel local_max={calib_list.MAX_PARALLEL_LOCAL} bsub=unlimited"
    )

    def spawn(run, opts):
        argv = opts_to_calib_argv(
            run, opts, cli_local=cli_local, cli_dry=cli_dry, cli_force=cli_force,
        )
        # Each child uses the normal launcher; the entry controls restart versus resume.
        return subprocess.Popen(
            argv,
            cwd=str(CALIB_SCRIPT.parents[2]),
            env=os.environ.copy(),
            stdin=subprocess.DEVNULL,
        )

    agg_rc, finished = calib_list.run_caliblist_pool(
        entries,
        cli_local=cli_local,
        spawn=spawn,
    )
    runner.info(
        f"caliblist done: {path} "
        f"({len(finished)} finished, agg_rc={agg_rc})"
    )
    sys.exit(int(agg_rc) if agg_rc else runner.EXIT_OK)


#_____________________________________________________________________________
def main():
    a = parse_args()
    if a.list is not None:
        if a.init:
            runner.fail("--list cannot be combined with --init")
        if a.run is not None:
            runner.fail("--list does not take a run argument")
        if a.list_entry is not None:
            runner.fail("--list cannot be combined with --list-entry")
        if a.pion or a.kaon or a.hodo or a.bcout or a.bcin or a.d5 or a.auto or a.momentum is not None:
            runner.fail("with --list, set particle/stages/auto/momentum in the YAML")
        do_list(a.list, cli_local=a.local, cli_dry=a.dry_run, cli_force=a.force)
        return

    if a.list_entry is not None:
        apply_list_entry(a)
        if a.dry_run and getattr(a, "write_table", False):
            runner.info("[dry-run] beam table update requested; table is left unchanged")
    max_iter = int(getattr(a, "max_iter", 10))
    hdphc_iter = int(getattr(a, "hdphc_iter", 10))

    if a.run is None:
        runner.fail("run number required (or pass --list PATH)")
    if a.run < 1:
        runner.fail("run must be >= 1")

    if a.init:
        do_init(a.run, local=a.local, dry_run=a.dry_run)

    calib_particle = cs.resolve_calib_particle(pion=a.pion, kaon=a.kaon)
    steps, target = cs.steps_from_flags(
        hodo=a.hodo, bcout=a.bcout, bcin=a.bcin, d5=a.d5,
    )
    announce_chain(a.run, calib_particle, steps)
    want_lanes = calib_lanes.requested_lanes(
        cs.expand_coarse_steps(
            steps, auto=a.auto, hdphc_iter=hdphc_iter, max_iter=max_iter,
        )
    )

    st = cd.read_state(a.run, lane="")
    lane_status = cd.collect_lane_statuses(a.run)
    has_prior = st is not None or bool(lane_status)
    interactive = (not a.dry_run)  # dry-run: no TTY ask; treat as continue/start

    if has_prior:
        st_eff = st or {
            "status": "running" if calib_lanes.busy_lanes(lane_status) else "stopped",
            "lanes": dict(lane_status),
            "mode": "lane_coordinator",
        }
        if a.list_entry is not None and not a.force and calib_approval.is_process_complete(st_eff) and not calib_lanes.busy_lanes(lane_status):
            runner.info(f"run{a.run:05d}: already complete; preserving state (set force: true to restart)")
            return
        choice = choose_action(
            a.run, st_eff,
            force=a.force,
            interactive=interactive and sys.stdin.isatty(),
            want_lanes=want_lanes,
        )
        if choice == CHOICE_CONTINUE and a.list_entry is not None:
            held = calib_lanes.overlap_lanes(want_lanes, calib_lanes.held_lanes(lane_status))
            if held or st_eff.get("status") == "held":
                runner.info(f"run{a.run:05d}: confirmation still pending; preserving hold (use Slack/ops to decide)")
                return
        if choice == CHOICE_RESTART:
            restart_from_beginning(a, calib_particle, steps, target, max_iter=max_iter, hdphc_iter=hdphc_iter)
            return
        if choice == CHOICE_ATTACH:
            attach_complementary_lanes(a, calib_particle, steps, target, max_iter=max_iter, hdphc_iter=hdphc_iter)
            return

        # continue: keep existing session unless missing
        path, _doc = cs.ensure_session(
            a.run,
            calib_particle=calib_particle,
            steps=steps,
            target=target,
            momentum=a.momentum,
            beam_particle=getattr(a, "beam_particle", None),
            charge=getattr(a, "charge", None),
            write_table=bool(getattr(a, "write_table", False)) and not a.dry_run,
            max_iter=max_iter,
            hdphc_iter=hdphc_iter,
            local=a.local,
            auto=a.auto,
            overwrite=False,
            prompt_if_missing=True,
        )
        if a.local:
            cs.set_session_local(path, True)
        _, cfg = cd.load_config(path)
        for i, raw in enumerate(cfg["steps"]):
            cd.normalize_step(raw, i)
        # Re-read in case ensure_session side effects (none expected)
        st = cd.read_state(a.run, lane="") or st_eff
        _maybe_ensure_slack_bot(a, cfg)
        continue_from_state(cfg, st, a.dry_run, want_lanes=want_lanes)
        return

    # No prior state → fresh start
    path, _doc = cs.ensure_session(
        a.run,
        calib_particle=calib_particle,
        steps=steps,
        target=target,
        momentum=a.momentum,
        beam_particle=getattr(a, "beam_particle", None),
        charge=getattr(a, "charge", None),
        write_table=bool(getattr(a, "write_table", False)) and not a.dry_run,
        max_iter=max_iter,
        hdphc_iter=hdphc_iter,
        local=a.local,
        auto=a.auto,
        overwrite=bool(a.force or a.list_entry is not None),
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
    _maybe_ensure_slack_bot(a, cfg)
    cd.do_start(cfg, a.dry_run)


if __name__ == "__main__":
    main()
