#!/usr/bin/env python3
"""Daily calibration launcher: ensure session and start (or resume if stopped).

Usage:
  ./calib.py <run> [--pion|--kaon] [--hodo] [--bcout] [--bcin] [--d5]
             [--local] [--dry] [--force]

State ops: ./ops.py status|stop|… <run>  (or Slack /calib …).
Prefer this over calib_session.py / calib_dispatch.py.
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import calib_dispatch as cd  # noqa: E402
import calib_session as cs  # noqa: E402
from common.lib import calib_approval, runner  # noqa: E402


#_____________________________________________________________________________
def parse_args():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int, help="run number")
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
                   help="rewrite session YAML (not while running)")
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
            "Hodo: if TDC windows still match the USER template, the chain "
            "stops (EXIT_CHECK). Set windows (tdcwin / set_tdc_window), then "
            f"./ops.py decide {run} refit && ./ops.py resume {run} "
            f"(or ./calib.py {run})."
        )


#_____________________________________________________________________________
def main():
    a = parse_args()
    if a.run < 1:
        runner.fail("run must be >= 1")
    calib_particle = cs.resolve_calib_particle(pion=a.pion, kaon=a.kaon)
    steps, target = cs.steps_from_flags(
        hodo=a.hodo, bcout=a.bcout, bcin=a.bcin, d5=a.d5,
    )
    announce_chain(a.run, calib_particle, steps)

    st = cd.read_state(a.run)
    if st is not None and a.force:
        status = calib_approval.effective_status(st)
        if st.get("status") == "running" or status == "running":
            runner.fail(
                f"run{a.run:05d}: cannot --force rewrite session while running; "
                f"./ops.py stop {a.run} first"
            )

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

    if st is None:
        cd.do_start(cfg, a.dry_run)

    status = calib_approval.effective_status(st)
    if st.get("status") == "stopped" or status == "stopped":
        cd.do_resume(cfg, a.dry_run)

    if st.get("status") == "running" or status == "running":
        runner.fail(
            f"run{a.run:05d}: already running; "
            f"use ./ops.py status|stop {a.run} or Slack /calib"
        )
    if st.get("status") == "held" or status == "held":
        runner.fail(
            f"run{a.run:05d}: held; decide via Slack /calib decide "
            f"or ./ops.py decide {a.run} refit|continue"
        )
    if calib_approval.is_process_complete(st):
        runner.fail(
            f"run{a.run:05d}: process complete ({status}); "
            f"approve/revoke via Slack or ./ops.py"
        )
    runner.fail(
        f"run{a.run:05d}: status={st.get('status')!r} ({status}); "
        f"use ./ops.py status {a.run}"
    )


if __name__ == "__main__":
    main()
