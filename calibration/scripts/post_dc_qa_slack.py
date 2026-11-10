#!/usr/bin/env python3
"""Post existing DC (BLC) QA PDF(s) to the run's Slack thread (resend without re-calib).

Usage: post_dc_qa_slack.py <run> {t0|drift|resi|reso} [--bcout|--bcin|--all] [--kaon] [--dry-run]

Requires calibration/config/slack_local.yml. Does not re-run fit, decode, or param update.
No PDF merge: posts runXXXXX_BLC{1|2}_{TDC|DRIFT|residual|pull}_{S}.pdf as-is.
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import dc_qa_slack, runner  # noqa: E402


#_____________________________________________________________________________
def main():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int)
    p.add_argument("mode", choices=dc_qa_slack.MODES)
    group = p.add_mutually_exclusive_group()
    group.add_argument("--bcout", action="store_true", help="BcOut / BLC2 (default)")
    group.add_argument("--bcin", action="store_true", help="BcIn / BLC1")
    group.add_argument("--all", action="store_true",
                       help="reso only: post BLC1 and BLC2 pull PDFs")
    p.add_argument("--kaon", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()
    if a.run < 1:
        runner.fail("run must be >= 1")
    if a.all and a.mode != "reso":
        runner.fail("--all is for mode reso only")

    ok, message = dc_qa_slack.post_dc_qa(
        a.run,
        a.mode,
        kaon=a.kaon,
        bcin=a.bcin,
        all_chambers=a.all,
        dry_run=a.dry_run,
    )
    if ok:
        runner.info(message)
        sys.exit(runner.EXIT_OK)
    runner.fail(message)


if __name__ == "__main__":
    main()
