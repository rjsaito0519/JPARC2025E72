#!/usr/bin/env python3
"""Post existing D5 QA PDF to the run's Slack thread (resend without re-calib).

Usage: post_d5_qa_slack.py <run> {pos_review|plane_resi_qa} [--kaon] [--dry-run]

Requires calibration/config/slack_local.yml. Does not re-run fit, decode, or param update.
No PDF merge: posts d5_*_runXXXXX_{S}.pdf as-is.
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import d5_qa_slack, runner  # noqa: E402


#_____________________________________________________________________________
def main():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int)
    p.add_argument("kind", choices=d5_qa_slack.KINDS)
    p.add_argument("--kaon", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()
    if a.run < 1:
        runner.fail("run must be >= 1")

    ok, message = d5_qa_slack.post_d5_qa(
        a.run, a.kind, kaon=a.kaon, dry_run=a.dry_run,
    )
    if ok:
        runner.info(message)
        sys.exit(runner.EXIT_OK)
    runner.fail(message)


if __name__ == "__main__":
    main()
