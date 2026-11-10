#!/usr/bin/env python3
"""Post existing D5 QA PDF to the run's Slack thread (resend without re-calib).

Usage: post_d5_qa_slack.py <run> {pos_review|plane_resi_qa} (--pdf PATH | --manifest PATH)
       [--kaon] [--dry-run]

Requires calibration/config/slack_local.yml. Does not re-run fit, decode, or param update.
No PDF merge: posts d5_*_runXXXXX_{S}.pdf as-is.
Manual resend must name --pdf or --manifest (no implicit latest).
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import d5_qa_slack, qa_session, runner  # noqa: E402


#_____________________________________________________________________________
def main():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int)
    p.add_argument("kind", choices=d5_qa_slack.KINDS)
    p.add_argument("--kaon", action="store_true")
    p.add_argument("--pdf", type=Path, default=None,
                   help="PDF path to post (required unless --manifest)")
    p.add_argument("--manifest", type=Path, default=None,
                   help="qa_manifests/<exec_id>.json from this fit")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()
    if a.run < 1:
        runner.fail("run must be >= 1")
    if (a.pdf is None) == (a.manifest is None):
        runner.fail("specify exactly one of --pdf or --manifest")

    manifest = None
    pdf_path = a.pdf
    if a.manifest is not None:
        try:
            manifest = qa_session.load_manifest(a.manifest)
        except (OSError, ValueError, TypeError) as e:
            runner.fail(f"manifest: {e}")
        ok_r, msg_r, pdf_path = qa_session.resolve_sealed_output(
            manifest,
            run=a.run,
            mode=a.kind,
            particle="kaon" if a.kaon else "pion",
        )
        if not ok_r:
            runner.fail(msg_r)

    ok, message = d5_qa_slack.post_d5_qa(
        a.run, a.kind, kaon=a.kaon, dry_run=a.dry_run,
        manifest=manifest, pdf_path=pdf_path,
    )
    if ok:
        runner.info(message)
        sys.exit(runner.EXIT_OK)
    runner.fail(message)


if __name__ == "__main__":
    main()
