#!/usr/bin/env python3
"""Post existing DC (BLC) QA PDF(s) to the run's Slack thread (resend without re-calib).

Usage: post_dc_qa_slack.py <run> {t0|drift|resi|reso} (--pdf PATH | --manifest PATH)
       [--bcout|--bcin|--all] [--kaon] [--dry-run]

Requires calibration/config/slack_local.yml. Does not re-run fit, decode, or param update.
No PDF merge: posts runXXXXX_BLC{1|2}_{TDC|DRIFT|residual|pull}_{S}.pdf as-is.
Manual resend must name --pdf or --manifest (no implicit latest).
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import dc_qa_slack, qa_session, runner  # noqa: E402


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
    p.add_argument("--pdf", type=Path, default=None,
                   help="PDF path to post (required unless --manifest)")
    p.add_argument("--manifest", type=Path, default=None,
                   help="qa_manifests/<exec_id>.json from this fit")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()
    if a.run < 1:
        runner.fail("run must be >= 1")
    if a.all and a.mode != "reso":
        runner.fail("--all is for mode reso only")
    if (a.pdf is None) == (a.manifest is None):
        runner.fail("specify exactly one of --pdf or --manifest")

    manifest = None
    pdf_paths = None
    if a.manifest is not None:
        try:
            manifest = qa_session.load_manifest(a.manifest)
        except (OSError, ValueError, TypeError) as e:
            runner.fail(f"manifest: {e}")
        ok_r, msg_r, pdf_paths = qa_session.resolve_sealed_outputs(
            manifest,
            run=a.run,
            mode=a.mode,
            particle="kaon" if a.kaon else "pion",
        )
        if not ok_r:
            runner.fail(msg_r)
    else:
        pdf_paths = [a.pdf]

    ok, message = dc_qa_slack.post_dc_qa(
        a.run,
        a.mode,
        kaon=a.kaon,
        bcin=a.bcin,
        all_chambers=a.all,
        dry_run=a.dry_run,
        manifest=manifest,
        pdf_paths=pdf_paths,
    )
    if ok:
        runner.info(message)
        sys.exit(runner.EXIT_OK)
    runner.fail(message)


if __name__ == "__main__":
    main()
