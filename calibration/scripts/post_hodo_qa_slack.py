#!/usr/bin/env python3
"""Post a merged Hodo QA PDF to the run's Slack thread (resend without re-calib).

Usage: post_hodo_qa_slack.py <run> {hdprm|t0|hdphc} (--pdf PATH | --manifest PATH)
       [--kaon] [--ftof] [--dry-run]

Requires calibration/config/slack_local.yml (copy from slack_local.example.yml).
Does not re-run fit, decode, or param update. Manual resend must name the PDF
or a qa_manifest explicitly (no implicit latest).
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import hodo_qa_pdf, qa_session, runner, slack_notify  # noqa: E402


#_____________________________________________________________________________
def main():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int)
    p.add_argument("mode", choices=hodo_qa_pdf.MERGE_MODES)
    p.add_argument("--kaon", action="store_true")
    p.add_argument("--ftof", action="store_true",
                   help="hdprm/hdphc: FTOF merged PDF name (ignored for t0)")
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
            mode=a.mode,
            particle="kaon" if a.kaon else "pion",
        )
        if not ok_r:
            runner.fail(msg_r)

    ok, message = slack_notify.post_hodo_qa(
        a.run,
        a.mode,
        kaon=a.kaon,
        ftof=a.ftof,
        pdf_path=pdf_path,
        dry_run=a.dry_run,
        manifest=manifest,
    )
    if ok:
        runner.info(message)
        sys.exit(runner.EXIT_OK)
    runner.fail(message)


if __name__ == "__main__":
    main()
