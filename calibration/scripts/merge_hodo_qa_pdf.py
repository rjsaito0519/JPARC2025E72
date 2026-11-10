#!/usr/bin/env python3
"""Merge per-detector Hodo QA PDFs into one book for a stage.

Usage: merge_hodo_qa_pdf.py <run> {hdprm|t0|hdphc} [--kaon] [--ftof] [--dry-run]

Uses pdfunite (page concat). Does not glob; detector order is config.HODO_DETECTORS.
Originals are kept. Output: results/img/runXXXXX/runXXXXX_hodo_<mode>[_ftof]_<S>.pdf
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from common.lib import hodo_qa_pdf, runner  # noqa: E402


#_____________________________________________________________________________
def main():
    p = runner.Parser(__doc__)
    p.add_argument("run", type=int)
    p.add_argument("mode", choices=hodo_qa_pdf.MERGE_MODES)
    p.add_argument("--kaon", action="store_true")
    p.add_argument("--ftof", action="store_true",
                   help="hdprm/hdphc: FTOF detector set (ignored for t0)")
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()
    if a.run < 1:
        runner.fail("run must be >= 1")

    ok, message, out = hodo_qa_pdf.merge(
        a.run, a.mode, kaon=a.kaon, ftof=a.ftof, dry_run=a.dry_run,
    )
    if ok:
        runner.info(message)
        if out:
            print(out)
        sys.exit(runner.EXIT_OK)
    runner.fail(message)


if __name__ == "__main__":
    main()
