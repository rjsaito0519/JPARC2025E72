#!/usr/bin/env python3
"""Compare frozen PHC candidate ROOT files without changing live calibration."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np

from common.lib import phc_peak_study as study


#_____________________________________________________________________________
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path, help='JSON list: run, label, root, hdprm, particle (Pi/K)')
    parser.add_argument('--out', required=True, type=Path, help='New scratch output directory')
    args = parser.parse_args(argv)
    entries = json.loads(args.manifest.read_text())
    groups = {}
    for entry in entries:
        if entry['particle'] not in ('Pi', 'K'):
            parser.error('particle must be Pi or K')
        if not Path(entry['root']).is_file() or not Path(entry['hdprm']).is_file():
            parser.error('root and hdprm must be existing files (absolute paths recommended)')
        groups.setdefault((int(entry['run']), entry['particle']), []).append(entry)
    if not groups:
        parser.error('manifest is empty')
    for group in groups.values():
        if len(group) < 2 or len({entry['label'] for entry in group}) != len(group):
            parser.error('each run/particle requires at least two uniquely labeled candidates')
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out/'manifest.json').write_text(json.dumps(entries, indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    reports = []
    with PdfPages(args.out/'btof_comparison.pdf') as pdf:
        for (run, particle), group in groups.items():
            candidates = []
            provenance = []
            for entry in group:
                root = Path(entry['root'])
                before = root.stat()
                hdprm_bytes = Path(entry['hdprm']).read_bytes()
                pairs, metadata = study.skim_pairs(root, entry['hdprm'], particle)
                after = root.stat()
                if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                    raise RuntimeError('ROOT changed while reading')
                if hdprm_bytes != Path(entry['hdprm']).read_bytes():
                    raise RuntimeError('HDPRM changed while reading')
                candidates.append(pairs)
                provenance.append(dict(entry, size=before.st_size, mtime_ns=before.st_mtime_ns,
                                       hdprm_sha256=hashlib.sha256(hdprm_bytes).hexdigest(), **metadata))
            if len({row["hdprm_sha256"] for row in provenance}) != 1:
                raise RuntimeError("HDPRM differs between candidates")
            candidates = study.align_pairs(candidates)
            train = study.event_partition(candidates[0]['event'])
            values = [study.btof_values(pairs) for pairs in candidates]
            fits = []
            for sample in values:
                fits.append({name: {str(window): study.gaussian_core(sample[mask], window, event_ids=candidates[0]["event"][mask])
                                     for window in (0.3, 0.5, 0.75, 1.0, 1.5)}
                             for name, mask in [('train', train), ('validation', ~train)]})
            mixture, choices = study.channel_mixture(candidates, train)
            report = {'run': run, 'particle': particle, 'provenance': provenance, 'fits': fits,
                      'n_pairs': len(train), 'n_events': len(np.unique(candidates[0]['event'])),
                      'mixture_choices': choices,
                      'mixture_validation': study.gaussian_core(mixture[~train], event_ids=candidates[0]["event"][~train]),
                      'mixture_status': 'predicted only; native decode still required',
                      'reference': group[0]['label'], 'parameters_written': False}
            reports.append(report)
            np.savez_compressed(args.out/f'run{run:05d}_{particle}_pairs.npz',
                                event=candidates[0]['event'], train=train, mixture=mixture,
                                **{f'candidate_{index}': sample for index, sample in enumerate(values)})
            figure, axes = plt.subplots(1, 2, figsize=(11, 4))
            edges = np.linspace(-1.5, 1.5, 121)
            for axis, name, mask in zip(axes, ('train', 'validation'), (train, ~train)):
                denominator = max(1, np.sum(mask))*(edges[1]-edges[0])
                for entry, sample in zip(group, values):
                    counts = np.histogram(sample[mask & np.isfinite(sample)], edges)[0]
                    axis.stairs(counts/denominator, edges, label=entry['label'])
                counts = np.histogram(mixture[mask & np.isfinite(mixture)], edges)[0]
                axis.stairs(counts/denominator, edges, label='mixture (predicted)', color='black', linestyle='--')
                axis.set(title=f'run{run:05d} {particle}: {name}', xlabel='Corrected BTOF (ns)',
                         ylabel='Counts / common pair / ns')
                axis.axvline(0, color='gray', linewidth=0.7)
                axis.legend(fontsize=8)
                axis.grid(alpha=0.2)
            figure.tight_layout()
            pdf.savefig(figure)
            plt.close(figure)
            (args.out/'comparison.json').write_text(json.dumps(reports, indent=2)+'\n')
            print(f'run{run:05d} {particle}: {len(train)} identical hit pairs; no parameters written')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
