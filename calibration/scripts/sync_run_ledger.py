#!/usr/bin/env python3
"""Manually sync source run summaries and read-only calibration status to a private sheet."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import yaml

from common.lib import run_ledger
from common.lib.calib_state_io import locked_path, write_json
from common.lib.google_sheets_ledger import SheetsClient, plan_updates, read_table


#_____________________________________________________________________________
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path.home()/'.config/e72/run-ledger.yml')
    parser.add_argument('--write', action='store_true', help='Write machine-owned cells, then read back by run key')
    parser.add_argument('--local-only', action='store_true', help='Fetch sources and export locally, without Google authentication')
    parser.add_argument('--tab', help='Override destination tab; use a test tab for the first write')
    parser.add_argument('--runs', help='Comma-separated runs or ranges, e.g. 2489,2977-2978')
    args = parser.parse_args(argv)
    if args.write and args.local_only:
        parser.error('--write and --local-only are mutually exclusive')
    settings = yaml.safe_load(args.config.expanduser().read_text())
    if not isinstance(settings, dict):
        raise ValueError('Config must be a mapping')
    unknown = set(settings)-{'destination', 'credential_file', 'tab', 'cache_dir', 'sources', 'junk_runs'}
    if unknown:
        raise ValueError('Unknown config keys: '+', '.join(sorted(unknown)))
    sources = settings.get('sources')
    if not isinstance(sources, list) or not sources:
        raise ValueError('Configure read-only sources in the private YAML; source IDs are not stored in Git')
    junk_runs = settings.get('junk_runs') or []
    if not isinstance(junk_runs, list) or any(type(run) is not int or run <= 0 for run in junk_runs):
        raise ValueError('junk_runs must be a list of positive run numbers')
    if len({row['name'] for row in sources}) != len(sources):
        raise ValueError('Source names must be unique')
    destination = settings.get('destination') or ''
    title = args.tab or settings.get('tab') or 'E72_Test'
    if not title or len(title) > 100:
        raise ValueError('Invalid destination tab')
    cache = Path(settings.get('cache_dir') or '~/.cache/e72/run-ledger').expanduser()
    cache.mkdir(parents=True, exist_ok=True, mode=0o700)
    runs = None
    if args.runs:
        runs = sorted({run for value in args.runs.split(',') for run in run_ledger.parse_run_cell(value)})
    client = None
    # Every config for this user/destination serializes through the same shared-home lock.
    import hashlib
    lock_key = hashlib.sha256((run_ledger.spreadsheet_id(destination) if destination else 'local-only').encode()).hexdigest()[:24]
    lock = Path.home()/'.config/e72/ledger-locks'/lock_key
    with locked_path(cache/'transaction'), locked_path(lock):
        old_path = cache/'ledger.json'
        previous = json.loads(old_path.read_text()) if old_path.exists() else None
        refreshed = run_ledger.refresh_sources(sources, cache)
        ledger = run_ledger.build_ledger(refreshed, previous, runs=runs, junk_runs=junk_runs)
        # Source/snapshot progress remains usable even when Google authentication or publication fails.
        write_json(old_path, ledger)
        records = [run_ledger.sheet_row(row) for number, row in sorted(ledger['runs'].items(), key=lambda item: int(item[0]))
                   if runs is None or int(number) in runs]
        write_json(cache/'export.json', {'rows': records, 'source_errors': ledger['source_errors'], 'generated_at': ledger['updated_at']})
        print(f"Sources: {len(refreshed)-len(ledger['source_errors'])}/{len(refreshed)} successful; runs: {len(records)}")
        for result in refreshed:
            print(f"  {result['source']['name']}: {'OK' if result['ok'] else result['error']}; last success {result['last_success'] or '?'}")
        if args.local_only:
            print(f'Local export: {cache/"export.json"}')
            return 2 if ledger['source_errors'] else 0
        try:
            credential_file = Path(settings['credential_file']).expanduser()
            client = SheetsClient(destination, credential_file, sources=sources)
            metadata = client.metadata()
            tabs = [row['properties'] for row in metadata.get('sheets', []) if row['properties']['title'] == title]
            grid = client.read(title) if tabs else []
            observed = read_table(grid)
            if runs is None:
                extra_runs = {int(key[3:]) for key in observed['rows']}-{int(number) for number in ledger['runs']}
                if extra_runs:
                    extra = run_ledger.build_ledger(refreshed, ledger, runs=extra_runs, now=ledger['updated_at'], junk_runs=junk_runs)
                    ledger = extra
                    records = [run_ledger.sheet_row(row) for _, row in sorted(ledger['runs'].items(), key=lambda item: int(item[0]))]
                    write_json(old_path, ledger)
                    write_json(cache/'export.json', {'rows': records, 'source_errors': ledger['source_errors'], 'generated_at': ledger['updated_at']})
            plan = plan_updates(grid, records, title)
            write_json(cache/'plan.json', {key: value for key, value in plan.items() if key != 'records'})
            print(f"{'WRITE' if args.write else 'DRY RUN'} tab={title}: new {plan['inserted']}, changed {plan['changed']}; human columns preserved")
            if not args.write:
                return 2 if ledger['source_errors'] else 0
            properties = tabs[0] if tabs else client.ensure_tab(title)
            client.grow(properties, plan)
            receipt = client.publish(title, plan)
            from common.lib.run_ledger_style import apply_style
            receipt['style_changed'] = apply_style(client, title)
            receipt['requests'] = client.requests
            receipt.update(ok=True, at=run_ledger.utc_now(), tab=title)
            write_json(cache/'last_publish.json', receipt)
            write_json(cache/'last_publish_success.json', receipt)
            print(f"Verified {len(records)} runs by key; API calls: {client.requests}")
            return 2 if ledger['source_errors'] else 0
        except Exception as error:
            if args.write:
                write_json(cache/'last_publish.json', {'ok': False, 'at': run_ledger.utc_now(),
                           'error': type(error).__name__, 'message': str(error)[:300] if isinstance(error, RuntimeError) else '', 'tab': title})
            raise


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, KeyError) as error:
        print(f'[Error] {error}', file=sys.stderr)
        raise SystemExit(1)
