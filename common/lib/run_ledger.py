"""Run-summary normalization, successful-source caches, and ledger projection."""
from __future__ import annotations

import csv
import copy
import hashlib
import io
import json
import os
import tempfile
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from . import calib_snapshot, config
from .calib_state_io import write_json

MACHINE_HEADERS = ['run_key', 'run番号', '元情報', '日付', '開始', '終了', '所要時間', 'タイトル',
                   '粒子', '運動量(原文)', '元シートメモ', '較正状態', 'Hodo', 'BcIn', 'BcOut', 'D5',
                   '解析状態', '承認照合', '対象版', '前回成功', 'ローカル最終確認', '元情報最終成功', '出典', '要確認', '用途', '分類根拠', 'ローカル最終試行', '工程別最終確認']
OPTIONAL_HEADERS = ['用途', '分類根拠', 'ローカル最終試行', '工程別最終確認']
HUMAN_HEADERS = ['メモ', 'タグ']
FIELDS = ['date', 'start', 'end', 'duration', 'title', 'particle', 'momentum', 'source_memo']
LAYOUTS = {
    '2025-11': {'positions': dict(zip(FIELDS, [0, 2, 3, 4, 5, 9, 10, 19])),
                'headers': {0: 'Date', 1: 'Run Number', 2: 'Start Time', 3: 'End Time', 5: 'Title', 9: 'Particle', 10: 'Momentum (MeV/c)', 19: 'memo'}},
    '2026-04': {'positions': dict(zip(FIELDS, [0, 2, 3, 4, 6, 10, 11, 21])),
                'headers': {0: 'Date', 1: 'Run Number', 2: 'Start Time', 3: 'End Time', 4: 'Duration', 5: 'Experiment', 6: 'Title', 10: 'Particle', 12: 'Vth', 13: 'Vth', 14: 'Vth', 21: 'memo'}},
}


#_____________________________________________________________________________
def utc_now():
    return datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


#_____________________________________________________________________________
def spreadsheet_id(value):
    value = str(value).strip()
    if value.startswith('https://docs.google.com/spreadsheets/d/'):
        value = value.split('/d/', 1)[1].split('/', 1)[0]
    if not re.fullmatch(r'[A-Za-z0-9_-]{20,150}', value):
        raise ValueError('Invalid spreadsheet ID')
    return value


#_____________________________________________________________________________
def run_key(run):
    return f'run{int(run):05d}'


#_____________________________________________________________________________
def parse_run_cell(value, *, max_span=1000):
    value = value.strip()
    match = re.fullmatch(r'(\d+)(?:\s*[-–—~〜]\s*(\d+))?', value)
    if not match:
        raise ValueError(f'Unrecognized run number/range: {value!r}')
    first = int(match[1])
    last = int(match[2]) if match[2] else first
    if not 1 <= first <= last <= 10**9 or last-first+1 > max_span:
        raise ValueError('Invalid or excessive run range')
    return list(range(first, last+1))


#_____________________________________________________________________________
def parse_source(body, source):
    text = body.decode('utf-8-sig') if isinstance(body, bytes) else body
    if '<html' in text[:1000].lower() or '<!doctype' in text[:1000].lower():
        raise ValueError('Received HTML instead of CSV')
    rows = list(csv.reader(io.StringIO(text)))
    layout = LAYOUTS[source['layout']]
    header_index = next((index for index, row in enumerate(rows[:30])
                         if len(row) > 1 and row[1].strip().casefold() == 'run number'), None)
    if header_index is None:
        raise ValueError('Run header not found')
    header = rows[header_index]
    for index, expected in layout['headers'].items():
        if index >= len(header) or header[index].strip().casefold() != expected.casefold():
            raise ValueError(f"{source['name']}: header changed at column {index+1}")
    records = []
    for row_number, row in enumerate(rows[header_index+1:], header_index+2):
        if len(row) < 2 or not row[1].strip():
            continue
        runs = parse_run_cell(row[1])
        padded = row + ['']*max(0, len(header)-len(row))
        fields = {name: padded[index].strip() for name, index in layout['positions'].items()}
        aggregate = dict(fields)
        if len(runs) > 1:
            for name in ('date', 'start', 'end', 'duration'):
                fields[name] = ''
        for run in runs:
            records.append({'run': run, 'run_key': run_key(run), 'source': source['name'], 'row': row_number,
                            'run_cell': row[1].strip(), 'is_range': len(runs) > 1, 'fields': dict(fields),
                            'aggregate_fields': aggregate if len(runs) > 1 else None,
                            'raw_header': header, 'raw_row': row,
                            'url': f"https://docs.google.com/spreadsheets/d/{source['spreadsheet_id']}/edit#gid={source['gid']}&range=A{row_number}"})
    if not records:
        raise ValueError('Unexpected empty run list; previous successful data retained')
    return records


#_____________________________________________________________________________
def fetch_csv(source):
    identifier = spreadsheet_id(source['spreadsheet_id'])
    gid = str(source['gid'])
    if not gid.isdigit():
        raise ValueError('Invalid source gid')
    url = f'https://docs.google.com/spreadsheets/d/{identifier}/export?format=csv&gid={gid}'
    request = urllib.request.Request(url, headers={'User-Agent': 'E72-run-ledger/1'})
    with urllib.request.urlopen(request, timeout=25) as response:
        body = response.read(20*1024*1024+1)
        if len(body) > 20*1024*1024:
            raise ValueError('CSV exceeds size limit')
        if 'text/csv' not in response.headers.get('Content-Type', '').lower():
            raise ValueError('Source did not return text/csv')
        return body


#_____________________________________________________________________________
def refresh_sources(sources, cache_dir, *, fetcher=fetch_csv, now=None):
    now = now or utc_now()
    cache_dir = Path(cache_dir)
    results = []
    for source in sources:
        name = source['name']
        if not re.fullmatch(r'[A-Za-z0-9_-]+', name):
            raise ValueError('Invalid source name')
        path = cache_dir/'sources'/f'{name}.json'
        old = json.loads(path.read_text()) if path.exists() else None
        if old:
            prior = old['source']
            identity = lambda item: (spreadsheet_id(item['spreadsheet_id']), str(item['gid']), item['layout'])
            if identity(prior) != identity(source):
                raise ValueError('Source identity changed; use a new source name/cache')
        try:
            body = fetcher(source)
            records = parse_source(body, source)
            digest = hashlib.sha256(body).hexdigest()
            raw = cache_dir/'raw'/f'{name}_{digest}.csv'
            raw.parent.mkdir(parents=True, exist_ok=True)
            if not raw.exists() or hashlib.sha256(raw.read_bytes()).hexdigest() != digest:
                staged = None
                try:
                    with tempfile.NamedTemporaryFile(dir=raw.parent, prefix='.'+name, delete=False) as stream:
                        staged = Path(stream.name)
                        stream.write(body)
                    os.replace(staged, raw)
                finally:
                    if staged is not None and staged.exists():
                        staged.unlink()
            cached = {'source': source, 'last_success': now, 'records': records,
                      'sha256': digest, 'raw_path': str(raw)}
            write_json(path, cached)
            results.append(dict(cached, ok=True, error=''))
        except (OSError, ValueError) as error:
            cached = old or {'source': source, 'last_success': '', 'records': []}
            # HTTP errors are identified by class/status, never by a response body or credentials.
            reason = (str(error)[:300] if isinstance(error, ValueError) else type(error).__name__) + (f' HTTP {error.code}' if hasattr(error, 'code') else '')
            results.append(dict(cached, ok=False, error=reason))
    return results


#_____________________________________________________________________________
def _equivalent(value):
    value = ' '.join(value.split())
    match = re.fullmatch(r'(\d+):(\d\d):(\d\d)', value)
    if match:
        return ':'.join(str(int(part)) for part in match.groups())
    return value


#_____________________________________________________________________________
def _source_fields(observations, previous_fields=None):
    present = {name: item for name, item in observations.items() if item.get('presence', 'present') == 'present'}
    stale = {name: item for name, item in observations.items() if item.get('presence') == 'stale'}
    selected = present or stale
    if not selected and previous_fields is not None:
        return {name: previous_fields.get(name, '') for name in FIELDS}, []
    fields = {}
    conflicts = []
    for name in FIELDS:
        candidates = []
        for source, observation in (selected or observations).items():
            for record in observation['records']:
                value = record['fields'].get(name, '')
                if value:
                    candidates.append((source, value))
        unique = {_equivalent(value) for _, value in candidates}
        if len(unique) > 1:
            fields[name] = ' | '.join(f'[{source}] {value}' for source, value in dict.fromkeys(candidates))
            conflicts.append(name)
        else:
            fields[name] = candidates[0][1] if candidates else ''
    return fields, conflicts



#_____________________________________________________________________________
def snapshot_view(snapshot, old, now):
    def states(document):
        result = dict(document.get('lanes') or {})
        if document.get('parent'):
            result['parent'] = document['parent']
        return result

    prior = copy.deepcopy(old.get('state_cache') or {})
    confirmed = dict(old.get('state_checked_at') or {})
    if not prior:
        previous = old.get('snapshot') or {}
        retained = previous.get('previous_snapshot') or previous
        prior = {lane: copy.deepcopy(state) for lane, state in states(retained).items()
                 if not state.get('error') and state.get('status') != 'unknown'}
        if previous.get('consistent') and not previous.get('previous_snapshot'):
            confirmed = {lane: old.get('checked_at', '') for lane in prior}
    current = states(snapshot) if snapshot['consistent'] else {}
    keys = set(current) if snapshot['consistent'] else set(prior) | {'hodo', 'bcin', 'bcout', 'd5'}
    cached = {lane: copy.deepcopy(value) for lane, value in prior.items() if lane in keys}
    confirmed = {lane: at for lane, at in confirmed.items() if lane in keys}
    display, stale = {}, []
    for lane in sorted(keys):
        state, previous = current.get(lane), prior.get(lane)
        readable = state is not None and state.get('status') != 'unknown'
        if readable:
            value = copy.deepcopy(state)
            if state.get('error'):
                stale.append(lane)
                if previous:
                    for field in ('version', 'parameters'):
                        if field in previous:
                            value[field] = copy.deepcopy(previous[field])
                else:
                    value['version'] = '未確認'
            else:
                confirmed[lane] = now
            cached[lane] = copy.deepcopy(value)
        else:
            stale.append(lane)
            value = copy.deepcopy(previous or state or {
                'status': 'unknown', 'label': '未確認（取得失敗）', 'version': '未確認', 'approval': '未確認'})
        if lane in stale:
            if readable and value.get('status') in ('failed', 'error'):
                value['label'] = value.get('label', '失敗')+'（対象照合未確認）'
            else:
                prefix = '対象版未確認' if readable else '前回取得値' if previous else '取得失敗'
                value['label'] = f'未確認（{prefix}）: '+value.get('label', '未確認')
            value['approval'] = '今回未確認'+(f"（前回照合: {previous.get('approval', '未確認')}）" if previous else '')
            value['version'] = str(value.get('version') or '未確認')+('（前回取得値・今回未確認）' if previous else '')
        display[lane] = value
    view = dict(snapshot, parent=display.get('parent'), lanes={lane: value for lane, value in display.items() if lane != 'parent'})
    if stale:
        view['label'] = (snapshot['label']+' / 一部更新失敗（前回値あり）'
                         if snapshot['consistent'] and snapshot['label'].startswith('失敗')
                         else '未確認（更新失敗・前回値あり）' if prior else '未確認（取得失敗）')
    return cached, confirmed, view, stale


#_____________________________________________________________________________
def build_ledger(sources, previous=None, *, runs=None, snapshotter=calib_snapshot.snapshot_run, now=None, junk_runs=()):
    now = now or utc_now()
    previous = previous or {'runs': {}}
    old_runs = previous.get('runs') or {}
    grouped = {}
    for result in sources:
        by_run = {}
        for record in result['records']:
            by_run.setdefault(str(record['run']), []).append(record)
        grouped[result['source']['name']] = by_run
    all_runs = {int(run) for run in old_runs}
    all_runs.update(int(run) for group in grouped.values() for run in group)
    local = Path(config.OUTPUT_DIR)/'calib'
    all_runs.update(int(path.name[3:]) for path in local.glob('run*') if path.is_dir() and path.name[3:].isdigit())
    if runs is not None:
        all_runs = set(map(int, runs))
    output = dict(old_runs)
    partial_errors = [result['source']['name'] for result in sources if not result['ok']]
    calibration_candidates = prephysics_runs(sources)
    for run in sorted(all_runs):
        old = old_runs.get(str(run)) or {}
        observations = dict(old.get('sources') or {})
        for result in sources:
            name = result['source']['name']
            records = grouped[name].get(str(run), [])
            prior = observations.get(name)
            if records:
                observations[name] = {'records': records, 'presence': 'present' if result['ok'] else 'stale',
                                      'last_success': result['last_success'], 'last_present_at': result['last_success']}
            elif prior:
                observations[name] = dict(prior,
                    presence='absent' if result['ok'] or prior.get('presence') == 'absent' else 'stale',
                    last_present_at=prior.get('last_present_at') or prior.get('last_success'),
                    last_success=result['last_success'] if result['ok'] else prior.get('last_success', ''),
                    fetch_failed=not result['ok'])
        fields, conflicts = _source_fields(observations, old.get('fields'))
        source_value_state = ('current' if any(row['presence'] == 'present' for row in observations.values()) else
                              'stale' if any(row['presence'] == 'stale' for row in observations.values()) else
                              'absent' if observations else 'unknown')
        presence = [row['presence'] for row in observations.values()]
        status = ('取得失敗（前回取得値）' if partial_errors and source_value_state == 'stale' else
                  '一部取得失敗（前回確認では元情報不在）' if partial_errors and source_value_state == 'absent' else
                  '一部取得失敗' if partial_errors else
                  '元シートに不在' if presence and all(value == 'absent' for value in presence) else
                  '範囲行を展開' if any(record['is_range'] for row in observations.values() for record in row['records']) else
                  '取得済み' if presence else '元情報未確認')
        snapshot = snapshotter(run)
        state_cache, state_checked_at, display_snapshot, stale_lanes = snapshot_view(snapshot, old, now)
        last_success = dict(old.get('last_success') or {})
        states = dict(snapshot.get('lanes') or {})
        if snapshot.get('parent'):
            states['parent'] = snapshot['parent']
        if snapshot['consistent']:
            for lane, state in states.items():
                usable = (not state.get('error') and not state.get('missing') and not state.get('unverified')
                          and (state.get('status') != 'approved' or state.get('approval') == '照合済み'))
                if usable and state.get('status') in ('done', 'awaiting_final_approval', 'approved') and state.get('rc') in (None, 0):
                    last_success[lane] = {'at': state.get('updated_at'), 'label': state['label'],
                                          'version': state.get('version'), 'particle': state.get('particle'),
                                          'scope': state.get('scope') or []}
        if snapshot.get('previous_approval'):
            last_success['historical_approval'] = snapshot['previous_approval']
        conflict_label = '元情報競合（前回取得値）' if source_value_state == 'stale' else '元情報競合'
        issues = [f'{conflict_label}: {", ".join(conflicts)}'] if conflicts else []
        if partial_errors:
            issues.append('取得失敗: '+', '.join(partial_errors))
        if not snapshot['consistent']:
            issues.append('較正snapshot未確認（前回値はcacheに保持）')
        for lane, state in states.items():
            if '要確認' in state.get('approval', '') or state.get('error'):
                issues.append(f'{lane}: {state.get("approval") or state.get("error")}')
            if state.get('process_complete') and (state.get('missing') or state.get('unverified')):
                issues.append(f'{lane}: QA/対象の未確認あり')
        if not snapshot['consistent'] or any(state.get('error') for state in states.values()):
            retained = old.get('snapshot') or {}
            snapshot['previous_snapshot'] = retained.get('previous_snapshot') or retained
        if stale_lanes:
            issues.append('較正表示の更新失敗（前回値/確認時刻を保持）: '+', '.join(stale_lanes))
        legacy_snapshot = old.get('snapshot') or {}
        legacy_states = list((legacy_snapshot.get('lanes') or {}).values())
        if legacy_snapshot.get('parent'):
            legacy_states.append(legacy_snapshot['parent'])
        legacy_valid = (legacy_snapshot.get('consistent') and not legacy_snapshot.get('previous_snapshot')
                        and not any(state.get('error') or state.get('status') == 'unknown' for state in legacy_states))
        prior_confirmed = (old.get('checked_at') or '') if 'attempted_at' in old or legacy_valid else ''
        checked_at = now if snapshot['consistent'] and not stale_lanes else prior_confirmed
        output[str(run)] = {'run': run, 'run_key': run_key(run), 'source_status': status,
                            'sources': observations, 'fields': fields, 'conflicts': conflicts, 'source_value_state': source_value_state,
                            'snapshot': snapshot, 'display_snapshot': display_snapshot,
                            'state_cache': state_cache, 'state_checked_at': state_checked_at,
                            'stale_lanes': stale_lanes, 'last_success': last_success,
                            'checked_at': checked_at, 'attempted_at': now,
                            'issues': issues}
        usage, reason = classify_usage(run, fields['title'], calibration_candidates, junk_runs)
        output[str(run)]['usage'] = usage
        output[str(run)]['usage_reason'] = reason
    return {'schema': 1, 'updated_at': now, 'runs': output, 'source_errors': partial_errors}


#_____________________________________________________________________________
def sheet_row(record):
    snapshot = record.get('display_snapshot') or record['snapshot']
    states = dict(snapshot.get('lanes') or {})
    if snapshot.get('parent'):
        states['parent'] = snapshot['parent']
    fields = record['fields']
    presence_labels = {'present': '現在', 'stale': '前回取得', 'absent': '履歴・元表に不在'}
    source_rows = [f"{source} [{presence_labels.get(observation['presence'], '未確認')}]: {item['url']} ({item['run_cell']})"
                   for source, observation in record['sources'].items() for item in observation['records']]
    result = dict(zip(MACHINE_HEADERS, ['']*len(MACHINE_HEADERS)))
    result.update({'run_key': record['run_key'], 'run番号': str(record['run']), '元情報': record['source_status'],
                   '日付': fields['date'], '開始': fields['start'], '終了': fields['end'], '所要時間': fields['duration'],
                   'タイトル': fields['title'], '粒子': fields['particle'], '運動量(原文)': fields['momentum'],
                   '元シートメモ': fields['source_memo'], '較正状態': snapshot['label'], '解析状態': '未連携',
                   '承認照合': '; '.join(f"{lane}: {state.get('approval', '未確認')}" for lane, state in states.items()),
                   '対象版': '; '.join(f"{lane}: {state.get('version', '?')} ({','.join(state.get('scope', []))})" for lane, state in states.items()),
                   '前回成功': json.dumps(record['last_success'], ensure_ascii=False, sort_keys=True) if record['last_success'] else '記録なし',
                   'ローカル最終確認': record['checked_at'],
                   'ローカル最終試行': record.get('attempted_at', record['checked_at']),
                   '工程別最終確認': '; '.join(f"{lane}: {record.get('state_checked_at', {}).get(lane) or '未確認'}"
                       +('（今回更新失敗）' if lane in record.get('stale_lanes', []) else '') for lane in states),
                   '元情報最終成功': '; '.join(f"{name}: {row['last_success']} ({row['presence']})" for name, row in record['sources'].items()),
                   '出典': '\n'.join(source_rows), '要確認': '; '.join(record['issues']),
                   '用途': record.get('usage', '不明'), '分類根拠': record.get('usage_reason', '未分類')})
    for lane, column in [('hodo', 'Hodo'), ('bcin', 'BcIn'), ('bcout', 'BcOut'), ('d5', 'D5')]:
        state = states.get(lane)
        result[column] = (state['label'] + f" [{state.get('particle', '?')}] {state.get('step') or ''}") if state else '未連携'
    return result


#_____________________________________________________________________________
def is_physics(title):
    return bool(re.search(r'physics\s*run|e72\s*physics', title, re.I))


#_____________________________________________________________________________
def is_plain_beam(title):
    if re.search(r'physics|through|study|scan|tune|tuning|junk|dummy|off|calibration|cosmic', title, re.I):
        return False
    return bool(re.search(r'\bk?beam\s*(?:\([^)]*\)\s*)?(?:ps\s*=|[x×]\s*1)', title, re.I))


#_____________________________________________________________________________
def prephysics_runs(sources):
    candidates = set()
    for source in sources:
        if not source['ok']:
            continue
        titles = {}
        for record in source['records']:
            titles.setdefault(record['run'], set()).add(record['fields']['title'])
        previous = None
        pending = []
        for run, variants in sorted(titles.items()):
            if previous is not None and run != previous+1:
                pending = []
            previous = run
            title = next(iter(variants)) if len(variants) == 1 else ''
            if is_physics(title):
                candidates.update(pending)
                pending = []
            elif is_plain_beam(title):
                pending.append(run)
            else:
                pending = []
    return candidates


#_____________________________________________________________________________
def classify_usage(run, title, calibration_candidates=(), junk_runs=()):
    if run in junk_runs:
        return 'その他', '確認済みjunk'
    if ' | ' in title:
        return '不明', '元タイトル競合・要確認'
    if re.search(r'\bjunk\b|\bdummy\b', title, re.I):
        return 'その他', '元タイトル: junk/dummy'
    if re.search(r'cosmic', title, re.I):
        return 'その他', '元タイトル: cosmic'
    if re.search(r'\b(?:chamber|blc\w*)\s+off\b', title, re.I):
        return '不明', '検出器OFF・個別確認が必要'
    if is_physics(title):
        return 'Physics', '元タイトル: Physics Run/E72Physics'
    if re.search(r'pedestal|calibration|beam\s*(?:through|thru)|\bftof\b', title, re.I):
        return '較正', '元タイトル: pedestal/calibration/beam through/FTOF'
    if re.search(r"\bthreshold\s*scan\b|\b(?:k\s*)?beam\s*tun(?:e|ing)\b|\btrigger\s*timing\s*adjustment\b|\bdaq\s*(?:test|check)\b", title, re.I):
        return 'その他', '元タイトル: 調整/threshold scan/DAQ試験'
    if run in calibration_candidates:
        return '較正', '連続するPhysics直前の単純Beam/KBeam trigger（タイトルによる用途分類）'
    if not title.strip() or title.strip() == '0' or re.search('beam', title, re.I):
        return '不明', '用途を個別確認'
    return 'その他', '較正/Physicsを示す明確な記載なし'
