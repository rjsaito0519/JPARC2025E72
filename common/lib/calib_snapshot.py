"""Read-only calibration snapshots for external ledgers; never reconcile state."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from . import calib_approval, config

LABELS = {
    'running': '実行中', 'held': '保留', 'failed': '失敗', 'error': '失敗',
    'stopped': '停止', 'done': '対象工程の処理完走（承認は別）',
    'awaiting_final_approval': '最終確認待ち', 'finished': '最終確認待ち',
    'approved': '対象範囲の承認済み', 'pending': '待機中', 'queued': '待機中',
}


#_____________________________________________________________________________
def _tokens(folder):
    paths = sorted(path for path in folder.rglob('*')
                   if path.is_file() and (path.suffix == '.json' or path.name in ('approval_journal.jsonl', 'channel_flags.jsonl')))
    return [(str(path), path.stat().st_ino, path.stat().st_size, path.stat().st_mtime_ns) for path in paths]


#_____________________________________________________________________________
def _entry(run, path, subject_builder):
    try:
        text = path.read_text(encoding='utf-8')
        state = json.loads(text)
        if not isinstance(state, dict) or not isinstance(state.get('status'), str):
            raise ValueError('missing state status')
    except (OSError, ValueError) as error:
        return {'status': 'unknown', 'label': '未確認（状態を読めない）', 'error': type(error).__name__}
    status = calib_approval.effective_status(state)
    result = {'status': status, 'raw_status': state['status'], 'label': LABELS.get(status, f'未確認 ({status})'),
              'step': state.get('step'), 'updated_at': state.get('updated_at'), 'rc': state.get('rc'),
              'process_complete': calib_approval.is_process_complete(state) or status == 'done',
              'particle': calib_approval.particle_suffix(state)[0],
              'scope': sorted(calib_approval.plan_step_names(state)),
              'version': 'state:'+hashlib.sha256(text.encode()).hexdigest()[:16],
              'approval': '未承認', 'unverified': state.get('unverified') or [], 'missing': []}
    if state.get('plan'):
        try:
            first = subject_builder(run, state)
            second = subject_builder(run, state)
            first_sha = calib_approval.subject_sha(first)
            if first_sha != calib_approval.subject_sha(second):
                raise RuntimeError('calibration subject changed while reading')
            result['version'] = first_sha
            result['missing'] = calib_approval.subject_missing(first)
            result['parameters'] = [{key: row.get(key) for key in ('key', 'path', 'sha', 'missing')} for row in first.get('params', [])]
            if status == 'approved':
                approval = state.get('approval') or {}
                pinned = approval.get('subject_sha') if isinstance(approval, dict) else None
                if not pinned or pinned != first_sha:
                    result['approval'] = '版不一致・要確認'
                    result['label'] = '承認記録あり（版不一致）'
                elif (result['missing'] or first.get('unverified') or first.get('open_channel_flags')
                      or not first.get('families')):
                    result['approval'] = '検証不能・要確認'
                    result['label'] = '承認記録あり（検証不能）'
                else:
                    result['approval'] = '照合済み'
        except RuntimeError:
            raise
        except (OSError, ValueError, KeyError, TypeError) as error:
            result['approval'] = '検証不能・要確認'
            result['error'] = f'subject:{type(error).__name__}'
            if status == 'approved':
                result['label'] = '承認記録あり（検証不能）'
    elif status == 'approved':
        result['approval'] = '対象不明・要確認'
        result['label'] = '承認記録あり（対象不明）'
    if result['process_complete'] and result['rc'] not in (None, 0):
        result['label'] = '状態不整合・要確認'
        result['approval'] = '検証不能・要確認'
    return result


#_____________________________________________________________________________
def snapshot_run(run, *, folder=None, subject_builder=None, retries=3):
    """Double-collect state and subjects; no status CLI, locks, or state writes."""
    run = int(run)
    folder = Path(folder) if folder is not None else Path(config.OUTPUT_DIR)/'calib'/f'run{run:05d}'
    subject_builder = subject_builder or calib_approval.build_subject
    for _ in range(retries):
        try:
            before = _tokens(folder)
            parent_path = folder/'calib_dispatch_state.json'
            parent = _entry(run, parent_path, subject_builder) if parent_path.is_file() else None
            lanes = {path.parent.name: _entry(run, path, subject_builder)
                     for path in sorted((folder/'lanes').glob('*/calib_dispatch_state.json')) if path.is_file()}
            # A coordinator map without its lane file is evidence of a lane, not a verified result.
            if parent_path.is_file():
                try:
                    raw_parent = json.loads(parent_path.read_text())
                    for name in (raw_parent.get('lanes') or {}):
                        lanes.setdefault(name, {'status': 'unknown', 'label': '未確認（lane記録なし）', 'approval': '検証不能'})
                except (ValueError, AttributeError):
                    pass
            previous_approval = None
            journal = folder/'approval_journal.jsonl'
            if journal.is_file():
                for line in journal.read_text().splitlines():
                    try:
                        event = json.loads(line)
                    except ValueError:
                        continue
                    if isinstance(event, dict) and event.get('event') == 'approve':
                        previous_approval = {key: event.get(key) for key in ('at', 'lane', 'subject_sha')}
            if before != _tokens(folder):
                continue
            current = list(lanes.values()) + ([parent] if parent else [])
            statuses = {row.get('status') for row in current}
            labels = []
            if statuses & {'failed', 'error'}:
                labels.append('失敗あり')
            if statuses & {'running', 'pending', 'queued'}:
                labels.append('実行中/待機中')
            if statuses & {'held', 'stopped'}:
                labels.append('保留/停止あり')
            if 'unknown' in statuses:
                labels.append('未確認あり')
            label = ' / '.join(labels) if labels else (parent['label'] if parent else ('工程別の記録あり' if lanes else '未連携'))
            return {'run': run, 'consistent': True, 'parent': parent, 'lanes': lanes,
                    'previous_approval': previous_approval, 'label': label,
                    'analysis': '未連携'}
        except (OSError, RuntimeError):
            continue
    return {'run': run, 'consistent': False, 'parent': None, 'lanes': {},
            'label': '未確認（同時更新・読取失敗）', 'analysis': '未連携'}
