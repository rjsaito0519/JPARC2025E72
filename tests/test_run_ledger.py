"""Run-key ownership, failed-source retention, and read-only calibration snapshots."""
import copy
import csv
import io
import json
import re
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from common.lib import calib_snapshot, run_ledger
from common.lib.google_sheets_ledger import HEADERS, SheetsClient, plan_updates, read_table


TEST_SOURCES = [
    {'name': layout, 'spreadsheet_id': 'test_source_' + layout.replace('-', '_') + '_identifier',
     'gid': index, 'layout': layout}
    for index, layout in enumerate(('2025-11', '2026-04'))
]

def source_csv(layout='2025-11', run='1916', title='beam', duration='0:23:49'):
    header = ['']*30
    for index, value in run_ledger.LAYOUTS[layout]['headers'].items():
        header[index] = value
    row = ['']*30
    row[1] = run
    positions = run_ledger.LAYOUTS[layout]['positions']
    for name, value in {'date': '2025/11/10', 'start': '12:08:00', 'end': '12:31:49',
                        'duration': duration, 'title': title, 'particle': 'pion', 'momentum': '-735'}.items():
        row[positions[name]] = value
    stream = io.StringIO()
    csv.writer(stream).writerows([['Total Run Time'], header, row])
    return stream.getvalue().encode()


def empty_snapshot(run):
    return {'consistent': True, 'label': '未連携', 'parent': None, 'lanes': {}, 'analysis': '未連携'}


class SourceTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.sources = copy.deepcopy(TEST_SOURCES)
        self.patch = patch.object(run_ledger.config, 'OUTPUT_DIR', self.root/'output')
        self.patch.start()
        self.addCleanup(self.patch.stop)

    def refresh(self, fetcher, sources=None, now='first'):
        return run_ledger.refresh_sources(sources or self.sources, self.root/'cache', fetcher=fetcher, now=now)

    def test_both_layouts_use_correct_title_and_momentum_columns(self):
        for source in self.sources:
            rows = run_ledger.parse_source(source_csv(source['layout']), source)
            self.assertEqual('beam', rows[0]['fields']['title'])
            self.assertEqual('-735', rows[0]['fields']['momentum'])
            self.assertEqual('1916', rows[0]['raw_row'][1])

    def test_range_does_not_copy_aggregate_time_to_each_run(self):
        rows = run_ledger.parse_source(source_csv(run='2189 - 2201'), self.sources[0])
        self.assertEqual(13, len(rows))
        self.assertEqual(2201, rows[-1]['run'])
        for row in rows:
            for key in ('date', 'start', 'end', 'duration'):
                self.assertEqual('', row['fields'][key])
            self.assertEqual('0:23:49', row['aggregate_fields']['duration'])

    def test_html_empty_and_changed_header_are_rejected(self):
        for body in (b'<html>login</html>', b'Run Number,Title\n', source_csv().replace(b'Title', b'Moved title')):
            with self.assertRaises(ValueError):
                run_ledger.parse_source(body, self.sources[0])

    def test_invalid_run_ranges_are_rejected(self):
        for value in ('2-1', '0', '1-100000', 'abc', '1; rm x', '1.5'):
            with self.assertRaises(ValueError):
                run_ledger.parse_run_cell(value)

    def test_source_failure_retains_last_success_and_marks_stale(self):
        first = self.refresh(lambda source: source_csv(source['layout']))
        ledger = run_ledger.build_ledger(first, snapshotter=empty_snapshot, now='first')
        second = self.refresh(lambda source: b'<html>login</html>', now='second')
        self.assertTrue(all(not row['ok'] for row in second))
        self.assertTrue(all(row['last_success'] == 'first' for row in second))
        current = run_ledger.build_ledger(second, ledger, snapshotter=empty_snapshot, now='second')
        self.assertEqual('beam', current['runs']['1916']['fields']['title'])
        self.assertEqual('一部取得失敗', current['runs']['1916']['source_status'])

    def test_one_failed_source_does_not_block_other_source_updates(self):
        self.refresh(lambda source: source_csv(source['layout']))
        current = self.refresh(lambda source: b'<html>no</html>' if source['name'] == '2025-11'
                               else source_csv(source['layout'], run='2000'), now='second')
        ledger = run_ledger.build_ledger(current, snapshotter=empty_snapshot)
        self.assertIn('1916', ledger['runs'])
        self.assertIn('2000', ledger['runs'])

    def test_deleted_run_is_retained_with_source_absence(self):
        first = self.refresh(lambda source: source_csv(source['layout']))
        old = run_ledger.build_ledger(first, snapshotter=empty_snapshot)
        current = self.refresh(lambda source: source_csv(source['layout'], run='2000'))
        ledger = run_ledger.build_ledger(current, old, snapshotter=empty_snapshot)
        self.assertEqual('元シートに不在', ledger['runs']['1916']['source_status'])
        self.assertEqual('beam', ledger['runs']['1916']['fields']['title'])

    def test_conflicting_sources_keep_both_values_and_request_review(self):
        current = self.refresh(lambda source: source_csv(source['layout'], title=source['name']))
        ledger = run_ledger.build_ledger(current, snapshotter=empty_snapshot)
        row = ledger['runs']['1916']
        self.assertIn('2025-11', row['fields']['title'])
        self.assertIn('2026-04', row['fields']['title'])
        self.assertIn('title', row['conflicts'])

    def test_equivalent_durations_are_not_conflicts(self):
        current = self.refresh(lambda source: source_csv(source['layout'], duration='00:23:49' if source['name'] == '2026-04' else '0:23:49'))
        row = run_ledger.build_ledger(current, snapshotter=empty_snapshot)['runs']['1916']
        self.assertNotIn('duration', row['conflicts'])

    def test_current_failure_keeps_previous_process_success(self):
        sources = self.refresh(lambda source: source_csv(source['layout']))
        state = {'consistent': True, 'label': '処理完走', 'parent': {'status': 'done', 'label': '処理完走', 'rc': 0,
                 'updated_at': 'yesterday', 'version': 'old', 'particle': 'pion', 'scope': ['hodo_t0']}, 'lanes': {}}
        old = run_ledger.build_ledger(sources, snapshotter=lambda run: copy.deepcopy(state))
        state['parent'].update(status='failed', label='失敗', rc=1, particle='kaon', scope=['hodo_hdphc'])
        state['label'] = '失敗'
        new = run_ledger.build_ledger(sources, old, snapshotter=lambda run: copy.deepcopy(state))
        self.assertEqual('失敗', new['runs']['1916']['snapshot']['label'])
        self.assertEqual('old', new['runs']['1916']['last_success']['parent']['version'])
        self.assertEqual('pion', new['runs']['1916']['last_success']['parent']['particle'])
        self.assertEqual(['hodo_t0'], new['runs']['1916']['last_success']['parent']['scope'])

    def test_changed_source_identity_does_not_reuse_another_sheet_cache(self):
        self.refresh(lambda source: source_csv(source['layout']))
        changed = copy.deepcopy(self.sources)
        changed[0]['gid'] += 1
        with self.assertRaisesRegex(ValueError, 'Source identity changed'):
            self.refresh(lambda source: source_csv(source['layout']), sources=changed)

    def test_stale_approval_does_not_replace_verified_previous_success(self):
        sources = self.refresh(lambda source: source_csv(source['layout']))
        state = {'consistent': True, 'label': '承認済み', 'parent': {'status': 'approved', 'label': '承認済み',
                 'approval': '照合済み', 'rc': 0, 'updated_at': 'first', 'version': 'verified'}, 'lanes': {}}
        old = run_ledger.build_ledger(sources, snapshotter=lambda run: copy.deepcopy(state))
        state['parent'].update(approval='版不一致・要確認', version='unreviewed')
        new = run_ledger.build_ledger(sources, old, snapshotter=lambda run: copy.deepcopy(state))
        self.assertEqual('verified', new['runs']['1916']['last_success']['parent']['version'])

    def test_unlabeled_momentum_column_is_guarded_by_neighbor_headers(self):
        source = self.sources[1]
        body = source_csv(source['layout']).replace(b'Vth,Vth,Vth', b',Vth,Vth,Vth')
        with self.assertRaises(ValueError):
            run_ledger.parse_source(body, source)

    def test_incomplete_raw_cache_file_is_repaired_before_reuse(self):
        first = self.refresh(lambda source: source_csv(source['layout']))
        raw = Path(first[0]['raw_path'])
        raw.write_bytes(b'partial file')
        second = self.refresh(lambda source: source_csv(source['layout']))
        self.assertTrue(second[0]['ok'])
        self.assertEqual(source_csv(), raw.read_bytes())


class SnapshotTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def write(self, path, value):
        target = self.root/path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(value))
        return target

    def subject(self, run, state):
        return {'families': ['hodo'], 'missing': [], 'params': [], 'unverified': [], 'open_channel_flags': []}

    def test_stale_approval_is_reported_without_modifying_state(self):
        path = self.write('calib_dispatch_state.json', {'status': 'approved', 'plan': {'steps': [{'name': 'hodo_t0'}]}, 'approval': {'subject_sha': 'old'}})
        before = path.read_bytes(), path.stat().st_mtime_ns
        result = calib_snapshot.snapshot_run(1, folder=self.root, subject_builder=self.subject)
        self.assertIn('版不一致', result['parent']['approval'])
        self.assertEqual(before, (path.read_bytes(), path.stat().st_mtime_ns))

    def test_matching_approval_is_verified(self):
        digest = calib_snapshot.calib_approval.subject_sha(self.subject(1, {}))
        self.write('calib_dispatch_state.json', {'status': 'approved', 'plan': {'steps': [{'name': 'hodo_t0'}]}, 'approval': {'subject_sha': digest}})
        result = calib_snapshot.snapshot_run(1, folder=self.root, subject_builder=self.subject)
        self.assertEqual('照合済み', result['parent']['approval'])

    def test_lane_state_wins_over_stale_coordinator(self):
        self.write('calib_dispatch_state.json', {'status': 'running', 'lanes': {'hodo': 'done', 'd5': 'done'}})
        self.write('lanes/hodo/calib_dispatch_state.json', {'status': 'held'})
        result = calib_snapshot.snapshot_run(1, folder=self.root, subject_builder=self.subject)
        self.assertEqual('held', result['lanes']['hodo']['status'])
        self.assertEqual('unknown', result['lanes']['d5']['status'])

    def test_changed_subject_is_not_exported_as_consistent(self):
        self.write('calib_dispatch_state.json', {'status': 'approved', 'plan': {'steps': [{'name': 'hodo_t0'}]}})
        counter = [0]
        def changing(run, state):
            counter[0] += 1
            return {'version': counter[0]}
        self.assertFalse(calib_snapshot.snapshot_run(1, folder=self.root, subject_builder=changing, retries=2)['consistent'])

    def test_partial_completion_does_not_become_global_approval(self):
        self.write('lanes/hodo/calib_dispatch_state.json', {'status': 'done', 'rc': 0})
        self.write('lanes/d5/calib_dispatch_state.json', {'status': 'failed', 'rc': 1})
        result = calib_snapshot.snapshot_run(1, folder=self.root, subject_builder=self.subject)
        self.assertNotIn('承認済み', result['label'])
        self.assertEqual('failed', result['lanes']['d5']['status'])

    def test_older_parent_approval_does_not_hide_current_running_lane(self):
        digest = calib_snapshot.calib_approval.subject_sha(self.subject(1, {}))
        self.write('calib_dispatch_state.json', {'status': 'approved', 'plan': {'steps': [{'name': 'hodo_t0'}]}, 'approval': {'subject_sha': digest}})
        self.write('lanes/hodo/calib_dispatch_state.json', {'status': 'running'})
        result = calib_snapshot.snapshot_run(1, folder=self.root, subject_builder=self.subject)
        self.assertIn('実行中', result['label'])
        self.assertNotIn('承認済み', result['label'])
        self.assertEqual('照合済み', result['parent']['approval'])


class FakeSheets(SheetsClient):
    def __init__(self, grid):
        self.grid = copy.deepcopy(grid)
        self.requests = 0
        self.posts = 0
        self.lost_response = False
        self.reject_write = False
        self.before_write = None

    def read(self, title):
        self.requests += 1
        return copy.deepcopy(self.grid)

    def request(self, method, suffix, body=None):
        self.requests += 1
        self.posts += 1
        if self.reject_write:
            raise RuntimeError('offline')
        if self.before_write:
            self.before_write(self)
        for update in body['data']:
            region = update['range'].rsplit('!', 1)[1]
            match = re.fullmatch(r'([A-Z]+)(\d+):([A-Z]+)(\d+)', region)
            column = 0
            for char in match[1]:
                column = column*26+ord(char)-64
            column -= 1
            row = int(match[2])-1
            while len(self.grid) <= row:
                self.grid.append([])
            values = update['values'][0]
            self.grid[row] += ['']*max(0, column+len(values)-len(self.grid[row]))
            self.grid[row][column:column+len(values)] = values
        if self.lost_response:
            raise RuntimeError('response lost after commit')
        return {}


class SheetOwnershipTests(unittest.TestCase):
    def record(self, run, status='実行中'):
        row = dict.fromkeys(run_ledger.MACHINE_HEADERS, '')
        row.update(run_key=run_ledger.run_key(run), **{'run番号': str(run), 'Hodo': status})
        return row

    def grid(self):
        rows = [HEADERS+['自由列']]
        for run, note in [(2, 'run2 memo'), (1, 'run1 memo')]:
            record = self.record(run)
            rows.append([note if header == 'メモ' else 'tag' if header == 'タグ' else 'personal' if header == '自由列'
                         else record.get(header, '') for header in rows[0]])
        return rows

    def test_reordered_rows_keep_memos_tags_and_unknown_columns(self):
        sheet = FakeSheets(self.grid())
        plan = plan_updates(sheet.grid, [self.record(1, '完走'), self.record(2, '保留')], 'Test')
        sheet.publish('Test', plan)
        table = read_table(sheet.grid)
        self.assertEqual('run1 memo', table['rows']['run00001']['values'][HEADERS.index('メモ')])
        self.assertEqual('run2 memo', table['rows']['run00002']['values'][HEADERS.index('メモ')])
        self.assertTrue(all(row[-1] == 'personal' for row in sheet.grid[1:]))

    def test_direct_note_edit_after_plan_is_preserved(self):
        sheet = FakeSheets(self.grid())
        plan = plan_updates(sheet.grid, [self.record(1, '完走')], 'Test')
        sheet.grid[2][HEADERS.index('メモ')] = 'edited while syncing'
        sheet.publish('Test', plan)
        self.assertEqual('edited while syncing', sheet.grid[2][HEADERS.index('メモ')])

    def test_structural_change_before_write_aborts_without_posts(self):
        sheet = FakeSheets(self.grid())
        plan = plan_updates(sheet.grid, [self.record(1, '完走')], 'Test')
        sheet.grid[1:] = list(reversed(sheet.grid[1:]))
        with self.assertRaisesRegex(RuntimeError, 'changed since planning'):
            sheet.publish('Test', plan)
        self.assertEqual(0, sheet.posts)

    def test_duplicate_destination_keys_stop_before_write(self):
        grid = self.grid()
        grid[2][0] = grid[1][0]
        with self.assertRaisesRegex(ValueError, 'Duplicate run_key'):
            plan_updates(grid, [self.record(1)], 'Test')

    def test_note_only_row_is_reserved_when_appending(self):
        grid = [HEADERS, ['', '', 'unkeyed note']]
        plan = plan_updates(grid, [self.record(1)], 'Test')
        self.assertTrue(all('3' in item['range'] for item in plan['updates']))

    def test_repeat_sync_is_idempotent_even_after_lost_write_response(self):
        sheet = FakeSheets([])
        records = [self.record(1), self.record(2)]
        sheet.lost_response = True
        self.assertTrue(sheet.publish('Test', plan_updates(sheet.grid, records, 'Test'))['verified'])
        next_plan = plan_updates(sheet.grid, records, 'Test')
        self.assertEqual([], next_plan['updates'])
        self.assertEqual(2, len(read_table(sheet.grid)['rows']))

    def test_failed_write_is_not_blindly_retried(self):
        sheet = FakeSheets(self.grid())
        sheet.reject_write = True
        with self.assertRaises(RuntimeError):
            sheet.publish('Test', plan_updates(sheet.grid, [self.record(1, 'new')], 'Test'))
        self.assertEqual(1, sheet.posts)

    def test_mid_write_row_sort_is_detected_not_claimed_successful(self):
        sheet = FakeSheets(self.grid())
        sheet.before_write = lambda current: current.grid.__setitem__(slice(1, None), list(reversed(current.grid[1:])))
        with self.assertRaises(RuntimeError):
            sheet.publish('Test', plan_updates(sheet.grid, [self.record(1, 'new')], 'Test'))

    def test_source_sheet_cannot_be_a_destination(self):
        for source in TEST_SOURCES:
            with self.assertRaisesRegex(ValueError, 'original source'):
                SheetsClient(source['spreadsheet_id'], '/not/read', sources=TEST_SOURCES)

    def test_partial_batch_retry_resolves_existing_keys_without_duplicate_append(self):
        class FailSecond(FakeSheets):
            def request(self, method, suffix, body=None):
                if self.posts == 1:
                    self.posts += 1
                    raise RuntimeError('second batch lost')
                return super().request(method, suffix, body)
        sheet = FailSecond([])
        records = [dict(self.record(run), **{'タイトル': 'x'*40000}) for run in range(1, 21)]
        with self.assertRaises(RuntimeError):
            sheet.publish('Test', plan_updates(sheet.grid, records, 'Test'))
        self.assertTrue(sheet.publish('Test', plan_updates(sheet.grid, records, 'Test'))['verified'])
        self.assertEqual(20, len(read_table(sheet.grid)['rows']))


class UsageAndStyleTests(unittest.TestCase):
    def test_usage_junk_and_cosmic_take_priority(self):
        for title in ('HTOF cosmic', 'Physics Run -> junk', 'dummy'):
            self.assertEqual('その他', run_ledger.classify_usage(1, title)[0])
        self.assertEqual('その他', run_ledger.classify_usage(3828, 'BEAM x1/90 chamber off', junk_runs=[3828])[0])
        self.assertEqual('不明', run_ledger.classify_usage(3828, 'BEAM x1/90 chamber off')[0])

    def test_physics_and_explicit_calibration(self):
        for title in ('Physics Run : Mp2 + Beam x1/3000', 'E72Physics : KBEAM x1/600'):
            self.assertEqual('Physics', run_ledger.classify_usage(1, title)[0])
        for title in ('pedestal run', 'HTOF calibration run', 'beam through proton', 'FTOF'):
            self.assertEqual('較正', run_ledger.classify_usage(1, title)[0])
        self.assertEqual('不明', run_ledger.classify_usage(1, '[a] Physics Run | [b] pedestal run')[0])
        self.assertEqual('不明', run_ledger.classify_usage(1, 'beam check')[0])

    def test_prephysics_requires_adjacent_simple_trigger_and_successful_source(self):
        def sources(titles, ok=True):
            return [{'ok': ok, 'records': [{'run': run, 'fields': {'title': title}} for run, title in titles]}]
        good = [(1, 'BEAM(4-10) x1/90 check'), (2, 'KBEAM x1/12'), (3, 'E72Physics')]
        self.assertEqual({1, 2}, run_ledger.prephysics_runs(sources(good)))
        self.assertEqual(set(), run_ledger.prephysics_runs(sources(good, ok=False)))
        for titles in ([(1, 'BEAM x1/90'), (3, 'Physics Run')],
                       [(1, 'BEAM x1/90'), (2, 'beam tuning'), (3, 'Physics Run')],
                       [(1, 'KBeamx /KVCxHTOF-Mp2'), (2, 'Physics Run')]):
            self.assertEqual(set(), run_ledger.prephysics_runs(sources(titles)))

    def test_old_headers_upgrade_without_changing_personal_columns(self):
        old_headers = [header for header in HEADERS if header not in run_ledger.OPTIONAL_HEADERS]+['自由列']
        old_row = ['run00001' if header == 'run_key' else 'my memo' if header == 'メモ'
                   else 'my tag' if header == 'タグ' else 'private' if header == '自由列' else '' for header in old_headers]
        sheet = FakeSheets([old_headers, old_row])
        record = dict.fromkeys(run_ledger.MACHINE_HEADERS, '')
        record.update(run_key='run00001', **{'用途': '較正'})
        sheet.publish('Test', plan_updates(sheet.grid, [record], 'Test'))
        headers = sheet.grid[0]
        self.assertEqual('my memo', sheet.grid[1][headers.index('メモ')])
        self.assertEqual('my tag', sheet.grid[1][headers.index('タグ')])
        self.assertEqual('private', sheet.grid[1][headers.index('自由列')])
        self.assertEqual('較正', sheet.grid[1][headers.index('用途')])
        self.assertEqual([], plan_updates(sheet.grid, [record], 'Test')['updates'])

    def test_style_rules_target_only_status_and_run_cells(self):
        from common.lib import run_ledger_style
        headers = list(reversed(HEADERS))
        rules = run_ledger_style.rules_for(5, headers)
        allowed = set(run_ledger_style.STATUS_COLUMNS) | {'run番号'}
        for rule in rules:
            self.assertTrue(run_ledger_style.owned(rule))
            for region in rule['ranges']:
                self.assertIn(headers[region['startColumnIndex']], allowed)
                self.assertEqual(1, region['startRowIndex'])
                self.assertEqual(region['startColumnIndex']+1, region['endColumnIndex'])
        self.assertFalse(run_ledger_style.owned({'booleanRule': {'condition': {'type': 'TEXT_EQ', 'values': [{'userEnteredValue': 'my rule'}]}}}))

    def test_style_signature_accepts_default_sheet_id_and_quantized_colors(self):
        from common.lib import run_ledger_style
        expected = run_ledger_style.rules_for(0, HEADERS)[0]
        observed = copy.deepcopy(expected)
        observed['ranges'][0].pop('sheetId')
        rgb = observed['booleanRule']['format']['backgroundColor']
        for key in rgb:
            rgb[key] = int(rgb[key]*255)/255
        self.assertEqual(run_ledger_style.signature(expected), run_ledger_style.signature(observed))


if __name__ == '__main__':
    unittest.main()
