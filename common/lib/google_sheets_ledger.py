"""Bounded Sheets REST access and machine-column-only ledger updates."""
from __future__ import annotations

import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request

from .run_ledger import HUMAN_HEADERS, MACHINE_HEADERS, OPTIONAL_HEADERS, run_key, spreadsheet_id

# Hashes guard the original sources without publishing their share identifiers.
PROTECTED_SOURCE_HASHES = {
    'e90d1de2b467bc4e4975447814c210837ec107e7c3c8fbddbe961df481fd91c6',
    'a997e6c3361d1df098d1b3676d8eec318a6b41cbd2bf0ae062ec4da78d5a1f38',
}

HEADERS = ['run_key', 'run番号', 'タイトル', *HUMAN_HEADERS, '較正状態', 'Hodo', 'BcIn', 'BcOut', 'D5', '解析状態']
HEADERS += [name for name in MACHINE_HEADERS if name not in HEADERS]


#_____________________________________________________________________________
def column_name(index):
    name = ''
    while index >= 0:
        index, remainder = divmod(index, 26)
        name = chr(65+remainder)+name
        index -= 1
    return name


#_____________________________________________________________________________
def quoted_tab(title):
    return "'"+title.replace("'", "''")+"'"


#_____________________________________________________________________________
def read_table(grid):
    if not grid or not any(grid[0]):
        if any(any(str(value) for value in row) for row in grid):
            raise ValueError('Nonempty sheet without ledger header')
        return {'headers': [], 'rows': {}, 'last_row': 1, 'structure': []}
    headers = list(grid[0])
    if len(headers) != len(set(headers)) or any(header not in headers for header in HEADERS if header not in OPTIONAL_HEADERS):
        raise ValueError('Missing/duplicate ledger headers; existing columns are not rewritten')
    key_column = headers.index('run_key')
    rows = {}
    structure = []
    last_row = 1
    for number, row in enumerate(grid[1:], 2):
        if not any(str(value) for value in row):
            continue
        last_row = number
        key = str(row[key_column]).strip() if key_column < len(row) else ''
        structure.append((number, key))
        if not key:
            continue  # A note-only row is reserved, never reused for an append.
        import re
        if not re.fullmatch(r'run\d{5,10}', key) or key != run_key(int(key[3:])):
            raise ValueError(f'Invalid run_key at sheet row {number}')
        if key in rows:
            raise ValueError(f'Duplicate run_key {key}; no writes performed')
        rows[key] = {'number': number, 'values': list(row)}
    return {'headers': headers, 'rows': rows, 'last_row': last_row, 'structure': structure}


#_____________________________________________________________________________
def structure_token(table):
    return hashlib.sha256(json.dumps([table['headers'], table['structure']], ensure_ascii=False).encode()).hexdigest()


#_____________________________________________________________________________
def plan_updates(grid, records, title):
    table = read_table(grid)
    headers = list(table['headers'] or HEADERS)
    updates = []
    for header in OPTIONAL_HEADERS:
        if header not in headers:
            headers.append(header)
            column = column_name(len(headers)-1)
            updates.append({'range': f'{quoted_tab(title)}!{column}1:{column}1', 'values': [[header]]})
    if not table['headers']:
        updates.append({'range': quoted_tab(title)+'!A1:'+column_name(len(headers)-1)+'1', 'values': [headers]})
    occupied = set()
    last_row = table['last_row']
    inserted = changed = 0
    for record in records:
        key = record['run_key']
        if key in occupied:
            raise ValueError(f'Duplicate planned run {key}')
        occupied.add(key)
        existing = table['rows'].get(key)
        if existing:
            number, old_values = existing['number'], existing['values']
        else:
            last_row += 1
            number, old_values = last_row, []
            inserted += 1
        cells = []
        for column, header in enumerate(headers):
            if header not in MACHINE_HEADERS:
                continue
            value = str(record.get(header, ''))
            if len(value) > 45000:
                raise ValueError(f'Cell too long for {key}/{header}')
            old_value = str(old_values[column]) if column < len(old_values) else ''
            if old_value != value:
                cells.append((column, value))
        if cells:
            changed += 1
        groups = []
        for column, value in cells:
            if groups and groups[-1][-1][0]+1 == column:
                groups[-1].append((column, value))
            else:
                groups.append([(column, value)])
        for group in groups:
            updates.append({'range': f'{quoted_tab(title)}!{column_name(group[0][0])}{number}:{column_name(group[-1][0])}{number}',
                            'values': [[value for _, value in group]]})
    return {'before': structure_token(table), 'updates': updates, 'inserted': inserted, 'changed': changed,
            'required_rows': last_row, 'required_columns': len(headers), 'records': records}


#_____________________________________________________________________________
class SheetsClient:
    def __init__(self, destination, credential_file, *, sources=(), max_requests=30):
        self.destination = spreadsheet_id(destination)
        protected = {spreadsheet_id(row['spreadsheet_id']) for row in sources}
        if self.destination in protected or hashlib.sha256(self.destination.encode()).hexdigest() in PROTECTED_SOURCE_HASHES:
            raise ValueError('Destination is an original source sheet; write path refused')
        self.credential_file = str(credential_file)
        self.max_requests = max_requests
        self.requests = 0
        self.credentials = None
        self.base = 'https://sheets.googleapis.com/v4/spreadsheets/'+self.destination

    def _token(self):
        if self.credentials is None:
            try:
                from google.oauth2 import service_account
            except ImportError as error:
                raise RuntimeError('google-auth is missing: install via conda and refresh envfs') from error
            self.credentials = service_account.Credentials.from_service_account_file(
                self.credential_file, scopes=['https://www.googleapis.com/auth/spreadsheets'])
        if not self.credentials.valid:
            from google.auth.transport.requests import Request
            transport = Request()
            def limited(*args, **kwargs):
                kwargs['timeout'] = 25
                return transport(*args, **kwargs)
            try:
                self.credentials.refresh(limited)
            except Exception as error:
                raise RuntimeError('Google authentication failed ('+type(error).__name__+')') from None
        return self.credentials.token

    def request(self, method, suffix, body=None):
        reading = method == 'GET' or suffix in ('/developerMetadata:search', '/values:batchGetByDataFilter')
        attempts = 3 if reading else 1
        for attempt in range(attempts):
            if self.requests >= self.max_requests:
                raise RuntimeError('Manual sync request budget exceeded; retry later')
            self.requests += 1
            headers = {'Authorization': 'Bearer '+self._token(), 'Content-Type': 'application/json'}
            payload = json.dumps(body, ensure_ascii=False).encode() if body is not None else None
            request = urllib.request.Request(self.base+suffix, headers=headers, data=payload, method=method)
            try:
                with urllib.request.urlopen(request, timeout=25) as response:
                    return json.load(response)
            except urllib.error.HTTPError as error:
                if not reading or error.code not in (429, 500, 502, 503, 504) or attempt+1 == attempts:
                    raise RuntimeError(f'Google Sheets HTTP {error.code}') from None
            except OSError:
                if attempt+1 == attempts:
                    raise RuntimeError('Google Sheets connection failed') from None
            time.sleep(2**attempt)

    def metadata(self):
        return self.request('GET', '?fields=sheets(properties,conditionalFormats)')

    def ensure_tab(self, title):
        document = self.metadata()
        found = [row['properties'] for row in document.get('sheets', []) if row['properties']['title'] == title]
        if found:
            return found[0]
        try:
            response = self.request('POST', ':batchUpdate', {'requests': [{'addSheet': {'properties': {'title': title}}}]})
            return response['replies'][0]['addSheet']['properties']
        except RuntimeError:
            # A lost response is not evidence that creation failed; never blindly add again.
            found = [row['properties'] for row in self.metadata().get('sheets', []) if row['properties']['title'] == title]
            if len(found) == 1:
                return found[0]
            raise

    def read(self, title):
        region = quoted_tab(title)
        suffix = '/values/'+urllib.parse.quote(region, safe='')+'?valueRenderOption=FORMATTED_VALUE'
        return self.request('GET', suffix).get('values', [])

    def grow(self, properties, plan):
        current = properties.get('gridProperties') or {}
        wanted = {'rowCount': max(current.get('rowCount', 1000), plan['required_rows']+20),
                  'columnCount': max(current.get('columnCount', 26), plan['required_columns']),
                  'frozenRowCount': max(current.get('frozenRowCount', 0), 1), 'frozenColumnCount': max(current.get('frozenColumnCount', 0), 2)}
        if all(current.get(key) == value for key, value in wanted.items()):
            return
        self.request('POST', ':batchUpdate', {'requests': [{'updateSheetProperties': {
            'properties': {'sheetId': properties['sheetId'], 'gridProperties': wanted},
            'fields': ','.join('gridProperties.'+key for key in wanted)}}]})

    def append_rows(self, sheet_id, rows):
        pending = []
        size = 0
        for values in rows:
            row = {'values': [{'userEnteredValue': {'stringValue': value}} if value is not None else {}
                              for value in values]}
            length = len(json.dumps(row, ensure_ascii=False).encode())
            if pending and size+length > 700000:
                self.request('POST', ':batchUpdate', {'requests': [{'appendCells': {
                    'sheetId': sheet_id, 'rows': pending, 'fields': 'userEnteredValue'}}]})
                pending, size = [], 0
            pending.append(row)
            size += length
        if pending:
            self.request('POST', ':batchUpdate', {'requests': [{'appendCells': {
                'sheetId': sheet_id, 'rows': pending, 'fields': 'userEnteredValue'}}]})

    def row_ids(self, title, keys, sheet_id):
        # Metadata follows rows during a sort/insert. Verify the initial association before using it.
        for attempt in range(2):
            document = self.request('POST', '/developerMetadata:search', {'dataFilters': [
                {'developerMetadataLookup': {'metadataKey': 'e72_run_key_v1', 'visibility': 'DOCUMENT'}}]})
            table = read_table(self.read(title))
            found = {}
            for match in document.get('matchedDeveloperMetadata', []):
                item = match['developerMetadata']
                location = item.get('location', {}).get('dimensionRange', {})
                if item.get('metadataKey') != 'e72_run_key_v1' or location.get('sheetId', 0) != sheet_id:
                    continue
                if location.get('dimension') != 'ROWS' or location.get('endIndex', 0)-location.get('startIndex', 0) != 1:
                    raise RuntimeError('Invalid ledger row metadata')
                key = item.get('metadataValue')
                if key in found:
                    raise RuntimeError('Duplicate ledger row metadata; no value updates performed')
                row = table['rows'].get(key)
                if row is None or row['number'] != location.get('startIndex', 0)+1:
                    self.request('POST', ':batchUpdate', {'requests': [{'deleteDeveloperMetadata': {
                        'dataFilter': {'developerMetadataLookup': {'metadataId': item['metadataId']}}}}]})
                    raise RuntimeError('Row metadata/key association changed; tag removed without value updates; rerun sync')
                found[key] = item['metadataId']
            missing = set(keys)-set(found)
            if not missing:
                return found, table
            if attempt:
                raise RuntimeError('Row metadata creation unconfirmed; retry sync')
            requests = []
            for key in sorted(missing):
                row = table['rows'].get(key)
                if row is None:
                    raise RuntimeError('Run disappeared before metadata association; retry sync')
                requests.append({'createDeveloperMetadata': {'developerMetadata': {
                    'metadataKey': 'e72_run_key_v1', 'metadataValue': key, 'visibility': 'DOCUMENT',
                    'location': {'dimensionRange': {'sheetId': sheet_id, 'dimension': 'ROWS',
                                                   'startIndex': row['number']-1, 'endIndex': row['number']}}}}})
            try:
                self.request('POST', ':batchUpdate', {'requests': requests})
            except RuntimeError:
                # The following read resolves lost responses without creating duplicate metadata.
                pass
        raise RuntimeError('Could not resolve ledger row metadata')

    def publish(self, title, plan):
        table = read_table(self.read(title))
        if structure_token(table) != plan['before']:
            raise RuntimeError('Sheet keys/header changed since planning; rerun sync')
        properties = next(row['properties'] for row in self.metadata()['sheets']
                          if row['properties']['title'] == title)
        sheet_id = properties['sheetId']
        headers = list(table['headers'] or HEADERS)
        if not table['headers']:
            self.append_rows(sheet_id, [headers])
        else:
            # Column/header edits remain prohibited during a sync; row movement is handled by metadata.
            extra = [name for name in OPTIONAL_HEADERS if name not in headers]
            if extra:
                updates = []
                for name in extra:
                    column = column_name(len(headers))
                    headers.append(name)
                    updates.append({'range': f'{quoted_tab(title)}!{column}1:{column}1', 'values': [[name]]})
                self.request('POST', '/values:batchUpdate', {'valueInputOption': 'RAW', 'data': updates})
        table = read_table(self.read(title))
        if table['headers'] != headers:
            raise RuntimeError('Sheet header changed; no run values updated')
        missing = [record for record in plan['records'] if record['run_key'] not in table['rows']]
        pending = None
        if missing:
            try:
                self.append_rows(sheet_id, [[str(record.get(name, '')) if name in MACHINE_HEADERS else None
                                            for name in headers] for record in missing])
            except RuntimeError as error:
                pending = error
            table = read_table(self.read(title))
            if not all(record['run_key'] in table['rows'] for record in missing):
                raise RuntimeError('Append unconfirmed; manual retry re-resolves run keys') from pending
        keys = [record['run_key'] for record in plan['records']]
        if not keys:
            return {'verified': True, 'requests': self.requests, 'changed': plan['changed'], 'inserted': plan['inserted']}
        identifiers, table = self.row_ids(title, keys, sheet_id)
        if table['headers'] != headers:
            raise RuntimeError('Sheet columns changed before value update')
        batches, batch, size = [], [], 0
        for record in plan['records']:
            old_values = table['rows'][record['run_key']]['values']
            values = [None]*len(headers)
            for column, name in enumerate(headers):
                if name not in MACHINE_HEADERS:
                    continue
                value = str(record.get(name, ''))
                old_value = str(old_values[column]) if column < len(old_values) else ''
                if value != old_value:
                    values[column] = value
            if all(value is None for value in values):
                continue
            update = {'dataFilter': {'developerMetadataLookup': {'metadataId': identifiers[record['run_key']]}},
                      'majorDimension': 'ROWS', 'values': [values]}
            length = len(json.dumps(update, ensure_ascii=False).encode())
            if batch and size+length > 700000:
                batches.append(batch)
                batch, size = [], 0
            batch.append(update)
            size += length
        if batch:
            batches.append(batch)
        for updates in batches:
            try:
                self.request('POST', '/values:batchUpdateByDataFilter', {'valueInputOption': 'RAW', 'data': updates})
            except RuntimeError as error:
                pending = error
                break
        observed = read_table(self.read(title))
        if observed['headers'] != headers:
            raise RuntimeError('Sheet header changed during update; stop column edits and retry')
        for record in plan['records']:
            row = observed['rows'].get(record['run_key'])
            if row is None:
                raise RuntimeError('Sheet write unconfirmed; manual retry re-resolves keys') from pending
            for header in MACHINE_HEADERS:
                column = headers.index(header)
                value = str(row['values'][column]) if column < len(row['values']) else ''
                if value != str(record.get(header, '')):
                    raise RuntimeError('Sheet changed or write incomplete; retry sync') from pending
        return {'verified': True, 'requests': self.requests, 'changed': plan['changed'], 'inserted': plan['inserted']}
