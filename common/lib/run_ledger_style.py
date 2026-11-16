"""Presentation rules restricted to machine-owned ledger cells."""
from .google_sheets_ledger import column_name

MARKER = 'e72_ledger_style_v1'
STATUS_COLUMNS = ('較正状態', 'Hodo', 'BcIn', 'BcOut', 'D5', '解析状態')
PALETTE = {
    'red': (252, 220, 220), 'blue': (219, 234, 254), 'yellow': (255, 239, 190),
    'green': (215, 239, 216), 'gray': (238, 238, 238),
    'physics': (230, 221, 248), 'calibration': (206, 239, 236), 'white': (255, 255, 255),
}


#_____________________________________________________________________________
def color(name):
    return dict(zip(('red', 'green', 'blue'), [min(1.0, (value+0.5)/255) for value in PALETTE[name]]))


#_____________________________________________________________________________
def rules_for(sheet_id, headers, row_count=1000):
    rules = []
    def add(column, predicate, shade):
        rules.append({'ranges': [{'sheetId': sheet_id, 'startRowIndex': 1, 'endRowIndex': row_count,
                                 'startColumnIndex': column, 'endColumnIndex': column+1}],
                      'booleanRule': {'condition': {'type': 'CUSTOM_FORMULA', 'values': [
                          {'userEnteredValue': f'=AND(N("{MARKER}")=0,{predicate})'}]},
                          'format': {'backgroundColor': color(shade)}}})
    for header in STATUS_COLUMNS:
        column = headers.index(header)
        cell = column_name(column)+'2'
        for pattern, shade in [
            ('^(失敗|状態不整合)', 'red'),
            ('^(実行中|待機中)', 'blue'),
            ('^(保留|停止|最終確認待ち|承認記録あり|未確認|対象工程の処理完走)', 'yellow'),
            ('^対象範囲の承認済み', 'green'),
            ('^(未連携|工程別の記録あり)', 'gray'),
        ]:
            add(column, f'REGEXMATCH(TO_TEXT({cell}),"{pattern}")', shade)
    usage_cell = '$'+column_name(headers.index('用途'))+'2'
    for usage, shade in [('Physics', 'physics'), ('較正', 'calibration'), ('その他', 'white'), ('不明', 'gray')]:
        add(headers.index('run番号'), f'{usage_cell}="{usage}"', shade)
    return rules


#_____________________________________________________________________________
def owned(rule):
    return any(MARKER in str(value.get('userEnteredValue', ''))
               for value in rule.get('booleanRule', {}).get('condition', {}).get('values', []))


#_____________________________________________________________________________
def signature(rule):
    boolean = rule['booleanRule']
    rgb = boolean['format'].get('backgroundColor', {})
    ranges = [dict(sheetId=0, startRowIndex=0, startColumnIndex=0) | region for region in rule['ranges']]
    return (ranges, boolean['condition'], tuple(int(rgb.get(key, 0)*255+0.001) for key in ('red', 'green', 'blue')))


#_____________________________________________________________________________
def apply_style(client, title):
    document = client.metadata()
    sheet = next(row for row in document['sheets'] if row['properties']['title'] == title)
    from .google_sheets_ledger import read_table
    table = read_table(client.read(title))
    headers = table['headers']
    sheet_id = sheet['properties']['sheetId']
    desired = rules_for(sheet_id, headers, sheet['properties']['gridProperties']['rowCount'])
    existing = sheet.get('conditionalFormats', [])
    if [signature(rule) for rule in existing if owned(rule)] == [signature(rule) for rule in desired]:
        return False
    requests = [{'deleteConditionalFormatRule': {'sheetId': sheet_id, 'index': index}}
                for index in reversed(range(len(existing))) if owned(existing[index])]
    requests += [{'addConditionalFormatRule': {'rule': rule, 'index': index}}
                 for index, rule in enumerate(desired)]
    for header, width in [('run番号', 95), ('タイトル', 560), *((header, 185) for header in STATUS_COLUMNS)]:
        column = headers.index(header)
        requests += [
            {'updateDimensionProperties': {
                'range': {'sheetId': sheet_id, 'dimension': 'COLUMNS', 'startIndex': column, 'endIndex': column+1},
                'properties': {'pixelSize': width}, 'fields': 'pixelSize'}},
            {'repeatCell': {'range': {'sheetId': sheet_id, 'startColumnIndex': column, 'endColumnIndex': column+1},
                            'cell': {'userEnteredFormat': {'wrapStrategy': 'WRAP', 'verticalAlignment': 'TOP'}},
                            'fields': 'userEnteredFormat.wrapStrategy,userEnteredFormat.verticalAlignment'}},
        ]
    requests.append({'autoResizeDimensions': {'dimensions': {
        'sheetId': sheet_id, 'dimension': 'ROWS', 'startIndex': 0, 'endIndex': table['last_row']}}})
    client.request('POST', ':batchUpdate', {'requests': requests})
    observed = next(row for row in client.metadata()['sheets'] if row['properties']['title'] == title)
    if [signature(rule) for rule in observed.get('conditionalFormats', []) if owned(rule)] != [signature(rule) for rule in desired]:
        raise RuntimeError('Ledger formatting was not confirmed; retry sync')
    return True
