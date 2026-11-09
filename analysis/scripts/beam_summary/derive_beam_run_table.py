#!/usr/bin/env python3
"""Derive the BeamRunCondition run table (src/BeamRunCondition.cc, kBeamRunRanges) from the two
run-summary Google sheets. Read-only on the sheets; writes a JSON of per-run decisions and the C++
table body. Rules: docs/beam_summary/beam-run-condition-handoff.md (data rules, owner decisions),
beam-run-condition-new-data-after-2544.md and the user decisions of 2026-09-28: every run whose
beam particle and momentum are known is listed (calibration, study, check and daq-test runs too);
dummy/junk/clock/tune/cosmic/no-beam runs and empty titles are not (pedestal runs are listed
when the particle and momentum are known); in a mixed-beam title
the beam is the non-pion particle (pions contaminate every beam); "slit close" only narrows the
beam; the two sheets cover different periods (a run in both is held for review).
See docs/beam_summary/beam-run-condition-table-20260928.md.

Download the sheets as CSV (link-shared, no login needed):
  curl -sL -o s1.csv 'https://docs.google.com/spreadsheets/d/1InLDUJ0r0K-O6rJH3NRi5tuxHurfRHQ6BXPYSFhVyuA/export?format=csv&gid=587812810'
  curl -sL -o s2.csv 'https://docs.google.com/spreadsheets/d/1aRQxyxYVs1FBuR2Pa83zUQRgdYaYYOJujOiqYGNCFHI/export?format=csv&gid=1541909500'
Then:
  python3 derive_beam_run_table.py --s1 s1.csv --s2 s2.csv --json derived.json --table table_body.txt
and paste table_body.txt between "kBeamRunRanges[] = {" and "};". The script also prints the review
list and runs the boundary tests of the handoff.
"""
import argparse, csv, re, json, collections
ap = argparse.ArgumentParser()
ap.add_argument('--s1', required=True, help='2025-11 run summary CSV (gid 587812810)')
ap.add_argument('--s2', required=True, help='2026-04 run summary CSV (gid 1541909500)')
ap.add_argument('--json', default='derived_beam.json')
ap.add_argument('--table', default='table_body.txt')
args = ap.parse_args()
OWNER_KAON = {2610, 2647, 2680, 2774, 2775, 2815, 2856, 2884, 2940, 2977, 2993}
OWNER_UNKNOWN = {2320, 2353, 2358, 2422, 2745, 2746, 2773, 2798, 2855, *range(2876, 2883), 2885, 2886, 2991, 2994, 3004}
OWNER_JUNK = {2394, 2433}
def rng(*spec):
    out = set()
    for x in spec:
        out |= set(range(x[0], x[1] + 1)) if isinstance(x, tuple) else {x}
    return out
EXCL_RE = re.compile(r'dummy|junk|clock|clcok|pedestal|daq test|cosmic|no beam|\btest\b')
# "pure" exclusion for rows that have an explicit particle column (the column is filled even on
# dummy/junk rows): title starts with dummy/junk/clock/cosmic (a pedestal run with a known beam
# particle and momentum is listed), "-> dummy/junk", junk word,
# or "no beam". Trigger mentions like "Clock1kHz" and "slit close" (narrowed beam) do not exclude.
PURE_EXCL = re.compile(r'^(dummy|junk|clock|clcok|cosmic)\b|->\s*(dummy|junk)|\bjunk\b|no beam')
def norm_particle(s):
    s = s.strip().lower()
    if s in ('pion', 'pi', 'pi+', 'pi-'): return 'Pion'
    if s in ('kaon', 'k', 'k+', 'k-'): return 'Kaon'
    if s in ('proton', 'pbar', 'p-bar', 'antiproton', 'p'): return 'Proton'
    return None
def mom_col(s):
    s = s.strip()
    return abs(float(s)) / 1000. if re.fullmatch(r'[-+]?\d+(\.\d+)?', s) else None
def mom_title(t):
    m = re.search(r'([-+]?\d\.\d+)\s*(?:gev/c)?', t.lower())
    m2 = re.search(r'([-+]?\d\.\d+)\s*gev/c', t.lower())
    if m2: return abs(float(m2.group(1)))
    if m and re.match(r'^\s*[-+]?\d\.\d+\s', t): return abs(float(m.group(1)))  # e.g. "-0.790 Beam"
    return None
def title_particle(t):
    tl = t.lower()
    pats = [('Proton', r'p-bar|pbar|antiproton|\bproton\b'), ('Kaon', r'\bk-|\bk\+'), ('Pion', r'\bpi-|\bpi\+|\bpi\b|\bpion')]
    return [p for p, r in pats if re.search(r, tl)]
out = {}; review = []
def add(run, sheet, particle, mom, why, title, flag=None):
    out[run] = dict(run=run, sheet=sheet, particle=particle, mom=mom, why=why, title=title[:120])
    if flag: review.append((run, flag, title[:120]))
# sheet 1: 2025-11 (B run, F title, J particle, K momentum)
for r in list(csv.reader(open(args.s1)))[2:]:
    if len(r) <= 10 or not r[1].strip().isdigit(): continue
    run = int(r[1]); t = ' '.join(r[5].split()); tl = t.lower(); pc = r[9].strip(); mc = mom_col(r[10]); mt = mom_title(t)
    if run in OWNER_JUNK or not t: add(run, 1, 'EXCL', None, 'junk/empty title', t); continue
    has_col = norm_particle(pc) is not None
    if PURE_EXCL.search(tl): add(run, 1, 'EXCL', None, 'excluded title (dummy/junk/clock/pedestal/cosmic/no beam)', t); continue
    if not has_col and EXCL_RE.search(tl) and 'physics' not in tl: add(run, 1, 'EXCL', None, 'excluded title', t); continue
    mom = mc if mc is not None else mt
    flag = None
    if mc is not None and mt is not None and abs(mc - mt) > 1e-6: flag = f'momentum column {mc} != title {mt}'
    p = norm_particle(pc)
    if pc and p is None: add(run, 1, 'EXCL', None, f'particle column "{pc}"', t); continue
    if run in OWNER_UNKNOWN: add(run, 1, 'Unknown', mom, 'owner: Unknown', t); continue
    if p: add(run, 1, p, mom, 'particle column', t, flag); continue
    if run in OWNER_KAON: add(run, 1, 'Kaon', mom, 'owner: kaon', t, flag); continue
    if 'physics run' in tl: add(run, 1, 'Kaon', mom, 'Physics Run -> kaon', t, flag); continue
    if 'beam tune' in tl or re.search(r'\btune\b', tl) or 'check' in tl: add(run, 1, 'Unknown', mom, 'check/tune', t); continue
    if 'e45 study' in tl: add(run, 1, 'Pion', mom, 'E45 study -> pion', t, flag); continue
    if 'htof calibration' in tl:
        tp = title_particle(t)
        add(run, 1, 'Pion' if tp == ['Pion'] else 'REVIEW', mom, 'HTOF calibration (title particle)', t, flag if tp == ['Pion'] else f'calib title particles {tp}'); continue
    if 'beam through' in tl or 'beam-through' in tl or 'p-bar study' in tl or 'ftof study' in tl:
        tp = title_particle(t)
        if len(tp) == 2 and 'Pion' in tp: tp = [x for x in tp if x != 'Pion']  # pions contaminate any beam
        if len(tp) == 1: add(run, 1, tp[0], mom, 'title particle', t, flag)
        else: add(run, 1, 'REVIEW', mom, 'mixed/unclear title particle', t, f'title particles {tp}')
        continue
    if 'study' in tl: add(run, 1, 'Unknown', mom, 'study (gate etc.)', t); continue
    if re.search(r'k:\s*\d+k', tl) and 'beam' in tl: add(run, 1, 'Kaon', mom, 'normal Beam run (K/pi intensity) -> kaon', t, flag); continue
    add(run, 1, 'REVIEW', mom, 'unclassified', t, 'unclassified')
# sheet 2: 2026-04 (B run, F experiment, G title, K particle, L momentum)
for r in list(csv.reader(open(args.s2)))[2:]:
    if len(r) <= 11 or not r[1].strip().isdigit(): continue
    run = int(r[1]); exp = r[5].strip(); t = ' '.join(r[6].split()); pc = r[10].strip(); mc = mom_col(r[11]); mt = mom_title(t)
    if run in out:
        prev = out[run]
        if run == 1916 and prev['particle'] == 'Pion' and abs((prev['mom'] or 0) - 0.735) < 1e-6 and norm_particle(pc) == 'Pion' and mc and abs(mc - 0.735) < 1e-6:
            continue  # known consistent duplicate (handoff)
        review.append((run, 'in both sheets', t)); continue
    if exp.lower() == 'junk': add(run, 2, 'EXCL', None, 'experiment junk', t); continue
    if PURE_EXCL.search(t.lower()): add(run, 2, 'EXCL', None, 'excluded title (dummy/junk/clock/pedestal/cosmic/no beam)', t); continue
    p = norm_particle(pc)
    if p is None: add(run, 2, 'EXCL', None, f'particle column "{pc}"', t); continue
    # owner/audit decision (beam-run-condition-new-data-after-2544.md l.164-165): DAQ/TPC/trigger
    # check and daq test rows are Unknown/Excluded even with a particle column
    mom = mc if mc is not None else mt
    flag = None
    if mc is not None and mt is not None and abs(mc - mt) > 1e-6: flag = f'momentum column {mc} != title {mt}'
    if 3688 <= run <= 3703: mom = 1.144; flag = None  # owner-confirmed override
    add(run, 2, p, mom, 'particle column', t, flag)
# owner assumption: 2059-2062 pion 1.144 (assumed)
for run in range(2059, 2063):
    if run in out and out[run]['particle'] == 'Pion' and out[run]['mom'] is None:
        out[run]['mom'] = 1.144; out[run]['why'] += ' + ASSUMED 1.144 (handoff)'
json.dump({'rows': out, 'review': review}, open(args.json, 'w'), indent=0)
print(collections.Counter(v['particle'] for v in out.values()))
print('no momentum but particle known:', [k for k, v in sorted(out.items()) if v['particle'] in ('Pion', 'Kaon', 'Proton') and v['mom'] is None])
print('review:'); [print('  ', x) for x in review]

# ---- C++ table body ----
d = json.load(open(args.json)); rows = {int(k): v for k, v in d['rows'].items()}
known = {n: (v['particle'], round(v['mom'], 3), v) for n, v in rows.items() if v['particle'] in ('Pion', 'Kaon', 'Proton')}
runs = sorted(known)
ranges = []
for n in runs:
    p, m, v = known[n]
    tag = ''
    if 2059 <= n <= 2062: tag = 'assumed'
    elif 3688 <= n <= 3703: tag = 'override'
    if ranges and ranges[-1]['last'] == n - 1 and ranges[-1]['p'] == p and ranges[-1]['m'] == m and ranges[-1]['tag'] == tag:
        ranges[-1]['last'] = n
    else:
        ranges.append(dict(first=n, last=n, p=p, m=m, tag=tag))
notes = {'assumed': '  // momentum not in the summary: 1.144 GeV/c assumed from run 2058',
         'override': '  // 1.144 GeV/c (the momentum column had 1.114, a typo)',
}
lines = []; prev_m = None
for r in ranges:
    if r['m'] != prev_m:
        if lines: lines.append('')
        lines.append(f"  // --- {r['m']:.3f} GeV/c ---"); prev_m = r['m']
    if r['tag']: lines.append(notes[r['tag']])
    lines.append(f"  {{{r['first']}, {r['last']}, BeamSpecies::{r['p'] + ',':<8} {r['m']:.3f}}},")
open(args.table, 'w').write('\n'.join(lines) + '\n')
print(len(ranges), 'ranges,', len(runs), 'runs')
# boundary tests from the handoff (+ user decisions)
exp = {1916: ('Pion', .735), 1933: ('Pion', .735), 1934: ('Kaon', .735), 2058: ('Pion', 1.144), 2105: ('Kaon', .645),
       3650: ('Pion', .964), 3704: ('Proton', .964), 3730: ('Pion', .964), 3743: ('Pion', 1.144), 3762: ('Kaon', .735),
       3782: ('Kaon', .902), 3848: ('Kaon', .933), 3868: ('Kaon', .814), 2208: ('Kaon', .735), 2518: ('Proton', 1.0),
       3774: ('Kaon', .735), 2059: ('Pion', 1.144), 3688: ('Proton', 1.144), 3689: ('Proton', 1.144), 3690: ('Proton', 1.144),
       2601: ('Pion', 1.0), 2607: ('Pion', 1.0), 2296: ('Kaon', .735), 2544: ('Proton', .3)}
nomatch = [2394, 2433, 1948, 2141, 3667, 2605, 2600, 99999, 2202, 2245]
bad = 0
for n, e in exp.items():
    g = known.get(n)
    if not g or g[0] != e[0] or abs(g[1] - e[1]) > 1e-9: print('FAIL', n, e, g and g[:2]); bad += 1
for n in nomatch:
    if n in known: print('FAIL no-match expected', n, known[n][:2]); bad += 1
print('boundary tests:', 'all passed' if not bad else f'{bad} failures')
