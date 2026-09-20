# -*- coding: utf-8 -*-
"""run_08: compare the official cytoHubba export vs the in-house implementation, Top10 per method.

Usage: python run_08_compare_cytohubba.py [official export csv]
Default: 04_analysis/cytohubba_official_scores.csv
Outputs: console comparison + 04_analysis/run08_official_vs_python.csv
"""
import csv, os, sys
from collections import Counter

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

P = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY_CSV = os.path.join(P, '04_analysis', 'run06_hub_scores.csv')
OFF_CSV = sys.argv[1] if len(sys.argv) > 1 else os.path.join(P, '04_analysis', 'cytohubba_official_scores.csv')
OUT = os.path.join(P, '04_analysis', 'run08_official_vs_python.csv')
# NOTE (2026-09-20): see run_09_epc_official_replication.py for the final EPC replication and comparison outputs (rho = 0.9965).
METHODS = ['MCC', 'MNC', 'Degree', 'Closeness', 'Radiality', 'EPC']


def read_rows(p):
    with open(p, encoding='utf-8-sig', errors='replace') as f:
        rd = csv.reader(f)
        return [r for r in rd if r and any(c.strip() for c in r)]


def exact_or_sub(hdr, name):
    lo = [h.strip().lower() for h in hdr]
    if name.lower() in lo:
        return lo.index(name.lower())
    for i, h in enumerate(lo):
        if name.lower() in h:
            return i
    return None


# ---------- in-house side ----------
print('== in-house:', os.path.basename(PY_CSV))
py_rows = read_rows(PY_CSV)
py_hdr = py_rows[0]
py_top = {}
for m in METHODS:
    ri = exact_or_sub(py_hdr, 'rank_' + m)
    rows_m = []
    for r in py_rows[1:]:
        try:
            rows_m.append((float(r[ri]), r[0].strip()))
        except Exception:
            pass
    rows_m.sort(key=lambda t: (t[0], t[1]))
    py_top[m] = [g for _, g in rows_m[:10]]
    print('  in-house %-9s top10: %s' % (m, ', '.join(py_top[m])))

# ---------- official side ----------
of_rows = read_rows(OFF_CSV)
of_hdr = of_rows[0]
print()
print('== official:', os.path.basename(OFF_CSV))
print('  header:', of_hdr)
ni = exact_or_sub(of_hdr, 'node_name')
if ni is None:
    ni = 0
of_top = {}
of_maps = {}
for m in METHODS:
    ci = exact_or_sub(of_hdr, m)
    vals = []
    bad = 0
    for r in of_rows[1:]:
        try:
            vals.append((float(r[ci]), r[ni].strip()))
        except Exception:
            bad += 1
    vals.sort(key=lambda t: (-t[0], t[1]))
    of_top[m] = [g for _, g in vals[:10]]
    of_maps[m] = dict((g, v) for v, g in vals)
    c = Counter(v for v, _ in vals)
    topv = vals[0][0]
    v10 = vals[9][0] if len(vals) >= 10 else None
    print('  official %-9s top10: %s' % (m, ', '.join(of_top[m])))
    print('      [diagnostic] parse failures=%-3d top score=%s (%d-way tie) 10th value=%s (%d-way tie)' % (
        bad, topv, c.get(topv, 0), v10, c.get(v10, 0) if v10 is not None else 0))

# ---------- comparison ----------
print()
print('== comparison (Top10 overlap) ==')
with open(OUT, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['Method', 'Overlap', 'PythonTop10', 'OfficialTop10'])
    for m in METHODS:
        a, b = set(py_top[m]), set(of_top[m])
        ov = len(a & b) if a and b else 0
        w.writerow([m, ov, ';'.join(py_top[m]), ';'.join(of_top[m])])
        print('  %-9s overlap %2d/10 | in-house only: %s | official only: %s' % (
            m, ov, ','.join(sorted(a - b)) or '-', ','.join(sorted(b - a)) or '-'))
print('written:', OUT)

# ---------- consistency check: ranks of in-house top10 nodes in the official table ----------
print()
print('== consistency check (rank ranges of in-house top10 nodes in the official table; ranges given for ties) ==')
for m in METHODS:
    print('----', m)
    allvals = sorted(of_maps[m].values(), reverse=True)
    for g in py_top[m]:
        v = of_maps[m].get(g)
        if v is None:
            print('     %-9s missing from official table' % g)
            continue
        gt = sum(1 for x in allvals if x > v + 1e-9)
        ge = sum(1 for x in allvals if x >= v - 1e-9)
        print('     %-9s official=%-22s rank range %d-%d' % (g, v, gt + 1, ge))
