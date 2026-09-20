# -*- coding: utf-8 -*-
"""B1: specificity controls for the 479-gene overlap (hypergeometric, permutation, stratified)."""
import csv
import json
import os
import sys
from collections import Counter

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

import numpy as np
from scipy import stats

BASE = r'C:\Users\user\workspace\PFOA_NAFLD_project'
ANA = os.path.join(BASE, 'data')
SPEC = os.path.join(BASE, 'data', 'specificity')

# universe = symbols measurable on GPL14951 (from GSE89632 top table)
uni = set()
for r in csv.DictReader(open(os.path.join(SPEC, 'GSE89632_NASH_vs_HC_topTable.csv'), encoding='utf-8-sig')):
    s = (r.get('Symbol') or '').strip().upper()
    if s:
        uni.add(s)
print('universe:', len(uni))

# pools
comp = {}
for r in csv.DictReader(open(os.path.join(SPEC, 'PFOA_target_pool_v1.tsv'), encoding='utf-8-sig'), delimiter='\t'):
    comp[r['Symbol'].strip().upper()] = (r['Sources'] or '').upper()
dis = {}
for r in csv.DictReader(open(os.path.join(SPEC, 'NAFLD_disease_pool_v1.tsv'), encoding='utf-8-sig'), delimiter='\t'):
    dis[r['Symbol'].strip().upper()] = (r['Sources'] or '').upper()
print('compound pool:', len(comp), '| disease pool:', len(dis))

K = len(comp)
n = len(dis)
N = len(uni)
inter = set(comp) & set(dis)
obs = len(inter)
exp = K * n / N
fold = obs / exp
p_hyper = stats.hypergeom.sf(obs - 1, N, K, n)
print('observed =', obs, '| expected = %.1f' % exp, '| fold = %.3f' % fold, '| hypergeom P =', '%.3e' % p_hyper)

# permutation: draw K genes from universe, overlap with disease (restricted to universe)
rng = np.random.default_rng(42)
uni_list = sorted(uni)
d_mask = np.array([1 if g in dis else 0 for g in uni_list], dtype=np.int32)
uni_idx = np.arange(N)
counts = np.empty(10000, dtype=np.int32)
for i in range(10000):
    pick = rng.choice(uni_idx, size=K, replace=False)
    counts[i] = d_mask[pick].sum()
z = (obs - counts.mean()) / counts.std()
p_perm = float((counts >= obs).mean())
print('permutation: mean=%.1f sd=%.1f z=%.1f p=%.4g (10k draws)' % (counts.mean(), counts.std(), z, p_perm))

# stratified analyses
ctd_only = {g for g, s in comp.items() if 'CTD' in s}
pred_only = {g for g, s in comp.items() if any(k in s for k in ('TARGETNET', 'PHARMMAPPER', 'SWISSTARGET', 'SEA'))}
pred_not_ctd = pred_only - ctd_only
print('CTD-only compound:', len(ctd_only), '| predicted (any):', len(pred_only), '| predicted-not-CTD:', len(pred_not_ctd))

for lab, cset in [('CTD_only', ctd_only), ('predicted_only_not_CTD', pred_not_ctd), ('predicted_any', pred_only)]:
    Kx = len(cset)
    obsx = len(cset & set(dis))
    expx = Kx * n / N
    ph = stats.hypergeom.sf(obsx - 1, N, Kx, n)
    print('%s: K=%d obs=%d exp=%.1f fold=%.2f P=%.3e' % (lab, Kx, obsx, expx, obsx / expx, ph))

core = ['TNF', 'IL6', 'TP53', 'ALB', 'IL1B', 'INS']
ext = ['AKT1', 'PPARG', 'EGFR', 'BCL2']
sig5 = ['EPHA2', 'FMO1', 'MYC', 'PHLDA1', 'WNT5A']
ctd_inter = ctd_only & set(dis)
print('CTD-only ∩ disease =', len(ctd_inter))
print('  core in CTD-only∩:', [g for g in core if g in ctd_inter])
print('  ext in CTD-only∩:', [g for g in ext if g in ctd_inter])
print('  sig5 in CTD-only∩:', [g for g in sig5 if g in ctd_inter])

out = {
    'universe': N, 'compound_K': K, 'disease_n': n, 'disease_in_universe': len(dis & uni), 'observed': obs,
    'expected': round(float(exp), 1), 'fold': round(float(fold), 3),
    'hypergeom_p': float(p_hyper),
    'permutation': {'mean': round(float(counts.mean()), 1), 'sd': round(float(counts.std()), 1),
                    'z': round(float(z), 1), 'p': p_perm, 'n': 10000},
    'stratified': {
        'CTD_only': {'K': len(ctd_only), 'observed': len(ctd_only & set(dis)),
                     'fold': round(len(ctd_only & set(dis)) / (len(ctd_only) * n / N), 2),
                     'p': float(stats.hypergeom.sf(len(ctd_only & set(dis)) - 1, N, len(ctd_only), n))},
        'predicted_not_CTD': {'K': len(pred_not_ctd), 'observed': len(pred_not_ctd & set(dis)),
                              'fold': round(len(pred_not_ctd & set(dis)) / (len(pred_not_ctd) * n / N), 2)},
    },
    'ctd_only_retains': {'core': [g for g in core if g in ctd_inter], 'ext': [g for g in ext if g in ctd_inter],
                         'sig5': [g for g in sig5 if g in ctd_inter]},
}
json.dump(out, open(os.path.join(BASE, 'data', 'specificity_controls.json'), 'w', encoding='utf-8'), indent=1)
print('saved run10_specificity_controls.json')
print('B1 DONE')
