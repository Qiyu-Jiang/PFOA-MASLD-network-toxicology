# -*- coding: utf-8 -*-
"""run_02b: merge the PFOA target pool (6 sources) and disease pool (5 sources) and compute their intersection"""
import csv, os, re
from collections import defaultdict, Counter

BASE = r'C:\Users\user\workspace\PFOA_NAFLD_project'
RAW = os.path.join(BASE, '02_data_raw')
KEEP = os.path.join(BASE, '03_data_processed')
ANA = os.path.join(BASE, '04_analysis')

def up(s):
    return (s or '').strip().upper()

# ================= compound side (PFOA target pool) =================
comp = defaultdict(lambda: {'sources': set(), 'notes': {}})

# 1) CTD-human
for r in csv.DictReader(open(os.path.join(KEEP, 'PFOA_targets_CTD_human.tsv'), encoding='utf-8-sig'), delimiter='\t'):
    s = up(r['GeneSymbol'])
    if s:
        comp[s]['sources'].add('CTD')
        comp[s]['notes']['CTD_rows'] = int(r['EvidenceRows'])
        comp[s]['notes']['CTD_inc'] = int(r['IncExpression']); comp[s]['notes']['CTD_dec'] = int(r['DecExpression'])

# 2) TargetNet
for r in csv.DictReader(open(os.path.join(KEEP, 'TargetNet_targets_symbols.tsv'), encoding='utf-8-sig'), delimiter='\t'):
    s = up(r['Symbol'])
    if s:
        comp[s]['sources'].add('TargetNet')
        try: comp[s]['notes']['TargetNet_prob'] = float(r['Prob'])
        except ValueError: pass

# 3) SwissTarget
for r in csv.DictReader(open(os.path.join(KEEP, 'PFOA_targets_SwissTarget.tsv'), encoding='utf-8-sig'), delimiter='\t'):
    s = up(r['GeneSymbol'])
    if s:
        comp[s]['sources'].add('SwissTarget')

# 4) PharmMapper
for r in csv.DictReader(open(os.path.join(KEEP, 'PharmMapper_targets_symbols.tsv'), encoding='utf-8-sig'), delimiter='\t'):
    s = up(r['Symbol'])
    if s:
        comp[s]['sources'].add('PharmMapper')
        try:
            nf = float(r['Norm_Fit'])
            if nf > comp[s]['notes'].get('PharmMapper_NormFit', -1):
                comp[s]['notes']['PharmMapper_NormFit'] = nf
        except ValueError: pass

# 5) SEA
comp['TP53']['sources'].add('SEA'); comp['TP53']['notes']['SEA'] = 'p53 Z=-0.41(ns)'
# SuperPred: not retrieved (network unreachable); ProTox: no gene list -- neither is counted

comp_out = os.path.join(KEEP, 'PFOA_target_pool_v1.tsv')
with open(comp_out, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['Symbol', 'N_Sources', 'Sources', 'CTD_rows', 'TargetNet_prob', 'PharmMapper_NormFit', 'CTD_IncExpr', 'CTD_DecExpr'])
    for s, d in sorted(comp.items(), key=lambda kv: -len(kv[1]['sources'])):
        n = d['notes']
        w.writerow([s, len(d['sources']), '+'.join(sorted(d['sources'])),
                    n.get('CTD_rows', ''), n.get('TargetNet_prob', ''), n.get('PharmMapper_NormFit', ''),
                    n.get('CTD_inc', ''), n.get('CTD_dec', '')])
print('compound pool:', len(comp), '->', comp_out)
print('  source sizes:', {k: sum(1 for d in comp.values() if k in d['sources']) for k in ['CTD','TargetNet','SwissTarget','PharmMapper','SEA']})

# ================= disease side (NAFLD disease pool) =================
dis = defaultdict(lambda: {'sources': set(), 'notes': {}})

# 1) gene-level DEG from GSE89632 NASH_vs_HC
deg_best = {}
for r in csv.DictReader(open(os.path.join(ANA, 'GSE89632_NASH_vs_HC_topTable.csv'), encoding='utf-8-sig')):
    s = up(r.get('Symbol'))
    if not s: continue
    try:
        lf, ap = float(r['logFC']), float(r['adj.P.Val'])
    except (ValueError, KeyError): continue
    if abs(lf) > 1 and ap < 0.05:
        if s not in deg_best or abs(lf) > abs(deg_best[s][0]):
            deg_best[s] = (lf, ap)
for s, (lf, ap) in deg_best.items():
    dis[s]['sources'].add('DEG')
    dis[s]['notes']['DEG_logFC'] = round(lf, 3); dis[s]['notes']['DEG_adjP'] = '%.2g' % ap
print('disease DEG set:', len(deg_best))

# 2) CTD-NAFLD (50-row parsed version)
ctd_d = os.path.join(RAW, 'CTD_NAFLD_genes.csv')
for r in csv.DictReader(open(ctd_d, encoding='utf-8-sig')):
    vals = list(r.values())
    s = up(vals[1]) if len(vals) > 1 else ''
    if s and s not in ('GENE',):
        dis[s]['sources'].add('CTD-NAFLD')

# 3) DisGeNET top30
dg = os.path.join(RAW, 'DisGeNET_NAFLD_top30.tsv')
rows = list(csv.reader(open(dg, encoding='utf-8-sig'), delimiter='\t'))
hdr_i = next(i for i, r in enumerate(rows) if r and 'Gene' in r)
gi = rows[hdr_i].index('Gene')
for r in rows[hdr_i+1:]:
    if len(r) > gi:
        s = up(r[gi])
        if s:
            dis[s]['sources'].add('DisGeNET')

# 4/5) GeneCards NAFLD + MASLD (Relevance > 5)
for tag in ['NAFLD', 'MASLD']:
    p = os.path.join(RAW, f'GeneCards_{tag}.csv')
    for r in csv.reader(open(p, encoding='utf-8-sig')):
        if len(r) >= 5 and r[0] and r[0] != 'Symbol':
            try: rel = float(r[3])
            except ValueError: continue
            if rel > 5:
                s = up(r[0])
                dis[s]['sources'].add('GeneCards')
                dis[s]['notes'][f'GC_{tag}_rel'] = round(rel, 2)

dis_out = os.path.join(KEEP, 'NAFLD_disease_pool_v1.tsv')
with open(dis_out, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['Symbol', 'N_Sources', 'Sources', 'DEG_logFC', 'DEG_adjP', 'GC_NAFLD_rel', 'GC_MASLD_rel', 'DisGeNET', 'CTD_NAFLD'])
    for s, d in sorted(dis.items(), key=lambda kv: -len(kv[1]['sources'])):
        n = d['notes']
        w.writerow([s, len(d['sources']), '+'.join(sorted(d['sources'])),
                    n.get('DEG_logFC', ''), n.get('DEG_adjP', ''), n.get('GC_NAFLD_rel', ''), n.get('GC_MASLD_rel', ''),
                    'Y' if 'DisGeNET' in d['sources'] else '', 'Y' if 'CTD-NAFLD' in d['sources'] else ''])
print('disease pool:', len(dis), '->', dis_out)
print('  source sizes:', {k: sum(1 for d in dis.values() if k in d['sources']) for k in ['DEG','GeneCards','DisGeNET','CTD-NAFLD']})

# ================= intersection =================
inter = sorted(set(comp) & set(dis))
print()
print('=' * 70)
print('PFOA×NAFLD intersection:', len(inter))

rows_out = []
for s in inter:
    c, d = comp[s], dis[s]
    ncs, nds = len(c['sources']), len(d['sources'])
    rows_out.append({
        'Symbol': s, 'N_comp': ncs, 'Comp_sources': '+'.join(sorted(c['sources'])),
        'N_dis': nds, 'Dis_sources': '+'.join(sorted(d['sources'])),
        'DEG_logFC': c and d['notes'].get('DEG_logFC', ''), 'DEG_adjP': d['notes'].get('DEG_adjP', ''),
        'CTD_rows': c['notes'].get('CTD_rows', ''), 'TargetNet_prob': c['notes'].get('TargetNet_prob', ''),
        'PharmMapper_NormFit': c['notes'].get('PharmMapper_NormFit', ''),
        'GC_NAFLD_rel': d['notes'].get('GC_NAFLD_rel', ''), 'GC_MASLD_rel': d['notes'].get('GC_MASLD_rel', ''),
    })
# high confidence: >=2 compound sources and >=2 disease sources
hiconf = [r for r in rows_out if r['N_comp'] >= 2 and r['N_dis'] >= 2]
print('high-confidence (>=2 compound & >=2 disease sources):', len(hiconf))

out = os.path.join(ANA, 'PFOAxNAFLD_intersection_v2.csv')
with open(out, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()) if rows_out else [])
    w.writeheader()
    for r in sorted(rows_out, key=lambda x: (-(x['N_comp'] + x['N_dis']), x['Symbol'])):
        w.writerow(r)
print('saved:', out)

print()
print('--- top 30 by combined evidence ---')
for r in sorted(rows_out, key=lambda x: (-(x['N_comp'] + x['N_dis']), x['Symbol']))[:30]:
    extra = r['DEG_logFC'] if r['DEG_logFC'] != '' else ''
    print(f"  {r['Symbol']:10s} comp={r['N_comp']}({r['Comp_sources'][:28]:28s}) dis={r['N_dis']}({r['Dis_sources'][:30]:30s}) DEG={extra}")

print()
print('--- high-confidence list ---')
print(', '.join(r['Symbol'] for r in hiconf))
