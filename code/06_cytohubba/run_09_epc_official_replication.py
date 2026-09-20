# -*- coding: utf-8 -*-
"""run_09 (save edition): official cytoHubba EPC adaptive-threshold replication + comparison outputs."""
import csv
import json
import os
import sys
import time

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

import numpy as np

BASE = r'C:\Users\user\workspace\PFOA_NAFLD_project'
ANA = os.path.join(BASE, '04_analysis')
OFF = os.path.join(ANA, 'cytohubba_official_scores.csv')
EDGES = os.path.join(ANA, 'ppi_edge.csv')

off = {}
with open(OFF, encoding='utf-8-sig') as f:
    rd = csv.reader(f)
    hdr = next(rd)
    ei = hdr.index('EPC')
    for r in rd:
        if r and r[0].strip():
            try:
                off[r[0].strip()] = float(r[ei])
            except Exception:
                pass

edges = []
with open(EDGES, encoding='utf-8-sig') as f:
    rd = csv.reader(f)
    next(rd)
    for r in rd:
        if len(r) >= 2 and r[0].strip() and r[1].strip():
            edges.append((r[0].strip(), r[1].strip()))

nodes = sorted({x for e in edges for x in e})
idx = {v: i for i, v in enumerate(nodes)}
N = len(nodes)
E = len(edges)
eu = np.array([idx[a] for a, b in edges], dtype=np.int32)
ev = np.array([idx[b] for a, b in edges], dtype=np.int32)
print('N=%d E=%d' % (N, E), flush=True)


def comp_sizes(keep):
    parent = np.arange(N, dtype=np.int32)
    size = np.ones(N, dtype=np.int32)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in zip(eu[keep], ev[keep]):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
            size[ra] += size[rb]
    out = np.empty(N, dtype=np.int32)
    for v in range(N):
        out[v] = size[find(v)]
    return out


rng = np.random.default_rng(20260916)

t0 = time.time()
best_vi, best_std = None, -1.0
curve = []
for i in range(50):
    vi_t = 0.5 + 0.01 * i
    mx = np.empty(30)
    for it in range(30):
        keep = rng.random(E) > vi_t
        mx[it] = comp_sizes(keep).max()
    mean = mx.mean()
    std = float(((mx - mean) ** 2).sum())
    curve.append((vi_t, mean, std))
    if std > best_std:
        best_std = std
        best_vi = vi_t
print('chosen vi_threshold = %.4f  [%.0fs]' % (best_vi, time.time() - t0), flush=True)

t1 = time.time()
acc = np.zeros(N, dtype=np.float64)
for k in range(1000):
    keep = rng.random(E) > best_vi
    acc += comp_sizes(keep)
epc = acc / 1000.0
print('EPC 1000 trials done [%.0fs]' % (time.time() - t1), flush=True)

sim = {nodes[i]: float(epc[i]) for i in range(N)}
top_sim = sorted(nodes, key=lambda g: -sim[g])[:10]
top_off = [g for g, _ in sorted(off.items(), key=lambda kv: (-kv[1], kv[0]))[:10]]
overlap = len(set(top_sim) & set(top_off))
common = [g for g in nodes if g in off]
rs = {g: i for i, g in enumerate(sorted(common, key=lambda g: -sim[g]))}
ro = {g: i for i, g in enumerate(sorted(common, key=lambda g: -off[g]))}
n = len(common)
d2 = sum((rs[g] - ro[g]) ** 2 for g in common)
rho = 1 - 6 * d2 / (n * (n * n - 1))
print('sim top10 :', ', '.join(top_sim))
print('off top10 :', ', '.join(top_off))
print('overlap   : %d/10' % overlap)
print('spearman rho over %d common nodes: %.4f' % (n, rho))

with open(os.path.join(ANA, 'run09_epc_replication_scores.csv'), 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['Node', 'EPC_replication', 'EPC_official'])
    for g in nodes:
        w.writerow([g, round(sim[g], 4), off.get(g, '')])
json.dump({'vi_threshold': best_vi, 'spearman_rho': round(rho, 4), 'top10_overlap': overlap,
           'sim_top10': top_sim, 'official_top10': top_off, 'common_nodes': n},
          open(os.path.join(ANA, 'run09_epc_comparison.json'), 'w', encoding='utf-8'), indent=1)
print('saved run09 outputs')
print('RUN09 DONE')
