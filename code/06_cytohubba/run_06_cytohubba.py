# -*- coding: utf-8 -*-
"""run_06: full in-house implementation of six cytoHubba algorithms (MCC/EPC/MNC/Degree/Closeness/Radiality)
Implemented from the original equations in Chin et al. 2014 (BMC Syst Biol).
Outputs: hub-score table, per-algorithm Top10, intersection core targets, six-panel figure.
"""
import csv, os, time, math
from collections import Counter, defaultdict
import numpy as np
import networkx as nx

BASE = r'C:\Users\user\workspace\PFOA_NAFLD_project'
ANA = os.path.join(BASE, '04_analysis'); FIG = os.path.join(BASE, '05_figures')

# ---------- load network ----------
edges = []
with open(os.path.join(ANA, 'ppi_edge.csv'), encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        edges.append((row['Source'], row['Target'], float(row['Score'])))
G = nx.Graph()
for a, b, s in edges:
    G.add_edge(a, b)
nodes = sorted(G.nodes())
N = len(nodes)
print(f'network: {N} nodes, {G.number_of_edges()} edges')

deg = dict(G.degree())

t0 = time.time()
# ---------- 1) Degree ----------
print('degree done')

# ---------- 2) MNC ----------
mnc = {}
for v in nodes:
    nbrs = list(G.neighbors(v))
    if not nbrs:
        mnc[v] = 0; continue
    sub = G.subgraph(nbrs)
    comps = nx.connected_components(sub)
    mnc[v] = max(len(c) for c in comps)
print(f'MNC done ({time.time()-t0:.1f}s)')

# ---------- 3) Closeness (enhanced: component-fraction-weighted sum of harmonic distances within component) ----------
def comp_map():
    cm = {}
    for i, comp in enumerate(nx.connected_components(G)):
        for v in comp:
            cm[v] = i
    return cm
cmap = comp_map()
comp_sizes = Counter(cmap.values())
clo = {}
for v in nodes:
    # BFS distances within component
    dist = nx.single_source_shortest_path_length(G, v)
    s = sum(1.0 / d for w, d in dist.items() if d > 0)
    clo[v] = s * comp_sizes[cmap[v]] / N
print(f'Closeness done ({time.time()-t0:.1f}s)')

# ---------- 4) Radiality (enhanced form from the paper) ----------
rad = {}
# within-component diameter
comp_diam = {}
for cid in comp_sizes:
    comp_nodes = [v for v in nodes if cmap[v] == cid]
    if len(comp_nodes) <= 1:
        comp_diam[cid] = 0
    else:
        # diameter via eccentricity
        sub = G.subgraph(comp_nodes)
        comp_diam[cid] = nx.diameter(sub) if len(comp_nodes) > 1 else 0
print('comp diameters:', {k: v for k, v in list(comp_diam.items())[:5]}, '... total comps:', len(comp_diam))
for v in nodes:
    cid = cmap[v]
    comp_nodes = [x for x in nodes if cmap[x] == cid]
    if len(comp_nodes) <= 1:
        rad[v] = 0.0; continue
    dist = nx.single_source_shortest_path_length(G, v)
    D = comp_diam[cid]
    s = sum(D + 1 - d for w, d in dist.items())
    maxd = max(dist.values())
    rad[v] = (len(comp_nodes) / N) * s / maxd if maxd > 0 else 0.0
print(f'Radiality done ({time.time()-t0:.1f}s)')

# ---------- 5) MCC (maximum clique weighted enumeration; igraph preferred) ----------
mcc = defaultdict(float)
backend = None
try:
    import igraph as ig
    backend = 'igraph'
    name2idx = {v: i for i, v in enumerate(nodes)}
    g = ig.Graph(n=N, edges=[(name2idx[a], name2idx[b]) for a, b, _ in edges])
    t1 = time.time()
    cliques = g.maximal_cliques(min=2)
    print(f'igraph maximal cliques (>=2): {len(cliques)} in {time.time()-t1:.1f}s')
    for cl in cliques:
        w = math.factorial(len(cl) - 1)
        for idx in cl:
            mcc[nodes[idx]] += w
except Exception as e:
    print('igraph unavailable/failed:', repr(e)[:200], '-> fallback networkx (capped)')
    backend = 'networkx'
    t1 = time.time()
    cnt = 0
    for cl in nx.find_cliques(G):
        cnt += 1
        if len(cl) >= 2:
            w = math.factorial(len(cl) - 1)
            for v in cl:
                mcc[v] += w
        if cnt % 500000 == 0:
            print(f'  cliques so far: {cnt} ({time.time()-t1:.0f}s)')
        if time.time() - t1 > 600:
            print('  !! clique enumeration exceeded 600s, stopping early (partial)')
            break
    print(f'networkx cliques enumerated: {cnt} in {time.time()-t1:.1f}s')
mcc = {v: mcc.get(v, 0.0) for v in nodes}
print(f'MCC done ({time.time()-t0:.1f}s), backend={backend}')

# ---------- 6) EPC (threshold percolation simulation) ----------
def epc_scores(threshold, reps=1000, seed=42):
    rng = np.random.default_rng(seed)
    idx = {v: i for i, v in enumerate(nodes)}
    eu = np.array([idx[a] for a, b, _ in edges], dtype=np.int32)
    ev = np.array([idx[b] for a, b, _ in edges], dtype=np.int32)
    acc = np.zeros(N, dtype=np.float64)
    parent = np.arange(N, dtype=np.int32)
    size = np.zeros(N, dtype=np.int32)

    def find(x, parent):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for k in range(reps):
        r = rng.random(len(edges))
        keep = r >= threshold
        parent[:] = np.arange(N)
        size[:] = 1
        for a, b in zip(eu[keep], ev[keep]):
            ra, rb = find(a, parent), find(b, parent)
            if ra != rb:
                parent[rb] = ra
                size[ra] += size[rb]
        # component size of each node
        for v in range(N):
            acc[v] += size[find(v, parent)]
    return {nodes[i]: acc[i] / reps / N for i in range(N)}

t2 = time.time()
epc = epc_scores(0.5, reps=1000, seed=42)
print(f'EPC(t=0.5) done ({time.time()-t2:.1f}s)')

# sensitivity: alternative thresholds (fewer reps)
epc_sens = {}
for t in [0.2, 0.3, 0.4, 0.6]:
    epc_sens[t] = epc_scores(t, reps=300, seed=42)

# ---------- summary ----------
algos = {'MCC': mcc, 'EPC': epc, 'MNC': mnc, 'Degree': deg, 'Closeness': clo, 'Radiality': rad}
ranked = {}
for name, sc in algos.items():
    ranked[name] = sorted(nodes, key=lambda v: -sc.get(v, 0))

top10 = {name: ranked[name][:10] for name in ranked}
print()
for name, lst in top10.items():
    print(f'Top10 {name}:', ', '.join(lst))

# intersection (appearing in the top10 of >= k algorithms)
in_top = Counter()
for name in top10:
    for v in top10[name]:
        in_top[v] += 1
print()
print('frequency distribution:', Counter(in_top.values()))
core = [v for v, c in in_top.items() if c >= 3]
core_sorted = sorted(core, key=lambda v: -in_top[v])
print(f'core targets (>=3 algorithms): {len(core_sorted)}')
print(', '.join(core_sorted))
strict = [v for v, c in in_top.items() if c >= 6]
print('in ALL 6:', ', '.join(strict) if strict else '(none)')

# save
with open(os.path.join(ANA, 'run06_hub_scores.csv'), 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['Node', 'MCC', 'EPC', 'MNC', 'Degree', 'Closeness', 'Radiality',
                'Rank_MCC', 'Rank_EPC', 'Rank_MNC', 'Rank_Degree', 'Rank_Closeness', 'Rank_Radiality'])
    for v in nodes:
        rk = {name: ranked[name].index(v) + 1 for name in ranked}
        w.writerow([v, round(mcc[v], 4), round(epc[v], 6), mnc[v], deg[v], round(clo[v], 5), round(rad[v], 4),
                    rk['MCC'], rk['EPC'], rk['MNC'], rk['Degree'], rk['Closeness'], rk['Radiality']])
print('saved run06_hub_scores.csv')

with open(os.path.join(ANA, 'run06_top10_per_algorithm.csv'), 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['Algorithm', 'Rank', 'Node', 'Score'])
    for name, lst in top10.items():
        for i, v in enumerate(lst, 1):
            w.writerow([name, i, v, round(algos[name][v], 6)])
print('saved run06_top10_per_algorithm.csv')

with open(os.path.join(ANA, 'run06_core_targets.csv'), 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['Node', 'In_Top10_Count', 'Algorithms'])
    for v in core_sorted:
        where = [name for name in top10 if v in top10[name]]
        w.writerow([v, in_top[v], '+'.join(where)])
print('saved run06_core_targets.csv')

# EPC sensitivity comparison (top10 stability)
print()
print('EPC top10 stability across thresholds:')
for t in [0.2, 0.3, 0.4, 0.5, 0.6]:
    sc = epc if t == 0.5 else epc_sens[t]
    lst = sorted(nodes, key=lambda v: -sc.get(v, 0))[:10]
    overlap = len(set(lst) & set(top10['EPC']))
    print(f'  t={t}: overlap-with-main {overlap}/10 | top3: {lst[:3]}')

# ---------- figure: 2x3 panels ----------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'sans-serif']
pos = nx.spring_layout(G, k=0.06, iterations=60, seed=42, weight=None)
fig, axes = plt.subplots(2, 3, figsize=(19, 12.5), dpi=110)
fig.patch.set_facecolor('#f7f8f9')
import matplotlib.colors as mcolors
cmap = plt.get_cmap('YlOrRd')
order = ['MCC', 'EPC', 'MNC', 'Degree', 'Closeness', 'Radiality']
for ax, name in zip(axes.ravel(), order):
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.08, width=0.4, edge_color='#8a949c')
    base_sizes = [4 + 0.5 * math.sqrt(deg.get(g, 0)) for g in G.nodes()]
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=base_sizes, node_color='#c9d2d6', linewidths=0)
    lst = top10[name]
    size2 = [10 + 2.2 * math.sqrt(deg.get(g, 0)) for g in lst]
    cols = [cmap(0.15 + 0.8 * i / 9) for i in range(len(lst))][::-1]  # higher rank -> darker colour
    nx.draw_networkx_nodes(G, pos, nodelist=lst, ax=ax, node_size=size2, node_color=cols, linewidths=0)
    for i, g in enumerate(lst):
        x, y = pos[g]
        ax.annotate(g, (x, y), fontsize=8, color='#15191d', xytext=(3, 3), textcoords='offset points')
    ax.set_title(f'{name} — Top 10', fontsize=12)
    ax.axis('off')
fig.suptitle('PFOA × NAFLD intersection network — cytoHubba six algorithms (top 10 highlighted)', fontsize=15, y=0.995)
fig.tight_layout()
p = os.path.join(FIG, 'run06_cytohubba_panels.png')
fig.savefig(p, facecolor=fig.get_facecolor())
print('saved:', p, os.path.getsize(p), 'bytes')

# core intersection figure
fig2, ax2 = plt.subplots(figsize=(9, 5.5), dpi=200)
cc = Counter(in_top.values())
bars = sorted(cc.items())
ax2.bar([f'≥{k}' for k, _ in bars], [v for _, v in bars], color='#145c4d')
ax2.set_xlabel('appears in top10 of ≥k algorithms'); ax2.set_ylabel('nodes')
ax2.set_title('Consensus strength of hub candidates', fontsize=12)
for s in ['top', 'right']:
    ax2.spines[s].set_visible(False)
fig2.tight_layout()
p2 = os.path.join(FIG, 'run06_consensus_bar.png')
fig2.savefig(p2, facecolor='white')
print('saved:', p2)
print('RUN06-DONE')
