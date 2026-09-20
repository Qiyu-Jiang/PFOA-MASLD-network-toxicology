
NOTE (2026-09-20): superseded intermediate version retained for transparency. The fixed-threshold percolation variant of EPC degenerates on this dense network (mass ties) and is NOT part of the final consensus; the final analysis uses the official cytoHubba plugin EPC, replicated in run_09_epc_official_replication.py (Spearman rho = 0.9965; top-10 overlap 8/10). The five-algorithm subset below reproduces the same strict/extended core (see data/hub_consensus_top10_six_algorithms.csv).
# -*- coding: utf-8 -*-
"""run_06 wrap-up: recompute the five-algorithm consensus + corrected figure (EPC demoted to an explanatory panel)"""
import csv, os, math
from collections import Counter
import numpy as np
import networkx as nx

BASE = r'C:\Users\user\workspace\PFOA_NAFLD_project'
ANA = os.path.join(BASE, '04_analysis'); FIG = os.path.join(BASE, '05_figures')

# read scores for all six algorithms
rows = list(csv.DictReader(open(os.path.join(ANA, 'run06_hub_scores.csv'), encoding='utf-8-sig')))
nodes = [r['Node'] for r in rows]
sc = {}
for alg in ['MCC', 'EPC', 'MNC', 'Degree', 'Closeness', 'Radiality']:
    sc[alg] = {r['Node']: float(r[alg]) for r in rows}

FIVE = ['MCC', 'MNC', 'Degree', 'Closeness', 'Radiality']
top10_5 = {a: sorted(nodes, key=lambda v: -sc[a][v])[:10] for a in FIVE}

in_top = Counter()
for a in FIVE:
    for v in top10_5[a]:
        in_top[v] += 1

core = sorted([v for v, c in in_top.items() if c >= 4], key=lambda v: -in_top[v])
five_of_five = [v for v, c in in_top.items() if c == 5]
print('5/5 core:', five_of_five)
print('>=4/5 core:', core)

# save core targets (final version)
with open(os.path.join(ANA, 'run06_core_targets.csv'), 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['Node', 'Votes_of_5', 'Tier', 'Algorithms'])
    for v in core:
        where = [a for a in FIVE if v in top10_5[a]]
        tier = 'strict-core (5/5)' if in_top[v] == 5 else 'extended-core (>=4/5)'
        w.writerow([v, in_top[v], tier, '+'.join(where)])
print('saved run06_core_targets.csv (final)')

# network and layout
edges = []
with open(os.path.join(ANA, 'ppi_edge.csv'), encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        edges.append((row['Source'], row['Target']))
G = nx.Graph(); G.add_edges_from(edges)
pos = nx.spring_layout(G, k=0.06, iterations=60, seed=42)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'sans-serif']
cmap = plt.get_cmap('YlOrRd')

fig, axes = plt.subplots(2, 3, figsize=(19, 12.5), dpi=110)
fig.patch.set_facecolor('#f7f8f9')
order = ['MCC', 'MNC', 'Degree', 'Closeness', 'Radiality']
for ax, name in zip(axes.ravel()[:5], order):
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.08, width=0.4, edge_color='#8a949c')
    base_sizes = [4 + 0.5 * math.sqrt(G.degree(g)) for g in G.nodes()]
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=base_sizes, node_color='#c9d2d6', linewidths=0)
    lst = top10_5[name]
    size2 = [10 + 2.2 * math.sqrt(G.degree(g)) for g in lst]
    cols = [cmap(0.15 + 0.8 * i / 9) for i in range(len(lst))][::-1]
    nx.draw_networkx_nodes(G, pos, nodelist=lst, ax=ax, node_size=size2, node_color=cols, linewidths=0)
    for g in lst:
        x, y = pos[g]
        ax.annotate(g, (x, y), fontsize=8, color='#15191d', xytext=(3, 3), textcoords='offset points')
    ax.set_title(f'{name} — Top 10', fontsize=12)
    ax.axis('off')

# sixth panel: EPC explanation + consensus bars
ax = axes.ravel()[5]
ax.axis('off')
txt = ("EPC method note\n\n"
       "On this network (457 nodes / 9,896 edges, mean degree 43),\n"
       "edge-percolation simulations produce mass ties:\n"
       "  t=0.5: 359/457 nodes tie at the top score\n"
       "  t=0.7: 279/457; t=0.8: 214/457\n"
       "-> node ranks cannot be distinguished; excluded from the consensus.\n\n"
       "(The published method works well on sparse networks;\n"
       "run the official EPC in Cytoscape if an independent check is required.)\n\n"
       "Core target consensus (five-algorithm Top10):\n"
       "  5/5 = " + ', '.join(five_of_five) + "\n"
       "  ≥4/5 = " + ', '.join(core))
ax.text(0.02, 0.98, txt, va='top', fontsize=10.5, family='Microsoft YaHei')
fig.suptitle('PFOA × NAFLD intersection network — hub screening (five robust algorithms, top 10 each)', fontsize=15, y=0.995)
fig.tight_layout()
p1 = os.path.join(FIG, 'run06_hub_panels.png')
fig.savefig(p1, facecolor=fig.get_facecolor())
print('saved:', p1, os.path.getsize(p1), 'bytes')

# consensus bar chart (five-algorithm definition)
fig2, ax2 = plt.subplots(figsize=(9, 5.5), dpi=200)
cc = Counter(in_top.values())
bars = sorted(cc.items(), reverse=True)
ax2.bar([f'{k}/5 algorithms' for k, _ in bars], [v for _, v in bars], color='#145c4d')
ax2.set_ylabel('number of genes')
ax2.set_title('Consensus strength of hub candidates (five-algorithm votes)', fontsize=12)
for s in ['top', 'right']:
    ax2.spines[s].set_visible(False)
fig2.tight_layout()
p2 = os.path.join(FIG, 'run06_consensus_bar.png')
fig2.savefig(p2, facecolor='white')
print('saved:', p2)
print('FINALIZE-DONE')
