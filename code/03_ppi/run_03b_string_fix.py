# -*- coding: utf-8 -*-
"""run_03b: STRING PPI corrected version (score not divided by 1000) + metrics + preview plots + hub preview"""
import csv, os, requests, math
from collections import Counter, defaultdict
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = r'C:\Users\user\workspace\PFOA_NAFLD_project'
ANA = os.path.join(BASE, '04_analysis')
KEEP = os.path.join(BASE, '03_data_processed')

inter = list(csv.DictReader(open(os.path.join(ANA, 'PFOAxNAFLD_intersection_v2.csv'), encoding='utf-8-sig')))
genes = [r['Symbol'] for r in inter if r['Symbol']]

r = requests.post('https://string-db.org/api/tsv/network',
    data={'identifiers': '\r'.join(genes), 'species': 9606, 'required_score': 400,
          'caller_identity': 'pfoa_nafld_study'}, timeout=180)
lines = r.text.strip().splitlines()
hdr = lines[0].split('\t')
iA, iB, iS = hdr.index('preferredName_A'), hdr.index('preferredName_B'), hdr.index('score')
edges = []
for ln in lines[1:]:
    p = ln.split('\t')
    try:
        s = float(p[iS])
    except (ValueError, IndexError):
        continue
    edges.append((p[iA], p[iB], s))   # fix: API scores are already on 0-1 scale
print('edges:', len(edges), '| score range:', round(min(e[2] for e in edges), 3), '-', round(max(e[2] for e in edges), 3))

deg, wdeg = Counter(), defaultdict(float)
for a, b, s in edges:
    deg[a] += 1; deg[b] += 1; wdeg[a] += s; wdeg[b] += s

edge_out = os.path.join(ANA, 'ppi_edge.csv')
with open(edge_out, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f); w.writerow(['Source', 'Target', 'Score'])
    for a, b, s in edges: w.writerow([a, b, round(s, 3)])
node_out = os.path.join(ANA, 'ppi_node.csv')
with open(node_out, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f); w.writerow(['Node', 'Degree', 'WeightedDegree'])
    for g in sorted(set(genes) | set(deg), key=lambda x: -deg.get(x, 0)):
        w.writerow([g, deg.get(g, 0), round(wdeg.get(g, 0), 2)])

E7 = [(a, b, s) for a, b, s in edges if s >= 0.7]
d7 = Counter()
for a, b, s in E7: d7[a] += 1; d7[b] += 1
print('edges >=0.7:', len(E7), '| nodes:', len(d7))
E5 = [(a, b, s) for a, b, s in edges if s >= 0.5]
print('edges >=0.5:', len(E5))

G = nx.Graph(); G.add_nodes_from(genes); G.add_weighted_edges_from(edges)
print('components:', len(list(nx.connected_components(G))))

# hub metrics (betweenness uses k-sampling for speed)
btw = nx.betweenness_centrality(G, k=min(150, len(G)), seed=42)
clo = nx.closeness_centrality(G)
eig = nx.eigenvector_centrality_numpy(G)
hub_rows = []
for g in set(genes) | set(deg):
    hub_rows.append((g, deg.get(g, 0), round(btw.get(g, 0), 4), round(clo.get(g, 0), 4), round(eig.get(g, 0), 4)))
hub_out = os.path.join(ANA, 'ppi_hub_preview.csv')
with open(hub_out, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f); w.writerow(['Node', 'Degree', 'Betweenness', 'Closeness', 'Eigenvector'])
    for row in sorted(hub_rows, key=lambda x: -x[1]):
        w.writerow(row)
print('saved:', hub_out)

print()
print('--- top 15 degree ---')
for g, d in deg.most_common(15):
    print(f'  {g:10s} deg={d:4d} btw={btw.get(g,0):.4f} clo={clo.get(g,0):.4f} eig={eig.get(g,0):.3f}')
print('--- top 10 betweenness ---')
for g, v in sorted(btw.items(), key=lambda kv: -kv[1])[:10]:
    print(f'  {g:10s} btw={v:.4f} deg={deg.get(g,0)}')

# ---- preview: full network; edge alpha/width scale with score, node size with degree; label top 20 ----
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'sans-serif']
fig, ax = plt.subplots(figsize=(16, 9), dpi=110)
fig.patch.set_facecolor('#f7f8f9')
pos = nx.spring_layout(G, k=0.06, iterations=60, seed=42, weight='weight')
nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.12, width=0.5, edge_color='#7d8790')
sizes = [8 + 2.2 * math.sqrt(deg.get(g, 0)) for g in G.nodes()]
cols = ['#145c4d' if deg.get(g, 0) >= 130 else '#9db3ac' for g in G.nodes()]
nx.draw_networkx_nodes(G, pos, ax=ax, node_size=sizes, node_color=cols, linewidths=0)
top20 = [g for g, _ in deg.most_common(20)]
for g in top20:
    x, y = pos[g]
    ax.annotate(g, (x, y), fontsize=8.5, color='#15191d', xytext=(3, 3), textcoords='offset points')
ax.set_title('PFOA × NAFLD intersection PPI network (STRING, 457 nodes / 9,896 edges, confidence ≥ 0.4)',
             fontsize=13, color='#15191d')
ax.axis('off')
fig.tight_layout()
out_png = os.path.join(BASE, '05_figures', 'ppi_network_preview.png')
fig.savefig(out_png, dpi=110, facecolor=fig.get_facecolor())
print('saved:', out_png)
