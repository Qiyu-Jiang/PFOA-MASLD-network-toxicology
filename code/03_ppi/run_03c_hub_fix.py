# -*- coding: utf-8 -*-
"""run_03c: compute additional hub metrics (eigenvector on the largest connected component) + preview plot"""
import csv, os, math
from collections import Counter, defaultdict
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = r'C:\Users\user\workspace\PFOA_NAFLD_project'
ANA = os.path.join(BASE, '04_analysis')

edges = []
deg = Counter(); wdeg = defaultdict(float)
with open(os.path.join(ANA, 'ppi_edge.csv'), encoding='utf-8-sig') as f:
    for row in csv.DictReader(f):
        a, b, s = row['Source'], row['Target'], float(row['Score'])
        edges.append((a, b, s)); deg[a] += 1; deg[b] += 1; wdeg[a] += s; wdeg[b] += s

genes = [r['Node'] for r in csv.DictReader(open(os.path.join(ANA, 'ppi_node.csv'), encoding='utf-8-sig'))]

G = nx.Graph(); G.add_nodes_from(genes); G.add_weighted_edges_from(edges)
comps = sorted(nx.connected_components(G), key=len, reverse=True)
GC = G.subgraph(comps[0]).copy()
print('largest component:', len(GC), 'nodes')

btw = nx.betweenness_centrality(G, k=min(150, len(G)), seed=42)
clo = nx.closeness_centrality(G)
eig_lcc = nx.eigenvector_centrality_numpy(GC)
eig = {g: eig_lcc.get(g, 0.0) for g in G.nodes()}

hub_rows = [(g, deg.get(g, 0), round(btw.get(g, 0), 4), round(clo.get(g, 0), 4), round(eig.get(g, 0), 4)) for g in G.nodes()]
hub_out = os.path.join(ANA, 'ppi_hub_preview.csv')
with open(hub_out, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f); w.writerow(['Node', 'Degree', 'Betweenness', 'Closeness', 'Eigenvector'])
    for row in sorted(hub_rows, key=lambda x: -x[1]):
        w.writerow(row)
print('saved:', hub_out)
print('--- top 15 degree ---')
for g, d in deg.most_common(15):
    print(f'  {g:10s} deg={d:4d} btw={btw.get(g,0):.4f} clo={clo.get(g,0):.4f} eig={eig.get(g,0):.4f}')
print('--- top 10 betweenness ---')
for g, v in sorted(btw.items(), key=lambda kv: -kv[1])[:10]:
    print(f'  {g:10s} btw={v:.4f} deg={deg.get(g,0)}')

# preview plot
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'sans-serif']
fig, ax = plt.subplots(figsize=(16, 9), dpi=110)
fig.patch.set_facecolor('#f7f8f9')
pos = nx.spring_layout(G, k=0.06, iterations=60, seed=42, weight='weight')
nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.12, width=0.5, edge_color='#7d8790')
sizes = [8 + 2.4 * math.sqrt(deg.get(g, 0)) for g in G.nodes()]
cols = ['#145c4d' if deg.get(g, 0) >= 130 else '#9db3ac' for g in G.nodes()]
nx.draw_networkx_nodes(G, pos, ax=ax, node_size=sizes, node_color=cols, linewidths=0)
for g in [x for x, _ in deg.most_common(20)]:
    x, y = pos[g]
    ax.annotate(g, (x, y), fontsize=8.5, color='#15191d', xytext=(3, 3), textcoords='offset points')
ax.set_title('PFOA × NAFLD intersection PPI network (STRING, 457 nodes / 9,896 edges, confidence ≥ 0.4; top 20 labeled)',
             fontsize=13, color='#15191d')
ax.axis('off'); fig.tight_layout()
out_png = os.path.join(BASE, '05_figures', 'ppi_network_preview.png')
fig.savefig(out_png, dpi=110, facecolor=fig.get_facecolor())
print('saved:', out_png, os.path.getsize(out_png), 'bytes')
