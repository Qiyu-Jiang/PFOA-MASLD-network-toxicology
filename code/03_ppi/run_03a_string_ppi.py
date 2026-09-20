# -*- coding: utf-8 -*-
"""run_03a: STRING PPI network (479 intersection genes) + Cytoscape node/edge export + metrics"""
import csv, os, io, requests
from collections import Counter, defaultdict

BASE = r'C:\Users\user\workspace\PFOA_NAFLD_project'
ANA = os.path.join(BASE, '04_analysis')
KEEP = os.path.join(BASE, '03_data_processed')

inter = list(csv.DictReader(open(os.path.join(ANA, 'PFOAxNAFLD_intersection_v2.csv'), encoding='utf-8-sig')))
genes = [r['Symbol'] for r in inter if r['Symbol']]
print('genes submitted:', len(genes))

ids = '\r'.join(genes)
r = requests.post(
    'https://string-db.org/api/tsv/network',
    data={'identifiers': ids, 'species': 9606, 'required_score': 400,
          'caller_identity': 'pfoa_nafld_study'},
    timeout=180)
print('STRING status:', r.status_code, '| response len:', len(r.text))
lines = r.text.strip().splitlines()
print('header:', lines[0] if lines else '(empty)')
print('edge rows:', len(lines) - 1)

if len(lines) < 2:
    raise SystemExit('no edges returned - check input')

hdr = lines[0].split('\t')
iA = hdr.index('preferredName_A'); iB = hdr.index('preferredName_B'); iS = hdr.index('score')
edges = []
for ln in lines[1:]:
    p = ln.split('\t')
    if len(p) > iS:
        try:
            s = float(p[iS])
        except ValueError:
            continue
        edges.append((p[iA], p[iB], s / 1000.0))
print('parsed edges:', len(edges))

deg, wdeg = Counter(), defaultdict(float)
for a, b, s in edges:
    deg[a] += 1; deg[b] += 1
    wdeg[a] += s; wdeg[b] += s

all_nodes = sorted(set(deg) | set(genes))
iso = [g for g in genes if g not in deg]
print('nodes with >=1 edge:', len(deg), '| isolated (no edge):', len(iso), iso[:20])

# export Cytoscape files
edge_out = os.path.join(ANA, 'ppi_edge.csv')
with open(edge_out, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['Source', 'Target', 'Score'])
    for a, b, s in edges:
        w.writerow([a, b, round(s, 3)])

node_out = os.path.join(ANA, 'ppi_node.csv')
with open(node_out, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['Node', 'Degree', 'WeightedDegree'])
    for g in sorted(all_nodes, key=lambda x: -deg.get(x, 0)):
        w.writerow([g, deg.get(g, 0), round(wdeg.get(g, 0), 2)])
print('saved:', edge_out)
print('saved:', node_out)

# connectivity
import networkx as nx
G = nx.Graph()
G.add_nodes_from(genes)
G.add_weighted_edges_from(edges)
comps = sorted(nx.connected_components(G), key=len, reverse=True)
print('connected components:', len(comps), '| largest:', len(comps[0]), '| 2nd:', len(comps[1]) if len(comps) > 1 else '-')

print()
print('--- top 25 by degree ---')
for g, d in deg.most_common(25):
    print(f'  {g:10s} degree={d:4d}  wdeg={wdeg[g]:.1f}')

# overview of the high-confidence-edge (>= 0.7) subnetwork
E7 = [(a, b, s) for a, b, s in edges if s >= 0.7]
d7 = Counter()
for a, b, s in E7:
    d7[a] += 1; d7[b] += 1
print()
print('edges >=0.7:', len(E7), '| nodes >=0.7:', len(d7))
print('top 15 by degree(>=0.7):', ' '.join(f'{g}({d})' for g, d in d7.most_common(15)))
