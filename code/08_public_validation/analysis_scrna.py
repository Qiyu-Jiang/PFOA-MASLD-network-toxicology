# -*- coding: utf-8 -*-
"""GSE344087 scRNA analysis: cell-type expression of hub / signature genes.

Reads the 10x h5 files (samples S001-S020 subset), assigns coarse cell types by
marker-score, and computes per-cell-type detection rate and mean expression for
the 14 target genes (hub + extended + signature). Outputs CSV + dotplot PNG.
"""
import os
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

import h5py
import scipy.sparse as sp

BASE = r'C:\Users\user\workspace\PFOA_NAFLD_project'
SC = os.path.join(BASE, '02_data_raw', 'scRNA', 'GSE344087')
OUT = os.path.join(BASE, '04_analysis', 'scrna_GSE344087')
os.makedirs(OUT, exist_ok=True)

SAMPLES = [
    dict(id='S001', fname='GSM9970276_MASLD_S001_LIVER_filtered_feature_bc_matrix.h5', expected=8703180, group='MASLD_fibrosis'),
    dict(id='S003', fname='GSM9970278_MASLD_S003_LIVER_filtered_feature_bc_matrix.h5', expected=3460300, group='MASLD_nofib'),
    dict(id='S005', fname='GSM9970280_MASLD_S005_LIVER_filtered_feature_bc_matrix.h5', expected=9122611, group='MASLD_nofib'),
    dict(id='S009', fname='GSM9970284_MASLD_S009_LIVER_filtered_feature_bc_matrix.h5', expected=1887436, group='MASLD_fibrosis'),
    dict(id='S017', fname='GSM9970292_MASLD_S017_LIVER_filtered_feature_bc_matrix.h5', expected=9227468, group='MASLD_fibrosis'),
    dict(id='S020', fname='GSM9970295_HC_S020_LIVER_filtered_feature_bc_matrix.h5', expected=13601709, group='Healthy'),
]

TARGETS = ['EPHA2', 'FMO1', 'MYC', 'PHLDA1', 'WNT5A',
           'TNF', 'IL6', 'IL1B', 'TP53', 'ALB', 'INS', 'AKT1', 'PPARG', 'EGFR', 'BCL2']  # signature (5) + hub (10)

MARKERS = {
    'Hepatocyte': ['ALB', 'APOA1', 'APOB', 'TTR', 'TF'],
    'Cholangiocyte': ['KRT19', 'EPCAM', 'SOX9'],
    'T/NK': ['CD3D', 'CD3E', 'NKG7', 'GNLY'],
    'B cell': ['MS4A1', 'CD79A'],
    'Myeloid': ['LYZ', 'CD68', 'CD14'],
    'Endothelial': ['PECAM1', 'VWF', 'CLDN5'],
    'Mesenchymal': ['COL1A1', 'COL1A2', 'ACTA2', 'PDGFRB'],
}

CATS = list(MARKERS.keys()) + ['Unassigned']


def load_sample(path):
    with h5py.File(path, 'r') as f:
        m = f['matrix']
        shape = tuple(int(x) for x in m['shape'][:])
        data = m['data'][:].astype(np.float32)
        indices = m['indices'][:].astype(np.int64)
        indptr = m['indptr'][:].astype(np.int64)
        names = [x.decode('utf-8', 'replace') for x in m['features']['name'][:]]
        ftypes = [x.decode('utf-8', 'replace') for x in m['features']['feature_type'][:]]
    mtx = sp.csc_matrix((data, indices, indptr), shape=shape).tocsr()
    keep = np.array([i for i, t in enumerate(ftypes) if t == 'Gene Expression'], dtype=np.int64)
    mtx = mtx[keep, :]
    genes = [names[i] for i in keep]
    return mtx, genes


def zscores(mtx, totals, gi):
    vals = np.asarray(mtx[gi, :].todense()).ravel()
    norm = np.log1p(vals / np.maximum(totals, 1) * 10000.0)
    return (norm - norm.mean()) / (norm.std() + 1e-9)


loaded = []
for s in SAMPLES:
    path = os.path.join(SC, s['fname'])
    if not os.path.exists(path):
        print('SKIP (missing):', s['id'], flush=True)
        continue
    size = os.path.getsize(path)
    if size < s['expected'] * 0.999:
        print('SKIP (incomplete %.2f MB):' % (size / 1e6), s['id'], flush=True)
        continue

    try:
        mtx, genes = load_sample(path)
    except Exception as e:
        print('SKIP (read error: %s):' % str(e)[:90], s['id'], flush=True)
        continue
    ncells = mtx.shape[1]
    totals = np.asarray(mtx.sum(axis=0)).ravel()
    gidx = {}
    wanted = TARGETS + [g for lst in MARKERS.values() for g in lst]
    for want in wanted:
        wu = want.upper()
        if wu in gidx:
            continue
        for i, g in enumerate(genes):
            if g.upper() == wu:
                gidx[wu] = i
                break

    T = np.zeros((len(TARGETS), ncells), dtype=np.float32)
    for ti, g in enumerate(TARGETS):
        gi = gidx.get(g.upper())
        if gi is None:
            continue
        vals = np.asarray(mtx[gi, :].todense()).ravel()
        T[ti] = np.log1p(vals / np.maximum(totals, 1) * 10000.0)

    Z = {}
    for mname, lst in MARKERS.items():
        arrs = []
        for g in lst:
            gi = gidx.get(g.upper())
            if gi is None:
                continue
            arrs.append(zscores(mtx, totals, gi))
        if arrs:
            Z[mname] = np.mean(arrs, axis=0)

    labels = ['Unassigned'] * ncells
    if Z:
        keys = list(Z.keys())
        stack = np.vstack([Z[k] for k in keys])
        best = stack.argmax(axis=0)
        bestv = stack.max(axis=0)
        for c in range(ncells):
            if bestv[c] > 0.3:
                labels[c] = keys[best[c]]

    from collections import Counter
    print('%s (%s): %d cells | %s' % (s['id'], s['group'], ncells, dict(Counter(labels))), flush=True)
    loaded.append((s['id'], s['group'], np.array(labels), T))

# aggregate
rows = []
comp_rows = []
for cat in CATS:
    n_tot = 0
    per_gene = {g: [0, 0.0] for g in TARGETS}  # count positive, sum
    for sid, grp, labels, T in loaded:
        sel = labels == cat
        n = int(sel.sum())
        n_tot += n
        if n == 0:
            continue
        sub = T[:, sel]
        for ti, g in enumerate(TARGETS):
            pos = int((sub[ti] > 0).sum())
            per_gene[g][0] += pos
            per_gene[g][1] += float(sub[ti].sum())
    for g in TARGETS:
        pos, ssum = per_gene[g]
        pct = pos / n_tot if n_tot else 0.0
        mean = ssum / n_tot if n_tot else 0.0
        rows.append((cat, g, n_tot, round(pct, 4), round(mean, 4)))
    comp_rows.append((cat, n_tot))

with open(os.path.join(OUT, 'scrna_celltype_gene_stats.csv'), 'w', newline='', encoding='utf-8-sig') as f:
    f.write('cell_type,gene,n_cells,pct_expressing,mean_logCP10K\n')
    for r in rows:
        f.write('%s,%s,%d,%.4f,%.4f\n' % r)

with open(os.path.join(OUT, 'scrna_composition.csv'), 'w', newline='', encoding='utf-8-sig') as f:
    f.write('cell_type,n_cells\n')
    for r in comp_rows:
        f.write('%s,%d\n' % r)

# dotplot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

cats = [c for c in CATS if any(r[0] == c and r[2] > 0 for r in rows)]
stat = {}
for cat, g, n, pct, mean in rows:
    stat[(cat, g)] = (n, pct, mean)

vmax = max((v[2] for v in stat.values()), default=1.0)
fig, ax = plt.subplots(figsize=(10.6, 7.4))
for xi, cat in enumerate(cats):
    for yi, g in enumerate(TARGETS):
        n, pct, mean = stat[(cat, g)]
        if n == 0 or pct == 0:
            continue
        ax.scatter(xi, yi, s=20 + pct * 420, c=[mean], cmap='Reds', vmin=0, vmax=vmax,
                   edgecolors='#8a5a12', linewidths=0.4, zorder=3)
ax.set_xticks(range(len(cats)))
ax.set_xticklabels(cats, rotation=25, ha='right')
ax.set_yticks(range(len(TARGETS)))
ax.set_yticklabels(TARGETS)
ax.set_ylim(-0.7, len(TARGETS) - 0.3)
ax.set_xlim(-0.7, len(cats) - 0.3)
ax.axhline(4.5, color='#cccccc', lw=0.8, zorder=1)
ax.invert_yaxis()
ax.set_title('Signature (top) and hub (bottom) gene expression across liver cell types (GSE344087)',
             fontsize=11)
sm = plt.cm.ScalarMappable(cmap='Reds', norm=plt.Normalize(vmin=0, vmax=vmax))
cb = plt.colorbar(sm, ax=ax, shrink=0.8, pad=0.02)
cb.set_label('mean expression (log1p CP10K)', fontsize=9)
ax.text(1.02, -0.06, 'dot size = % cells expressing', transform=ax.transAxes, fontsize=8.5,
        color='#444444')
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'scrna_dotplot.png'), dpi=250)
print('saved dotplot + CSVs to', OUT, flush=True)
print('cells total:', sum(x[2].size for x in loaded), flush=True)
print('DONE', flush=True)
