# -*- coding: utf-8 -*-
"""run_04c: Hallmark enrichment (Enrichr route) + plot"""
import gseapy as gp
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

BASE = r'C:\Users\user\workspace\PFOA_NAFLD_project'
ANA = os.path.join(BASE, '04_analysis'); FIG = os.path.join(BASE, '05_figures')
genes = [g.strip() for g in open(os.path.join(ANA, 'intersection_genes_479.txt'), encoding='utf-8') if g.strip()]
print('genes:', len(genes))

res = gp.enrichr(gene_list=genes, gene_sets=['MSigDB_Hallmark_2020'],
                 organism='human', outdir=None, cutoff=1.0)
df = res.results
df_sig = df[df['Adjusted P-value'] < 0.05].copy()
df_sig = df_sig.sort_values('Adjusted P-value')
out = os.path.join(ANA, 'run04_Hallmark.csv')
df_sig.to_csv(out, index=False, encoding='utf-8-sig')
print('Hallmark significant sets:', len(df_sig))
for _, r in df_sig.head(10).iterrows():
    print(' ', r['Term'], '| adjP = %.2g' % r['Adjusted P-value'], '| overlap', r['Overlap'])

# plot (top 12, horizontal bars, -log10 adjP)
top = df_sig.head(12).iloc[::-1]
fig, ax = plt.subplots(figsize=(9.5, 6), dpi=220)
vals = -__import__('numpy').log10(top['Adjusted P-value'].values)
labels = [t.replace('MSigDB_Hallmark_2020__', '').replace('_', ' ') for t in top['Term']]
ax.barh(labels, vals, color='#145c4d', height=0.62)
for i, (v, ov) in enumerate(zip(vals, top['Overlap'].values)):
    ax.text(v + 0.05, i, str(ov), va='center', fontsize=7.5, color='#4a545c')
ax.set_xlabel('-log10(adjusted P)')
ax.set_title('Hallmark enrichment (Enrichr, top 12) - PFOA x NAFLD intersection', fontsize=11.5)
ax.tick_params(labelsize=8.5)
for s in ['top', 'right']:
    ax.spines[s].set_visible(False)
fig.tight_layout()
png = os.path.join(FIG, 'run04_Hallmark_barplot.png')
fig.savefig(png, facecolor='white')
print('saved:', png, os.path.getsize(png), 'bytes')
print('RUN04C-DONE')
