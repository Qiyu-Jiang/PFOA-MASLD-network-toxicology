# -*- coding: utf-8 -*-
"""run_05c: ML diagnostic signature v2 (tightened feature selection + per-cohort standardization for external validation)"""
import gzip, csv, os, re
import numpy as np
from collections import defaultdict, Counter

BASE = r'C:\Users\user\workspace\PFOA_NAFLD_project'
RAW = os.path.join(BASE, '02_data_raw'); ANA = os.path.join(BASE, '04_analysis'); FIG = os.path.join(BASE, '05_figures')

def load_series(path):
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
        lines = f.read().splitlines()
    titles = None; ser_s = ser_e = None
    for i, ln in enumerate(lines):
        if ln.startswith('!Sample_title'): titles = [x.strip('"') for x in ln.split('\t')[1:]]
        if ln.startswith('!series_matrix_table_begin'): ser_s = i
        if ln.startswith('!series_matrix_table_end'): ser_e = i
    expr = {}
    for ln in lines[ser_s + 1:ser_e]:
        p = ln.split('\t'); pid = p[0].strip('"')
        try: expr[pid] = [float(x.strip('"')) for x in p[1:]]
        except ValueError: continue
    return titles, expr

def load_probe_map(top_table_csv):
    m = {}
    with open(top_table_csv, encoding='utf-8-sig') as f:
        r = csv.DictReader(f); pid_col = r.fieldnames[0]
        for row in r:
            sym = (row.get('Symbol') or '').strip().upper()
            if sym: m[row[pid_col]] = sym
    return m

t1, e1 = load_series(os.path.join(RAW, 'GSE89632_series_matrix.txt.gz'))
t2, e2 = load_series(os.path.join(RAW, 'GSE48452_series_matrix.txt.gz'))
m1 = load_probe_map(os.path.join(ANA, 'GSE89632_NASH_vs_HC_topTable.csv'))
m2 = {}
for row in csv.DictReader(open(os.path.join(ANA, 'GSE48452_probe2symbol.csv'), encoding='utf-8-sig')):
    pid = row['PROBEID']; sym = (row['SYMBOL'] or '').strip().upper()
    if pid not in m2 or (not m2[pid] and sym): m2[pid] = sym
m2 = {k: v for k, v in m2.items() if v}

def build(expr, pmap):
    by_sym = defaultdict(list)
    for pid, vals in expr.items():
        s = pmap.get(pid)
        if s: by_sym[s].append(vals)
    out = {}
    for s, rows in by_sym.items():
        arr = np.array(rows)
        out[s] = arr[np.argmax(arr.mean(axis=1))]
    return out

X1 = build(e1, m1); X2 = build(e2, m2)
grp1 = [str(t).split('_')[1] if '_' in str(t) else 'UNK' for t in t1]
grp2 = [str(t).split(',')[1].strip() if ',' in str(t) else 'UNK' for t in t2]
inter = [r['Symbol'].upper() for r in csv.DictReader(open(os.path.join(ANA, 'PFOAxNAFLD_intersection_v2.csv'), encoding='utf-8-sig')) if r['Symbol']]
feats = [s for s in inter if s in X1 and s in X2]
M1 = np.array([X1[s] for s in feats], dtype=float).T
M2 = np.array([X2[s] for s in feats], dtype=float).T
y1 = np.array([1 if g in ('SS', 'NASH') else 0 for g in grp1])
y2 = np.array([1 if g in ('S', 'N') else 0 for g in grp2])
print('feats:', len(feats), '| train', M1.shape, '| external', M2.shape)
np.savez(os.path.join(ANA, 'run05_matrices.npz'), M1=M1, M2=M2, y1=y1, y2=y2, feats=np.array(feats))

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve

sc1 = StandardScaler().fit(M1)
Z1 = sc1.transform(M1)
sc2 = StandardScaler().fit(M2)          # per-cohort standardization (standard practice for cross-platform validation)
Z2w = sc2.transform(M2)

import warnings
warnings.filterwarnings('ignore')
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

lasso = LogisticRegressionCV(Cs=np.logspace(-3, 1, 25), cv=cv, penalty='l1', solver='liblinear',
                             scoring='roc_auc', max_iter=8000, random_state=42)
lasso.fit(Z1, y1)
coefs = lasso.coef_[0]
lasso_set = [feats[i] for i in range(len(feats)) if abs(coefs[i]) > 1e-8]  # corrected 2026-09: non-zero coefficients only
print('LASSO nonzero:', int((np.abs(coefs) > 1e-8).sum()), '| top30 taken')

rf = RandomForestClassifier(n_estimators=1000, random_state=42, n_jobs=-1)
rf.fit(Z1, y1)
imp = rf.feature_importances_
rf_top20 = [feats[i] for i in np.argsort(-imp)[:20]]

svm = SVC(kernel='linear', C=1.0)
rfe = RFE(estimator=svm, n_features_to_select=20, step=0.1)
rfe.fit(Z1, y1)
svm_top20 = [feats[i] for i in np.where(rfe.support_)[0]]

votes = Counter(lasso_set + rf_top20 + svm_top20)
sig = [g for g, c in votes.items() if c >= 2]
sig = sorted(sig, key=lambda g: (-votes[g], g))
print('consensus sign (>=2):', len(sig))
print('SIGNATURE:', sig)

idx = [feats.index(g) for g in sig]
clf = LogisticRegression(max_iter=3000, random_state=42)
proba_cv = cross_val_predict(clf, Z1[:, idx], y1, cv=cv, method='predict_proba')[:, 1]
auc_cv = roc_auc_score(y1, proba_cv)
clf.fit(Z1[:, idx], y1)
p_A = clf.predict_proba(sc1.transform(M2)[:, idx])[:, 1]   # training parameters transferred directly (previous convention)
p_B = clf.predict_proba(Z2w[:, idx])[:, 1]                 # per-cohort standardization (current convention)
auc_A = roc_auc_score(y2, p_A)
auc_B = roc_auc_score(y2, p_B)
print('train CV AUC: %.3f' % auc_cv)
print('external AUC (train-scaler): %.3f' % auc_A)
print('external AUC (per-cohort scaler): %.3f' % auc_B)
print('p_B stats: min=%.3f max=%.3f mean=%.3f std=%.3f' % (p_B.min(), p_B.max(), p_B.mean(), p_B.std()))

# save
with open(os.path.join(ANA, 'run05_signature_v2.csv'), 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['Gene', 'Votes', 'LASSO_coef', 'RF_importance', 'SVM_RFE'])
    for g in sig:
        i = feats.index(g)
        w.writerow([g, votes[g], round(coefs[i], 4), round(imp[i], 5), 'Y' if g in svm_top20 else ''])
with open(os.path.join(ANA, 'run05_method_selections_v2.csv'), 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['Method', 'Genes'])
    w.writerow(['LASSO_nonzero', ';'.join(lasso_set)])
    w.writerow(['RF_top20', ';'.join(rf_top20)])
    w.writerow(['SVM_RFE_top20', ';'.join(svm_top20)])

# plot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'sans-serif']
fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=200)
for ax, (yv, pv, ttl) in zip(axes, [
        (y1, proba_cv, 'Training (GSE89632, 5-fold CV) AUC = %.3f' % auc_cv),
        (y2, p_B, 'External (GSE48452, per-cohort scaling) AUC = %.3f' % auc_B)]):
    fpr, tpr, _ = roc_curve(yv, pv)
    ax.plot(fpr, tpr, color='#145c4d', lw=2)
    ax.plot([0, 1], [0, 1], ls='--', color='#9aa5ad', lw=1)
    ax.set_xlabel('False positive rate'); ax.set_ylabel('True positive rate')
    ax.set_title(ttl, fontsize=11); ax.set_aspect('equal')
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(FIG, 'run05_ROC_v2.png'), facecolor='white')

nz = [feats[i] for i in np.argsort(-np.abs(coefs)) if abs(coefs[i]) > 1e-8]
top15 = np.argsort(-imp)[:15]
fig2, axes2 = plt.subplots(1, 2, figsize=(13, 5.5), dpi=200)
axes2[0].barh([feats[i] for i in top15][::-1], [imp[i] for i in top15][::-1], color='#145c4d', height=0.62)
axes2[0].set_title('Random Forest importance (top 15)', fontsize=11); axes2[0].tick_params(labelsize=8)
axes2[1].barh(nz[::-1], [coefs[feats.index(g)] for g in nz][::-1], color='#0f3f6e', height=0.62)
axes2[1].set_title('LASSO coefficients (non-zero features)', fontsize=11); axes2[1].tick_params(labelsize=8)
for ax in axes2:
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
fig2.tight_layout()
fig2.savefig(os.path.join(FIG, 'run05_feature_selection_v2.png'), facecolor='white')
print('RUN05C-DONE')
