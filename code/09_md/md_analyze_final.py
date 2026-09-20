# -*- coding: utf-8 -*-
"""Final manuscript-grade analysis of the 100-ns PFOA complex MD trajectories.

Reads extracted trajectories in 07_md_results, computes:
 - receptor backbone RMSD (vs crystal start & vs first saved frame)
 - ligand RMSD (vs crystal pose; + internal ligand RMSD)
 - per-CA RMSF
 - protein-ligand contacts (all <4.5 A; polar <3.5 A lig O ... N/O)
 - hydrogen bonds (donor-H...acceptor geometry, protein donors -> ligand O)
 - per-residue contact occupancy
Saves CSVs + publication figures (300 dpi) to 07_md_results/analysis_final.
"""
import json
import os
import sys

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mdtraj as md
from scipy.spatial import cKDTree

RES = r'C:\md_final_work'
RUNS = r'C:\md_final_work'
OUT = os.path.join(RES, 'analysis_final')
os.makedirs(OUT, exist_ok=True)
KEYS = ['AKT1', 'PPARA', 'PPARG']
COLORS = {'AKT1': '#145c4d', 'PPARA': '#b3541e', 'PPARG': '#5b6e8c'}
FRAME_NS = 0.05


def kabsch_rmsd(P, Q):
    """RMSD of P onto Q with optimal rotation (rows = atoms)."""
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    Pr = (R @ Pc.T).T
    return float(np.sqrt((((Pr - Qc) ** 2).sum(1)).mean()))


def analyze(key):
    dcd = os.path.join(RES, key + '_v2', 'prod', 'traj.dcd')
    top_path = os.path.join(RUNS, key + '_v2', 'complex_start.pdb')
    traj = md.load(dcd, top=top_path)
    n = traj.n_frames
    top = traj.topology
    sel = top.select
    bb = sel('name C or name CA or name N or name O')
    lig = sel('resname LIG')
    ca = sel('name CA')
    lig_o_list = sel('resname LIG and element O')
    prot_all = np.array([i for i in range(top.n_atoms) if top.atom(i).residue.name != 'LIG'])
    prot_NO = np.array([i for i in prot_all if top.atom(i).element.symbol in ('N', 'O')])
    ref = md.load(top_path)
    ref_lig = ref.xyz[0][lig]

    # ---- stage 1: align on crystal, compute crystal-referenced metrics ----
    traj.superpose(ref, atom_indices=bb)
    x = traj.xyz
    xr = ref.xyz[0]
    rec_rmsd_c = np.sqrt(((x[:, bb, :] - xr[bb]) ** 2).sum(-1).mean(-1)) * 10.0
    lig_rmsd_c = np.sqrt(((x[:, lig, :] - xr[lig]) ** 2).sum(-1).mean(-1)) * 10.0
    # ligand internal rmsd (align ligand only)
    lig_int = np.array([kabsch_rmsd(x[fi, lig, :], ref_lig) for fi in range(n)]) * 10.0
    # contacts
    nc = np.zeros(n)
    npol = np.zeros(n)
    min_d = np.zeros(n)
    occ = {}      # protein atom idx -> count of frames in contact
    occ_frames_res = {}  # (chain, resSeq, name) -> count of frames in contact
    hb_pairs = {}  # (donor, lig_o) -> count
    hb_frames = np.zeros(n)
    # donor-H map
    bh = {}
    for b in top.bonds:
        a1, a2 = b[0], b[1]
        if a1.element.symbol == 'H':
            bh.setdefault(a2.index, []).append(a1.index)
        elif a2.element.symbol == 'H':
            bh.setdefault(a1.index, []).append(a2.index)
    donors = [i for i in prot_NO if i in bh]
    donor_pos = {i: k for k, i in enumerate(donors)}
    for fi in range(n):
        P = x[fi]
        tree = cKDTree(P[prot_all])
        pairs = tree.query_ball_point(P[lig], 0.45)
        cnt = 0
        hit = set()
        hitres = set()
        for li, plist in enumerate(pairs):
            cnt += len(plist)
            for pi in plist:
                hit.add(int(pi))
        for pi in hit:
            occ[prot_all[pi]] = occ.get(prot_all[pi], 0) + 1
            r = top.atom(prot_all[pi]).residue
            hitres.add((r.chain.index, r.resSeq, r.name))
        for kk in hitres:
            occ_frames_res[kk] = occ_frames_res.get(kk, 0) + 1
        nc[fi] = cnt
        # min distance ligand-protein
        dm = tree.query(P[lig], k=1)[0]
        min_d[fi] = dm.min() * 10.0
        # polar contacts + hbonds
        t2 = cKDTree(P[prot_NO])
        pp = t2.query_ball_point(P[lig_o_list], 0.35)
        npl = 0
        seen = set()
        for oi, dlist in enumerate(pp):
            npl += len(dlist)
            o_idx = lig_o_list[oi]
            for pi in dlist:
                a_idx = prot_NO[pi]
                if a_idx in bh:
                    for h in bh[a_idx]:
                        v1 = P[a_idx] - P[h]
                        v2 = P[o_idx] - P[h]
                        cang = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9)
                        ang = np.degrees(np.arccos(np.clip(cang, -1, 1)))
                        if ang > 120.0:
                            k2 = (a_idx, o_idx)
                            hb_pairs[k2] = hb_pairs.get(k2, 0) + 1
                            seen.add(k2)
        npol[fi] = npl
        hb_frames[fi] = len(seen)

    # ---- stage 2: align on first frame, compute first-frame-referenced rmsd ----
    traj.superpose(traj, frame=0, atom_indices=bb)
    x2 = traj.xyz
    x0 = x2[0]
    rec_rmsd_f = np.sqrt(((x2[:, bb, :] - x0[bb]) ** 2).sum(-1).mean(-1)) * 10.0
    lig_rmsd_f = np.sqrt(((x2[:, lig, :] - x0[lig]) ** 2).sum(-1).mean(-1)) * 10.0
    # rmsf (from stage-1 aligned coords x)
    xc = x[:, ca, :]
    rmsf = np.sqrt(((xc - xc.mean(0)) ** 2).sum(-1).mean(0)) * 10.0
    resseq = np.array([top.atom(i).residue.resSeq for i in ca])
    resname = [top.atom(i).residue.name for i in ca]

    t = np.arange(n) * FRAME_NS
    # residue-level occupancy
    occ_res = {}
    for ai, c in occ.items():
        r = top.atom(ai).residue
        kk = (r.chain.index, r.resSeq, r.name)
        occ_res[kk] = max(occ_res.get(kk, 0), c)
    occ_rows = sorted([(v / n * 100.0, k[2], k[1]) for k, v in occ_frames_res.items()], reverse=True)

    # hbond summary
    hb_rows = sorted([(v / n * 100.0, top.atom(k[0]).residue.name, top.atom(k[0]).residue.resSeq,
                       top.atom(k[0]).name, top.atom(k[1]).name) for k, v in hb_pairs.items()], reverse=True)

    # last-20ns stats
    m = t >= (t[-1] - 20.0)
    summ = dict(
        key=key, frames=int(n), ns=round(n * FRAME_NS, 2),
        rec_rmsd_c_plateau=round(float(rec_rmsd_c[m].mean()), 2), rec_rmsd_c_sd=round(float(rec_rmsd_c[m].std()), 2),
        rec_rmsd_f_plateau=round(float(rec_rmsd_f[m].mean()), 2), rec_rmsd_f_sd=round(float(rec_rmsd_f[m].std()), 2),
        lig_rmsd_c_plateau=round(float(lig_rmsd_c[m].mean()), 2), lig_rmsd_c_max=round(float(lig_rmsd_c.max()), 2),
        lig_rmsd_f_plateau=round(float(lig_rmsd_f[m].mean()), 2),
        lig_internal_plateau=round(float(lig_int[m].mean()), 2),
        contacts_plateau=round(float(nc[m].mean()), 0), contacts_min_after2ns=int(nc[t > 2].min()),
        polar_plateau=round(float(npol[m].mean()), 1),
        hb_frames_pct=round(float((hb_frames > 0).mean() * 100.0), 1),
        min_dist_plateau=round(float(min_d[m].mean()), 2),
    )
    print(json.dumps(summ, ensure_ascii=False), flush=True)
    print('   top contact residues: %s' % ', '.join('%s%s(%.0f%%)' % (r[1], r[2], r[0]) for r in occ_rows[:12]), flush=True)
    print('   persistent H-bonds: %s' % (', '.join('%s%s:%s-%s(%.0f%%)' % (h[1], h[2], h[3], h[4], h[0]) for h in hb_rows[:6]) or 'none >10%'), flush=True)

    np.savetxt(os.path.join(OUT, key + '_series.csv'),
               np.column_stack([t, rec_rmsd_c, rec_rmsd_f, lig_rmsd_c, lig_rmsd_f, lig_int, nc, npol, hb_frames, min_d]),
               delimiter=',', comments='',
               header='time_ns,rec_rmsd_vs_crystal_A,rec_rmsd_vs_f0_A,lig_rmsd_vs_crystal_A,lig_rmsd_vs_f0_A,lig_internal_rmsd_A,n_contacts,n_polar,hbond_count,min_lig_prot_dist_A')
    with open(os.path.join(OUT, key + '_contacts.csv'), 'w', encoding='utf-8-sig') as f:
        f.write('pct_frames,resname,resseq\n')
        for pct, rn, ri in occ_rows:
            f.write('%.1f,%s,%d\n' % (pct, rn, ri))
    with open(os.path.join(OUT, key + '_hbonds.csv'), 'w', encoding='utf-8-sig') as f:
        f.write('pct_frames,donor_res,donor_resseq,donor_atom,ligand_atom\n')
        for pct, rn, ri, da, la in hb_rows:
            f.write('%.1f,%s,%d,%s,%s\n' % (pct, rn, ri, da, la))
    np.savetxt(os.path.join(OUT, key + '_rmsf.csv'), np.column_stack([resseq, rmsf]), delimiter=',',
               header='resSeq,rmsf_Ca_A', comments='')
    return dict(summ=summ, t=t, rec_c=rec_rmsd_c, rec_f=rec_rmsd_f, lig_c=lig_rmsd_c,
                lig_f=lig_rmsd_f, nc=nc, npol=npol, rmsf=rmsf, resseq=resseq, occ=occ_rows,
                hb=hb_rows, hbf=hb_frames)


res = {}
for k in KEYS:
    print('== %s ==' % k, flush=True)
    res[k] = analyze(k)

# ---- summary csv ----
with open(os.path.join(OUT, 'md_final_summary.csv'), 'w', encoding='utf-8-sig') as f:
    ks = list(res['AKT1']['summ'].keys())
    f.write(','.join(ks) + '\n')
    for k in KEYS:
        f.write(','.join(str(res[k]['summ'][x]) for x in ks) + '\n')

# ---- figures ----
plt.rcParams.update({'font.size': 9})
fig, axes = plt.subplots(2, 2, figsize=(10.2, 6.8), dpi=300)
ax = axes[0][0]
for k in KEYS:
    ax.plot(res[k]['t'], res[k]['rec_c'], lw=1.1, color=COLORS[k], label=k)
ax.set_title('A  Receptor backbone RMSD (vs initial structure)')
ax.set_xlabel('time (ns)'); ax.set_ylabel('RMSD (A)'); ax.legend(frameon=False)
ax = axes[0][1]
for k in KEYS:
    ax.plot(res[k]['t'], res[k]['lig_c'], lw=1.1, color=COLORS[k], label=k)
ax.set_title('B  Ligand RMSD (vs docked pose)')
ax.set_xlabel('time (ns)'); ax.set_ylabel('RMSD (A)'); ax.legend(frameon=False)
ax = axes[1][0]
for k in KEYS:
    ax.plot(res[k]['t'], res[k]['nc'], lw=0.9, color=COLORS[k], label=k)
ax.set_title('C  Protein-ligand contacts (<4.5 A)')
ax.set_xlabel('time (ns)'); ax.set_ylabel('contact count'); ax.legend(frameon=False)
ax = axes[1][1]
for k in KEYS:
    ax.plot(res[k]['resseq'], res[k]['rmsf'], lw=0.9, color=COLORS[k], label=k)
ax.set_title('D  Per-residue RMSF (CA)')
ax.set_xlabel('residue'); ax.set_ylabel('RMSF (A)'); ax.legend(frameon=False)
for a in axes.flat:
    a.spines['top'].set_visible(False); a.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'fig_md_stability.png'))
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.9), dpi=300)
for ax, k in zip(axes, KEYS):
    rows = [r for r in res[k]['occ'][:12]][::-1]
    ys = np.arange(len(rows))
    ax.barh(ys, [r[0] for r in rows], color=COLORS[k], alpha=0.9)
    ax.set_yticks(ys)
    ax.set_yticklabels(['%s%s' % (r[1], r[2]) for r in rows], fontsize=8)
    ax.set_xlabel('% frames in contact')
    ax.set_title('%s - top contact residues' % k)
    ax.set_xlim(0, 100)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUT, 'fig_md_contacts.png'))
plt.close(fig)

# per-system 4-panel
for k in KEYS:
    r = res[k]
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.2), dpi=300)
    ax = axes[0][0]
    ax.plot(r['t'], r['rec_c'], lw=1.1, color=COLORS[k], label='vs initial structure')
    ax.plot(r['t'], r['rec_f'], lw=0.9, color='#888888', label='vs first frame')
    ax.set_title('%s - receptor backbone RMSD' % k); ax.set_xlabel('time (ns)'); ax.set_ylabel('RMSD (A)'); ax.legend(frameon=False)
    ax = axes[0][1]
    ax.plot(r['t'], r['lig_c'], lw=1.1, color='#b3541e', label='vs docked pose')
    ax.plot(r['t'], r['lig_f'], lw=0.9, color='#888888', label='vs first frame')
    ax.set_title('%s - ligand RMSD' % k); ax.set_xlabel('time (ns)'); ax.set_ylabel('RMSD (A)'); ax.legend(frameon=False)
    ax = axes[1][0]
    ax.plot(r['t'], r['nc'], lw=1.0, color='#444444', label='all contacts <4.5 A')
    ax.plot(r['t'], r['npol'], lw=1.1, color=COLORS[k], label='polar contacts <3.5 A')
    ax.plot(r['t'], r['hbf'], lw=1.1, color='#8a5a12', label='H-bonds (geometric)')
    ax.set_title('%s - protein-ligand interactions' % k); ax.set_xlabel('time (ns)'); ax.set_ylabel('count'); ax.legend(frameon=False, fontsize=7)
    ax = axes[1][1]
    ax.plot(r['resseq'], r['rmsf'], lw=0.9, color=COLORS[k])
    ax.set_title('%s - per-residue RMSF (CA)' % k); ax.set_xlabel('residue'); ax.set_ylabel('RMSF (A)')
    for a in axes.flat:
        a.spines['top'].set_visible(False); a.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, 'figS_%s_panels.png' % k))
    plt.close(fig)

print('FIGURES DONE')
print('ALL DONE')
