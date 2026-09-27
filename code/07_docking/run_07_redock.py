# -*- coding: utf-8 -*-
"""run_07_redock v2：8 靶点 native-ligand redocking 验证（AutoDock Vina 1.2.5）
与主对接（run_07_docking.py）完全一致的受体预处理/盒子/搜索参数；
配体 = 晶体自带 native ligand，多拷贝时取质心离对接盒最近的一组；
RMSD = 重原子 Kabsch + 元素分组内匈牙利最优指派（免 scipy）。
产出：07_docking/redock/redock_summary.csv + per-target logs
"""
import os, subprocess, re, csv, time, math
import numpy as np

BASE = r'C:\Users\Administrator\.openclaw-autoclaw\workspace\PFOA_NAFLD_研究方案'
DK = os.path.join(BASE, '07_docking')
REC = os.path.join(DK, 'receptors')
OUT = os.path.join(DK, 'redock')
os.makedirs(OUT, exist_ok=True)
VINA = os.path.join(DK, 'vina.exe')
OBEX = r'D:\Program Files\Autoclaw\resources\python\Lib\site-packages\openbabel\bin\obabel.exe'

TARGETS = [
    ('PPARA', '1K7L', (-17.6, -14.4, -5.4), '544'),
    ('PPARG', '2PRG', (49.7, -37.0, 19.3), 'BRL'),
    ('FXR_NR1H4', '4OIV', (12.9, 14.1, -20.8), 'XX9'),
    ('EGFR', '1M17', (22.0, 0.3, 52.8), 'AQ4'),
    ('CASP3', '1GFW', (37.6, 33.6, 28.0), 'MSI'),
    ('MMP9', '1GKC', (53.3, 22.5, 129.7), 'NFH'),
    ('TNF', '2AZ5', (-19.4, 74.7, 33.8), '307'),
    ('AKT1', '3O96', (8.4, -6.8, 12.6), 'IQO'),
]

AD2ELEM = {'A': 'C', 'C': 'C', 'N': 'N', 'NA': 'N', 'NS': 'N', 'OA': 'O', 'OS': 'O',
           'SA': 'S', 'S': 'S', 'P': 'P', 'F': 'F', 'Cl': 'Cl', 'CL': 'Cl',
           'Br': 'Br', 'BR': 'Br', 'I': 'I'}

def hungarian(cost):
    """e-maxx O(n^3) 匈牙利；返回 ans[i] = cost 矩阵第 i 行分配的列。"""
    n = cost.shape[0]
    INF = 1e18
    u = np.zeros(n + 1); v = np.zeros(n + 1)
    p = np.zeros(n + 1, dtype=int); way = np.zeros(n + 1, dtype=int)
    for i in range(1, n + 1):
        p[0] = i; j0 = 0
        minv = np.full(n + 1, INF); used = np.zeros(n + 1, dtype=bool)
        while True:
            used[j0] = True
            i0 = p[j0]; delta = INF; j1 = -1
            for j in range(1, n + 1):
                if not used[j]:
                    cur = cost[i0 - 1, j - 1] - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j] = cur; way[j] = j0
                    if minv[j] < delta:
                        delta = minv[j]; j1 = j
            for j in range(n + 1):
                if used[j]:
                    u[p[j]] += delta; v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]; p[j0] = p[j1]; j0 = j1
            if j0 == 0:
                break
    ans = np.zeros(n, dtype=int)
    for j in range(1, n + 1):
        if p[j] > 0:
            ans[p[j] - 1] = j - 1
    return ans

def clean_receptor(pid, keep_zn=False):
    src = os.path.join(REC, f'{pid}.pdb')
    dst = os.path.join(REC, f'{pid}_clean.pdb')
    out = []
    for line in open(src, encoding='utf-8', errors='replace'):
        if line.startswith('ATOM'):
            if line[16] not in (' ', 'A'):
                continue
            out.append(line)
        elif keep_zn and line.startswith('HETATM') and line[17:20].strip() == 'ZN':
            out.append(line)
        elif line.startswith('TER'):
            out.append(line)
    open(dst, 'w', newline='\n').write(''.join(out) + 'END\n')
    return dst

def pdbqt_receptor(clean_pdb, pid):
    dst = os.path.join(REC, f'{pid}.pdbqt')
    if not (os.path.exists(dst) and os.path.getsize(dst) > 1000):
        subprocess.run([OBEX, clean_pdb, '-O', dst, '-xr', '--partialcharge', 'gasteiger'],
                       capture_output=True, text=True, timeout=600)
    return dst

def extract_native_copies(pid, resname):
    """全部拷贝行组：[(chain, resseq), [lines]]；altloc 仅 ' '/'A'"""
    groups = {}
    for line in open(os.path.join(REC, f'{pid}.pdb'), encoding='utf-8', errors='replace'):
        if line.startswith('HETATM') and line[17:20].strip() == resname:
            if line[16] not in (' ', 'A'):
                continue
            groups.setdefault((line[21], line[22:27]), []).append(line)
    return groups

def pick_copy_nearest(groups, center):
    best_k, best_d = None, None
    for k, lines in groups.items():
        xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines])
        d = float(np.linalg.norm(xyz.mean(0) - np.array(center)))
        if best_d is None or d < best_d:
            best_k, best_d = k, d
    return best_k, groups[best_k], best_d

def pdbqt_ligand(lig_lines, pid, resname):
    lig_pdb = os.path.join(OUT, f'{pid}_{resname}_lig.pdb')
    open(lig_pdb, 'w', newline='\n').write(''.join(lig_lines) + 'END\n')
    lig_pdbqt = os.path.join(OUT, f'{pid}_{resname}_lig.pdbqt')
    r = subprocess.run([OBEX, lig_pdb, '-O', lig_pdbqt, '--partialcharge', 'gasteiger'],
                       capture_output=True, text=True, timeout=300)
    ok = os.path.exists(lig_pdbqt) and os.path.getsize(lig_pdbqt) > 100
    return lig_pdb, lig_pdbqt, ok, ((r.stderr or r.stdout).strip()[-100:])

def crystal_heavy(lig_lines):
    els, xyz = [], []
    for line in lig_lines:
        el = line[76:78].strip().upper()
        if not el:
            el = line[12:16].strip()[0].upper()
        if el in ('H', 'D'):
            continue
        els.append(el)
        xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return els, xyz

def docked_heavy(pdbqt_path):
    txt = open(pdbqt_path, encoding='utf-8', errors='replace').read()
    m = re.search(r'MODEL\s+(\d+)\s*\n(.*?)ENDMDL', txt, re.S)
    block = m.group(2) if m else txt
    els, xyz = [], []
    for line in block.splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            tok = line.split()
            adt = tok[-1] if tok else ''
            if not adt or adt.upper().startswith('H'):
                continue
            el = AD2ELEM.get(adt)
            if el is None:
                el = adt[0].upper()
                if el == 'H':
                    continue
            els.append(el)
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return els, xyz

def rmsd_assign(cels, cxyz, dels, dxyz):
    if sorted(cels) != sorted(dels):
        return None, f'elem-multiset mismatch cry={len(cels)} dock={len(dels)}'
    P = np.zeros((len(cels), 3)); Q = np.zeros((len(dels), 3))
    # 元素分组：晶体侧按出现序编号，docked 侧对每个晶体原子找最优 docked 原子
    idx_by_el_c = {}
    for i, e in enumerate(cels):
        idx_by_el_c.setdefault(e, []).append(i)
    idx_by_el_d = {}
    for i, e in enumerate(dels):
        idx_by_el_d.setdefault(e, []).append(i)
    pairs = []
    for e in idx_by_el_c:
        ci = idx_by_el_c[e]; di = idx_by_el_d[e]
        sub = np.array([[math.dist(cxyz[i], dxyz[j]) for j in di] for i in ci])
        assign = hungarian(sub)
        for a, b in zip(range(len(ci)), assign):
            pairs.append((ci[a], di[b]))
    for a, b in pairs:
        P[a] = cxyz[a]; Q[a] = dxyz[b]
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    Pr = (R @ Pc.T).T
    return float(math.sqrt(((Pr - Qc) ** 2).sum() / len(P))), ''

summary = []
for name, pid, (cx, cy, cz), lig in TARGETS:
    print('=' * 70)
    print(f'--- {name} ({pid}) native ligand {lig} ---')
    t_start = time.time()
    keep_zn = (pid == '1GKC')
    clean = clean_receptor(pid, keep_zn)
    rec_pdbqt = pdbqt_receptor(clean, pid)
    if os.path.getsize(rec_pdbqt) < 1000:
        alt = os.path.join(REC, f'{pid}_nozn.pdbqt')
        if os.path.exists(alt) and os.path.getsize(alt) > 1000:
            rec_pdbqt = alt
            print('!! ZN-conversion empty, fallback to', os.path.basename(alt))
    print('receptor pdbqt:', os.path.getsize(rec_pdbqt), 'bytes')

    groups = extract_native_copies(pid, lig)
    if not groups:
        print('!! no native ligand atoms')
        summary.append({'Target': name, 'PDB': pid, 'Native_Lig': lig, 'Lig_Atoms': 0,
                        'RMSD_A': 'FAIL', 'Redock_kcal': 'FAIL', 'Note': 'no ligand atoms'})
        continue
    copy_key, lig_lines, dist = pick_copy_nearest(groups, (cx, cy, cz))
    lig_pdb, lig_pdbqt, ok, msg = pdbqt_ligand(lig_lines, pid, lig)
    if not ok:
        print('!! obabel ligand fail:', msg)
        summary.append({'Target': name, 'PDB': pid, 'Native_Lig': lig, 'Lig_Atoms': len(lig_lines),
                        'RMSD_A': 'FAIL', 'Redock_kcal': 'FAIL', 'Note': 'obabel fail'})
        continue
    print(f'ligand copy {copy_key} ({len(lig_lines)} atoms, centroid {dist:.1f} A from box center)')

    out_pdbqt = os.path.join(OUT, f'{name}_{pid}_redock_out.pdbqt')
    log = os.path.join(OUT, f'{name}_{pid}_redock.log')
    cmd = [VINA, '--receptor', rec_pdbqt, '--ligand', lig_pdbqt,
           '--center_x', str(cx), '--center_y', str(cy), '--center_z', str(cz),
           '--size_x', '25', '--size_y', '25', '--size_z', '25',
           '--exhaustiveness', '32', '--num_modes', '10', '--seed', '42', '--cpu', '8',
           '--out', out_pdbqt]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=2400)
    dt = time.time() - t0
    open(log, 'w', encoding='utf-8').write(r.stdout + '\n' + (r.stderr or ''))
    modes = re.findall(r'^\s+(\d+)\s+(-?\d+\.\d+)\s+', r.stdout, re.M)
    if not modes:
        print('!! vina no modes:', (r.stderr or '')[:150])
        summary.append({'Target': name, 'PDB': pid, 'Native_Lig': lig, 'Lig_Atoms': len(lig_lines),
                        'RMSD_A': 'FAIL', 'Redock_kcal': 'FAIL', 'Note': 'vina fail'})
        continue
    best = float(modes[0][1])
    cels, cxyz = crystal_heavy(lig_lines)
    dels, dxyz = docked_heavy(out_pdbqt)
    rmsd_val, note = rmsd_assign(cels, cxyz, dels, dxyz)
    if rmsd_val is None:
        print('!!', note)
        summary.append({'Target': name, 'PDB': pid, 'Native_Lig': lig, 'Lig_Atoms': len(cels),
                        'RMSD_A': 'SEQ_MISMATCH', 'Redock_kcal': f'{best:.1f}', 'Note': note})
    else:
        flag = 'PASS' if rmsd_val <= 2.0 else 'over'
        print(f'redocked {dt:.0f}s | mode1 {best:.1f} kcal/mol | heavy-atom RMSD {rmsd_val:.2f} A [{flag}]')
        summary.append({'Target': name, 'PDB': pid, 'Native_Lig': lig, 'Lig_Atoms': len(cels),
                        'RMSD_A': f'{rmsd_val:.2f}', 'Redock_kcal': f'{best:.1f}',
                        'Note': 'within 2.0 A' if rmsd_val <= 2.0 else 'above 2.0 A'})
    print(f'(target wall time {time.time() - t_start:.0f}s)')

csv_path = os.path.join(OUT, 'redock_summary.csv')
with open(csv_path, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.DictWriter(f, fieldnames=['Target', 'PDB', 'Native_Lig', 'Lig_Atoms', 'RMSD_A', 'Redock_kcal', 'Note'])
    w.writeheader(); w.writerows(summary)
print()
print('=' * 70)
print('REDOCK SUMMARY:')
npass = sum(1 for s in summary if isinstance(s['RMSD_A'], str) and s['RMSD_A'].replace('.', '').isdigit() and float(s['RMSD_A']) <= 2.0)
for row in summary:
    print(f"  {row['Target']:12s} {row['PDB']} {row['Native_Lig']:4s} RMSD={row['RMSD_A']} kcal={row['Redock_kcal']} {row['Note']}")
print(f'PASS (<=2.0 A): {npass}/8')
print('saved:', csv_path)
print('REDOCK-DONE')
