# -*- coding: utf-8 -*-
"""run_07_redock_pass2: compact-box pose-recovery validation (standard scheme).

Box center = centroid of the largest native-ligand copy; box edge lengths = per-axis
ligand extent + 8 A (minimum 10 A). All other parameters match the primary docking
run (exhaustiveness 32, seed 42, num_modes 10). Element normalization fix: both
sides are upper-cased. Outputs redock_summary_compact.csv.
"""
import os, subprocess, re, csv, time, math
import numpy as np

BASE = os.environ.get('PFOA_PROJECT_ROOT', r'C:\Users\user\workspace\PFOA_NAFLD_project')
DK = os.path.join(BASE, '07_docking')
REC = os.path.join(DK, 'receptors')
OUT = os.path.join(DK, 'redock')
os.makedirs(OUT, exist_ok=True)
VINA = os.path.join(DK, 'vina.exe')
OBEX = os.environ.get('OBABEL_EXE', 'obabel')

TARGETS = [
    ('PPARA', '1K7L', '544'),
    ('PPARG', '2PRG', 'BRL'),
    ('FXR_NR1H4', '4OIV', 'XX9'),
    ('EGFR', '1M17', 'AQ4'),
    ('CASP3', '1GFW', 'MSI'),
    ('MMP9', '1GKC', 'NFH'),
    ('TNF', '2AZ5', '307'),
    ('AKT1', '3O96', 'IQO'),
]

AD2ELEM = {'A': 'C', 'C': 'C', 'N': 'N', 'NA': 'N', 'NS': 'N', 'OA': 'O', 'OS': 'O',
           'SA': 'S', 'S': 'S', 'P': 'P', 'F': 'F', 'CL': 'Cl', 'BR': 'Br', 'I': 'I'}

def hungarian(cost):
    """e-maxx O(n^3) Hungarian algorithm; returns ans[i] = column assigned to row i of cost."""
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

def crystal_heavy(lines):
    els, xyz = [], []
    for line in lines:
        el = line[76:78].strip().upper()
        if not el:
            el = line[12:16].strip()[0].upper()
        if el in ('H', 'D'):
            continue
        els.append(el.upper())
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
            el = AD2ELEM.get(adt, adt[0].upper())
            els.append(el.upper())
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return els, xyz

def rmsd_assign(cels, cxyz, dels, dxyz):
    if sorted(cels) != sorted(dels):
        return None, f'elem-multiset mismatch cry={len(cels)} dock={len(dels)} cryset={sorted(set(cels))} dockset={sorted(set(dels))}'
    P = np.zeros((len(cels), 3)); Q = np.zeros((len(dels), 3))
    idx_c = {}; idx_d = {}
    for i, e in enumerate(cels): idx_c.setdefault(e, []).append(i)
    for i, e in enumerate(dels): idx_d.setdefault(e, []).append(i)
    for e in idx_c:
        ci = idx_c[e]; di = idx_d[e]
        sub = np.array([[math.dist(cxyz[i], dxyz[j]) for j in di] for i in ci])
        assign = hungarian(sub)
        for a, b in zip(range(len(ci)), assign):
            P[ci[a]] = cxyz[ci[a]]; Q[ci[a]] = dxyz[di[b]]
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    Pr = (R @ Pc.T).T
    return float(math.sqrt(((Pr - Qc) ** 2).sum() / len(P))), ''

summary = []
for name, pid, lig in TARGETS:
    print('=' * 70)
    print(f'--- {name} ({pid}) {lig} compact-box ---')
    t0 = time.time()
    # Ligand copy: largest HETATM copy (same copy pass1 selected as centroid-nearest
    # to the 25 A screening box); reuse that copy directly here.
    groups = {}
    for line in open(os.path.join(REC, f'{pid}.pdb'), encoding='utf-8', errors='replace'):
        if line.startswith('HETATM') and line[17:20].strip() == lig and line[16] in (' ', 'A'):
            groups.setdefault((line[21], line[22:27]), []).append(line)
    copy_key, lines = max(groups.items(), key=lambda kv: len(kv[1]))
    xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines])
    centroid = xyz.mean(0)
    extent = xyz.max(0) - xyz.min(0)
    size = np.maximum(extent + 8.0, 10.0)
    size = np.ceil(size).astype(int)
    print(f'copy {copy_key} ({len(lines)} atoms) | centroid {centroid.round(2)} | compact box {size.tolist()} A')

    lig_pdb = os.path.join(OUT, f'{pid}_{lig}_lig.pdb')
    open(lig_pdb, 'w', newline='\n').write(''.join(lines) + 'END\n')
    lig_pdbqt = os.path.join(OUT, f'{pid}_{lig}_lig.pdbqt')
    subprocess.run([OBEX, lig_pdb, '-O', lig_pdbqt, '--partialcharge', 'gasteiger'],
                   capture_output=True, text=True, timeout=300)

    # Receptor: pdbqt already prepared by pass1.
    rec_pdbqt = os.path.join(REC, f'{pid}.pdbqt')
    if os.path.getsize(rec_pdbqt) < 1000:
        rec_pdbqt = os.path.join(REC, f'{pid}_nozn.pdbqt')

    out_pdbqt = os.path.join(OUT, f'{name}_{pid}_compact_out.pdbqt')
    log = os.path.join(OUT, f'{name}_{pid}_compact.log')
    cmd = [VINA, '--receptor', rec_pdbqt, '--ligand', lig_pdbqt,
           '--center_x', f'{centroid[0]:.3f}', '--center_y', f'{centroid[1]:.3f}', '--center_z', f'{centroid[2]:.3f}',
           '--size_x', str(size[0]), '--size_y', str(size[1]), '--size_z', str(size[2]),
           '--exhaustiveness', '32', '--num_modes', '10', '--seed', '42', '--cpu', '8',
           '--out', out_pdbqt]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=2400)
    open(log, 'w', encoding='utf-8').write(r.stdout + '\n' + (r.stderr or ''))
    modes = re.findall(r'^\s+(\d+)\s+(-?\d+\.\d+)\s+', r.stdout, re.M)
    if not modes:
        print('!! vina no modes')
        summary.append({'Target': name, 'PDB': pid, 'Native_Lig': lig,
                        'RMSD_compact_A': 'FAIL', 'kcal_compact': 'FAIL', 'Box': 'compact'})
        continue
    best = float(modes[0][1])
    cels, cxyz = crystal_heavy(lines)
    dels, dxyz = docked_heavy(out_pdbqt)
    rv, note = rmsd_assign(cels, cxyz, dels, dxyz)
    if rv is None:
        print('!!', note)
        summary.append({'Target': name, 'PDB': pid, 'Native_Lig': lig,
                        'RMSD_compact_A': 'MISMATCH', 'kcal_compact': f'{best:.1f}', 'Box': f'{size[0]}x{size[1]}x{size[2]}'})
    else:
        flag = 'PASS' if rv <= 2.0 else 'over'
        print(f'redocked {time.time()-t0:.0f}s | mode1 {best:.1f} kcal/mol | RMSD {rv:.2f} A [{flag}]')
        summary.append({'Target': name, 'PDB': pid, 'Native_Lig': lig,
                        'RMSD_compact_A': f'{rv:.2f}', 'kcal_compact': f'{best:.1f}', 'Box': f'{size[0]}x{size[1]}x{size[2]}'})

csv_path = os.path.join(OUT, 'redock_summary_compact.csv')
with open(csv_path, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.DictWriter(f, fieldnames=['Target', 'PDB', 'Native_Lig', 'RMSD_compact_A', 'kcal_compact', 'Box'])
    w.writeheader(); w.writerows(summary)
print()
print('COMPACT-BOX SUMMARY:')
npass = 0
for row in summary:
    try:
        ok = float(row['RMSD_compact_A']) <= 2.0
    except ValueError:
        ok = False
    npass += ok
    print(f"  {row['Target']:12s} {row['PDB']} {row['Native_Lig']:4s} RMSD={row['RMSD_compact_A']} kcal={row['kcal_compact']} box={row['Box']}")
print(f'PASS (<=2.0 A): {npass}/8')
print('PASS2-DONE')
