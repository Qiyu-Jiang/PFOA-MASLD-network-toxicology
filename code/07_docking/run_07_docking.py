# -*- coding: utf-8 -*-
"""run_07: batch molecular docking of 8 targets (AutoDock Vina 1.2.5) + result aggregation"""
import os, subprocess, re, csv, time

BASE = r'C:\Users\user\workspace\PFOA_NAFLD_project'
DK = os.path.join(BASE, '07_docking')
REC = os.path.join(DK, 'receptors')
OUT = os.path.join(DK, 'results')
TMP = r'C:\Users\user\workspace\tmp'
os.makedirs(OUT, exist_ok=True)
VINA = os.path.join(DK, 'vina.exe')
LIG = os.path.join(TMP, 'PFOA_lig_v4.pdbqt')
DATADIR = r'D:\tools\python\site-packages'
ENV = dict(os.environ)

TARGETS = [
    ('PPARA', '1K7L', (-17.6, -14.4, -5.4), '544 (agonists)'),
    ('PPARG', '2PRG', (49.7, -37.0, 19.3), 'BRL (rosiglitazone)'),
    ('FXR_NR1H4', '4OIV', (12.9, 14.1, -20.8), 'XX9 (agonist)'),
    ('EGFR', '1M17', (22.0, 0.3, 52.8), 'AQ4 (erlotinib)'),
    ('CASP3', '1GFW', (37.6, 33.6, 28.0), 'MSI (inhibitor)'),
    ('MMP9', '1GKC', (53.3, 22.5, 129.7), 'NFH (hydroxamate)'),
    ('TNF', '2AZ5', (-19.4, 74.7, 33.8), '307 (SPD-304)'),
    ('AKT1', '3O96', (8.4, -6.8, 12.6), 'IQO (inhibitor)'),
]

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
    r = subprocess.run(['obabel', clean_pdb, '-O', dst, '-xr', '--partialcharge', 'gasteiger'],
                       capture_output=True, text=True, timeout=600)
    return dst, (r.stderr or r.stdout).strip()

summary = []
for name, pid, (cx, cy, cz), pocket in TARGETS:
    print('=' * 70)
    print(f'--- {name} ({pid}) | pocket: {pocket} ---')
    # NOTE: MMP9 (1GKC) is docked WITHOUT its catalytic zinc (removed during receptor preparation);
    # this matches the manuscript (a metal-site limitation). Set keep_zn=True to dock with the zinc.
    keep_zn = (pid == '1GKC' and False)
    clean = clean_receptor(pid, keep_zn)
    rec_pdbqt, conv_msg = pdbqt_receptor(clean, pid)
    print('receptor pdbqt:', os.path.getsize(rec_pdbqt), 'bytes |', conv_msg[-80:])

    out_pdbqt = os.path.join(OUT, f'{name}_{pid}_out.pdbqt')
    log = os.path.join(OUT, f'{name}_{pid}_vina.log')
    cmd = [VINA,
           '--receptor', rec_pdbqt, '--ligand', LIG,
           '--center_x', str(cx), '--center_y', str(cy), '--center_z', str(cz),
           '--size_x', '25', '--size_y', '25', '--size_z', '25',
           '--exhaustiveness', '32', '--num_modes', '10', '--seed', '42', '--cpu', '8',
           '--out', out_pdbqt]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    dt = time.time() - t0
    stdout = r.stdout
    open(log, 'w', encoding='utf-8').write(stdout + '\n' + (r.stderr or ''))
    modes = re.findall(r'^\s+(\d+)\s+(-?\d+\.\d+)\s+', stdout, re.M)
    affs = [float(m[1]) for m in modes]
    print(f'docked in {dt:.0f}s | modes: {len(affs)} | best: {min(affs) if affs else "FAIL"} kcal/mol')
    if affs:
        top3 = sorted(affs)[:3]
        print('top3:', top3)
        summary.append({'Target': name, 'PDB': pid, 'Pocket': pocket,
                        'Best_kcal_mol': top3[0],
                        'Top2': top3[1] if len(top3) > 1 else '',
                        'Top3': top3[2] if len(top3) > 2 else '',
                        'Secs': round(dt)})
    else:
        print('no modes parsed! stderr:', (r.stderr or '')[:200])
        summary.append({'Target': name, 'PDB': pid, 'Pocket': pocket,
                        'Best_kcal_mol': 'FAIL', 'Top2': '', 'Top3': '', 'Secs': round(dt)})

# aggregate
csv_path = os.path.join(OUT, 'docking_summary.csv')
with open(csv_path, 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.DictWriter(f, fieldnames=['Target', 'PDB', 'Pocket', 'Best_kcal_mol', 'Top2', 'Top3', 'Secs'])
    w.writeheader(); w.writerows(summary)
print()
print('=' * 70)
print('SUMMARY:')
for row in summary:
    print(f"  {row['Target']:12s} {row['PDB']}  best={row['Best_kcal_mol']} kcal/mol  (top2={row['Top2']}, top3={row['Top3']})")
print('saved:', csv_path)

# plot (binding-energy bar chart)
ok_rows = [r for r in summary if r['Best_kcal_mol'] != 'FAIL']
if ok_rows:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = ['Microsoft YaHei', 'SimHei', 'sans-serif']
    rows = sorted(ok_rows, key=lambda r: r['Best_kcal_mol'])
    fig, ax = plt.subplots(figsize=(9.5, 5.5), dpi=220)
    vals = [r['Best_kcal_mol'] for r in rows]
    labels = [f"{r['Target']} ({r['PDB']})" for r in rows]
    colors = ['#145c4d' if v <= -7 else ('#3d7a68' if v <= -5.5 else '#9db3ac') for v in vals]
    ax.barh(labels[::-1], [-v for v in vals[::-1]], color=colors[::-1], height=0.6)
    for i, v in enumerate(vals[::-1]):
        ax.text(-v + 0.05, i, f'{v:.1f}', va='center', fontsize=9, color='#333')
    ax.set_xlabel('binding affinity, -kcal/mol (longer = stronger)')
    ax.set_title('PFOA (deprotonated) docking against 8 targets — AutoDock Vina best mode', fontsize=11.5)
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    pp = os.path.join(BASE, '05_figures', 'run07_docking_barplot.png')
    fig.savefig(pp, facecolor='white')
    print('saved figure:', pp)
print('RUN07-DONE')
