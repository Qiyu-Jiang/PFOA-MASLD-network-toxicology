# -*- coding: utf-8 -*-
"""MM-GBSA single-trajectory binding free-energy estimation for the three complexes."""
import json
import os
import sys
import time

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

import numpy as np
import openmm
from openmm import app, unit
from openmmforcefields.generators import SystemGenerator
from openff.toolkit import Molecule
import mdtraj as md

RUNS = r'C:\Users\user\workspace\md_runs'
WORK = r'C:\md_final_work'
OUT = os.path.join(WORK, 'mmgbsa')
os.makedirs(OUT, exist_ok=True)
KEYS = ['AKT1', 'PPARA', 'PPARG']
KJ2KCAL = 1.0 / 4.184
FRAME_IDX = list(range(1000, 1961, 40))  # 25 frames, 50-98 ns


def set_groups(system):
    for f in system.getForces():
        n = f.__class__.__name__
        if n == 'NonbondedForce':
            f.setForceGroup(1)
        elif n == 'CustomGBForce':
            f.setForceGroup(2)
        else:
            f.setForceGroup(0)


def get_platform():
    for cand in ['OpenCL', 'CPU']:
        try:
            p = openmm.Platform.getPlatformByName(cand)
            props = {'Precision': 'mixed'} if cand == 'OpenCL' else {'Threads': '16'}
            return p, props
        except Exception:
            continue
    raise RuntimeError('no platform')


def energies(ctx, positions):
    ctx.setPositions(positions)
    tot = ctx.getState(getEnergy=True).getPotentialEnergy() / unit.kilojoule_per_mole
    mm = ctx.getState(getEnergy=True, groups={0, 1}).getPotentialEnergy() / unit.kilojoule_per_mole
    gb = ctx.getState(getEnergy=True, groups={2}).getPotentialEnergy() / unit.kilojoule_per_mole
    return tot, mm, gb


plat, props = get_platform()
print('platform:', plat.getName(), flush=True)

summary = {}
for k in KEYS:
    t0 = time.time()
    base = os.path.join(RUNS, k + '_v2')
    cplx = openmm.XmlSerializer.deserialize(open(os.path.join(base, 'system.xml'), encoding='utf-8').read())
    set_groups(cplx)
    pdb = app.PDBFile(os.path.join(base, 'complex_start.pdb'))
    rec_idx = [a.index for a in pdb.topology.atoms() if a.residue.name != 'LIG']
    lig_idx = [a.index for a in pdb.topology.atoms() if a.residue.name == 'LIG']

    mol = Molecule.from_file(os.path.join(base, 'lig_fixed.sdf'), allow_undefined_stereo=True)
    mol.name = 'LIG'
    mol.assign_partial_charges(partial_charge_method='mmff94')

    # receptor-only topology
    m1 = app.Modeller(pdb.topology, pdb.positions)
    m1.delete([a for a in m1.topology.atoms() if a.residue.name == 'LIG'])
    # ligand-only topology
    m2 = app.Modeller(pdb.topology, pdb.positions)
    m2.delete([a for a in m2.topology.atoms() if a.residue.name != 'LIG'])

    sg = SystemGenerator(forcefields=['amber14-all.xml', 'implicit/obc2.xml'],
                         small_molecule_forcefield='openff-2.2.0', molecules=[mol], cache=None,
                         forcefield_kwargs={'constraints': app.HBonds, 'hydrogenMass': 3.0 * unit.amu})
    rec_sys = sg.create_system(m1.topology)
    lig_sys = sg.create_system(m2.topology)
    set_groups(rec_sys)
    set_groups(lig_sys)
    integ = openmm.VerletIntegrator(0.001)
    c_c = openmm.Context(cplx, openmm.VerletIntegrator(0.001), plat, props)
    c_r = openmm.Context(rec_sys, openmm.VerletIntegrator(0.001), plat, props)
    c_l = openmm.Context(lig_sys, openmm.VerletIntegrator(0.001), plat, props)

    traj = md.load(os.path.join(WORK, k + '_v2', 'prod', 'traj.dcd'),
                   top=os.path.join(WORK, k + '_v2', 'complex_start.pdb'))
    xyz = traj.xyz

    rows = []
    for fi in FRAME_IDX:
        pos = xyz[fi]
        e_ct, e_cm, e_cg = energies(c_c, unit.Quantity(pos, unit.nanometer))
        e_rt, e_rm, e_rg = energies(c_r, unit.Quantity(pos[rec_idx], unit.nanometer))
        e_lt, e_lm, e_lg = energies(c_l, unit.Quantity(pos[lig_idx], unit.nanometer))
        dg = (e_ct - e_rt - e_lt) * KJ2KCAL
        dg_mm = (e_cm - e_rm - e_lm) * KJ2KCAL
        dg_gb = (e_cg - e_rg - e_lg) * KJ2KCAL
        rows.append((fi, fi * 0.05, dg, dg_mm, dg_gb))
        print('   %s frame %d (%.1f ns): dG=%.1f kcal/mol (MM=%.1f, GB=%.1f)'
              % (k, fi, fi * 0.05, dg, dg_mm, dg_gb), flush=True)
    arr = np.array(rows)[:, 2:]
    summary[k] = dict(mean=float(arr[:, 0].mean()), sd=float(arr[:, 0].std()),
                      mm_mean=float(arr[:, 1].mean()), gb_mean=float(arr[:, 2].mean()),
                      n=len(rows))
    print('%s: dG_MM-GBSA = %.1f ± %.1f kcal/mol (MM %.1f, GB %.1f) [%.0fs]'
          % (k, arr[:, 0].mean(), arr[:, 0].std(), arr[:, 1].mean(), arr[:, 2].mean(), time.time() - t0), flush=True)
    with open(os.path.join(OUT, k + '_mmgbsa_frames.csv'), 'w', encoding='utf-8-sig') as f:
        f.write('frame,time_ns,dG_kcal,dG_MM_kcal,dG_GB_kcal\n')
        for fi, t, dg, dgm, dgg in rows:
            f.write('%d,%.2f,%.2f,%.2f,%.2f\n' % (fi, t, dg, dgm, dgg))

with open(os.path.join(OUT, 'mmgbsa_summary.json'), 'w', encoding='utf-8') as f:
    json.dump(summary, f, indent=2)
print('SUMMARY:', json.dumps(summary), flush=True)
print('MMGBSA DONE', flush=True)
