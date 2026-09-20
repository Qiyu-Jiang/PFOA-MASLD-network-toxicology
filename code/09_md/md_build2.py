# -*- coding: utf-8 -*-
"""Build fixed implicit-solvent MD systems for PFOA complexes (v2).

Fixes vs v1:
  - Receptor chain-break splitting: consecutive residues whose C(i)-N(i+1)
    distance > 2.0 A are placed in separate chains (TER + new chain letters),
    preventing spurious cross-gap peptide bonds that caused huge bond energy.
  - keepIds=True on all PDB writes (preserve original residue numbering).
  - No NonbondedForce cutoff patch: both NonbondedForce and GB force use
    NoCutoff (consistent; OpenCL-compatible).
  - Progress-printing minimization + speed probe.
Targets: AKT1 (3O96), PPARA (1K7L), PPARG (2PRG). Ligand = PFOA (docked pose 1).

Usage: python md_build2.py AKT1|PPARA|PPARG|ALL
"""
import os
import subprocess
import sys
import time
import traceback

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass

DJ = r'C:\Users\user\workspace\PFOA_NAFLD_project\07_docking'
RUNS = r'C:\Users\user\workspace\md_runs'
OBABEL = r'D:\tools\python\Scripts\obabel.exe'
LOGDIR = r'C:\Users\user\workspace\tmp\md2'

SYS = {
    'AKT1': dict(rec=r'receptors\3O96_clean.pdb', chains=['A'],
                 pose=r'results\AKT1_3O96_out.pdbqt'),
    'PPARA': dict(rec=r'receptors\1K7L_clean.pdb', chains=['A'],
                  pose=r'results\PPARA_1K7L_out.pdbqt'),
    'PPARG': dict(rec=r'receptors\2PRG_clean.pdb', chains=['A'],
                  pose=r'results\PPARG_2PRG_out.pdbqt'),
}
PFOA_SMILES = '[O-]C(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'
BREAK_CUTOFF = 2.0  # angstrom


class Tee:
    def __init__(self, path):
        self.f = open(path, 'a', encoding='utf-8')
        self.stdout = sys.stdout

    def write(self, s):
        self.stdout.write(s)
        self.f.write(s)

    def flush(self):
        self.stdout.flush()
        self.f.flush()


def log(msg):
    print(msg, flush=True)


def split_receptor_text(in_path, keep_chains, out_path):
    """Select chains; split at chain breaks (C-N > cutoff); rewrite chain letters; add TER."""
    residues = []  # list of dict(ch=, id=, icode=, name=, lines=[], atoms={name: (x,y,z)})
    order_ok = True
    for ln in open(in_path, encoding='utf-8', errors='replace'):
        if ln.startswith('ATOM'):
            ch = ln[21]
            if ch not in keep_chains:
                continue
            try:
                rid = int(ln[22:26])
            except ValueError:
                continue
            name = ln[12:16].strip()
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            if not residues or (residues[-1]['ch'], residues[-1]['id']) != (ch, rid):
                residues.append(dict(ch=ch, id=rid, icode=ln[26], name=ln[17:20].strip(),
                                     lines=[], atoms={}))
            residues[-1]['lines'].append(ln.rstrip('\n'))
            residues[-1]['atoms'][name] = xyz

    # walk & split
    out = []
    letter = [0]

    def next_letter():
        c = chr(ord('A') + letter[0])
        letter[0] += 1
        return c

    seq = []          # list of (segment_letter, residue)
    cur = next_letter()
    segments = []
    seg_start = residues[0]['id'] if residues else 0
    for i, res in enumerate(residues):
        if i > 0:
            prev = residues[i - 1]
            c = prev['atoms'].get('C')
            n = res['atoms'].get('N')
            brk = False
            if c is not None and n is not None:
                d = sum((a - b) ** 2 for a, b in zip(c, n)) ** 0.5
                brk = d > BREAK_CUTOFF
            if brk or prev['ch'] != res['ch']:
                segments.append((cur, seg_start, prev['id']))
                cur = next_letter()
                seg_start = res['id']
                out.append('TER')
        for ln in res['lines']:
            out.append(ln[:21] + cur + ln[22:])
        if len(res['lines']) and res['name']:
            last = res
    if residues:
        segments.append((cur, seg_start, residues[-1]['id']))
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(out) + '\n')
    return segments


def run(key):
    import numpy as np
    import openmm
    from openmm import app, unit
    from openmmforcefields.generators import SystemGenerator
    from pdbfixer import PDBFixer
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from openff.toolkit import Molecule

    cfg = SYS[key]
    out = os.path.join(RUNS, key + '_v2')
    os.makedirs(out, exist_ok=True)
    rec_p = os.path.join(DJ, cfg['rec'])
    pose_p = os.path.join(DJ, cfg['pose'])
    log('== %s: receptor=%s pose=%s' % (key, cfg['rec'], cfg['pose']))

    # 1) receptor: select chains + split at breaks
    split_p = os.path.join(out, 'receptor_split.pdb')
    segments = split_receptor_text(rec_p, cfg['chains'], split_p)
    log('segments (%d): %s' % (len(segments), segments))
    n_chains = len(segments)

    # 2) PDBFixer cleanup
    rec_fix = os.path.join(out, 'receptor_fixed.pdb')
    fixer = PDBFixer(filename=split_p)
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    log('missing residues: %d' % len(fixer.missingResidues))
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    with open(rec_fix, 'w') as fh:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    rs = [r.id for r in fixer.topology.residues()]
    log('receptor fixed: %d atoms, %d residues, first=%s last=%s (keepIds)'
        % (fixer.topology.getNumAtoms(), len(rs), rs[0], rs[-1]))

    # 3) ligand from docked pose (mode 1)
    lig_sdf_raw = os.path.join(out, 'lig_raw.sdf')
    subprocess.run([OBABEL, pose_p, '-O', lig_sdf_raw, '-f', '1', '-l', '1'], check=True)
    ref = Chem.MolFromSmiles(PFOA_SMILES)
    raw = Chem.MolFromMolFile(lig_sdf_raw, sanitize=False, removeHs=False)
    raw.UpdatePropertyCache(strict=False)
    fixed = AllChem.AssignBondOrdersFromTemplate(ref, raw)
    molH = Chem.AddHs(fixed, addCoords=True)
    lig_fixed = os.path.join(out, 'lig_fixed.sdf')
    Chem.MolToMolFile(molH, lig_fixed)
    log('ligand atoms: %d' % molH.GetNumAtoms())

    lig_pdb = os.path.join(out, 'lig_fixed.pdb')
    subprocess.run([OBABEL, lig_fixed, '-O', lig_pdb], check=True)
    _lines = []
    _cnt = {}
    for _ln in open(lig_pdb, encoding='utf-8', errors='replace').read().splitlines():
        if _ln.startswith(('ATOM', 'HETATM')) and len(_ln) >= 20:
            _name = _ln[12:16].strip()
            _el = ''.join(ch for ch in _name if ch.isalpha()) or 'X'
            _cnt[_el] = _cnt.get(_el, 0) + 1
            _new_name = '%s%d' % (_el, _cnt[_el])
            _ln = _ln[:12] + ('%-4s' % _new_name) + _ln[16:17] + 'LIG' + _ln[20:21] + 'L' + _ln[22:]
        _lines.append(_ln)
    with open(lig_pdb, 'w', encoding='utf-8') as _fh:
        _fh.write('\n'.join(_lines) + '\n')

    mol = Molecule.from_file(lig_fixed, allow_undefined_stereo=True)
    mol.name = 'LIG'
    try:
        mol.assign_partial_charges(partial_charge_method='mmff94')
        log('charges: mmff94, net=%.3f' % mol.total_charge.m)
    except Exception as e1:
        log('mmff94 failed (%s) -> gasteiger' % str(e1)[:80])
        mol.assign_partial_charges(partial_charge_method='gasteiger')

    # 4) complex
    rec2 = app.PDBFile(rec_fix)
    lig2 = app.PDBFile(lig_pdb)
    modeller = app.Modeller(rec2.topology, rec2.positions)
    modeller.add(lig2.topology, lig2.positions)
    complex_pdb = os.path.join(out, 'complex_start.pdb')
    with open(complex_pdb, 'w') as fh:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)
    log('complex atoms: %d' % modeller.topology.getNumAtoms())

    # 5) system (implicit solvent, NoCutoff everywhere; OBC2 GB; HMR for 4 fs)
    sysgen = SystemGenerator(forcefields=['amber14-all.xml', 'implicit/obc2.xml'],
                             small_molecule_forcefield='openff-2.2.0',
                             molecules=[mol], cache=None,
                             forcefield_kwargs={'constraints': app.HBonds,
                                                'hydrogenMass': 3.0 * unit.amu})
    system = sysgen.create_system(modeller.topology)
    for f in system.getForces():
        extra = ''
        if isinstance(f, openmm.NonbondedForce):
            extra = ' method=%s' % f.getNonbondedMethod()
        if f.__class__.__name__ == 'CustomGBForce':
            extra = ' method=%s' % getattr(f, 'getNonbondedMethod', lambda: '?')()
        log('  force: %-22s%s' % (f.__class__.__name__, extra))
    system_xml = os.path.join(out, 'system.xml')
    with open(system_xml, 'w') as fh:
        fh.write(openmm.XmlSerializer.serialize(system))
    log('system particles: %d' % system.getNumParticles())

    # 6) QC: bond deviations vs equilibrium lengths
    pos_arr = np.array([list(p.value_in_unit(unit.angstrom)) for p in modeller.positions])
    worst = []
    for f in system.getForces():
        if f.__class__.__name__ == 'HarmonicBondForce':
            for k in range(f.getNumBonds()):
                i, j, r0, kk = f.getBondParameters(k)
                d = float(np.linalg.norm(pos_arr[i] - pos_arr[j]))
                worst.append((abs(d - r0.value_in_unit(unit.angstrom)), d, r0.value_in_unit(unit.angstrom), i, j))
    worst.sort(key=lambda t: -t[0])
    atoms = list(modeller.topology.atoms())
    log('QC bond check: top-3 deviations:')
    for dev, d, r0a, i, j in worst[:3]:
        log('   dev=%.3f d=%.3f r0=%.3f | %s%s %s -- %s%s %s' % (
            dev, d, r0a, atoms[i].residue.name, atoms[i].residue.id, atoms[i].name,
            atoms[j].residue.name, atoms[j].residue.id, atoms[j].name))
    nbad = sum(1 for w in worst if w[0] > 0.5)
    log('QC bond check: bonds with dev>0.5A: %d %s' % (nbad, 'PASS' if nbad == 0 else 'FAIL'))

    # 7) minimize with progress + speed probe (OpenCL preferred: GB force is ~400x faster than CPU here)
    integrator = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                                 0.004 * unit.picoseconds)
    try:
        platform = openmm.Platform.getPlatformByName('OpenCL')
        sim = app.Simulation(modeller.topology, system, integrator, platform, {'Precision': 'mixed'})
        log('platform: OpenCL (mixed)')
    except Exception as e_plat:
        log('OpenCL failed (%s) -> CPU fallback' % str(e_plat)[:120])
        platform = openmm.Platform.getPlatformByName('CPU')
        sim = app.Simulation(modeller.topology, system, integrator, platform, {'Threads': '16'})
    sim.context.setPositions(modeller.positions)
    t0 = time.time()
    prev = None
    for rnd in range(10):
        openmm.LocalEnergyMinimizer.minimize(sim.context, maxIterations=250)
        st = sim.context.getState(getEnergy=True)
        e = st.getPotentialEnergy() / unit.kilojoule_per_mole
        log('minimize round %d: E = %.1f kJ/mol (%.1fs)' % (rnd + 1, e, time.time() - t0))
        if prev is not None and abs(prev - e) < 2.0:
            break
        prev = e
    state = sim.context.getState(getPositions=True)
    with open(os.path.join(out, 'state_min.xml'), 'w') as fh:
        fh.write(openmm.XmlSerializer.serialize(state))

    # speed probe (synced)
    sim.context.getState(getEnergy=True)
    t0 = time.time()
    sim.step(500)
    sim.context.getState(getEnergy=True)
    dt_a = time.time() - t0
    sim.context.getState(getEnergy=True)
    t0 = time.time()
    sim.step(1000)
    sim.context.getState(getEnergy=True)
    dt_b = time.time() - t0
    sps = dt_b / 1000.0
    log('speed: %.2f ms/step (500: %.2fs, 1000: %.2fs) => ~%.1f ns/day @4fs (synced)'
        % (sps * 1000, dt_a, dt_b, 0.3456 / sps))
    log('DONE2 %s' % key)


if __name__ == '__main__':
    os.makedirs(LOGDIR, exist_ok=True)
    arg = sys.argv[1] if len(sys.argv) > 1 else 'ALL'
    keys = list(SYS) if arg == 'ALL' else [arg]
    for k in keys:
        sys.stdout = Tee(os.path.join(LOGDIR, 'build2_%s.log' % k))
        try:
            run(k)
        except Exception:
            traceback.print_exc()
        sys.stdout.flush()
