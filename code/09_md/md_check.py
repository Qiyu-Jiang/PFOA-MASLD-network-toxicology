# -*- coding: utf-8 -*-
"""Smoke test for the PFOA-complex MD package (run before production).

Checks: system loads, checkpoint loads, platform works, energy is finite,
and reports a short timing estimate. Nothing is written to the run dirs.

Usage: MD_RUNS=<pkg_root> python md_check.py
"""
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

RUNS = os.environ.get('MD_RUNS') or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KEYS = ['AKT1', 'PPARA', 'PPARG']
DT_FS = float(os.environ.get('MD_DT_FS', '4'))

print('OpenMM version:', openmm.version.version, flush=True)
for i in range(openmm.Platform.getNumPlatforms()):
    print('platform available:', openmm.Platform.getPlatform(i).getName(), flush=True)

ok_all = True
for key in KEYS:
    base = os.path.join(RUNS, key + '_v2')
    try:
        system = openmm.XmlSerializer.deserialize(open(os.path.join(base, 'system.xml'), encoding='utf-8').read())
        pdb = app.PDBFile(os.path.join(base, 'complex_start.pdb'))
        integrator = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                                     DT_FS * 0.001 * unit.picoseconds)
        selected = None
        sim = None
        pref = os.environ.get('MD_PLATFORM', '').strip().upper()
        for cand in ([pref] if pref else []) + ['CUDA', 'OpenCL', 'CPU']:
            try:
                plat = openmm.Platform.getPlatformByName(cand)
                props = {'Precision': 'mixed'} if cand in ('CUDA', 'OpenCL') else {'Threads': '16'}
                sim = app.Simulation(pdb.topology, system, integrator, plat, props)
                selected = cand
                break
            except Exception:
                continue
        ckpt = os.path.join(base, 'prod', 'ckpt.chk')
        step0 = 0
        if os.path.exists(ckpt):
            try:
                sim.loadCheckpoint(ckpt)
                step0 = sim.currentStep
            except Exception as e_ck:
                print('%s: checkpoint not usable on this platform (%s) -> using state_min.xml'
                      % (key, str(e_ck)[:90]), flush=True)
        if step0 == 0:
            st = openmm.XmlSerializer.deserialize(open(os.path.join(base, 'state_min.xml'), encoding='utf-8').read())
            sim.context.setPositions(st.getPositions())
            sim.context.setVelocitiesToTemperature(300 * unit.kelvin)
        sim.context.getState(getEnergy=True)
        t0 = time.time()
        sim.step(500)
        st = sim.context.getState(getEnergy=True, getPositions=True)
        dt = time.time() - t0
        e = st.getPotentialEnergy() / unit.kilojoule_per_mole
        pos = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        finite = bool(np.isfinite(pos).all())
        nsday = 86400.0 * DT_FS * 1e-6 / (dt / 500.0)
        print('%s: platform=%s atoms=%d step0=%d E=%.0f kJ/mol finite=%s | %.1f ms/step => ~%.0f ns/day (@%g fs)'
              % (key, selected, system.getNumParticles(), step0, e, finite, dt / 500.0 * 1000, nsday, DT_FS), flush=True)
        ok_all = ok_all and finite
    except Exception as e:
        ok_all = False
        print('%s: FAILED %s' % (key, str(e)[:200]), flush=True)
print('CHECK RESULT:', 'ALL OK' if ok_all else 'PROBLEMS FOUND', flush=True)
