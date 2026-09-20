# -*- coding: utf-8 -*-
"""Production MD runner (OpenCL) for PFOA complexes — resumable.

Usage: MD_DT_FS=4 python md_prod2.py <AKT1|PPARA|PPARG> [total_ns]

- Loads md_runs/<KEY>_v2/{system.xml, state_min.xml}
- Fresh start: short NVT equilibration (50 ps, Langevin 300 K), then production.
- Resumable: re-run finds prod/ckpt.chk and continues from it.
- Saves: prod/traj.dcd (every 0.05 ns), prod/ckpt.chk (every 0.05 ns), prod/prod.log
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


def main():
    key = sys.argv[1]
    total_ns = float(sys.argv[2]) if len(sys.argv) > 2 else 50.0
    dt_fs = float(os.environ.get('MD_DT_FS', '4'))
    ns_per_step = dt_fs * 1e-6
    maxh = float(os.environ.get('MD_MAXH', '0'))  # 0 = unlimited

    base = os.path.join(RUNS, key + '_v2')
    prod = os.path.join(base, 'prod')
    os.makedirs(prod, exist_ok=True)
    sys.stdout = Tee(os.path.join(prod, 'prod.log'))

    system = openmm.XmlSerializer.deserialize(open(os.path.join(base, 'system.xml'), encoding='utf-8').read())
    pdb = app.PDBFile(os.path.join(base, 'complex_start.pdb'))
    integrator = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                                                 dt_fs * 0.001 * unit.picoseconds)
    pref = os.environ.get('MD_PLATFORM', '').strip().upper()
    order = ([pref] if pref else []) + ['CUDA', 'OpenCL', 'CPU']
    sim = None
    for cand in order:
        try:
            plat = openmm.Platform.getPlatformByName(cand)
            props = {'Precision': 'mixed'} if cand in ('CUDA', 'OpenCL') else {'Threads': '16'}
            sim = app.Simulation(pdb.topology, system, integrator, plat, props)
            print('[%s] platform: %s' % (key, cand), flush=True)
            break
        except Exception as e_plat:
            print('[%s] platform %s unavailable: %s' % (key, cand, str(e_plat)[:90]), flush=True)
    if sim is None:
        raise RuntimeError('no usable compute platform')

    ckpt = os.path.join(prod, 'ckpt.chk')
    dcd = os.path.join(prod, 'traj.dcd')
    stride = int(0.05 * 1e6 / dt_fs)       # 0.05 ns
    chunk = int(0.1 * 1e6 / dt_fs)         # 0.1 ns
    equil_steps = int(0.05 * 1e6 / dt_fs)  # 50 ps
    total_steps = int(total_ns * 1e6 / dt_fs)
    t_start = time.time()

    resumed = os.path.exists(ckpt)
    if resumed:
        try:
            sim.loadCheckpoint(ckpt)
            step0 = sim.currentStep
            print('[%s] resumed from checkpoint at step %d (%.3f ns) @ %.1f fs' % (key, step0, step0 * ns_per_step, dt_fs), flush=True)
        except Exception as e_ck:
            resumed = False
            print('[%s] WARNING: checkpoint unusable on this platform (%s) -> fresh start from state_min.xml'
                  % (key, str(e_ck)[:120]), flush=True)
    if not resumed:
        st = openmm.XmlSerializer.deserialize(open(os.path.join(base, 'state_min.xml'), encoding='utf-8').read())
        sim.context.setPositions(st.getPositions())
        sim.context.setVelocitiesToTemperature(300 * unit.kelvin)
        step0 = 0
        if os.path.exists(dcd):
            try:
                os.replace(dcd, dcd + '.localpilot')
                print('[%s] previous trajectory kept as %s' % (key, os.path.basename(dcd) + '.localpilot'), flush=True)
            except Exception:
                pass
        t0 = time.time()
        sim.step(equil_steps)
        st2 = sim.context.getState(getEnergy=True)
        e = st2.getPotentialEnergy() / unit.kilojoule_per_mole
        print('[%s] equilibrated 50 ps (E=%.0f kJ/mol, %.0fs) @ %.1f fs' % (key, e, time.time() - t0, dt_fs), flush=True)
        sim.saveCheckpoint(ckpt)

    sim.reporters.append(app.DCDReporter(dcd, stride, append=resumed))
    sim.reporters.append(app.CheckpointReporter(ckpt, stride))

    wall0 = time.time()
    step_now = sim.currentStep
    while step_now < total_steps:
        if maxh > 0 and (time.time() - t_start) > maxh * 3600:
            print('[%s] wall-time cap reached (%.1fh) — pausing at %.3f ns (resumable)' % (key, maxh, step_now * ns_per_step), flush=True)
            break
        n_this = min(chunk, total_steps - step_now)
        sim.step(n_this)
        step_now = sim.currentStep
        el = time.time() - wall0
        done_ns = (step_now - step0) * ns_per_step
        rate = done_ns / (el / 86400.0) if el > 0 else 0.0
        st = sim.context.getState(getEnergy=True, getPositions=True)
        e = st.getPotentialEnergy() / unit.kilojoule_per_mole
        pos = st.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        finite = bool(np.isfinite(pos).all())
        print('[%s] step %d | %.2f/%.1f ns | E=%.0f kJ/mol | %.1f ns/day | wall %.2fh | finite=%s'
              % (key, step_now, step_now * ns_per_step, total_ns, e, rate, el / 3600.0, finite), flush=True)
        if not finite:
            print('[%s] WARNING: non-finite positions — stopping for inspection' % key, flush=True)
            break
    print('[%s] production loop end at %d steps (%.2f ns)' % (key, step_now, step_now * ns_per_step), flush=True)


if __name__ == '__main__':
    main()
