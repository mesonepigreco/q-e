#!/usr/bin/env python3
"""Reproducible two-atom Si TB09-PAW diagnostic (CPU, nonmagnetic).

Requires an explicitly augmented PAW UPF from prepare_paw_tau.py. A new work
folder is mandatory. The gap is searched on 81 points along Gamma-X; this
one material/path is not a validation of general PAW completeness.
"""
import argparse, hashlib, json, os, re, shutil, subprocess
from pathlib import Path
import xml.etree.ElementTree as ET

HARTREE_EV = 27.211386245988

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pw', type=Path, required=True)
    p.add_argument('--pseudo', type=Path, required=True)
    p.add_argument('--workdir', type=Path, required=True)
    p.add_argument('--cutoff', type=float, default=40.)
    p.add_argument('--density-cutoff', type=float, default=320.)
    p.add_argument('--mesh', type=int, default=8)
    a = p.parse_args()
    exe, pseudo = a.pw.resolve(), a.pseudo.resolve()
    tree = ET.parse(pseudo).getroot()
    paw = tree.find('PP_PAW')
    if paw is None or paw.get('tb09_core_tau', '').strip().upper() not in ('T', 'TRUE', '.TRUE.'):
        p.error('PAW UPF must explicitly supply frozen-core tau')
    if a.cutoff <= 0 or a.density_cutoff < a.cutoff or a.mesh < 1:
        p.error('Invalid cutoff or mesh')
    work = a.workdir.resolve(); work.mkdir(parents=True, exist_ok=False)
    shutil.copy2(pseudo, work / 'si-paw.upf')
    metadata = {'executable': str(exe), 'executable_sha256': digest(exe),
                'pseudo': str(pseudo), 'pseudo_sha256': digest(pseudo),
                'cutoff_Ry': a.cutoff, 'density_cutoff_Ry': a.density_cutoff,
                'mesh': a.mesh, 'lattice_A': 5.431, 'nat': 2,
                'scope': 'experimental reconstruction; Si Gamma-X diagnostic'}
    (work / 'manifest.json').write_text(json.dumps(metadata, indent=2) + '\n')
    common = f'''&CONTROL
 calculation='CALC', prefix='si', outdir='./scratch', pseudo_dir='./', verbosity='high'
/
&SYSTEM
 ibrav=2, A=5.431, nat=2, ntyp=1, ecutwfc={a.cutoff}, ecutrho={a.density_cutoff},
 nbnd=8, occupations='fixed', input_dft='TB09', nspin=1
/
&ELECTRONS
 conv_thr=1d-12, mixing_beta=.2, electron_maxstep=160,
 diagonalization='cg', diago_cg_maxiter=200, diago_full_acc=.true., EXTRA
/
ATOMIC_SPECIES
Si 28.0855 si-paw.upf
ATOMIC_POSITIONS crystal
Si 0 0 0
Si .25 .25 .25
'''
    env = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'}
    for name, kpoints, extra in [
        ('scf', f'K_POINTS automatic\n{a.mesh} {a.mesh} {a.mesh} 0 0 0\n', ''),
        ('bands', 'K_POINTS crystal_b\n2\n0 0 0 80\n0 .5 .5 1\n', "startingpot='file', diago_thr_init=1d-10")]:
        text = common.replace('CALC', name).replace('EXTRA', extra) + kpoints
        inp = work / (name + '.in'); out = work / (name + '.out')
        inp.write_text(text)
        with out.open('w') as stream:
            run = subprocess.run([str(exe), '-in', inp.name], cwd=work, env=env,
                                 stdout=stream, stderr=subprocess.STDOUT)
        log = out.read_text()
        if run.returncode or 'JOB DONE.' not in log or (name == 'scf' and 'convergence has been achieved' not in log):
            raise RuntimeError(f'{name} failed; inspect {out}')
        if 'negative reconstructed tau:' in log:
            raise RuntimeError(f'Unphysical reconstructed kinetic density; inspect {out}')
        shutil.copy2(work / 'scratch/si.save/data-file-schema.xml', work / (name + '.xml'))
        metadata[name + '_input_sha256'] = digest(inp)
    xml = ET.parse(work / 'bands.xml').getroot()
    states = xml.findall('.//output/band_structure/ks_energies')
    bands = [[float(x) * HARTREE_EV for x in k.find('eigenvalues').text.split()] for k in states]
    vbm = max(row[3] for row in bands); cbm = min(row[4] for row in bands)
    scf_log = (work / 'scf.out').read_text()
    metadata.update(gap_eV=cbm-vbm, vbm_eV=vbm, cbm_eV=cbm,
                    c=float(re.findall(r'TB09 PAW c =\s*([\d.Ee+-]+)', scf_log)[-1]),
                    bands_eV=bands, kpoints=[k.find('k_point').text.strip() for k in states])
    (work / 'result.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(f'TB09-PAW Si: gap={cbm-vbm:.8f} eV, c={metadata["c"]:.8f}')

if __name__ == '__main__':
    main()
