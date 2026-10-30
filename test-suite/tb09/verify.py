#!/usr/bin/env python3
"""Reproduce Si TB09 gaps and optionally check MPI/spin/Gamma consistency.

Requires Python 3.8+, a CPU/Libxc-enabled pw.x and Si.pbe-rrkj.UPF.
Run in a fresh work directory; no automatic downloads or fitted parameters.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess
import xml.etree.ElementTree as ET

HA_EV = 27.211386245988
PSEUDO_SHA = "d42702eef81fd417a503bc33e1756081d012b398c46ad9c2955aad963a1760c7"


def local(tag):
    return tag.split("}")[-1]


def eigenvalues(xml):
    return [[float(x) for x in next(e for e in b if local(e.tag) == "eigenvalues").text.split()]
            for b in ET.parse(xml).iter() if local(b.tag) == "ks_energies"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pw", type=Path, required=True)
    ap.add_argument("--pseudo", type=Path, required=True)
    ap.add_argument("--workdir", type=Path, required=True)
    ap.add_argument("--resume", action="store_true", help="reuse completed jobs with identical inputs and binaries")
    ap.add_argument("--cutoff", type=int, default=40)
    ap.add_argument("--mesh", type=int, default=8)
    ap.add_argument("--functional", choices=["TB09", "PBE"], default="TB09")
    ap.add_argument("--mpi-pw", type=Path)
    ap.add_argument("--launcher", default="mpiexec -n 2")
    ap.add_argument("--checks", action="store_true", help="also check spin, Gamma, and unsupported forces")
    args = ap.parse_args()
    pw, pseudo = args.pw.resolve(), args.pseudo.resolve()
    if hashlib.sha256(pseudo.read_bytes()).hexdigest() != PSEUDO_SHA:
        ap.error("Use the specified Si.pbe-rrkj.UPF; its checksum is part of this benchmark")
    root = args.workdir.resolve()
    root.mkdir(parents=True, exist_ok=args.resume)
    signature = {"pw_sha256": hashlib.sha256(pw.read_bytes()).hexdigest(),
                 "mpi_pw_sha256": hashlib.sha256(args.mpi_pw.resolve().read_bytes()).hexdigest() if args.mpi_pw else None,
                 "pseudo_sha256": PSEUDO_SHA, "cutoff": args.cutoff, "mesh": args.mesh,
                 "functional": args.functional, "launcher": args.launcher, "checks": args.checks}
    config = root / "configuration.json"
    if args.resume:
        if not config.exists() or json.loads(config.read_text()) != signature:
            ap.error("Cannot resume: configuration or binary checksum differs")
    else:
        config.write_text(json.dumps(signature, indent=2) + "\n")
    records = {}

    def run(label, command, spin=1, gamma=False, bands=False, force=False):
        path = root / label
        path.mkdir(exist_ok=args.resume)
        inp = f"""&CONTROL
 calculation='scf', prefix='si', outdir='./scratch',
 pseudo_dir='{pseudo.parent}', verbosity='high', tprnfor={'.true.' if force else '.false.'}
/
&SYSTEM
 ibrav=2, A=5.431, nat=2, ntyp=1,
 ecutwfc={args.cutoff}, ecutrho={4*args.cutoff}, nbnd=8,
 occupations='fixed', input_dft='{args.functional}', nspin={spin}
 {('tot_magnetization=0,' if spin == 2 else '')}
/
&ELECTRONS
 conv_thr=1.0d-12, mixing_beta=0.2, electron_maxstep=160,
 diagonalization='david', diago_full_acc=.true.
/
ATOMIC_SPECIES
Si 28.0855 {pseudo.name}
ATOMIC_POSITIONS crystal
Si 0.0 0.0 0.0
Si 0.25 0.25 0.25
"""
        inp += "K_POINTS gamma\n" if gamma else f"K_POINTS automatic\n{args.mesh} {args.mesh} {args.mesh} 0 0 0\n"

        reused = set()

        def execute(name, text):
            input_path, output_path = path / (name + ".in"), path / (name + ".out")
            if (args.resume and input_path.exists() and output_path.exists()
                    and (path / (name + "-schema.xml")).exists() and input_path.read_text() == text):
                cached = output_path.read_text()
                if "JOB DONE." in cached and (name != "scf" or "convergence has been achieved" in cached):
                    print("Reusing", label, name, flush=True)
                    reused.add(name)
                    return cached
            input_path.write_text(text)
            print("Running", label, name, flush=True)
            with (path / (name + ".out")).open("w") as out:
                proc = subprocess.run(command + ["-in", name + ".in"], cwd=path,
                                      stdout=out, stderr=subprocess.STDOUT, timeout=1800)
            output = (path / (name + ".out")).read_text()
            if force:
                if proc.returncode == 0 or "total-energy forces and stress are not available" not in output:
                    raise RuntimeError("TB09 forces were not rejected correctly")
                return output
            if proc.returncode or "JOB DONE." not in output:
                raise RuntimeError(str(path) + ": " + output[-2000:])
            if name == "scf" and "convergence has been achieved" not in output:
                raise RuntimeError(str(path) + ": SCF did not converge")
            return output

        out = execute("scf", inp)
        if force:
            records[label] = {"forces_rejected": True}
            return None
        xml = path / "scratch/si.save/data-file-schema.xml"
        scf_xml = path / "scf-schema.xml"
        if "scf" not in reused:
            scf_xml.write_bytes(xml.read_bytes())
        ev = eigenvalues(scf_xml)
        cs = re.findall(r"TB09 c =\s*([0-9.+Ee-]+)", out)
        records[label] = {"c": float(cs[-1]) if cs else None}
        if bands:
            binput = inp.replace("calculation='scf'", "calculation='bands'")
            binput = binput.replace("conv_thr=", "startingpot='file', conv_thr=")
            binput = binput[:binput.index("K_POINTS")] + "K_POINTS crystal_b\n2\n0 0 0 80\n0 0.5 0.5 1\n"
            execute("bands", binput)
            band_xml = path / "bands-schema.xml"
            if "bands" not in reused:
                band_xml.write_bytes(xml.read_bytes())
            be = eigenvalues(band_xml)
            records[label]["gap_eV"] = (min(e[4] for e in be) - max(e[3] for e in be)) * HA_EV
        return ev

    def compare(label, ev, ref, tolerance=1.e-7):
        if len(ev) != len(ref) or any(len(a) != len(b) for a, b in zip(ev, ref)):
            raise RuntimeError(label + ": eigenvalue shape mismatch")
        err = max(abs(a-b) for ea, eb in zip(ev, ref) for a, b in zip(ea, eb))
        records[label]["max_eigenvalue_difference_Ha"] = err
        if err > tolerance:
            raise RuntimeError(label + ": eigenvalues differ by " + str(err) + " Ha")
        if args.functional == "TB09" and abs(records[label]["c"] - records["serial"]["c"]) > tolerance:
            raise RuntimeError(label + ": c differs from serial")

    ref = run("serial", [str(pw)], bands=True)
    if args.mpi_pw:
        mpi = shlex.split(args.launcher) + [str(args.mpi_pw.resolve())]
        for pools in [1, 2]:
            label = "mpi-pools-" + str(pools)
            compare(label, run(label, mpi + ["-npool", str(pools)]), ref)
    if args.checks:
        spin = run("spin-2", [str(pw)], spin=2)
        compare("spin-2", spin, [e + e for e in ref], tolerance=1.e-6)
        gamma = run("gamma", [str(pw)], gamma=True)
        if args.mpi_pw:
            ge = run("gamma-mpi", mpi, gamma=True)
            err = max(abs(a-b) for ea, eb in zip(ge, gamma) for a, b in zip(ea, eb))
            records["gamma-mpi"]["max_eigenvalue_difference_Ha"] = err
            if err > 1.e-7 or (args.functional == "TB09" and
                              abs(records["gamma-mpi"]["c"] - records["gamma"]["c"]) > 1.e-7):
                raise RuntimeError("Gamma serial/MPI mismatch")
        if args.functional == "TB09":
            run("forces", [str(pw)], force=True)
    report = {"functional": args.functional, "cutoff_Ry": args.cutoff, "kmesh": args.mesh,
              "pseudo_sha256": PSEUDO_SHA, "results": records}
    (root / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
