# Experimental reconstructed-density TB09–PAW path

This branch implements a dedicated multiplicative model-potential transform.
It does not implement a general variational meta-GGA PAW functional. ABINIT
is used as an independent reference executable; its TB09 algorithm is not
transplanted into QE. No gap is used to fit `c` or any other parameter.

Supported scope is currently CPU, nonmagnetic (`nspin=1`), all-PAW species,
without spin orbit. Non-overlapping augmentation spheres are the previously
tested scope; a new coreless-hydrogen overlap experiment is described below. Forces,
stress and relaxation remain disabled. TB09 has no parent exchange energy;
printed SCF energies must not be used as material total energies. Optical
matrix elements and PAW linear response have not been validated by this work.
Do not infer validated dielectric constants from a validated eigenvalue gap.

## Density and operator construction

At every potential update the code constructs:

```
n_AE   = n_FFT + n_core_PS + sum_a(n_val_AE^a - n_val_PS^a - n_hat^a
                                  + n_core_AE^a - n_core_PS^a)
tau_AE = tau_PW + tau_core_PS + sum_a(tau_val_AE^a - tau_val_PS^a
                                    + tau_core_AE^a - tau_core_PS^a)
```

`n_FFT` already includes compensation charges `n_hat`. There is no fictitious
kinetic density associated with that compensation charge. Kinetic density
means `tau = 0.5 sum |grad(phi)|^2` in Hartree/bohr^3. Occupation matrices
already contain twice the off-diagonal occupations. Angular derivatives use
Gaunt products and the angular-Laplacian product identity.

Density gradients and Laplacians are calculated from this same density
reconstruction. The smooth periodic crystalline environment is interpolated
onto atomic quadratures. The original nonlinear `|grad(n_AE)|/n_AE` is
integrated after reconstruction, using a partition of unity, and divided by
the full cell volume. The published coefficient is then applied once:

```
c = -0.012 + 1.023 sqrt(cell_average)
```

The atomic partition equals one below 0.7 times the augmentation radius and
tapers to zero with a quintic polynomial at the sphere boundary. Its
complement is used on the FFT grid. A partition-volume diagnostic is printed;
this diagnostic alone does not establish accuracy of the nonlinear integral.
The angular quadrature currently has maximum angular degree at least 16.

An arbitrary regular LDA reference replaces the singular potential inside
atomic partitions on the FFT grid. The one-center PAW matrices use the true
reconstructed AE TB09 potential and subtract the *actual interpolated FFT
reference* including its compensation-charge matrix elements. The latter
cancel the corresponding FFT augmentation terms in QE's operator assembly.
The smooth reference is not a second TB09 evaluation on an isolated pseudo
atom. All TB09 kinetic-operator derivatives are zero: tau is an input to a
multiplicative potential, not an energy derivative.

## Atomic data

Existing PAW UPFs generally do not supply both AE and pseudo frozen-core
kinetic densities. The prototype adds explicit optional `tb09_core_tau`
metadata and `PP_AE_TAUCORE` / `PP_PS_TAUCORE` fields to `PP_PAW`. It rejects
missing data rather than treating an unknown core as having zero tau.
These are experimental extensions, not an upstream UPF schema standard.

For a matching open PAW-XML/UPF pair:

```
python3 test-suite/tb09/prepare_paw_tau.py Si.xml Si.UPF Si-TB09.UPF
```

The converter checks species, valence, core normalization and AE/pseudo
partial waves, converts XML Y00 coefficients to scalar densities, interpolates
to the UPF mesh, and records input/output checksums. It does not change the
projectors or atomic potential. It requires NumPy and SciPy. Do not use a
similar-looking XML from a different atomic dataset.

The Si diagnostic uses the open JTH-PBE 2.0 XML/UPF pair, pinned to atomic-data
repository commit `a03e3f0b7ac8395550cca88880d8cb7dcd478bb3`.
Original XML SHA256:
`79075fcce1418817b2dce0ebe90d6a970776fc0832375a2ef99af0ac9dbe863a`.
Original UPF SHA256:
`379e58f289bb863c874b941811b33513c369a5ffae86e61108b132730cec8568`.
Core charge is 9.9999999973 electrons; partial-wave pairing error is below
3.6e-7 relative. These checks validate data conventions, not the full PAW model.

## Reproduce and assess

```
python3 test-suite/tb09/paw_math_test.py
python3 test-suite/tb09/verify_paw.py --pw build/bin/pw.x \
  --pseudo Si-TB09.UPF --workdir /tmp/si-tb09-paw-40 \
  --cutoff 40 --density-cutoff 320 --mesh 8
```

Each work directory must be new. The benchmark uses two-atom diamond Si,
a=5.431 Angstrom, eight bands, Gamma-centered SCF mesh and 81 Gamma-X band
points. It saves independent SCF/band XML, binary/input/pseudopotential
checksums, `c` and the extracted gap. An SCF failure or diagnosed negative
reconstructed kinetic density rejects the benchmark result.

Analytic tests compile the actual radial and kinetic-density helpers with
bounds checks. They cover polynomial differentiation, a hydrogenic cusp,
Hermite derivatives, a Cartesian mixed s-p orbital with analytic gradient,
and exact nearest-image search in a strongly skew lattice. They do not prove
PAW partial-wave completeness or validate every angular channel/material.

Required evidence before promoting this prototype to production includes
cutoff and k-mesh convergence, angular-quadrature and smooth-reference
sensitivity, checks of core tails / periodic images for each dataset, and
broader material coverage. The present frozen smooth-core tail treatment uses
the nearest image of each atom; tails from additional images are neglected.
This approximation must be bounded or replaced before claiming a general
periodic all-electron reconstruction. The previous validation covers only non-overlapping augmentation spheres;
this is especially relevant to high-pressure systems.

Borlido/Botti's Si references are 1.1994 eV for VASP PAW and 1.176 eV for
WIEN2k all-electron in the [2020 validation supplement](https://doi.org/10.1021/acs.jctc.0c00214.s004). They are distinct
numerical references. The current 5.431 Angstrom geometry is a diagnostic,
not a recovered copy of their ICSD 51688 CIF. Agreement with one Si number
does not reproduce their full benchmark or establish superiority.


## Si diagnostic measured on 2026-10-30

| Calculation | Gap (eV) |
| --- | ---: |
| QE reconstructed PAW, 40/320 Ry, 8^3, 81 Gamma-X points | 1.1324453 |
| QE reconstructed PAW, 60/640 Ry, 8^3, 81 Gamma-X points | 1.1323423 |
| Earlier corrected NC, 60 Ry, 12^3 | 0.9238103 |
| Borlido et al. 2020 WIEN2k mBJLDA reference | 1.176 |
| Borlido et al. VASP mBJ@LDA reference | 1.1994 |

The PAW cutoff refinement changes the gap by 0.000103 eV. A much coarser
4^3 SCF mesh gave 1.1763 eV; that apparent agreement with the AE reference
was removed by k-mesh refinement. No asymptotic k-mesh bound is established.
The remaining AE discrepancy is 0.04366 eV and the VASP discrepancy 0.06706 eV.
The PAW–NC difference is evidence of density-representation sensitivity;
these settings/datasets do not isolate its full cause.

ABINIT 10.8.3 / Libxc 7 passes its official low-cutoff Si PAW TB09 test (JTH1.0,
6/10 Ha, 2^3 mesh, converged in 19 iterations). Attempts with the paired JTH2
at 20/80 Ha and 4^3 mesh, including PBE restarts, two compensation conventions
and a kinetic-density floor diagnostic, did not yield an accepted converged
gap. Raising the official JTH1.0 case to these settings also failed to
converge. Consequently no matched quantitative QE–ABINIT PAW agreement is
claimed. The original unmodified ABINIT source supplies a reference executable,
not part of the QE implementation.

The numerical checks and explicit limits are recorded in
[paw-validation.json](paw-validation.json). No large hydrogen PW runs were
performed and no upstream merge request was opened.


## Hydrogen overlap experiment: checkpoint, 2026-10-30

The selected hydrogen dataset is the locally tested
`H.pbe-kjpaw_psl.1.0.0.UPF`, with original SHA256
`dd48be99c7a91c4163ed123cdac938af189a9f62f382c41e3f664d00a586839d`.
Its augmentation radius is approximately 1.2523 bohr, so molecular hydrogen
at the intended pressures has overlapping augmentation spheres. The choice
is based on the available local energetic tests, not a demonstrated QMC
ranking of several PAW datasets. No refractive-index result is available yet.

This checkpoint adds an experimental overlap partition:
`w_out = product(1-theta_a)` and
`w_a = (1-w_out)*theta_a/sum(theta_a)`. All atomic field corrections are
summed before evaluating the nonlinear model. Coreless neighbor filtering,
a bounding-box rejection and a spherical-channel fast path reduce its cost.
The overlap branch currently permits only coreless hydrogen, rejects spheres
containing another nucleus, and retains the self-image overlap guard.

The serial `qe_pw_exe` target compiles successfully. The new overlap path has
not yet been tested with an SCF calculation, analytic partition checks or
optical calculations. The Si numbers above predate this checkpoint. In
particular, partition normalization does not establish the validity of an
overlapping one-center PAW reconstruction: cross-atom reconstruction terms
and positivity of reconstructed fields require investigation. Treat this
commit as a work checkpoint, not a validated high-pressure PAW implementation.

Resume with small hydrogen checks, then the pristine C2/c-24 SSCHA centroids
at 150, 175, 200, 225 and 250 GPa. Compare historical TB09, BLYP plus a 0.8 eV
scissor, corrected NC TB09 and this PAW path at 532.5 nm using controlled
geometries and optical settings. No calculation on a 300+ atom cell should
be launched on the local computer. PAW optical matrix elements and the
separate epsilon normalization correction still require validation.

Local provenance and the fixed-geometry energetic analysis are retained in
the ignored `audit/hydrogen-paw-pristine/` directory, including
`selected-paw.json`, `analyze_local_pseudos.py` and
`local-pseudo-energetics.json`. These local research files, pseudopotentials
and build artifacts are not included in the Git commit.
