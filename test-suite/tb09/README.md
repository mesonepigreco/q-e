# TB09 / Tran–Blaha mBJLDA regression and silicon benchmark

Upstream accepts contributions through GitLab merge requests; GitHub is a
mirror (see [the submission template](../../.github/pull_request_template.md)).
No upstream submission is part of this benchmark.

This change repairs four discrepancies in the CPU PW implementation at
upstream commit `61569eb480231b649c47f631df9c5c83537df461` (QE 7.6):

* Compute the density-dependent cell coefficient `c` at each potential update.
* Pass the spin-resolved density Laplacian to Libxc's Becke–Roussel term.
* Combine TB09 exchange with Perdew–Wang LDA correlation, rather than TPSS.
* Preserve the local model potential at zero density gradient / kinetic-energy
  threshold points; retain the density cutoff.

The cell average uses the spin-summed density and its gradient, counts only
physical FFT points, and reduces over the FFT communicator within each pool.
A density at or below 1e-12 bohr^-3 contributes zero to this average. These
regularization and pseudodensity choices are explicit implementation limits.
The code prints `c` and the average. XClib requires the `lapl_` input for TB09;
standalone callers must also supply `c` through `set_libxc_ext_param`.

TB09 has no parent exchange-energy functional. PW rejects forces, stress,
structural relaxation and GPU execution for this implementation. Existing
norm-conserving / collinear meta-GGA restrictions remain. No meaningful TB09
total energy or structural optimization is claimed.

## Equation tests

`XClib/test_tb09.f90` checks signed Laplacians (both physical branches of the
BR equation), zero-gradient points, equal and unequal spin channels, the
TB09 alias, LDA-only correlation and rejection of a missing Laplacian.
`reference_equations.py` reproduces the fixed local reference potentials from
the BR/TB09 equations using bisection, with only Python's standard library.
The tolerance in the compiled test is 1e-9 Hartree. These are local-potential
fixtures, independent of any fitted material gap.

```
python3 test-suite/tb09/reference_equations.py
cmake --build build --target qe_xclib_tb09_test qe_xclib_test
ctest --test-dir build -R test_qe_xclib --output-on-failure
```

Configure with `QE_ENABLE_TEST=ON` and `QE_ENABLE_LIBXC=ON`. The legacy
all-terms/all-short fixtures explicitly use zero Laplacian and Libxc's local
`c=1`; their TB09 blocks test exchange only. The separate compiled test
checks the alias and LDA correlation. The legacy driver also avoids indexing
QE functional metadata with the translated Libxc ID (208).

## Reproduce the silicon benchmark

Use the unmodified norm-conserving PBE pseudopotential
[Si.pbe-rrkj.UPF](https://pseudopotentials.quantum-espresso.org/upf_files/Si.pbe-rrkj.UPF),
SHA256 `d42702eef81fd417a503bc33e1756081d012b398c46ad9c2955aad963a1760c7`.
It has no nonlinear core correction. Supply its local path; the driver
checks its checksum and does not download anything.

```
python3 test-suite/tb09/verify.py --pw build/bin/pw.x \
  --pseudo /path/to/Si.pbe-rrkj.UPF --workdir /tmp/si-tb09-40
python3 test-suite/tb09/verify.py --pw build/bin/pw.x \
  --pseudo /path/to/Si.pbe-rrkj.UPF --workdir /tmp/si-tb09-60 \
  --cutoff 60 --mesh 12
```

Each work directory must be new unless `--resume` is given. Resuming checks
the binary checksums, configuration and exact input before reusing a completed
job. The driver uses diamond Si, a=5.431 Angstrom,
2 atoms, 8 bands, fixed occupations, ecutrho=4*ecutwfc, no SOC, a Gamma-centered
mesh and 81 points along Gamma–X. It extracts `min(band 5)-max(band 4)` in eV
from the XML, with the VBM at Gamma and CBM on Delta. It saves inputs, outputs,
SCF and band XML, and a JSON summary. The driver uses conv_thr=1e-12 Ry and
`diago_full_acc` to make comparisons of empty states meaningful.
Use `--functional PBE` for the control.

For integration checks with a separate MPI build:

```
python3 test-suite/tb09/verify.py --pw build-serial/bin/pw.x \
  --mpi-pw build-mpi/bin/pw.x --launcher 'mpiexec -n 2' \
  --pseudo /path/to/Si.pbe-rrkj.UPF --workdir /tmp/si-tb09-parallel --checks
```

This compares SCF eigenvalues and `c` with one and two k-point pools, repeats
with two unpolarized spin channels, checks Gamma-only serial/MPI FFTs, and
requires an explicit failure on force requests. It tests numerical consistency,
not agreement with a different density representation. MPI comparisons use
1e-7 Hartree / 1e-7 in `c`; the spin comparison uses 1e-6 to allow independent
SCF stopping errors.

## Validation performed

With the versions below, all six XClib CTests pass in serial and MPI builds,
and XClib also builds without Libxc. The expanded equation fixtures pass in
both builds. In the integration run, the two MPI pool layouts differ from
serial by at most 4.82e-8 Hartree in eigenvalues and 5e-8 in `c`; the spin
comparison differs by 1.25e-7 Hartree. Gamma-only serial/MPI eigenvalues agree
to 8.03e-8 Hartree. Force requests fail with the intended diagnostic. The
resume path also preserves the separate SCF and band XML records.
Machine-readable summaries are in `validation.json`.

## Measured gaps and remaining discrepancy

GCC/gfortran 11.4, Libxc 7.0.0 (`7bd5bb41415968db94c499a2f093309c9a2dcf53`),
internal FFT, system BLAS/LAPACK. The table below uses conv_thr=1e-10 Ry;
PBE is unchanged by the patch. These two refinements are a diagnostic, not a
proof of an asymptotic convergence bound. Numerical records are in
`silicon-results.json`.

| Potential | 40 Ry / 8^3 | 60 Ry / 12^3 | Borlido/Botti 2019 |
| --- | ---: | ---: | ---: |
| PBE, stock and corrected | 0.598623 | 0.599249 | 0.6187 |
| TB09, stock QE | 2.345374 | 2.356616 | 1.1994 |
| TB09, corrected QE | 0.924960 | 0.923810 | 1.1994 |

The original 2009 all-electron paper reports 1.17 eV for Si. The requested
Borlido/Botti PAW reference is the comparison target here. The corrected
norm-conserving result remains 0.275590 eV below it at the finer settings.
That discrepancy is unresolved. A pseudo-valence density, even with a partial
core correction, does not reconstruct the all-electron density used to
calibrate mBJ. The different density representations are a plausible source
of residual error, but this calculation does not establish its full cause.
The patch corrects the demonstrated equation/interface defects; it does not
claim PAW reconstruction or reproduction of Botti's complete benchmark.
Do not fit `c` to Si to hide this difference.

The original 2019 supplementary workbook was also reanalyzed: 472 materials,
471 numeric mBJ gaps (NaYbP2S6 missing), mBJ MAE=0.496558 eV and
MAPE=29.5937%, reproducing the rounded published statistics. These are
statistics from published data, not new calculations for 472 materials.
Workbook SHA256:
`21fc5c597411093c98e6209f3a94e907f6e81b6e27ee87611bcd73bd59ac2db0`.
The 2020 follow-up uses a modified 473-material set and is not mixed into
these statistics. Recompute them with the optional `openpyxl` dependency:

```
python3 test-suite/tb09/botti_statistics.py /path/to/ct9b00322_si_002.xlsx
```

## Sources and historical scope

* Tran and Blaha, PRL 102, 226401 (2009):
  https://doi.org/10.1103/PhysRevLett.102.226401
* Becke and Roussel, PRA 39, 3761 (1989):
  https://doi.org/10.1103/PhysRevA.39.3761
* Borlido, Aull, Huran, Tran, Marques and Botti, JCTC 15, 5069–5079 (2019):
  https://doi.org/10.1021/acs.jctc.9b00322 ;
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6739738/
* Original benchmark workbook: https://doi.org/10.1021/acs.jctc.9b00322.s002
* Borlido et al., npj Computational Materials 6, 96 (2020):
  https://doi.org/10.1038/s41524-020-00360-0
* The requested Germaneau CPC (2013) implementation paper:
  https://www.sciencedirect.com/science/article/abs/pii/S0010465513000775
  Only the accessible preview was available; its complete numerical table
  and original submitted code were not reproduced.

The historical https://gitlab.com/mesonepigreco/q-e fork was inspected at
`8254f37a91a64fcda0b3fb76a09a2106b3afba8c` (2019). In its
`Modules/metagga.f90`, the Laplacian input is assigned the squared gradient,
`c` is not updated, TPSS correlation is used, and `sx=sx*rho` reads an
unassigned exchange energy after a potential-only call. The modern wrapper
already initializes that energy to zero. Those historical findings are not
an attribution of authorship, and this patch targets current QE rather than
backporting to the old Libxc interface.
