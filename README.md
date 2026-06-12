# Vector Meson Generator Tuning

Tuning of the DIFFRAD Monte Carlo generator cross-section parameters for
exclusive vector meson electroproduction at CLAS12, using neural-network
fast simulation ([clas12-fastmc](https://github.com/vkubarovsky/clas12-fastmc))
for detector acceptance.

> **One-program policy (since 2026-06-12).** The supported generator is
> `diffrad_vm.f90` in [MC_vector_mesons](https://github.com/vkubarovsky/MC_vector_mesons)
> — one repository, one Fortran source, all bug fixes included
> (tmin_k energy pairing, dipole normalization sign, 64-bit counters,
> sig_hard_fix RC scheme). **Any new or repeated tuning starts from that
> file.**
>
> The generators that produced the June 2026 results stay here as frozen
> reference — do not edit them, do not start new work from them:
> `phi/diffrad_vpk_exp.f90` and `jpsi/diffrad_jpsi_dipole.f90`.
> The combined version as of the tuning is in git history:
> `git show 16083a9:combined/diffrad_vm.f90`.

## Status (2026-06-12)

```
phi/      ep -> e' p phi,   phi  -> K+ K-     (FD electron)   — COMPLETE
jpsi/     ep -> e' p J/psi, J/psi-> e+ e-     (FT electron)   — tuned (alf2=4.122, mg2=3.112); RC study done (see jpsi/rc_test/SUMMARY.md)
rho/      ep -> e' p rho,   rho  -> pi+ pi-                   — planned
```

The tuned parameter sets for both channels are compiled into
`MC_vector_mesons/diffrad_vm.f90` as defaults, overridable via the
generator input file (keys: `alf1 alf2 alf3 nuT cR`, `bt` for phi,
`mg2` for J/psi).

## phi — final tuned parameters (June 2026)

dsigma/dt = alf1 * (1 - W_th^2/W^2)^alf2 * W^alf3 * (1 + Q^2/m_phi^2)^(-nuT)
            * bt * exp(bt*t) * (1 + cR*Q^2/m_phi^2)

| Parameter | Value | Meaning |
|---|---|---|
| alf2 | -1.245 +- 0.185 | threshold exponent |
| alf3 |  0.762 +- 0.467 | W power law |
| bt   |  1.284 +- 0.114 | t-slope [GeV^-2] |
| nuT  |  2.344 +- 0.119 | Q^2 suppression |
| cR   |  1.0 (fixed)    | sigma_L/sigma_T coefficient |

Average over 5 CLAS12 RGA datasets (Fall18 inb/outb, Spring18 inb/outb,
Spring19 inb); errors are the dataset-to-dataset spread.

## J/psi — tuned parameters (June 2026, dipole t-form)

| Parameter | Value | Meaning |
|---|---|---|
| alf2 | 4.122 | threshold exponent |
| alf3 | 0.32 (fixed) | W power law |
| nuT  | 3.0 (fixed) | Q^2 suppression |
| mg2  | 3.112 | dipole mass^2 [GeV^2] |
| cR   | 0.4 (fixed) | sigma_L/sigma_T coefficient |

Tuned against Mariana's 68-event FT skim; see `jpsi/tex/jpsi_tuning_report.pdf`.
RC self-consistency test: eta_total = 0.9273 +- 0.0093 (-7.3% correction),
details in `jpsi/rc_test/SUMMARY.md`.

## External dependencies (not in this repo)

| What | Where | Source |
|---|---|---|
| fastMC code | `~/fastmc/scripts_phi`, `~/fastmc/scripts_jpsi` | github.com/vkubarovsky/clas12-fastmc |
| NN models (.pt) | `~/Downloads/volatile/clas12/vpk/fastmc/<meson>/...` | mirror of ifarm `/volatile/clas12/vpk/fastmc/` |
| phi data skims | `~/Downloads/bhawani_phi_data/*.lund` | Bhawani |
| jpsi data | `Mariana_All_top2_jpsi.lund` (68 events) | Mariana |
| Generated MC | `~/Downloads/volatile/clas12/vpk/fastmc/<meson>/<config>/...` | produced by the chain |

## Running the tuning chain (e.g. phi)

```bash
# get the current production generator
cp ~/MC_vector_mesons/diffrad_vm.f90 phi/
cd phi
gfortran -O2 -o diffrad_gen_exp.exe diffrad_vm.f90
export KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 run_full_chain.py        # generate -> fastMC -> M(KK) fits -> acceptance fit -> LaTeX report
```

Outputs: `runs/<dataset>/` (fit results JSON, plots), `tex/phi_tuning_report.pdf`.

The chain self-documents one tuning iteration; to iterate, put the fitted
averages into the configs (`cross_section_model.parameters`) and rerun.
