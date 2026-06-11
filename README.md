# Vector Meson Generator Tuning

Tuning of the DIFFRAD Monte Carlo generator cross-section parameters for
exclusive vector meson electroproduction at CLAS12, using neural-network
fast simulation ([clas12-fastmc](https://github.com/vkubarovsky/clas12-fastmc))
for detector acceptance.

Each meson is a **fully independent package** — one particle, one generator,
one tuning chain. Code is deliberately duplicated rather than shared, so a
converged analysis can never be broken by work on another channel.

```
phi/      ep -> e' p phi,   phi  -> K+ K-     (FD electron)   — COMPLETE
jpsi/     ep -> e' p J/psi, J/psi-> e+ e-     (FT electron)   — in progress
rho/      ep -> e' p rho,   rho  -> pi+ pi-                   — planned
```

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
Spring19 inb); errors are the dataset-to-dataset spread. These values are
compiled into `phi/diffrad_vpk_exp.f90` as defaults and can be overridden
via the generator input file (keys: alf1 alf2 alf3 nuT bt cR).

## External dependencies (not in this repo)

| What | Where | Source |
|---|---|---|
| fastMC code | `~/fastmc/scripts_phi`, `~/fastmc/scripts_jpsi` | github.com/vkubarovsky/clas12-fastmc |
| NN models (.pt) | `~/Downloads/volatile/clas12/vpk/fastmc/<meson>/...` | mirror of ifarm `/volatile/clas12/vpk/fastmc/` |
| phi data skims | `~/Downloads/bhawani_phi_data/*.lund` | Bhawani |
| jpsi data | `Mariana_All_top2_jpsi.lund` (68 events) | Mariana |
| Generated MC | `~/Downloads/volatile/clas12/vpk/fastmc/<meson>/<config>/...` | produced by the chain |

## Running the phi chain

```bash
cd phi
gfortran -O2 -o diffrad_gen_exp.exe diffrad_vpk_exp.f90
export KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 run_full_chain.py        # generate -> fastMC -> M(KK) fits -> acceptance fit -> LaTeX report
```

Outputs: `runs/<dataset>/` (fit results JSON, plots), `tex/phi_tuning_report.pdf`.

The chain self-documents one tuning iteration; to iterate, put the fitted
averages into the configs (`cross_section_model.parameters`) and rerun.
