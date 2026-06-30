# phi Roberts W-shape retune (2026-06-30)

## What changed
The phi cross-section W-shape is now **fixed to photoproduction** instead of
fitted from CLAS12:
- `alf2 = 2`    (Donnachie-Landshoff threshold exponent; sigma_T vanishes at
  threshold, no divergence)
- `alf3 = 0.20` (Roberts et al. arXiv:2510.08845, phi photoproduction Pomeron
  exponent epsilon)

Only `bt` (t-slope) and `nuT` (Q^2 exponent) are refit; `alf1`, `cR` unchanged.

## Why
The old free fit gave `alf2 = -1.245` (negative threshold exponent), which is
unphysical (sigma diverges at threshold) and blew up the radiative-correction
integral. The CLAS12 W-range [2.0, 3.5] starts ~40 MeV above threshold and is
strongly correlated with Q^2, so it cannot constrain the W-shape -- the
per-dataset `alf3` scattered over 0.22-1.44. So we take the W-shape from where
it is cleanly measured (real-photon photoproduction) and only fit what CLAS12
determines.

## Result (5 RGA datasets, regenerated MC at the fixed W-shape, 2M ev each)
| | alf2 | alf3 | bt | nuT | chi2 |
|---|---|---|---|---|---|
| old free fit | -1.245 | 0.762 | 1.285 | 2.337 | 11.3 |
| **this** | **2.0 (fixed)** | **0.20 (fixed)** | **1.326 +-0.082** | **2.253 +-0.181** | **13.0** |

chi2 statistically identical (Delta=+1.7); bt/nuT consistent and tighter.
The W-distribution chi2 alone is ~1.6 -- the fixed Roberts W-shape describes
the data well.

## How it was run
- Per dataset: `cross_section_model` alf2=2, alf3=0.20 (these are the
  generation params AND the fixed fit values, so reweighting is consistent);
  `fit_parameters` free=[bt,nuT], fixed=[alf2,alf3,cR].
- Generator `diffrad_gen_exp.exe` recompiled from `diffrad_vpk_exp.f90`.
- `phi_pipeline.py <config>` regenerates MC (alf2=2) + fastMC, then
  `acceptance_correct.py <ds>` does the bt/nuT fit. nev=2M (ample for the
  binned chi2; bump to 10M for a final production pass).

## Production
Applied to MC_vector_mesons/diffrad_vm.f90 (branch phi-roberts-wshape):
alf2=2.0, alf3=0.20, nuT=2.253, bt=1.326. RC verified clean
(RC_vector_mesons_DIFFRAD/tuned_roberts_phi_run: eta=0.745-0.840, 0/60
blow-ups vs 19/60 before).

Original alf2=-1.245 tuning preserved on `main` and in
phi/retune_alf2pos_backup/.
