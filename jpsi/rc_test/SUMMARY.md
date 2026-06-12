# J/psi RC Test A — MC self-consistency (Born vs RC)

Date: 2026-06-11 (evening, run via SSH after the interactive session stalled
on a permission prompt at 18:18)

Generator: `diffrad_jpsi_dipole.exe` (compiled 18:10, includes the t-sign fix
and Bug #1 fix v = 2*Mp*omega). Tuned dipole model: alf1=400, alf2=4.105,
alf3=0.32, nuT=3.0, mg2=3.106, cR=0.4.

## Setup

Two runs, identical inputs except `iborn`, same seed (iy=778899), Ebeam=10.6,
FT kinematics from config.json (Q2 0.01–0.12, W 4.05–4.55, y 0.5–0.999,
theta_e 0–7.5 deg), 20000 events each.

- `input_born.txt` (iborn=1) -> `smoke_born.lund` — 0.8 s
- `input_rc.txt`   (iborn=0) -> `smoke_rc.lund`   — 37 min

## Results

| quantity        | Born          | RC            |
|-----------------|---------------|---------------|
| sigma [nb]      | 268.1 ± 1.9   | 248.6 ± 1.8   |
| guard failures  | 0 (all 8 counters) | 0 (all 8 counters) |
| hard ISR/FSR    | —             | 0 / 0 (hard_frac = 0) |

**eta_total = sigma_RC / sigma_Born = 0.9273 ± 0.0093** (−7.3% correction)

eta(W), 10 bins over 4.05–4.55:

| W bin | eta(W) |
|-------|--------|
| 4.05–4.10 | 0.62 ± 0.56 |
| 4.10–4.15 | 0.67 ± 0.22 |
| 4.15–4.20 | 0.80 ± 0.11 |
| 4.20–4.25 | 0.92 ± 0.07 |
| 4.25–4.30 | 0.86 ± 0.04 |
| 4.30–4.35 | 0.88 ± 0.03 |
| 4.35–4.40 | 0.91 ± 0.02 |
| 4.40–4.45 | 0.94 ± 0.02 |
| 4.45–4.50 | 0.95 ± 0.02 |
| 4.50–4.55 | 0.93 ± 0.02 |

Plot: `eta_w_20k.png`

## Assessment

- eta(W) behaves physically: ~0.95 at W = 4.5 falling toward the threshold
  (smaller eta near W_th = 4.035) — radiating a photon pushes W_true below
  threshold, so near-threshold events are suppressed the most. The low-W
  bins are statistics-limited (3–22 Born events) but the trend is monotone
  and smooth where statistics exist.
- All guard/failure counters are zero in both modes (only the normal
  accept/reject and t-kinematics rejections fire).
- Zero hard-photon events in RC mode: expected with the tuned steep
  threshold (alf2=4.105) — sigma_F (Bardin–Shumeiko finite hard remainder)
  is negative and clamped to 0. The RC to the rate is carried entirely by
  the soft factor. Known small bias O(0.1%) from the clamp; the v_cut split
  fix discussed in the main session would make sigma_hard positive-definite.
- Same-seed pairing verified: first events of both runs land on identical
  kinematics.

**Test A: PASS** (internal consistency; no pathologies).

## 100k parallel run (run_parallel.sh, 10 jobs x 10k, distinct seeds)

Output in `lund100k/`. sigma_Born = 266.44 ± 0.84 nb,
sigma_RC = 247.48 ± 0.78 nb, **eta_total = 0.9288 ± 0.0042** —
confirms the 20k result. nhard = 2/100000 (sigma_F marginally positive
at isolated points). eta(W) now cleanly monotone: 0.34 ± 0.20 in the
threshold bin (4.05–4.10) rising to 0.96 at W = 4.55 (`eta_w_100k.png`).
eta(Q2) from plot_eta.py: flat, 0.935 -> 0.90 over Q2 = 0.01–0.12
(`lund100k/eta_mc_vs_exact.png`). plot_compare.py validation plots in
`lund100k/` (t-slope fit b = 0.727 ± 0.011 GeV^-2).

## v_cut split quick test (`vcut_test/`)

Implemented the naive v_cut split in `diffrad_jpsi_vcut.f90`:
delinf integrates soft part to vcut = 1e-4 only; Born-subtraction
removed from podinl so qqt integrates the exact bremsstrahlung over
[vcut, vmax]. 2000-event RC run:

- Hard-photon events restored: 171/2000 (8.6%), real gamma (pid 22)
  written to LUND; v-spectrum ~1/v piling up above the cut as expected
  (mean v = 6.5e-4 GeV^2).
- **FAILS criterion (a):** sigma_total = 180.6 ± 4.0 nb vs 248.6 ± 1.8
  baseline (−27%, ~17 sigma).
- All 171 tagged FSR, 0 ISR — also suspicious.

Diagnosis: (1) without the subtraction the exact integrand has
m_e^2-width collinear peaks in ta that simptx / the 40-point log-v grid
do not resolve -> most of the leading log(Q2/m_e2) ~ 12 is lost;
(2) exponentiating the soft factor at vcut instead of vmax changes
higher orders by O((alpha*L)^2) ~ 10%.

Proposed correct scheme: keep the smooth subtracted qqt = sigma_F and
add the analytic soft-log piece:
sigma_hard = sigma_F + sigma_Born*(alpha/pi)*(log(Q2/m_e2)-1)*log(vmax^2/vcut^2)
(positive by construction, sigma_total unchanged from the validated
Test A rate). Sample v from the approximate kernel weighted by the
Born suppression at shifted W'. Verify against Test B exact integration.

## Test B: exact idiffrad vs MC (`exact/`)

Ported the tuned dipole sigma_T/L (hardwired params) into the original
Akushevich `idiffrad.f` -> `exact/idiffrad_jpsi.f`. Ran exact integration
at the 10 MC W-bin centers, Q2 = 0.05, |t| = 1.5, Ebeam = 10.6, no v cut
(4.9 s). eta_exact rises 0.708 (W=4.075) -> 0.931 (W=4.525).

Comparison (`testb_mc_vs_exact.png`): MC eta(W) tracks the exact curve in
shape and magnitude. MC sits ~0.02-0.03 above exact at high W — expected,
since exact is at fixed (Q2,t) while MC is integrated over Q2 [0.01,0.12],
t, y (eta falls with Q2: 0.935 -> 0.90 across the FT range). Threshold
bin consistent within MC statistics. **Test B: MC RC factor confirmed by
exact integration at the few-% level.** Rigorous version: recompute MC
eta in a narrow (Q2, t) slice around the exact points.

## phi overnight run (5 jobs x 10k, tuned phi source)

Output: `~/Downloads/DIFFRAD_lund/vpk_tuned/phi/`, log
`MC_vector_mesons/phi/run_parallel_overnight.log`. Completed in ~25 min.
sigma_Born = 592920 ± 2652 nb, sigma_RC = 581080 ± 2599 nb,
**eta_phi = 0.980**. **nhard = 563 (1.13%) — hard photons present for
phi** (sigma_F > 0 away from threshold), confirming the J/psi zero-hard
result is channel-specific threshold physics, not a code defect.
Note: RC template uses cutv = 1.2 GeV^2, which caps the radiated-photon
phase space; the hard fraction is defined within that cut.
UPDATE (morning): with the Born template's phase space matched to the RC
template (y 0.2-0.9, xB 0.02-0.6), eta_phi = 0.898 ± 0.006 — the
historical 0.90. The 0.980 was a template phase-space mismatch artifact
(Born y >= 0.25 vs RC y >= 0.2). Recommend harmonizing gen_input_born.txt.

## sig_hard fix — IMPLEMENTED AND VERIFIED (`sighard_fix/`, 2026-06-12)

`diffrad_jpsi_sighard.f90` (copy; production source untouched):
- sigma_total keeps the validated B-S rate with SIGNED sigma_F (clamp
  bias removed).
- Hard-event sampling cross section = sigma_F +
  sib*(alpha/pi)(dlm-1)ln(vmax^2/vcut^2), vcut = 1e-2 GeV^2 (omega ~
  5 MeV). Positive by construction; no numerical integration of the
  collinear peaks (the failure mode of the naive v_cut prototype).
- New `sample_v_soft`: v sampled from (1/v) x Born-suppression
  sigma(W^2-v)/sigma(W^2); replaces podinl-based sampling.

Verification (10x2k parallel):
- sigma_Born = 267.0 ± 1.9 (baseline 268.1/266.4) OK
- **sigma_RC = 247.4 ± 1.7 vs baseline 248.6 ± 1.8 / 247.5 ± 0.8 —
  criterion (a) PASS** (the naive prototype gave 180.6)
- eta = 0.9267 (baseline 0.9273/0.9288); eta(W) consistent within stats
- **nhard = 3575/20000 (17.9%) real photons in LUND** — bookkeeping
  fraction at vcut=1e-2 (mostly quasi-soft, median v = 0.054); spectrum
  = flat-in-log 1/v kernel cut off by threshold suppression at v ~ 1
  (`vspectrum.png`)
- ISR/FSR = 8/3567 — follows the existing p_ISR = E2'^2/(E1^2+E2'^2)
  model at high y; same formula as production code (model choice worth
  reviewing separately)
- counters clean (W-threshold fail = 10 per 2k job, post-photon guard)

Remaining before production: port the patch to the phi/combined
generator + phi stability check (criterion c); validate the sampled
v-spectrum against idiffrad's exact v-distribution mode (cutv < 0).

`plot_eta2.py` (new, in OneDrive jpsi/ and rc_test/): auto-binned
eta(W) and eta(Q2) for any channel; original plot_eta.py untouched.
