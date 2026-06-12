# diffrad_vm — combined vector meson generator

> **FROZEN (2026-06-12).** This copy is the record of the tuning result.
> The supported production version lives in
> [MC_vector_mesons](https://github.com/vkubarovsky/MC_vector_mesons) —
> all further development happens there.

One executable for both tuned channels (June 2026):

| ivec | Meson | t-form | Tuned defaults |
|---|---|---|---|
| 3 | phi -> K+K- | exponential `bt*exp(bt*t)` | alf2=-1.245, alf3=0.762, nuT=2.344, bt=1.284, cR=1.0 |
| 4 | J/psi -> e+e- | dipole `3(mg2-tmin)^3/(mg2-t)^4` | alf2=4.128, alf3=0.32, nuT=3.0, mg2=3.170, cR=0.4 |

Build:
```bash
gfortran -O2 -o diffrad_vm.exe diffrad_vm.f90
```

Run:
```bash
./diffrad_vm.exe -input gen_input.dat -lund output.lund
```

Model parameters are optional input-file keys: `alf1 alf2 alf3 nuT cR`
(shared), `bt` (phi), `mg2` (J/psi). Any key not present in the input file
falls back to the tuned default of the selected meson (ivec). The standalone
per-meson generators in `../phi/` and `../jpsi/` remain the reference
versions used for the tuning itself.

Includes the June 2026 fixes: dipole normalization sign
(`(mg2+|tmin|)^3`, was negative near threshold) and 64-bit attempt counters
(INTEGER*4 overflow killed generations with nev >= 2M).
