#!/usr/bin/env python3
"""Iterative phi retune, W-shape FIXED to Roberts Eq.23b P-dyn
(alf2=2, alf3=0.20). Free: bt, nuT. Iterate regenerate->refit until
bt,nuT converge. alf1 (normalization) stays free/uninvolved (shape fit)."""
import json, subprocess, sys, os, numpy as np
from pathlib import Path

BASE = Path("/Users/vpk/vector_mesons_generator_tuning/phi")
DS = ["rgafall18_inb", "rgafall18_outb", "rgasp18_inb", "rgasp18_outb", "rgasp19_inb"]
ALF2, ALF3 = 2.0, 0.20
NEV = 2000000
MAXIT, TOL = 6, 0.005
env = {**os.environ, "KMP_DUPLICATE_LIB_OK": "TRUE", "OMP_NUM_THREADS": "1",
       "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}

def cfgpath(d): return BASE / "runs" / d / "config.json"

for d in DS:
    c = json.load(open(cfgpath(d)))
    p = c["cross_section_model"]["parameters"]
    p["alf2"]["value"] = ALF2
    p["alf3"]["value"] = ALF3
    c["fit_parameters"]["free"] = ["bt", "nuT"]
    c["fit_parameters"]["fixed"] = ["alf2", "alf3", "cR"]
    c["generation"]["nev"] = NEV
    json.dump(c, open(cfgpath(d), "w"), indent=2)
print(f"init: alf2={ALF2}, alf3={ALF3} fixed (Roberts Eq.23b P-dyn); free=[bt,nuT]; nev={NEV}", flush=True)

hist = []
for it in range(1, MAXIT + 1):
    print(f"\n===== ITERATION {it} =====", flush=True)
    bt_new, nu_new, ch, bt_old, nu_old = {}, {}, {}, {}, {}
    for d in DS:
        c = json.load(open(cfgpath(d)))
        pp = c["cross_section_model"]["parameters"]
        bt_old[d], nu_old[d] = pp["bt"]["value"], pp["nuT"]["value"]
        subprocess.run([sys.executable, str(BASE / "phi_pipeline.py"), str(cfgpath(d))],
                       capture_output=True, text=True, env=env, cwd=str(BASE), timeout=1800)
        subprocess.run([sys.executable, str(BASE / "acceptance_correct.py"), d],
                       capture_output=True, text=True, env=env, cwd=str(BASE), timeout=600)
        r = json.load(open(BASE / "runs" / d / "acceptance_results.json"))
        bt_new[d], nu_new[d], ch[d] = r["fit_params"]["bt"], r["fit_params"]["nuT"], r["chi2_total"]
        c = json.load(open(cfgpath(d)))
        c["cross_section_model"]["parameters"]["bt"]["value"] = bt_new[d]
        c["cross_section_model"]["parameters"]["nuT"]["value"] = nu_new[d]
        json.dump(c, open(cfgpath(d), "w"), indent=2)
        print(f"  {d:15} bt {bt_old[d]:.3f}->{bt_new[d]:.3f}  nuT {nu_old[d]:.3f}->{nu_new[d]:.3f}  chi2={ch[d]:.1f}", flush=True)
    dbt = max(abs(bt_new[d]-bt_old[d])/max(bt_old[d],1e-6) for d in DS)
    dnu = max(abs(nu_new[d]-nu_old[d])/max(nu_old[d],1e-6) for d in DS)
    avbt, sbt = np.mean([bt_new[d] for d in DS]), np.std([bt_new[d] for d in DS])
    avnu, snu = np.mean([nu_new[d] for d in DS]), np.std([nu_new[d] for d in DS])
    avch = np.mean([ch[d] for d in DS])
    hist.append((it, avbt, sbt, avnu, snu, avch, dbt, dnu))
    print(f"  -> AVG bt={avbt:.3f}+-{sbt:.3f} nuT={avnu:.3f}+-{snu:.3f} chi2={avch:.1f} | maxdelta bt={dbt:.3%} nuT={dnu:.3%}", flush=True)
    if dbt < TOL and dnu < TOL:
        print(f"\nCONVERGED after {it} iterations (delta < {TOL:.1%})", flush=True); break

print("\n===== SUMMARY =====")
print(f"{'it':>2} {'bt':>15} {'nuT':>15} {'chi2':>7} {'d_bt':>7} {'d_nuT':>7}")
for it,avbt,sbt,avnu,snu,avch,dbt,dnu in hist:
    print(f"{it:2d} {avbt:6.3f}+-{sbt:.3f} {avnu:6.3f}+-{snu:.3f} {avch:7.1f} {dbt:7.2%} {dnu:7.2%}")
a=hist[-1]
print(f"\nFINAL (alf2={ALF2}, alf3={ALF3} fixed): bt={a[1]:.3f}+-{a[2]:.3f}, nuT={a[3]:.3f}+-{a[4]:.3f}, chi2={a[5]:.1f}")
