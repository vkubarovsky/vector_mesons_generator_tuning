#!/usr/bin/env python3
"""
jpsi_pipeline.py — Generate J/psi MC and run NN fastMC.

  Step 1: Generate MC with diffrad_jpsi_dipole per beam energy
  Step 2: FastMC — per-event acceptance weight = mean accept over the
          fastmc configs of that beam energy (e' in FT, e+/e-/p in FD)
  Output: npz cache per energy + generated-vs-accepted plots

Usage:  python3 jpsi_pipeline.py [config.json]
"""

import json
import subprocess
import sys
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path.home() / "fastmc/scripts_jpsi"))
from fast_mc import load_model_auto

BASE      = Path(__file__).resolve().parent
GENERATOR = BASE / "diffrad_jpsi_dipole.exe"
DATA_ROOT = Path.home() / "Downloads/volatile/clas12/vpk/fastmc"
LUND_ROOT = DATA_ROOT / "jpsi_mc"           # generated LUND files
CACHE_ROOT = DATA_ROOT / "jpsi_tuning_cache" # npz caches

MP  = 0.938272; MP2 = MP**2
MJP = 3.0969;   MJP2 = MJP**2
WTH2 = (MP + MJP)**2
FT_VZ = -3.0   # FT measures no vertex; models trained with vz_rec ~ -3 cm


def compute_tmin(q2, xb):
    """Kinematic |t|_min from Q2, xB."""
    w2 = MP2 + q2 * (1 - xb) / xb
    w = np.sqrt(np.maximum(w2, 0.01))
    ecm_i = (w2 + q2 + MP2) / (2 * w)
    pcm_i = np.sqrt(np.maximum(ecm_i**2 - MP2, 0))
    ecm_f = (w2 + MJP2 - MP2) / (2 * w)
    pcm_f = np.sqrt(np.maximum(ecm_f**2 - MJP2, 0))
    return (ecm_i - ecm_f)**2 - (pcm_i - pcm_f)**2


def step1_generate(cfg, ebeam_str):
    """Run diffrad_jpsi_dipole for one beam energy."""
    gen = cfg["generation"]
    mp = {k: v["value"] for k, v in cfg["cross_section_model"]["parameters"].items()}
    nev = gen["energies"][ebeam_str]["nev"]
    ebeam = float(ebeam_str)

    lund_dir = LUND_ROOT / f"e{ebeam_str}"
    lund_dir.mkdir(parents=True, exist_ok=True)
    lund_file = lund_dir / "generated.lund"
    input_file = lund_dir / "gen_input.txt"

    with open(input_file, "w") as f:
        f.write(f"! J/psi generation, Ebeam={ebeam}\n")
        f.write(f"bmom    {ebeam}\n")
        f.write(f"tmom    0.0\n")
        f.write(f"lepton  1\n")
        f.write(f"ivec    4\n")
        f.write(f"cutv    0.0\n")
        f.write(f"nev     {nev}\n")
        f.write(f"iy      {gen.get('seed', 778899)}\n")
        f.write(f"Q2      {gen['gen_q2min']}  {gen['gen_q2max']}\n")
        f.write(f"y       {gen['gen_ymin']}  {gen['gen_ymax']}\n")
        f.write(f"t       {gen['gen_tmin']}  {gen['gen_tmax']}\n")
        f.write(f"W       {gen['gen_wmin']}  {gen['gen_wmax']}\n")
        f.write(f"xB      {gen['gen_xbmin']}  {gen['gen_xbmax']}\n")
        f.write(f"tslope  {gen['tslope']}\n")
        f.write(f"iborn   1\n")
        f.write(f"alf1    {mp['alf1']}\n")
        f.write(f"alf2    {mp['alf2']}\n")
        f.write(f"alf3    {mp['alf3']}\n")
        f.write(f"nuT     {mp['nuT']}\n")
        f.write(f"mg2     {mp['mg2']}\n")
        f.write(f"cR      {mp['cR']}\n")
        f.write(f"theta_electron  {gen['gen_theta_e_min']}  {gen['gen_theta_e_max']}\n")

    print(f"  Generating {nev} events at {ebeam} GeV "
          f"(alf2={mp['alf2']}, alf3={mp['alf3']}, nuT={mp['nuT']}, "
          f"mg2={mp['mg2']}, cR={mp['cR']})...", flush=True)
    t0 = time.time()
    result = subprocess.run(
        [str(GENERATOR), "-input", str(input_file), "-lund", str(lund_file)],
        cwd=str(BASE), capture_output=True, text=True, timeout=3600)
    print(f"  Generator done ({time.time()-t0:.0f}s), exit={result.returncode}")
    if result.returncode != 0:
        print(result.stdout[-800:]); sys.exit(1)
    return lund_file


def read_lund_jpsi(fpath, ebeam):
    """Read 6-particle ivec=4 LUND: beam e-, e'(FT), p, J/psi, e+, e-(decay)."""
    n_re = 0
    cols = {k: [] for k in
            ["q2", "xb", "t", "w", "t_tmin",
             "p_eft", "th_eft", "ph_eft",
             "p_p", "th_p", "ph_p",
             "p_ep", "th_ep", "ph_ep",
             "p_ed", "th_ed", "ph_ed", "mee"]}

    def sph(px, py, pz):
        p = np.sqrt(px*px + py*py + pz*pz)
        th = np.degrees(np.arctan2(np.hypot(px, py), pz))
        ph = np.degrees(np.arctan2(py, px))
        return p, th, ph

    with open(fpath) as f:
        while True:
            hdr = f.readline()
            if not hdr:
                break
            c = hdr.split()
            if not c:
                continue
            npart = int(c[0])
            parts = []
            for _ in range(npart):
                pc = f.readline().split()
                parts.append((int(pc[3]),
                              float(pc[6]), float(pc[7]), float(pc[8]), float(pc[9])))
            if npart < 6:
                continue
            # layout: 0=beam, 1=e'(FT), 2=p, 3=jpsi, 4=e+, 5=e-(decay)
            _, ex, ey, ez, eE = parts[1]
            _, px_, py_, pz_, pE = parts[2]
            _, qx, qy, qz, qE = parts[4]
            _, dx, dy, dz, dE_ = parts[5]

            pe, the, phe = sph(ex, ey, ez)
            cos_th = ez / pe if pe > 0 else 1.0
            q2 = 2.0 * ebeam * eE * (1.0 - cos_th)
            nu = ebeam - eE
            if nu <= 0:
                continue
            xb = q2 / (2.0 * MP * nu)
            w2 = MP2 + 2.0 * MP * nu - q2
            if w2 <= 0:
                continue
            w = np.sqrt(w2)
            dEp = pE - MP
            t = abs(dEp**2 - px_**2 - py_**2 - pz_**2)
            mee2 = (qE + dE_)**2 - (qx+dx)**2 - (qy+dy)**2 - (qz+dz)**2
            mee = np.sqrt(max(mee2, 0))

            pp, thp, php = sph(px_, py_, pz_)
            pq, thq, phq = sph(qx, qy, qz)
            pd, thd, phd = sph(dx, dy, dz)

            cols["q2"].append(q2); cols["xb"].append(xb)
            cols["t"].append(t); cols["w"].append(w)
            cols["t_tmin"].append(0.0)  # filled vectorized below
            cols["p_eft"].append(pe); cols["th_eft"].append(the); cols["ph_eft"].append(phe)
            cols["p_p"].append(pp); cols["th_p"].append(thp); cols["ph_p"].append(php)
            cols["p_ep"].append(pq); cols["th_ep"].append(thq); cols["ph_ep"].append(phq)
            cols["p_ed"].append(pd); cols["th_ed"].append(thd); cols["ph_ed"].append(phd)
            cols["mee"].append(mee)
            n_re += 1

    mc = {k: np.asarray(v) for k, v in cols.items()}
    mc["t_tmin"] = mc["t"] - compute_tmin(mc["q2"], mc["xb"])
    print(f"  Parsed {n_re} MC events")
    return mc


def step2_fastmc(mc, fastmc_configs):
    """Per-event acceptance weight = mean accept over configs.

    e'(FT) uses the e- model with vz fixed at FT_VZ; e+/e-d/p use FD models
    with vz = FT_VZ (generator has no vertex spread; -3 cm = target center).
    """
    n = len(mc["q2"])
    vz = np.full(n, FT_VZ)
    weight = np.zeros(n)

    for config in fastmc_configs:
        model_dir = str(DATA_ROOT / "jpsi" / config / "v1" / "models")
        acc = np.ones(n, dtype=bool)
        for mname, pk, thk, phk in [("e-", "p_eft", "th_eft", "ph_eft"),
                                    ("p",  "p_p",  "th_p",  "ph_p"),
                                    ("e+", "p_ep", "th_ep", "ph_ep"),
                                    ("e-d","p_ed", "th_ed", "ph_ed")]:
            m = load_model_auto(model_dir, mname)
            if m is None:
                print(f"  ERROR: no model {mname} in {model_dir}"); sys.exit(1)
            acc &= m.accept(mc[pk], mc[thk], mc[phk], vz)
        weight += acc
        print(f"    {config}: acc(all4) = {100*acc.mean():.2f}%")

    weight /= len(fastmc_configs)
    print(f"  Mean acceptance weight: {weight.mean():.4f}")
    return weight


def step_plots(mc_all, out_png):
    """Generated vs accepted for both energies combined (acceptance-weighted)."""
    fig, axes = plt.subplots(2, 5, figsize=(25, 8))
    plot_vars = [("q2", r"$Q^2$ [GeV$^2$]", 0.0, 0.12),
                 ("w",  r"$W$ [GeV]", 4.0, 4.6),
                 ("t",  r"$|t|$ [GeV$^2$]", 0, 6),
                 ("t_tmin", r"$|t|-|t|_{min}$ [GeV$^2$]", 0, 4),
                 ("p_eft", r"$p_{e'}$ (FT) [GeV]", 0, 3)]
    for col, (var, label, vmin, vmax) in enumerate(plot_vars):
        for row, scale in enumerate(["linear", "log"]):
            ax = axes[row, col]
            for eb, mc, wgt in mc_all:
                ax.hist(mc[var], bins=50, range=(vmin, vmax),
                        histtype="step", lw=1.2, label=f"gen {eb}")
                ax.hist(mc[var], bins=50, range=(vmin, vmax), weights=wgt,
                        histtype="step", lw=2, label=f"acc {eb}")
            ax.set_xlabel(label); ax.grid(True, alpha=0.3)
            if scale == "log":
                ax.set_yscale("log")
            if col == 0 and row == 0:
                ax.legend(fontsize=8)
    fig.suptitle("J/psi MC: generated vs fastMC-accepted (weighted)", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    print(f"  Saved: {out_png}")


def main():
    cfg_path = BASE / "config.json" if len(sys.argv) < 2 else Path(sys.argv[1])
    with open(cfg_path) as f:
        cfg = json.load(f)

    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    mc_all = []
    for ebeam_str, espec in cfg["generation"]["energies"].items():
        print(f"\n{'='*60}\n  Ebeam = {ebeam_str} GeV\n{'='*60}")
        cache = CACHE_ROOT / f"mc_e{ebeam_str}.npz"
        if cache.exists():
            print(f"  Loading cache {cache.name}")
            z = np.load(cache)
            mc = {k: z[k] for k in z.files if k != "acc_weight"}
            wgt = z["acc_weight"]
        else:
            lund = step1_generate(cfg, ebeam_str)
            mc = read_lund_jpsi(lund, float(ebeam_str))
            wgt = step2_fastmc(mc, espec["fastmc_configs"])
            np.savez_compressed(cache, acc_weight=wgt, **mc)
            print(f"  Cached -> {cache} ({cache.stat().st_size/1e6:.0f} MB)")
        mc_all.append((ebeam_str, mc, wgt))

    step_plots(mc_all, BASE / "step_gen_vs_acc.png")
    print("\nPipeline done.")


if __name__ == "__main__":
    main()
