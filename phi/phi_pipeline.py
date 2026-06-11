#!/usr/bin/env python3
"""
phi_pipeline.py — Full automated pipeline for phi electroproduction analysis.

  Step 1: Generate MC events with diffrad (exp(bt) model)
  Step 2: Run NN-based fast MC (acceptance + smearing)
  Step 3: Compare MC+fastMC vs Bhawani data
  Step 4: Compute and plot acceptance
  Step 5: Acceptance-correct data, fit parameters

Usage:
    python3 phi_pipeline.py runs/rgafall18_inb/config.json
    python3 phi_pipeline.py runs/rgafall18_outb/config.json

All parameters come from the config.json file.
All output (plots, results) goes into the same directory as config.json.
"""

import json
import subprocess
import sys
import time
import shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
from pathlib import Path

sys.path.insert(0, str(Path.home() / "fastmc/scripts_phi"))
from fast_mc import load_model_auto

# ═══════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════

LUND_DIR  = Path.home() / "Downloads/bhawani_phi_data"
GEN_DIR   = Path(__file__).resolve().parent
DATA_ROOT = Path.home() / "Downloads/volatile/clas12/vpk/fastmc"
GENERATOR = GEN_DIR / "diffrad_gen_exp.exe"

MP  = 0.938272;  MP2 = MP**2
MPHI = 1.019412; MPHI2 = MPHI**2
WTH2 = (MP + MPHI)**2

MIN_GEN = 50   # minimum generated events per bin for acceptance


# ═══════════════════════════════════════════════════════════════════
# Physics
# ═══════════════════════════════════════════════════════════════════

def compute_tmin(q2, xb):
    w2 = MP2 + q2 * (1 - xb) / xb
    w = np.sqrt(np.maximum(w2, 0.01))
    ecm_i = (w2 + q2 + MP2) / (2 * w)
    pcm_i = np.sqrt(np.maximum(ecm_i**2 - MP2, 0))
    ecm_f = (w2 + MPHI2 - MP2) / (2 * w)
    pcm_f = np.sqrt(np.maximum(ecm_f**2 - MPHI2, 0))
    return (pcm_i - pcm_f)**2 - ((q2 + MPHI2) / (2 * w))**2   # |t|_min > 0


def sigma_model(q2, xb, t_neg, alf2, alf3, bt, nuT, cR):
    """
    dσ/dt ∝ (1 - W²_th/W²)^α₂ · W^α₃ · (1+Q²/m²φ)^(-νT)
            · bt · exp(bt·t) · (1 + cR·Q²/m²φ)
    t < 0, bt > 0
    """
    w2 = MP2 + q2 * (1 - xb) / xb
    w = np.sqrt(np.maximum(w2, 0.01))
    thresh = np.maximum(1.0 - WTH2 / w2, 1e-30)**alf2
    sig_T = thresh * w**alf3 * (1.0 + q2 / MPHI2)**(-nuT) \
            * bt * np.exp(bt * t_neg)
    sig_L = cR * (q2 / MPHI2) * sig_T
    return sig_T + sig_L


def chi2_shape(h_data, h_model):
    nd, nm = h_data.sum(), h_model.sum()
    if nd == 0 or nm == 0:
        return 1e6
    hd, hm = h_data / nd, h_model / nm
    err2 = h_data / nd**2
    mask = err2 > 0
    ndf = max(mask.sum() - 1, 1)
    return np.sum((hd[mask] - hm[mask])**2 / err2[mask]) / ndf


# ═══════════════════════════════════════════════════════════════════
# I/O helpers
# ═══════════════════════════════════════════════════════════════════

def read_lund_6p(fpath, ebeam):
    """Read 6-particle LUND from diffrad generator."""
    q2, xb, t, w = [], [], [], []
    pe, pp, pkp, pkm = [], [], [], []
    theta_e, theta_p, theta_kp, theta_km = [], [], [], []
    phi_e, phi_p, phi_kp, phi_km = [], [], [], []
    vz_e, vz_p, vz_kp, vz_km = [], [], [], []

    with open(fpath) as f:
        while True:
            hdr = f.readline()
            if not hdr:
                break
            cols = hdr.split()
            if not cols:
                continue
            npart = int(cols[0])
            if npart != 6:
                for _ in range(npart):
                    f.readline()
                continue

            parts = []
            for _ in range(6):
                line = f.readline()
                if not line:
                    break
                pc = line.split()
                parts.append({
                    "pid": int(pc[3]),
                    "px": float(pc[6]), "py": float(pc[7]),
                    "pz": float(pc[8]), "E": float(pc[9]),
                    "mass": float(pc[10]),
                    "vz": float(pc[13]) if len(pc) > 13 else 0.0,
                })

            # Layout: 0=beam_e, 1=scat_e, 2=recoil_p, 3=phi(333), 4=K+, 5=K-
            e = parts[1]; pro = parts[2]; kp = parts[4]; km = parts[5]

            pe_val = np.sqrt(e["px"]**2 + e["py"]**2 + e["pz"]**2)
            Ee = e["E"]; nu = ebeam - Ee
            cos_th = e["pz"] / pe_val if pe_val > 0 else 1.0
            q2_val = 2.0 * ebeam * Ee * (1.0 - cos_th)
            xb_val = q2_val / (2.0 * MP * nu) if nu > 0 else 0
            w2_val = MP2 + q2_val * (1 - xb_val) / xb_val if xb_val > 0 else 0
            w_val = np.sqrt(max(w2_val, 0))

            dpx, dpy, dpz = pro["px"], pro["py"], pro["pz"]
            dE = pro["E"] - MP
            t_val = abs(dE**2 - dpx**2 - dpy**2 - dpz**2)

            q2.append(q2_val); xb.append(xb_val)
            t.append(t_val); w.append(w_val)

            for part, p_list, th_list, ph_list, vz_list in [
                (e, pe, theta_e, phi_e, vz_e),
                (pro, pp, theta_p, phi_p, vz_p),
                (kp, pkp, theta_kp, phi_kp, vz_kp),
                (km, pkm, theta_km, phi_km, vz_km),
            ]:
                px, py, pz = part["px"], part["py"], part["pz"]
                pmag = np.sqrt(px**2 + py**2 + pz**2)
                th = np.degrees(np.arctan2(np.sqrt(px**2 + py**2), pz))
                ph = np.degrees(np.arctan2(py, px))
                p_list.append(pmag); th_list.append(th)
                ph_list.append(ph); vz_list.append(part["vz"])

    return {
        "q2": np.array(q2), "xb": np.array(xb), "t": np.array(t), "w": np.array(w),
        "pe": np.array(pe), "pp": np.array(pp),
        "pkp": np.array(pkp), "pkm": np.array(pkm),
        "theta_e": np.array(theta_e), "theta_p": np.array(theta_p),
        "theta_kp": np.array(theta_kp), "theta_km": np.array(theta_km),
        "phi_e": np.array(phi_e), "phi_p": np.array(phi_p),
        "phi_kp": np.array(phi_kp), "phi_km": np.array(phi_km),
        "vz_e": np.array(vz_e), "vz_p": np.array(vz_p),
        "vz_kp": np.array(vz_kp), "vz_km": np.array(vz_km),
    }


def read_data_lund(fpath, ebeam):
    """Read 4-particle Bhawani LUND file (e-, p, K+, K-)."""
    q2, xb, t, w = [], [], [], []
    pe, pp, pkp, pkm = [], [], [], []
    theta_e_list = []

    with open(fpath) as f:
        while True:
            hdr = f.readline()
            if not hdr:
                break
            cols = hdr.split()
            if not cols:
                continue
            npart = int(cols[0])

            particles = {}
            for _ in range(npart):
                line = f.readline()
                if not line:
                    break
                pc = line.split()
                pid = int(pc[3])
                px, py, pz = float(pc[6]), float(pc[7]), float(pc[8])
                E = float(pc[9])
                particles[pid] = {"px": px, "py": py, "pz": pz, "E": E}

            if 11 not in particles or 2212 not in particles:
                continue

            e = particles[11]
            pe_val = np.sqrt(e["px"]**2 + e["py"]**2 + e["pz"]**2)
            theta_e_val = np.degrees(np.arctan2(
                np.sqrt(e["px"]**2 + e["py"]**2), e["pz"])) if pe_val > 0 else 0
            Ee = e["E"]; nu = ebeam - Ee
            cos_th = e["pz"] / pe_val if pe_val > 0 else 1.0
            q2_val = 2.0 * ebeam * Ee * (1.0 - cos_th)
            xb_val = q2_val / (2.0 * MP * nu) if nu > 0 else 0
            w2_val = MP2 + q2_val * (1 - xb_val) / xb_val if xb_val > 0 else 0

            pro = particles[2212]
            dpx, dpy, dpz = pro["px"], pro["py"], pro["pz"]
            dE = pro["E"] - MP
            t_val = abs(dE**2 - dpx**2 - dpy**2 - dpz**2)

            pp_val = np.sqrt(pro["px"]**2 + pro["py"]**2 + pro["pz"]**2)
            pkp_val = np.sqrt(particles[321]["px"]**2 + particles[321]["py"]**2 +
                              particles[321]["pz"]**2) if 321 in particles else 0
            pkm_val = np.sqrt(particles[-321]["px"]**2 + particles[-321]["py"]**2 +
                              particles[-321]["pz"]**2) if -321 in particles else 0

            q2.append(q2_val); xb.append(xb_val)
            t.append(t_val); w.append(np.sqrt(max(w2_val, 0)))
            pe.append(pe_val); pp.append(pp_val)
            pkp.append(pkp_val); pkm.append(pkm_val)
            theta_e_list.append(theta_e_val)

    return {
        "q2": np.array(q2), "xb": np.array(xb), "t": np.array(t), "w": np.array(w),
        "pe": np.array(pe), "pp": np.array(pp),
        "pkp": np.array(pkp), "pkm": np.array(pkm),
        "theta_e": np.array(theta_e_list),
    }


def apply_cuts(d, cuts):
    mask = ((d["q2"] >= cuts["q2min"]) & (d["q2"] <= cuts["q2max"]) &
            (d["w"] >= cuts["wmin"]) & (d["pe"] >= cuts["pemin"]))
    if cuts.get("wmax", 100) < 100:
        mask &= d["w"] <= cuts["wmax"]
    if "pk" in d:
        mask &= d["pk"] <= cuts["pkmax"]
    if cuts.get("pemax", 100) < 100:
        mask &= d["pe"] <= cuts["pemax"]
    if cuts.get("tmin", 0) > 0:
        mask &= d["t"] >= cuts["tmin"]
    if cuts.get("tmax", 100) < 100:
        mask &= d["t"] <= cuts["tmax"]
    if cuts.get("xbmin", 0) > 0:
        mask &= d["xb"] >= cuts["xbmin"]
    if cuts.get("xbmax", 1) < 1:
        mask &= d["xb"] <= cuts["xbmax"]
    if cuts.get("theta_e_min", 0) > 0 and "theta_e" in d:
        mask &= d["theta_e"] >= cuts["theta_e_min"]
    if cuts.get("theta_e_max", 180) < 180 and "theta_e" in d:
        mask &= d["theta_e"] <= cuts["theta_e_max"]
    return mask


# ═══════════════════════════════════════════════════════════════════
# Pipeline steps
# ═══════════════════════════════════════════════════════════════════

def step1_generate(cfg):
    """Generate MC events with diffrad."""
    print(f"\n{'='*70}")
    print(f"  STEP 1: Generate MC events")
    print(f"{'='*70}")

    out_dir = cfg["out_dir"]
    gen = cfg["generation"]
    mp = cfg["model_parameters"]

    input_file = out_dir.resolve() / "gen_input.txt"
    # Write large LUND file to volatile mirror, not the repo
    lund_dir = DATA_ROOT / "phi" / cfg["model"] / "v11" / "lund"
    lund_dir.mkdir(parents=True, exist_ok=True)
    lund_file = lund_dir / "generated.lund"

    with open(input_file, "w") as f:
        f.write(f"! Auto-generated for {cfg['name']}\n")
        f.write(f"bmom    {cfg['ebeam']}\n")
        f.write(f"tmom    0.0\n")
        f.write(f"lepton  1\n")
        f.write(f"ivec    3\n")
        f.write(f"cutv    0.0\n")
        f.write(f"nev     {gen['nev']}\n")
        f.write(f"iy      {gen.get('seed', 778899)}\n")
        f.write(f"Q2      {gen['gen_q2min']}  {gen['gen_q2max']}\n")
        f.write(f"y       {gen.get('gen_ymin', 0.25)}  {gen.get('gen_ymax', 0.9)}\n")
        f.write(f"t       {gen.get('gen_tmin', 0.2)}   {gen.get('gen_tmax', 8.0)}\n")
        f.write(f"W       {gen['gen_wmin']}  {gen.get('gen_wmax', 100.0)}\n")
        f.write(f"xB      {gen.get('gen_xbmin', 0.001)}  {gen.get('gen_xbmax', 0.6)}\n")
        f.write(f"tslope  {mp['bt']}\n")
        f.write(f"iborn   1\n")
        # Phi cross-section model parameters
        f.write(f"alf1    {mp.get('alf1', 400.0)}\n")
        f.write(f"alf2    {mp['alf2']}\n")
        f.write(f"alf3    {mp['alf3']}\n")
        f.write(f"nuT     {mp['nuT']}\n")
        f.write(f"bt      {mp['bt']}\n")
        f.write(f"cR      {mp['cR']}\n")
        if "gen_pemin" in gen:
            f.write(f"momentum_electron  {gen['gen_pemin']}  {gen['gen_pemax']}\n")
        if "gen_theta_e_min" in gen:
            f.write(f"theta_electron     {gen['gen_theta_e_min']}  {gen['gen_theta_e_max']}\n")

    print(f"  Ebeam={cfg['ebeam']}, nev={gen['nev']}, bt={mp['bt']}")
    print(f"  Q2=[{gen['gen_q2min']}, {gen['gen_q2max']}], "
          f"W=[{gen['gen_wmin']}, {gen.get('gen_wmax', 100.0)}]")
    if "gen_pemin" in gen:
        print(f"  pe=[{gen['gen_pemin']}, {gen['gen_pemax']}], "
              f"θe=[{gen.get('gen_theta_e_min', 0)}, {gen.get('gen_theta_e_max', 180)}]")
    print(f"  Running diffrad...", flush=True)

    t0 = time.time()
    result = subprocess.run(
        [str(GENERATOR), "-input", str(input_file), "-lund", str(lund_file)],
        cwd=str(GEN_DIR), capture_output=True, text=True, timeout=1800)

    if result.returncode != 0:
        print(f"  ERROR:\n{result.stderr}")
        sys.exit(1)

    dt = time.time() - t0
    nlines = sum(1 for _ in open(lund_file))
    nev_out = nlines // 7
    print(f"  Generated {nev_out} events in {dt:.1f}s")

    mc = read_lund_6p(lund_file, cfg["ebeam"])
    mc["pk"] = np.maximum(mc["pkp"], mc["pkm"])
    mc["tmin"] = compute_tmin(mc["q2"], mc["xb"])
    mc["t_tmin"] = mc["t"] - mc["tmin"]
    print(f"  Parsed {len(mc['q2'])} events")

    # Plot
    fig, axes = plt.subplots(2, 5, figsize=(25, 8))
    for col, (var, label, vmin, vmax) in enumerate([
        ("q2", r"$Q^2$ [GeV$^2$]", 0, 10),
        ("xb", r"$x_B$", 0, 0.8),
        ("t",  r"$|t|$ [GeV$^2$]", 0, 5),
        ("w",  r"$W$ [GeV]", 1.5, 5),
        ("pe", r"$p_e$ [GeV]", 0, cfg["ebeam"]),
    ]):
        for row, scale in enumerate(["linear", "log"]):
            ax = axes[row, col]
            ax.hist(mc[var], bins=50, range=(vmin, vmax),
                    histtype="step", lw=2, color="blue")
            ax.set_xlabel(label); ax.grid(True, alpha=0.3)
            if scale == "log": ax.set_yscale("log")
            if row == 0: ax.set_title(label)
    fig.suptitle(f"Step 1: Generated ({len(mc['q2'])} events)\n"
                 f"{cfg['name']}, Ebeam={cfg['ebeam']}, bt={mp['bt']}",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(str(out_dir / "step1_generated.png"), dpi=150)
    plt.close(fig)
    print(f"  Saved: step1_generated.png")
    return mc


def step2_fastmc(mc, cfg):
    """Run NN fast MC on generated events."""
    print(f"\n{'='*70}")
    print(f"  STEP 2: Run NN fast MC")
    print(f"{'='*70}")

    model_dir = str(DATA_ROOT / "phi" / cfg["model"] / "v11" / "models")
    print(f"  Model: {cfg['model']}")

    particle_map = [
        ("e-",  "pe",  "theta_e",  "phi_e",  "vz_e"),
        ("p",   "pp",  "theta_p",  "phi_p",  "vz_p"),
        ("K+",  "pkp", "theta_kp", "phi_kp", "vz_kp"),
        ("K-",  "pkm", "theta_km", "phi_km", "vz_km"),
    ]

    models = {}
    for name, *_ in particle_map:
        m = load_model_auto(model_dir, name)
        if m is None:
            print(f"  ERROR: no model for {name}"); sys.exit(1)
        models[name] = m
        print(f"  Loaded {name}")

    n = len(mc["q2"])
    accepted = np.ones(n, dtype=bool)

    t0 = time.time()
    for name, p_key, th_key, ph_key, vz_key in particle_map:
        mask = models[name].accept(mc[p_key], mc[th_key], mc[ph_key], mc[vz_key])
        accepted &= mask
        print(f"  {name:3s}: {mask.sum()}/{n} ({100*mask.sum()/n:.1f}%)")

    dt = time.time() - t0
    n_acc = accepted.sum()
    print(f"\n  All 4 accepted: {n_acc}/{n} ({100*n_acc/n:.1f}%)")
    print(f"  FastMC time: {dt:.2f}s ({n/dt:.0f} ev/s)")

    out_dir = cfg["out_dir"]
    fig, axes = plt.subplots(2, 5, figsize=(25, 8))
    for col, (var, label, vmin, vmax) in enumerate([
        ("q2", r"$Q^2$ [GeV$^2$]", 0, 10),
        ("xb", r"$x_B$", 0, 0.8),
        ("t",  r"$|t|$ [GeV$^2$]", 0, 5),
        ("w",  r"$W$ [GeV]", 1.5, 5),
        ("pe", r"$p_e$ [GeV]", 0, cfg["ebeam"]),
    ]):
        for row, scale in enumerate(["linear", "log"]):
            ax = axes[row, col]
            ax.hist(mc[var], bins=50, range=(vmin, vmax),
                    histtype="step", lw=1.5, color="blue", label="Generated")
            ax.hist(mc[var][accepted], bins=50, range=(vmin, vmax),
                    histtype="step", lw=2, color="red", label="After fastMC")
            ax.set_xlabel(label); ax.grid(True, alpha=0.3)
            if scale == "log": ax.set_yscale("log")
            if row == 0: ax.set_title(label)
            if col == 0 and row == 0: ax.legend(fontsize=9)
    fig.suptitle(f"Step 2: Generated vs After fastMC\n"
                 f"Generated: {n}, Accepted: {n_acc} ({100*n_acc/n:.1f}%)",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(str(out_dir / "step2_fastmc.png"), dpi=150)
    plt.close(fig)
    print(f"  Saved: step2_fastmc.png")
    return accepted


def step3_compare(mc, accepted, data, cfg):
    """Compare MC+fastMC with data."""
    print(f"\n{'='*70}")
    print(f"  STEP 3: MC vs Data comparison")
    print(f"{'='*70}")

    cuts = cfg["cuts"]
    out_dir = cfg["out_dir"]

    mc_gen_cut = apply_cuts(mc, cuts)
    mc_acc_cut = mc_gen_cut & accepted
    d_cut = apply_cuts(data, cuts)

    print(f"  MC generated (cuts): {mc_gen_cut.sum()}")
    print(f"  MC accepted  (cuts): {mc_acc_cut.sum()}")
    print(f"  Data         (cuts): {d_cut.sum()}")

    if "tmin" not in data:
        data["tmin"] = compute_tmin(data["q2"], data["xb"])
        data["t_tmin"] = data["t"] - data["tmin"]

    nbins = 30
    plot_vars = [
        ("q2", r"$Q^2$ [GeV$^2$]", cuts["q2min"], cuts["q2max"]),
        ("xb", r"$x_B$", 0.05, 0.55),
        ("t",  r"$|t|$ [GeV$^2$]", 0, 4),
        ("w",  r"$W$ [GeV]", cuts["wmin"], 3.5),
        ("t_tmin", r"$|t|-|t_{min}|$ [GeV$^2$]", 0, 3),
    ]

    for scale in ["linear", "log"]:
        fig, axes = plt.subplots(2, len(plot_vars),
                                 figsize=(5.5 * len(plot_vars), 10),
                                 gridspec_kw={"height_ratios": [3, 1]})
        for col, (var, xlabel, vmin, vmax) in enumerate(plot_vars):
            bins = np.linspace(vmin, vmax, nbins + 1)
            bc = 0.5 * (bins[:-1] + bins[1:]); bw = np.diff(bins)[0]
            hd, _ = np.histogram(data[var][d_cut], bins=bins)
            hm, _ = np.histogram(mc[var][mc_acc_cut], bins=bins)
            hd_n = hd / (hd.sum() * bw) if hd.sum() > 0 else hd * 0.0
            hm_n = hm / (hm.sum() * bw) if hm.sum() > 0 else hm * 0.0
            hd_err = np.sqrt(hd) / (hd.sum() * bw) if hd.sum() > 0 else hd * 0.0

            ax = axes[0, col]
            ax.errorbar(bc, hd_n, yerr=hd_err, fmt="ko", ms=4, capsize=2,
                        label=f"Data ({d_cut.sum()})")
            ax.step(bc, hm_n, where="mid", lw=2, color="red",
                    label=f"MC ({mc_acc_cut.sum()})")
            ax.set_ylabel("Normalized"); ax.set_xlim(vmin, vmax)
            ax.set_title(xlabel); ax.legend(fontsize=9); ax.grid(True, alpha=0.3)
            if scale == "log": ax.set_yscale("log")

            ax = axes[1, col]
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where(hd_n > 0, hm_n / hd_n, np.nan)
            valid = np.isfinite(ratio)
            if valid.any(): ax.plot(bc[valid], ratio[valid], "ro", ms=4)
            ax.axhline(1, color="gray", ls="--")
            ax.set_xlabel(xlabel); ax.set_ylabel("MC/Data")
            ax.set_ylim(0, 2.5); ax.set_xlim(vmin, vmax); ax.grid(True, alpha=0.3)

        cut_str = (f"Q²∈({cuts['q2min']},{cuts['q2max']}), W>{cuts['wmin']}, "
                   f"pe>{cuts['pemin']}, pK<{cuts['pkmax']}")
        fig.suptitle(f"Step 3: MC+fastMC vs Data ({scale}) — {cfg['name']}\n"
                     f"{cut_str}", fontsize=12)
        fig.tight_layout()
        sfx = "" if scale == "linear" else "_log"
        fig.savefig(str(out_dir / f"step3_mc_vs_data{sfx}.png"), dpi=150)
        plt.close(fig)
        print(f"  Saved: step3_mc_vs_data{sfx}.png")

    return mc_gen_cut, mc_acc_cut, d_cut


def step4_acceptance(mc, mc_gen_cut, mc_acc_cut, cfg):
    """Compute and plot acceptance."""
    print(f"\n{'='*70}")
    print(f"  STEP 4: Acceptance")
    print(f"{'='*70}")

    out_dir = cfg["out_dir"]
    cuts = cfg["cuts"]
    nbins = 25

    acc_vars = [
        ("q2", r"$Q^2$ [GeV$^2$]", cuts["q2min"], cuts["q2max"]),
        ("xb", r"$x_B$", 0.05, 0.55),
        ("t",  r"$|t|$ [GeV$^2$]", 0.2, 4),
        ("w",  r"$W$ [GeV]", cuts["wmin"], 3.5),
        ("pe", r"$p_e$ [GeV]", cuts["pemin"], cfg["ebeam"]),
        ("t_tmin", r"$|t|-|t_{min}|$ [GeV$^2$]", 0, 3),
    ]

    fig, axes = plt.subplots(2, len(acc_vars),
                             figsize=(5 * len(acc_vars), 8),
                             gridspec_kw={"height_ratios": [3, 1]})

    acceptance_data = {}
    for col, (var, xlabel, vmin, vmax) in enumerate(acc_vars):
        bins = np.linspace(vmin, vmax, nbins + 1)
        bc = 0.5 * (bins[:-1] + bins[1:])
        h_gen, _ = np.histogram(mc[var][mc_gen_cut], bins=bins)
        h_acc, _ = np.histogram(mc[var][mc_acc_cut], bins=bins)
        with np.errstate(divide="ignore", invalid="ignore"):
            acc = np.where(h_gen >= MIN_GEN, h_acc / h_gen, 0.0)
            acc_err = np.where(h_gen >= MIN_GEN,
                               np.sqrt(acc * (1 - acc) / h_gen), 0.0)
        acceptance_data[var] = (bins, acc, acc_err, h_gen, h_acc)

        ax = axes[0, col]
        ax.step(bc, h_gen, where="mid", lw=1.5, color="blue", label="Generated")
        ax.step(bc, h_acc, where="mid", lw=2, color="red", label="Accepted")
        ax.set_xlabel(xlabel); ax.set_ylabel("Events")
        ax.set_title(xlabel); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
        ax.set_yscale("log")

        ax = axes[1, col]
        good = h_gen >= MIN_GEN
        ax.errorbar(bc[good], acc[good], yerr=acc_err[good], fmt="ro", ms=4, capsize=2)
        ax.set_xlabel(xlabel); ax.set_ylabel("Acceptance")
        ax.set_ylim(0, 1); ax.set_xlim(vmin, vmax); ax.grid(True, alpha=0.3)

    n_gen = mc_gen_cut.sum(); n_acc = mc_acc_cut.sum()
    fig.suptitle(f"Step 4: Acceptance — {cfg['name']}\n"
                 f"model={cfg['model']}, overall={n_acc}/{n_gen} "
                 f"({100*n_acc/max(n_gen,1):.1f}%)", fontsize=13)
    fig.tight_layout()
    fig.savefig(str(out_dir / "step4_acceptance.png"), dpi=150)
    plt.close(fig)
    print(f"  Saved: step4_acceptance.png")
    print(f"  Overall: {n_acc}/{n_gen} ({100*n_acc/max(n_gen,1):.1f}%)")
    return acceptance_data


def step5_fit(mc, mc_gen_cut, mc_acc_cut, data, d_cut, cfg):
    """Acceptance-correct data and fit free parameters."""
    print(f"\n{'='*70}")
    print(f"  STEP 5: Acceptance correction + Fit")
    print(f"{'='*70}")

    out_dir = cfg["out_dir"]
    cuts = cfg["cuts"]
    mp = cfg["model_parameters"]
    fp = cfg["fit_parameters"]
    nbins_fit = 25

    free_params = fp["free"]
    fixed_vals = {k: mp[k] for k in fp.get("fixed", [])}
    bounds = fp["bounds"]

    if "tmin" not in data:
        data["tmin"] = compute_tmin(data["q2"], data["xb"])
        data["t_tmin"] = data["t"] - data["tmin"]

    fit_vars = [
        ("q2", r"$Q^2$ [GeV$^2$]", cuts["q2min"], cuts["q2max"]),
        ("t",  r"$|t|$ [GeV$^2$]", 0.5, 4.0),
        ("w",  r"$W$ [GeV]", cuts["wmin"], 3.5),
    ]
    plot_vars = fit_vars + [
        ("t_tmin", r"$|t|-|t_{min}|$ [GeV$^2$]", 0.0, 3.0),
        ("xb", r"$x_B$", 0.1, 0.55),
    ]

    # Build acceptance-corrected histograms
    hist_data, hist_bins = {}, {}
    for var, xlabel, vmin, vmax in plot_vars:
        bins = np.linspace(vmin, vmax, nbins_fit + 1)
        h_gen, _ = np.histogram(mc[var][mc_gen_cut], bins=bins)
        h_acc, _ = np.histogram(mc[var][mc_acc_cut], bins=bins)
        h_dat, _ = np.histogram(data[var][d_cut], bins=bins)
        with np.errstate(divide="ignore", invalid="ignore"):
            acceptance = np.where(h_gen >= MIN_GEN, h_acc / h_gen, 0.0)
            h_corr = np.where(acceptance > 0, h_dat / acceptance, 0.0)
            h_corr_err = np.where(acceptance > 0,
                                  np.sqrt(np.maximum(h_dat, 1)) / acceptance, 0.0)
        good = acceptance > 0
        h_corr[~good] = 0; h_corr_err[~good] = 0
        hist_data[var] = (h_corr, h_corr_err, good)
        hist_bins[var] = bins

    # MC arrays for reweighting
    mc_q2 = mc["q2"][mc_gen_cut].astype(np.float64)
    mc_xb = mc["xb"][mc_gen_cut].astype(np.float64)
    mc_t_pos = mc["t"][mc_gen_cut].astype(np.float64)
    mc_t_neg = -mc_t_pos
    mc_w  = mc["w"][mc_gen_cut].astype(np.float64)
    mc_ttm = mc["t_tmin"][mc_gen_cut].astype(np.float64)
    mc_dict = {"q2": mc_q2, "xb": mc_xb, "t": mc_t_pos,
               "w": mc_w, "t_tmin": mc_ttm}

    # Generation cross section
    gen_alf2 = mp["alf2"]; gen_alf3 = mp["alf3"]
    gen_bt = mp["bt"]; gen_nuT = mp["nuT"]; gen_cR = mp["cR"]
    s_old = sigma_model(mc_q2, mc_xb, mc_t_neg,
                        gen_alf2, gen_alf3, gen_bt, gen_nuT, gen_cR)

    # Subsample for speed
    N_FIT = min(300000, len(mc_q2))
    rng = np.random.default_rng(42)
    idx = np.sort(rng.choice(len(mc_q2), N_FIT, replace=False))
    f_q2 = mc_q2[idx]; f_xb = mc_xb[idx]; f_t_neg = mc_t_neg[idx]
    f_old = s_old[idx]
    f_dict = {"q2": f_q2, "xb": f_xb, "t": mc_t_pos[idx],
              "w": mc_w[idx], "t_tmin": mc_ttm[idx]}

    f_bin = {}
    for var, _, vmin, vmax in fit_vars:
        f_bin[var] = np.digitize(f_dict[var], hist_bins[var]) - 1

    print(f"  Free parameters: {free_params}")
    print(f"  Fixed: {fixed_vals}")
    print(f"  Fit subsample: {N_FIT}")

    def make_params(x):
        """Map fit vector x to full parameter dict."""
        p = dict(fixed_vals)
        for i, name in enumerate(free_params):
            p[name] = x[i]
        # Fill any remaining from generation values
        for k in ["alf2", "alf3", "bt", "nuT", "cR"]:
            if k not in p:
                p[k] = mp[k]
        return p

    def total_chi2(x):
        p = make_params(x)
        if p["bt"] <= 0.01: return 1e8
        s_new = sigma_model(f_q2, f_xb, f_t_neg,
                            p["alf2"], p["alf3"], p["bt"], p["nuT"], p["cR"])
        with np.errstate(divide="ignore", invalid="ignore"):
            w = np.where(f_old > 0, s_new / f_old, 0.0)
        chi2 = 0.0
        for var, _, vmin, vmax in fit_vars:
            nb = len(hist_bins[var]) - 1
            h_corr, _, good = hist_data[var]
            h_model = np.zeros(nb)
            bidx = f_bin[var]
            for i in range(nb): h_model[i] = w[bidx == i].sum()
            chi2 += chi2_shape(h_corr * good, h_model * good)
        return chi2

    # Grid scan
    print("  Grid scan...", flush=True)
    t0 = time.time()
    grid_ranges = {}
    for name in free_params:
        lo, hi = bounds[name]
        grid_ranges[name] = np.linspace(lo, hi, 12)

    # Build grid (up to 3D)
    best_chi2 = 1e10
    best_x = [mp.get(n, 0) for n in free_params]

    if len(free_params) == 2:
        for v0 in grid_ranges[free_params[0]]:
            for v1 in grid_ranges[free_params[1]]:
                c = total_chi2([v0, v1])
                if c < best_chi2:
                    best_chi2 = c; best_x = [v0, v1]
    elif len(free_params) == 3:
        for v0 in grid_ranges[free_params[0]]:
            for v1 in grid_ranges[free_params[1]]:
                for v2 in grid_ranges[free_params[2]]:
                    c = total_chi2([v0, v1, v2])
                    if c < best_chi2:
                        best_chi2 = c; best_x = [v0, v1, v2]
    elif len(free_params) == 4:
        # Coarser grid for 4D: 8 points per dimension = 4096 evaluations
        grid_ranges_4d = {}
        for name in free_params:
            lo, hi = bounds[name]
            grid_ranges_4d[name] = np.linspace(lo, hi, 8)
        for v0 in grid_ranges_4d[free_params[0]]:
            for v1 in grid_ranges_4d[free_params[1]]:
                for v2 in grid_ranges_4d[free_params[2]]:
                    for v3 in grid_ranges_4d[free_params[3]]:
                        c = total_chi2([v0, v1, v2, v3])
                        if c < best_chi2:
                            best_chi2 = c; best_x = [v0, v1, v2, v3]
    else:
        # >4: skip grid, go straight to DE
        pass

    print(f"  Grid: chi2={best_chi2:.1f}, "
          + ", ".join(f"{n}={v:.3f}" for n, v in zip(free_params, best_x))
          + f" ({time.time()-t0:.0f}s)", flush=True)

    # DE refinement
    print("  DE refinement...", flush=True)
    t0 = time.time()
    de_bounds = [tuple(bounds[n]) for n in free_params]
    res = differential_evolution(total_chi2, de_bounds,
                                 seed=42, maxiter=300, tol=0.001,
                                 x0=best_x, popsize=15)
    best_x = list(res.x); best_chi2 = res.fun
    best_params = make_params(best_x)
    print(f"  DE done ({time.time()-t0:.0f}s)")

    # Print results
    print(f"\n  {'='*50}")
    print(f"  BEST FIT")
    print(f"  chi2/ndf = {best_chi2:.2f}")
    for name in ["alf2", "alf3", "bt", "nuT", "cR"]:
        val = best_params[name]
        gen_val = mp[name]
        status = "(free)" if name in free_params else "(fixed)"
        print(f"  {name:5s} = {val:8.3f}   gen={gen_val:8.3f}  {status}")
    print(f"  {'='*50}")

    # Best weights (full sample)
    s_best = sigma_model(mc_q2, mc_xb, mc_t_neg,
                         best_params["alf2"], best_params["alf3"],
                         best_params["bt"], best_params["nuT"], best_params["cR"])
    with np.errstate(divide="ignore", invalid="ignore"):
        w_best = np.where(s_old > 0, s_best / s_old, 0)

    # Per-variable chi2
    for var, xlabel, vmin, vmax in plot_vars:
        bins = hist_bins[var]
        h_corr, _, good = hist_data[var]
        h_model, _ = np.histogram(mc_dict[var], bins=bins, weights=w_best)
        c = chi2_shape(h_corr * good, h_model * good)
        fitted = "(fitted)" if var in [v for v, _, _, _ in fit_vars] else ""
        print(f"  chi2({var:>6s}) = {c:.2f}  {fitted}")

    # ─── Plots ───────────────────────────────────────────────
    ncols = len(plot_vars)
    param_str = ", ".join(f"{n}={best_params[n]:.2f}" for n in free_params)

    for scale in ["linear", "log"]:
        fig, axes = plt.subplots(2, ncols, figsize=(5.5 * ncols, 10),
                                 gridspec_kw={"height_ratios": [3, 1]})
        for col, (var, xlabel, vmin, vmax) in enumerate(plot_vars):
            bins = hist_bins[var]
            bc = 0.5 * (bins[:-1] + bins[1:]); bw = np.diff(bins)[0]
            h_corr, h_corr_err, good = hist_data[var]
            h_orig, _ = np.histogram(mc_dict[var], bins=bins)
            h_fit, _ = np.histogram(mc_dict[var], bins=bins, weights=w_best)

            def norm(h, g=good, bwidth=bw):
                s = np.sum(h[g] * bwidth)
                return h / s if s > 0 else h * 0
            h_corr_n = norm(h_corr)
            cs = np.sum(h_corr[good] * bw)
            h_corr_err_n = h_corr_err / cs if cs > 0 else h_corr_err * 0
            h_orig_n = norm(h_orig); h_fit_n = norm(h_fit)

            is_fitted = var in [v for v, _, _, _ in fit_vars]
            ax = axes[0, col]
            ax.errorbar(bc[good], h_corr_n[good], yerr=h_corr_err_n[good],
                        fmt="ko", ms=5, capsize=3, label="Data/Acceptance")
            ax.step(bc, h_orig_n, where="mid", lw=1.5, color="gray",
                    ls="--", label="Generation model")
            ax.step(bc, h_fit_n, where="mid", lw=2.5, color="red", label="Best fit")
            ax.set_ylabel("Normalized"); ax.set_xlim(vmin, vmax)
            ax.set_title(xlabel + ("" if is_fitted else " (not fitted)"))
            ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
            if scale == "log": ax.set_yscale("log")

            ax = axes[1, col]
            with np.errstate(divide="ignore", invalid="ignore"):
                ratio = np.where((h_corr_n > 0) & good, h_fit_n / h_corr_n, np.nan)
            valid = np.isfinite(ratio)
            if valid.any(): ax.plot(bc[valid], ratio[valid], "ro", ms=4)
            ax.axhline(1, color="gray", ls="--")
            ax.set_xlabel(xlabel); ax.set_ylabel("Fit/Data")
            ax.set_ylim(0, 2.5); ax.set_xlim(vmin, vmax); ax.grid(True, alpha=0.3)

        fig.suptitle(
            f"Step 5: Corrected data vs fit ({scale}) — {cfg['name']}\n"
            f"{param_str}   [χ²={best_chi2:.1f}]\n"
            f"Cuts: Q²∈({cuts['q2min']},{cuts['q2max']}), W>{cuts['wmin']}, "
            f"pe>{cuts['pemin']}, pK<{cuts['pkmax']}",
            fontsize=12)
        fig.tight_layout()
        sfx = "" if scale == "linear" else "_log"
        fig.savefig(str(out_dir / f"step5_fit{sfx}.png"), dpi=150)
        plt.close(fig)
        print(f"  Saved: step5_fit{sfx}.png")

    # Save results to JSON
    param_desc = {
        "alf2": "Threshold exponent: (1 - W2_th/W2)^alf2",
        "alf3": "W power law: W^alf3",
        "bt":   "t-slope [GeV^-2]: exp(bt * t)",
        "nuT":  "Q2 suppression: (1 + Q2/mphi2)^(-nuT)",
        "cR":   "sigma_L/sigma_T ratio: sigma_L = cR * Q2/mphi2 * sigma_T",
    }
    result = {
        "chi2": best_chi2,
        "formula": "dsig/dt = (1-W2th/W2)^alf2 * W^alf3 * (1+Q2/mphi2)^(-nuT) * bt*exp(bt*t) * (1+cR*Q2/mphi2)",
        "fit_result": {
            n: {"value": best_params[n],
                "status": "free" if n in free_params else "fixed",
                "gen_value": mp[n],
                "description": param_desc.get(n, "")}
            for n in ["alf2", "alf3", "bt", "nuT", "cR"]
        },
    }
    with open(out_dir / "results.json", "w") as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: results.json")

    return result


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print("Usage: python3 phi_pipeline.py <config.json>")
        print("  e.g.: python3 phi_pipeline.py runs/rgafall18_inb/config.json")
        sys.exit(0)

    config_path = Path(sys.argv[1])
    if not config_path.exists():
        print(f"ERROR: {config_path} not found"); sys.exit(1)

    with open(config_path) as f:
        raw = json.load(f)

    # Flatten config into a single dict
    ds = raw["dataset"]
    out_dir = config_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # Extract parameter values from cross_section_model (new format)
    # or model_parameters (old format) for backward compatibility
    if "cross_section_model" in raw:
        mp = {k: v["value"] for k, v in raw["cross_section_model"]["parameters"].items()}
    else:
        mp = raw["model_parameters"]

    cfg = {
        "name":             ds["name"],
        "model":            ds["model"],
        "ebeam":            ds["ebeam"],
        "out_dir":          out_dir,
        "generation":       raw["generation"],
        "model_parameters": mp,
        "fit_parameters":   raw["fit_parameters"],
        "cuts":             raw["cuts"],
    }
    gen = cfg["generation"]
    cuts = cfg["cuts"]

    print(f"\n{'#'*70}")
    print(f"  PHI PIPELINE — {cfg['name']}")
    print(f"  Config: {config_path}")
    print(f"  Ebeam={cfg['ebeam']}, nev={gen['nev']}")
    print(f"  Model params: alf2={mp['alf2']}, alf3={mp['alf3']}, bt={mp['bt']}, "
          f"nuT={mp['nuT']}, cR={mp['cR']}")
    print(f"  Free: {cfg['fit_parameters']['free']}")
    print(f"  Cuts: Q2=[{cuts['q2min']},{cuts['q2max']}], W>{cuts['wmin']}, "
          f"pe>{cuts['pemin']}, pK<{cuts['pkmax']}")
    print(f"  Output: {out_dir}")
    print(f"{'#'*70}")

    t_start = time.time()

    mc = step1_generate(cfg)
    accepted = step2_fastmc(mc, cfg)

    data_file = LUND_DIR / ds["data_file"]
    if not data_file.exists():
        print(f"  ERROR: {data_file} not found"); sys.exit(1)
    print(f"\n  Loading data: {data_file.name}...", flush=True)
    data = read_data_lund(data_file, cfg["ebeam"])
    data["pk"] = np.maximum(data["pkp"], data["pkm"])
    data["tmin"] = compute_tmin(data["q2"], data["xb"])
    data["t_tmin"] = data["t"] - data["tmin"]
    print(f"  Data events: {len(data['q2'])}")

    mc_gen_cut, mc_acc_cut, d_cut = step3_compare(mc, accepted, data, cfg)
    step4_acceptance(mc, mc_gen_cut, mc_acc_cut, cfg)
    result = step5_fit(mc, mc_gen_cut, mc_acc_cut, data, d_cut, cfg)

    dt_total = time.time() - t_start
    print(f"\n{'#'*70}")
    print(f"  PIPELINE COMPLETE — {dt_total:.0f}s")
    print(f"  {cfg['name']}")
    fr = result["fit_result"]
    print(f"  Result: " + ", ".join(f"{k}={fr[k]['value']:.3f}" for k in cfg["fit_parameters"]["free"]))
    print(f"  Chi2:   {result['chi2']:.2f}")
    print(f"  Plots:  {out_dir}/")
    print(f"{'#'*70}\n")


if __name__ == "__main__":
    main()
