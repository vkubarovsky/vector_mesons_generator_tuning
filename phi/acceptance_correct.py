#!/usr/bin/env python3
"""
acceptance_correct.py — Compute acceptance from MC+fastMC, correct N(φ) from
M(KK) fits, fit cross-section model, make single-page plot per dataset.

Uses existing MC LUND files (no regeneration).
Uses N(φ) from mkk_binned_results.json for signal yields.

Usage:
    python3 acceptance_correct.py                   # all datasets
    python3 acceptance_correct.py rgafall18_outb     # one dataset
"""

import json
import sys
import time
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
MP  = 0.938272;  MP2 = MP**2
MPHI = 1.019412; MPHI2 = MPHI**2
WTH2 = (MP + MPHI)**2

BASE = Path("/Users/vpk/vector_mesons_generator_tuning/phi")
DATA_ROOT = Path.home() / "Downloads/volatile/clas12/vpk/fastmc"
MIN_GEN = 30  # minimum generated events per bin for acceptance

DATASETS = ["rgafall18_inb", "rgafall18_outb", "rgasp18_inb",
            "rgasp18_outb", "rgasp19_inb"]

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
    return (ecm_i - ecm_f)**2 - (pcm_i - pcm_f)**2


def sigma_model(q2, xb, t_neg, alf2, alf3, bt, nuT, cR):
    """
    dσ/dt ∝ (1 - W²_th/W²)^α₂ · W^α₃ · (1+Q²/m²φ)^(-νT)
            · bt · exp(bt·t) · (1 + cR·Q²/m²φ)
    t_neg < 0, bt > 0
    """
    w2 = MP2 + q2 * (1 - xb) / xb
    w = np.sqrt(np.maximum(w2, 0.01))
    thresh = np.maximum(1.0 - WTH2 / w2, 1e-30)**alf2
    sig_T = thresh * w**alf3 * (1.0 + q2 / MPHI2)**(-nuT) \
            * bt * np.exp(bt * t_neg)
    sig_L = cR * (q2 / MPHI2) * sig_T
    return sig_T + sig_L


# ═══════════════════════════════════════════════════════════════════
# I/O
# ═══════════════════════════════════════════════════════════════════

def read_lund_6p(fpath, ebeam):
    """Read 6-particle LUND from diffrad generator."""
    q2, xb, t_arr, w_arr = [], [], [], []
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
                    "vz": float(pc[13]) if len(pc) > 13 else 0.0,
                })

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
            t_arr.append(t_val); w_arr.append(w_val)

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

    mc = {
        "q2": np.array(q2), "xb": np.array(xb),
        "t": np.array(t_arr), "w": np.array(w_arr),
        "pe": np.array(pe), "pp": np.array(pp),
        "pkp": np.array(pkp), "pkm": np.array(pkm),
        "theta_e": np.array(theta_e), "theta_p": np.array(theta_p),
        "theta_kp": np.array(theta_kp), "theta_km": np.array(theta_km),
        "phi_e": np.array(phi_e), "phi_p": np.array(phi_p),
        "phi_kp": np.array(phi_kp), "phi_km": np.array(phi_km),
        "vz_e": np.array(vz_e), "vz_p": np.array(vz_p),
        "vz_kp": np.array(vz_kp), "vz_km": np.array(vz_km),
    }
    mc["pk"] = np.maximum(mc["pkp"], mc["pkm"])
    mc["tmin"] = compute_tmin(mc["q2"], mc["xb"])
    mc["t_tmin"] = mc["t"] - mc["tmin"]
    return mc


def apply_cuts(d, cuts):
    mask = ((d["q2"] >= cuts["q2min"]) & (d["q2"] <= cuts["q2max"]) &
            (d["w"] >= cuts["wmin"]) & (d["pe"] >= cuts["pemin"]))
    if cuts.get("wmax", 100) < 100:
        mask &= d["w"] <= cuts["wmax"]
    if "pk" in d:
        mask &= d["pk"] <= cuts["pkmax"]
    if cuts.get("pemax", 100) < 100:
        mask &= d["pe"] <= cuts["pemax"]
    if cuts.get("theta_e_min", 0) > 0 and "theta_e" in d:
        mask &= d["theta_e"] >= cuts["theta_e_min"]
    if cuts.get("theta_e_max", 180) < 180 and "theta_e" in d:
        mask &= d["theta_e"] <= cuts["theta_e_max"]
    return mask


# ═══════════════════════════════════════════════════════════════════
# Main processing
# ═══════════════════════════════════════════════════════════════════

def process_dataset(dataset):
    print(f"\n{'#'*70}")
    print(f"  {dataset}")
    print(f"{'#'*70}")

    # Load config
    cfg_path = BASE / "runs" / dataset / "config.json"
    with open(cfg_path) as f:
        raw = json.load(f)
    ds = raw["dataset"]
    ebeam = ds["ebeam"]
    model_name = ds["model"]
    cuts = raw["cuts"]
    mkk_bins = raw["mkk_bins"]

    if "cross_section_model" in raw:
        mp = {k: v["value"] for k, v in raw["cross_section_model"]["parameters"].items()}
    else:
        mp = raw["model_parameters"]

    out_dir = BASE / "runs" / dataset

    # Load N(φ) from M(KK) fits
    mkk_json = out_dir / "mkk_binned_results.json"
    with open(mkk_json) as f:
        mkk_results = json.load(f)
    print(f"  Loaded N(φ) from {mkk_json.name}")

    # ─── Load or build MC arrays ────────────────────────────────
    # npz caches are large (~700 MB) — keep them in the Downloads mirror, not the repo
    cache_dir = Path.home() / "Downloads/volatile/clas12/vpk/fastmc/phi_tuning_cache" / dataset
    cache_dir.mkdir(parents=True, exist_ok=True)
    npz_path = cache_dir / "mc_arrays.npz"
    if npz_path.exists():
        print(f"  Loading cached MC: {npz_path.name}")
        t0 = time.time()
        npz = np.load(npz_path)
        mc = {k: npz[k] for k in npz.files}
        accepted = mc.pop("accepted")
        print(f"  Loaded {len(mc['q2'])} events ({time.time()-t0:.1f}s)")
    else:
        # Read LUND + run fastMC (slow, done once)
        lund_file = DATA_ROOT / "phi" / model_name / "v11" / "lund" / "generated.lund"
        print(f"  Reading MC LUND: {lund_file.name}...", flush=True)
        t0 = time.time()
        mc = read_lund_6p(lund_file, ebeam)
        print(f"  Parsed {len(mc['q2'])} MC events ({time.time()-t0:.1f}s)")

        model_dir = str(DATA_ROOT / "phi" / model_name / "v11" / "models")
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
                print(f"  ERROR: no model for {name}"); return None
            models[name] = m

        n = len(mc["q2"])
        accepted = np.ones(n, dtype=bool)
        t0 = time.time()
        for name, p_key, th_key, ph_key, vz_key in particle_map:
            mask = models[name].accept(mc[p_key], mc[th_key], mc[ph_key], mc[vz_key])
            accepted &= mask
            print(f"  {name:3s}: {mask.sum()}/{n} ({100*mask.sum()/n:.1f}%)")
        n_acc = accepted.sum()
        print(f"  All 4: {n_acc}/{n} ({100*n_acc/n:.1f}%) [{time.time()-t0:.1f}s]")

        # Save for future runs
        np.savez_compressed(str(npz_path), accepted=accepted, **mc)
        print(f"  Cached MC → {npz_path.name} ({npz_path.stat().st_size/1e6:.0f} MB)")

    # Apply kinematic cuts
    mc_gen_cut = apply_cuts(mc, cuts)
    mc_acc_cut = mc_gen_cut & accepted
    print(f"  After cuts: gen={mc_gen_cut.sum()}, acc={mc_acc_cut.sum()}")

    # ─── Compute acceptance & correct N(φ) ────────────────────
    # Variables to process
    var_info = {
        "q2":     {"label": "Q²",         "unit": "[GeV²]", "log_y": True},
        "t":      {"label": "|t|",         "unit": "[GeV²]", "log_y": True},
        "t_tmin": {"label": "|t|-t_min",   "unit": "[GeV²]", "log_y": True},
        "xb":     {"label": "x_B",         "unit": "",       "log_y": False},
        "w":      {"label": "W",           "unit": "[GeV]",  "log_y": False},
    }

    corrected_data = {}
    for var_name, vinfo in var_info.items():
        bin_edges = np.array(mkk_bins[var_name], dtype=float)
        n_bins = len(bin_edges) - 1
        bc = 0.5 * (bin_edges[:-1] + bin_edges[1:])
        bw = np.diff(bin_edges)

        # MC acceptance in these bins
        h_gen, _ = np.histogram(mc[var_name][mc_gen_cut], bins=bin_edges)
        h_acc, _ = np.histogram(mc[var_name][mc_acc_cut], bins=bin_edges)
        with np.errstate(divide="ignore", invalid="ignore"):
            acc = np.where(h_gen >= MIN_GEN, h_acc / h_gen, 0.0)
            acc_err = np.where(h_gen >= MIN_GEN,
                               np.sqrt(acc * (1 - acc) / h_gen), 0.0)

        # N(φ) from M(KK) fits
        if var_name in mkk_results:
            bins_data = mkk_results[var_name]["bins"]
            n_phi = np.array([b["N_phi"] for b in bins_data])
            dn_phi = np.array([b["dN_phi"] for b in bins_data])
        else:
            print(f"  WARNING: {var_name} not in mkk_results, skipping")
            continue

        # Acceptance-corrected yields
        with np.errstate(divide="ignore", invalid="ignore"):
            n_corr = np.where(acc > 0, n_phi / acc, 0.0)
            # Error propagation: δ(N/ε) = (N/ε) * sqrt((δN/N)² + (δε/ε)²)
            rel_err_phi = np.where(n_phi > 0, dn_phi / n_phi, 0)
            rel_err_acc = np.where(acc > 0, acc_err / acc, 0)
            dn_corr = n_corr * np.sqrt(rel_err_phi**2 + rel_err_acc**2)

        good = (acc > 0) & (n_phi > 0)

        corrected_data[var_name] = {
            "bin_edges": bin_edges, "bc": bc, "bw": bw,
            "n_phi": n_phi, "dn_phi": dn_phi,
            "acc": acc, "acc_err": acc_err,
            "n_corr": n_corr, "dn_corr": dn_corr,
            "good": good, "h_gen": h_gen, "h_acc": h_acc,
        }

        print(f"\n  {vinfo['label']}:")
        for i in range(n_bins):
            status = "" if good[i] else " (skip)"
            print(f"    [{bin_edges[i]:.1f},{bin_edges[i+1]:.1f}): "
                  f"N(φ)={n_phi[i]:.0f}, acc={acc[i]:.3f}, "
                  f"N_corr={n_corr[i]:.0f}{status}")

    # ─── Fit model to corrected distributions ─────────────────
    # Full MC arrays (for final model histograms)
    mc_q2_full = mc["q2"][mc_gen_cut].astype(np.float64)
    mc_xb_full = mc["xb"][mc_gen_cut].astype(np.float64)
    mc_t_neg_full = -mc["t"][mc_gen_cut].astype(np.float64)
    mc_t_pos_full = mc["t"][mc_gen_cut].astype(np.float64)
    mc_w_full = mc["w"][mc_gen_cut].astype(np.float64)
    mc_ttm_full = mc["t_tmin"][mc_gen_cut].astype(np.float64)

    mc_vars_full = {"q2": mc_q2_full, "t": mc_t_pos_full,
                    "t_tmin": mc_ttm_full, "xb": mc_xb_full, "w": mc_w_full}

    s_old_full = sigma_model(mc_q2_full, mc_xb_full, mc_t_neg_full,
                             mp["alf2"], mp["alf3"], mp["bt"], mp["nuT"], mp["cR"])

    # Subsample for fitting speed (300k is plenty for 5-bin shapes)
    N_FIT = min(300000, len(mc_q2_full))
    rng = np.random.default_rng(42)
    idx_sub = np.sort(rng.choice(len(mc_q2_full), N_FIT, replace=False))

    mc_q2 = mc_q2_full[idx_sub]
    mc_xb = mc_xb_full[idx_sub]
    mc_t_neg = mc_t_neg_full[idx_sub]
    mc_t_pos = mc_t_pos_full[idx_sub]
    mc_w = mc_w_full[idx_sub]
    mc_ttm = mc_ttm_full[idx_sub]
    s_old = s_old_full[idx_sub]

    mc_vars = {"q2": mc_q2, "t": mc_t_pos, "t_tmin": mc_ttm,
               "xb": mc_xb, "w": mc_w}
    print(f"  Fit subsample: {N_FIT} events")

    # Pre-bin the subsample for fast chi2 (no np.histogram in loop!)
    fit_vars = ["q2", "t", "w"]  # only vars with model parameters

    fit_bin_data = {}
    for var_name in fit_vars:
        if var_name not in corrected_data:
            continue
        cd = corrected_data[var_name]
        n_bins = len(cd["bin_edges"]) - 1
        idx = np.digitize(mc_vars[var_name], cd["bin_edges"]) - 1
        # Build list of event indices per bin
        bin_masks = []
        for b in range(n_bins):
            bin_masks.append(np.where(idx == b)[0])
        good = cd["good"]
        s_data = np.sum(cd["n_corr"][good])
        fit_bin_data[var_name] = {
            "bin_masks": bin_masks, "n_bins": n_bins,
            "good": good, "s_data": s_data,
            "hd_n": cd["n_corr"] / s_data if s_data > 0 else cd["n_corr"],
            "err_n": cd["dn_corr"] / s_data if s_data > 0 else cd["dn_corr"],
        }

    # Free and fixed parameters
    fp = raw["fit_parameters"]
    free_params = fp["free"]
    bounds_dict = fp["bounds"]

    def make_params(x):
        p = dict(mp)  # start from generation values
        for i, name in enumerate(free_params):
            p[name] = x[i]
        return p

    def total_chi2(x):
        p = make_params(x)
        if p["bt"] <= 0.01:
            return 1e8
        s_new = sigma_model(mc_q2, mc_xb, mc_t_neg,
                            p["alf2"], p["alf3"], p["bt"], p["nuT"], p["cR"])
        with np.errstate(divide="ignore", invalid="ignore"):
            w = np.where(s_old > 0, s_new / s_old, 0.0)

        chi2_total = 0.0
        for var_name, fbd in fit_bin_data.items():
            good = fbd["good"]
            s_data = fbd["s_data"]
            if s_data <= 0:
                continue
            # Fast binned sum using precomputed masks
            h_model = np.array([np.sum(w[m]) for m in fbd["bin_masks"]])
            s_model = np.sum(h_model[good])
            if s_model <= 0:
                continue
            hm_n = h_model / s_model
            for i in range(fbd["n_bins"]):
                if good[i] and fbd["err_n"][i] > 0:
                    chi2_total += ((fbd["hd_n"][i] - hm_n[i]) / fbd["err_n"][i])**2

        return chi2_total

    # Optimization with differential evolution (bounded)
    print(f"\n  Fitting: free={free_params}")
    x0 = [mp[n] for n in free_params]
    bounds_list = [tuple(bounds_dict[n]) for n in free_params]

    t0 = time.time()
    res = differential_evolution(total_chi2, bounds_list,
                                 seed=42, maxiter=300, tol=0.001,
                                 x0=x0, popsize=15)
    best_x = res.x
    best_chi2 = res.fun
    best_params = make_params(list(best_x))
    print(f"  Fit done ({time.time()-t0:.1f}s), chi2={best_chi2:.1f}")

    for name in ["alf2", "alf3", "bt", "nuT", "cR"]:
        status = "(free)" if name in free_params else "(fixed)"
        print(f"    {name:5s} = {best_params[name]:8.4f}  {status}")

    # Compute BOTH original model and best-fit model histograms
    # Original model (generation params) — use full MC
    orig_hists = {}
    for var_name in var_info:
        if var_name not in corrected_data:
            continue
        cd = corrected_data[var_name]
        # Original model: weights = 1 (generation model IS the MC)
        h_orig, _ = np.histogram(mc_vars_full[var_name], bins=cd["bin_edges"])
        orig_hists[var_name] = h_orig

    # Best-fit model — reweight full MC (s_old_full already computed above)
    s_best = sigma_model(mc_q2_full, mc_xb_full, mc_t_neg_full,
                         best_params["alf2"], best_params["alf3"],
                         best_params["bt"], best_params["nuT"], best_params["cR"])
    with np.errstate(divide="ignore", invalid="ignore"):
        w_best = np.where(s_old_full > 0, s_best / s_old_full, 0.0)

    fit_hists = {}
    chi2_per_var = {}
    for var_name in var_info:
        if var_name not in corrected_data:
            continue
        cd = corrected_data[var_name]
        h_fit, _ = np.histogram(mc_vars_full[var_name], bins=cd["bin_edges"],
                                weights=w_best)
        fit_hists[var_name] = h_fit

        # Compute per-variable chi2
        good = cd["good"]
        s_data = np.sum(cd["n_corr"][good])
        s_model = np.sum(h_fit[good])
        if s_data > 0 and s_model > 0:
            hd_n = cd["n_corr"] / s_data
            hm_n = h_fit / s_model
            err_n = cd["dn_corr"] / s_data
            c2 = 0
            npts = 0
            for i in range(len(good)):
                if good[i] and err_n[i] > 0:
                    c2 += ((hd_n[i] - hm_n[i]) / err_n[i])**2
                    npts += 1
            chi2_per_var[var_name] = (c2, npts)
        else:
            chi2_per_var[var_name] = (0, 0)

    # ─── PLOT: single page 3×2 — original model + fit ────────
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    fig.subplots_adjust(hspace=0.40, wspace=0.35)

    plot_order = ["q2", "t", "t_tmin", "xb", "w"]
    positions = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]

    for idx, var_name in enumerate(plot_order):
        if var_name not in corrected_data:
            continue
        row, col = positions[idx]
        ax = axes[row, col]
        cd = corrected_data[var_name]
        vinfo = var_info[var_name]
        good = cd["good"]

        bc = cd["bc"]
        bw = cd["bw"]
        n_corr = cd["n_corr"]
        dn_corr = cd["dn_corr"]
        h_orig = orig_hists[var_name]
        h_fit = fit_hists[var_name]

        s_data = np.sum(n_corr[good])

        # Original model (blue, dashed step)
        s_orig = np.sum(h_orig[good])
        if s_data > 0 and s_orig > 0:
            h_orig_s = h_orig * (s_data / s_orig)
            ax.step(cd["bin_edges"], np.append(h_orig_s, h_orig_s[-1]),
                    where="post", lw=2, color="blue", ls="--", label="Gen model")

        # Fitted model (red filled)
        s_fit = np.sum(h_fit[good])
        if s_data > 0 and s_fit > 0:
            h_fit_s = h_fit * (s_data / s_fit)
            ax.bar(bc, h_fit_s, width=bw, color="lightcoral",
                   edgecolor="none", alpha=0.5, label="Fit model")

        # Data points on top
        if good.any():
            ax.errorbar(bc[good], n_corr[good], yerr=dn_corr[good],
                        xerr=bw[good]/2, fmt="ko", ms=4, capsize=2,
                        lw=1.2, label="Data/ε", zorder=5)

        ax.set_xlabel(f"{vinfo['label']} {vinfo['unit']}", fontsize=10)
        ax.set_ylabel("N(φ) / ε", fontsize=10)
        ax.set_xlim(cd["bin_edges"][0], cd["bin_edges"][-1])
        ax.grid(True, alpha=0.4, lw=0.5)
        ax.legend(fontsize=7, loc="upper right")

        if vinfo["log_y"]:
            ax.set_yscale("log")
            if good.any():
                ymin = max(1, np.min(n_corr[good]) * 0.3)
                ax.set_ylim(bottom=ymin)
        else:
            ax.set_ylim(bottom=0)

    # ─── Combined results box (bottom right) ─────────────────
    ax_box = axes[1, 2]
    ax_box.axis("off")

    lines = []
    lines.append(f"  {dataset}")
    lines.append(f"  {'─'*32}")
    lines.append(f"  {'Param':6s} {'Gen':>8s} {'Fit':>8s}  {'Status'}")
    lines.append(f"  {'─'*32}")
    for name in ["alf2", "alf3", "bt", "nuT", "cR"]:
        status = "free" if name in free_params else "fixed"
        lines.append(f"  {name:5s}  {mp[name]:8.3f} {best_params[name]:8.3f}  ({status})")
    lines.append(f"  {'─'*32}")
    lines.append(f"  χ²(fit) = {best_chi2:.1f}")
    lines.append(f"")
    lines.append(f"  M(KK): μ=1019.5 σ=4.7 (fix)")

    box_text = "\n".join(lines)
    props = dict(boxstyle="round,pad=0.6", facecolor="lightyellow",
                 edgecolor="black", alpha=0.9)
    ax_box.text(0.5, 0.95, box_text, transform=ax_box.transAxes,
                fontsize=9, va="top", ha="center", family="monospace",
                bbox=props, multialignment="left")

    fig.suptitle(f"{dataset}: Acceptance-corrected N(φ) — gen model (blue) vs fit (red)",
                 fontsize=12, fontweight="bold")

    out_path = out_dir / "acceptance_corrected.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Saved: {out_path}")

    # Save results
    result = {
        "dataset": dataset,
        "chi2_total": float(best_chi2),
        "fit_params": {n: float(best_params[n]) for n in ["alf2", "alf3", "bt", "nuT", "cR"]},
        "free_params": free_params,
        "chi2_per_var": {v: {"chi2": float(c), "npts": int(n)}
                         for v, (c, n) in chi2_per_var.items()},
        "gen_params": dict(mp),
        "acceptance": {
            var_name: {
                "bin_edges": corrected_data[var_name]["bin_edges"].tolist(),
                "h_gen": corrected_data[var_name]["h_gen"].tolist(),
                "h_acc": corrected_data[var_name]["h_acc"].tolist(),
                "acc": corrected_data[var_name]["acc"].tolist(),
                "acc_err": corrected_data[var_name]["acc_err"].tolist(),
                "n_phi": corrected_data[var_name]["n_phi"].tolist(),
                "dn_phi": corrected_data[var_name]["dn_phi"].tolist(),
                "n_corr": corrected_data[var_name]["n_corr"].tolist(),
                "dn_corr": corrected_data[var_name]["dn_corr"].tolist(),
                "h_orig_model": orig_hists[var_name].tolist(),
                "h_fit_model": fit_hists[var_name].tolist(),
            }
            for var_name in corrected_data
        },
    }
    json_path = out_dir / "acceptance_results.json"
    with open(json_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {json_path}")

    return result


# ═══════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    if len(sys.argv) > 1:
        targets = [sys.argv[1]]
    else:
        targets = DATASETS

    results = []
    for ds in targets:
        r = process_dataset(ds)
        if r:
            results.append(r)

    if results:
        print(f"\n{'='*70}")
        print(f"  SUMMARY")
        print(f"{'='*70}")
        print(f"  {'Dataset':20s} {'chi2':>8s}  {'alf2':>7s} {'alf3':>7s} {'bt':>7s} {'nuT':>7s}")
        print(f"  {'-'*60}")
        for r in results:
            p = r["fit_params"]
            print(f"  {r['dataset']:20s} {r['chi2_total']:8.1f}  "
                  f"{p['alf2']:7.3f} {p['alf3']:7.3f} {p['bt']:7.3f} {p['nuT']:7.3f}")
