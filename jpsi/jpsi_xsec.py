#!/usr/bin/env python3
"""
jpsi_xsec.py — Cross-section shapes and acceptance comparison in Mariana's bins.

1) dsigma/dt (gamma* p -> J/psi p), flux-folded, in the three correlated
   (W, Q2) bins of Mariana's Table 16. Acceptance-corrected data vs the
   tuned model. Shapes only: absolute normalization requires luminosity.
2) Cross section vs Q2 and vs W integrated over t (same treatment).
3) Our NN-fastMC mixture acceptance in Mariana's bins vs her Table 17
   efficiencies (eta1, eta2a, eta2b).

Reads: MC caches (jpsi_tuning_cache), data LUND, fit params from
jpsi_results.json (falls back to config gen values).

Usage: python3 jpsi_xsec.py [config.json]
"""

import json
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

BASE = Path(__file__).resolve().parent
CACHE_ROOT = Path.home() / "Downloads/volatile/clas12/vpk/fastmc/jpsi_tuning_cache"

MP  = 0.938272; MP2 = MP**2
MJP = 3.0969;   MJP2 = MJP**2
WTH2 = (MP + MJP)**2

# Mariana, Table 16: correlated (Q2, W) bins
MARIANA_BINS = [
    {"id": 1, "q2": (0.029, 0.071), "w": (4.202, 4.256)},
    {"id": 2, "q2": (0.018, 0.047), "w": (4.320, 4.387)},
    {"id": 3, "q2": (0.016, 0.043), "w": (4.426, 4.487)},
]
# Mariana, Table 17: efficiencies
MARIANA_EFF = [
    {"id": 1, "eta1": 0.1463, "eta2a": 0.1184, "eta2b": 0.1406},
    {"id": 2, "eta1": 0.1101, "eta2a": 0.1096, "eta2b": 0.1386},
    {"id": 3, "eta1": 0.1090, "eta2a": 0.1127, "eta2b": 0.1313},
]


def compute_tmin(q2, xb):
    w2 = MP2 + q2 * (1 - xb) / xb
    w = np.sqrt(np.maximum(w2, 0.01))
    ecm_i = (w2 + q2 + MP2) / (2 * w)
    pcm_i = np.sqrt(np.maximum(ecm_i**2 - MP2, 0))
    ecm_f = (w2 + MJP2 - MP2) / (2 * w)
    pcm_f = np.sqrt(np.maximum(ecm_f**2 - MJP2, 0))
    return (pcm_i - pcm_f)**2 - (ecm_i - ecm_f)**2   # |t|_min > 0


def sigma_model(q2, xb, t_pos, tmin_pos, p):
    w2 = MP2 + q2 * (1 - xb) / xb
    w = np.sqrt(np.maximum(w2, 0.01))
    thresh = np.maximum(1.0 - WTH2 / w2, 1e-30)**p["alf2"]
    sig_T = thresh * w**p["alf3"] * (1.0 + q2 / MJP2)**(-p["nuT"]) \
            * 3.0 * (p["mg2"] + tmin_pos)**3 / (p["mg2"] + t_pos)**4
    return sig_T * (1.0 + p["cR"] * q2 / MJP2)


def read_data(cfg):
    rows = {k: [] for k in ["ebeam", "q2", "w", "xb", "t", "mee"]}
    with open(cfg["channel"]["data_file"]) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        hdr = lines[i].split()
        npart, ebeam = int(hdr[0]), float(hdr[6])
        parts = [lines[i+1+k].split() for k in range(npart)]
        i += 1 + npart
        p4 = lambda pc: np.array([float(pc[6]), float(pc[7]), float(pc[8]), float(pc[9])])
        eft, ed, ep, pr = p4(parts[0]), p4(parts[1]), p4(parts[2]), p4(parts[3])
        pmag = np.linalg.norm(eft[:3])
        cos_th = eft[2] / pmag
        q2 = 2 * ebeam * eft[3] * (1 - cos_th)
        nu = ebeam - eft[3]
        xb = q2 / (2 * MP * nu)
        w = np.sqrt(max(MP2 + 2 * MP * nu - q2, 0))
        dE = pr[3] - MP
        t = abs(dE**2 - pr[0]**2 - pr[1]**2 - pr[2]**2)
        ee = ed + ep
        mee = np.sqrt(max(ee[3]**2 - ee[0]**2 - ee[1]**2 - ee[2]**2, 0))
        for k, v in [("ebeam", ebeam), ("q2", q2), ("w", w), ("xb", xb),
                     ("t", t), ("mee", mee)]:
            rows[k].append(v)
    return {k: np.asarray(v) for k, v in rows.items()}


def main():
    cfg_path = BASE / "config.json" if len(sys.argv) < 2 else Path(sys.argv[1])
    with open(cfg_path) as f:
        cfg = json.load(f)

    # tuned parameters: prefer fit result, fall back to config
    res_path = BASE / "jpsi_results.json"
    if res_path.exists():
        with open(res_path) as f:
            pars = json.load(f)["fit_params"]
        print(f"Using fitted params: { {k: round(v,3) for k,v in pars.items()} }")
    else:
        pars = {k: v["value"] for k, v in cfg["cross_section_model"]["parameters"].items()}
        print(f"Using config gen params")

    d = read_data(cfg)
    mee_cut = (d["mee"] >= cfg["cuts"]["mee_min"]) & (d["mee"] < cfg["cuts"]["mee_max"])
    comp = cfg["channel"]["data_composition"]

    mcs = {}
    for eb in cfg["generation"]["energies"]:
        z = np.load(CACHE_ROOT / f"mc_e{eb}.npz")
        mc = {k: z[k] for k in z.files}
        n_acc = mc["acc_weight"].sum()
        lam = comp[eb] / n_acc
        mcs[eb] = (mc, lam)

    # ════════════════════════════════════════════════════════════════
    #  1) Acceptance in Mariana's (W,Q2) bins  +  comparison table
    # ════════════════════════════════════════════════════════════════
    print("\n=== Acceptance in Mariana's bins (Table 16/17) ===")
    print(f"{'bin':>3s} {'Q2 range':>16s} {'W range':>16s} "
          f"{'eps(our)':>9s} {'eta1':>7s} {'eta2a':>7s} {'eta2b':>7s} {'N_dat':>6s}")
    acc_rows = []
    for mb, me in zip(MARIANA_BINS, MARIANA_EFF):
        g_sum = a_sum = 0.0
        for eb, (mc, lam) in mcs.items():
            inbin = ((mc["q2"] >= mb["q2"][0]) & (mc["q2"] < mb["q2"][1]) &
                     (mc["w"] >= mb["w"][0]) & (mc["w"] < mb["w"][1]))
            g_sum += lam * inbin.sum()
            a_sum += lam * mc["acc_weight"][inbin].sum()
        eps = a_sum / g_sum if g_sum > 0 else 0.0
        din = (mee_cut & (d["q2"] >= mb["q2"][0]) & (d["q2"] < mb["q2"][1]) &
               (d["w"] >= mb["w"][0]) & (d["w"] < mb["w"][1]))
        nd = int(din.sum())
        print(f"{mb['id']:>3d} [{mb['q2'][0]:.3f},{mb['q2'][1]:.3f}] "
              f"[{mb['w'][0]:.3f},{mb['w'][1]:.3f}]  {eps:8.4f} "
              f"{me['eta1']:7.4f} {me['eta2a']:7.4f} {me['eta2b']:7.4f} {nd:6d}")
        acc_rows.append({"bin": mb["id"], "q2": mb["q2"], "w": mb["w"],
                         "eps_ours": eps, **{k: me[k] for k in ["eta1","eta2a","eta2b"]},
                         "n_data": nd})

    # ════════════════════════════════════════════════════════════════
    #  2) dsigma/dt (flux-folded) in the three bins — shape, data vs model
    # ════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    dsdt_out = []
    for iax, mb in enumerate(MARIANA_BINS):
        ax = axes[iax]
        # bin-specific t range: start at the 1st percentile of generated |t|
        t_mc_all, w_mc_all = [], []
        for eb, (mc, lam) in mcs.items():
            inbin = ((mc["q2"] >= mb["q2"][0]) & (mc["q2"] < mb["q2"][1]) &
                     (mc["w"] >= mb["w"][0]) & (mc["w"] < mb["w"][1]))
            t_mc_all.append(mc["t"][inbin])
        tlow = np.floor(np.percentile(np.concatenate(t_mc_all), 1) * 10) / 10
        edges = np.linspace(tlow, 3.0, 5)
        widths = np.diff(edges)
        ctr = 0.5 * (edges[:-1] + edges[1:])

        h_gen = np.zeros(4); h_acc = np.zeros(4); h_mod = np.zeros(4)
        for eb, (mc, lam) in mcs.items():
            inbin = ((mc["q2"] >= mb["q2"][0]) & (mc["q2"] < mb["q2"][1]) &
                     (mc["w"] >= mb["w"][0]) & (mc["w"] < mb["w"][1]))
            tmin_p = compute_tmin(mc["q2"][inbin], mc["xb"][inbin])
            s_gen = sigma_model(mc["q2"][inbin], mc["xb"][inbin], mc["t"][inbin], tmin_p,
                                {k: v["value"] for k, v in cfg["cross_section_model"]["parameters"].items()})
            s_new = sigma_model(mc["q2"][inbin], mc["xb"][inbin], mc["t"][inbin], tmin_p, pars)
            rw = np.where(s_gen > 0, s_new / np.maximum(s_gen, 1e-300), 0.0)
            hg, _ = np.histogram(mc["t"][inbin], bins=edges)
            ha, _ = np.histogram(mc["t"][inbin], bins=edges, weights=mc["acc_weight"][inbin])
            hm, _ = np.histogram(mc["t"][inbin], bins=edges, weights=rw)
            h_gen += lam * hg; h_acc += lam * ha; h_mod += lam * hm
        eps = np.where(h_gen > 0, h_acc / h_gen, 0.0)

        din = (mee_cut & (d["q2"] >= mb["q2"][0]) & (d["q2"] < mb["q2"][1]) &
               (d["w"] >= mb["w"][0]) & (d["w"] < mb["w"][1]))
        h_dat, _ = np.histogram(d["t"][din], bins=edges)
        good = eps > 0
        n_corr = np.where(good, h_dat / np.where(good, eps, 1), 0.0)
        dn_corr = np.where(good, np.sqrt(np.maximum(h_dat, 1)) / np.where(good, eps, 1), 0.0)

        # dsigma/dt shape: divide by bin width, normalize data integral to 1
        norm_d = (n_corr * widths)[good].sum()
        ds_d  = np.where(good, n_corr / widths, 0) / max(norm_d, 1e-30) * norm_d  # keep counts scale
        dds_d = np.where(good, dn_corr / widths, 0)
        # model normalized to data integral over good bins
        norm_m = (h_mod * 1.0)[good].sum()
        sc = (n_corr[good].sum() / norm_m) if norm_m > 0 else 1.0
        ds_m = h_mod * sc / widths

        ax.stairs(ds_m, edges, color="red", lw=2, fill=True, alpha=0.35, label="Tuned model")
        ax.errorbar(ctr[good], ds_d[good], yerr=dds_d[good],
                    xerr=widths[good]/2, fmt="ko", ms=5, capsize=2, label=r"Data/$\varepsilon$")
        ax.set_yscale("log")
        ax.set_xlabel(r"$|t|$ [GeV$^2$]")
        ax.set_ylabel(r"$dN^{corr}/dt$ [GeV$^{-2}$] (shape)")
        ax.set_title(f"Bin {mb['id']}: $Q^2$[{mb['q2'][0]},{mb['q2'][1]}], "
                     f"$W$[{mb['w'][0]},{mb['w'][1]}]", fontsize=10)
        ax.grid(alpha=0.3)
        if iax == 0:
            ax.legend(fontsize=9)
        dsdt_out.append({"bin": mb["id"], "t_edges": edges.tolist(),
                         "eps": eps.tolist(), "n_dat": h_dat.tolist(),
                         "dsdt_data": ds_d.tolist(), "dsdt_err": dds_d.tolist(),
                         "dsdt_model": ds_m.tolist()})
    fig.suptitle(r"$d\sigma/dt$ ($\gamma^*p\to J/\psi\, p$, flux-folded) — "
                 r"Mariana bins, shape comparison", fontsize=12)
    fig.tight_layout()
    fig.savefig(BASE / "jpsi_dsdt_mariana_bins.png", dpi=130)
    plt.close(fig)
    print(f"Saved: jpsi_dsdt_mariana_bins.png")

    # ════════════════════════════════════════════════════════════════
    #  3) sigma vs Q2 and vs W integrated over t (fine model curve)
    # ════════════════════════════════════════════════════════════════
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    cuts = cfg["cuts"]
    for iax, (var, lo, hi, nb, xlabel) in enumerate([
            ("q2", cuts["q2min"], cuts["q2max"], 3, r"$Q^2$ [GeV$^2$]"),
            ("w",  cuts["wmin"],  cuts["wmax"],  3, r"$W$ [GeV]")]):
        ax = axes[iax]
        edges = np.linspace(lo, hi, nb + 1)
        widths = np.diff(edges); ctr = 0.5*(edges[:-1]+edges[1:])
        fine = np.linspace(lo, hi, 13)
        h_gen = np.zeros(nb); h_acc = np.zeros(nb)
        h_mod_f = np.zeros(12)
        for eb, (mc, lam) in mcs.items():
            sel = ((mc["q2"] >= cuts["q2min"]) & (mc["q2"] < cuts["q2max"]) &
                   (mc["w"] >= cuts["wmin"]) & (mc["w"] < cuts["wmax"]) &
                   (mc["t"] < cuts["tmax"]))
            tmin_p = compute_tmin(mc["q2"][sel], mc["xb"][sel])
            gp = {k: v["value"] for k, v in cfg["cross_section_model"]["parameters"].items()}
            s_gen = sigma_model(mc["q2"][sel], mc["xb"][sel], mc["t"][sel], tmin_p, gp)
            s_new = sigma_model(mc["q2"][sel], mc["xb"][sel], mc["t"][sel], tmin_p, pars)
            rw = np.where(s_gen > 0, s_new / np.maximum(s_gen, 1e-300), 0.0)
            hg, _ = np.histogram(mc[var][sel], bins=edges)
            ha, _ = np.histogram(mc[var][sel], bins=edges, weights=mc["acc_weight"][sel])
            hf, _ = np.histogram(mc[var][sel], bins=fine, weights=rw)
            h_gen += lam * hg; h_acc += lam * ha; h_mod_f += lam * hf
        eps = np.where(h_gen > 0, h_acc / h_gen, 0.0)
        dsel = (mee_cut & (d["q2"] >= cuts["q2min"]) & (d["q2"] < cuts["q2max"]) &
                (d["w"] >= cuts["wmin"]) & (d["w"] < cuts["wmax"]) & (d["t"] < cuts["tmax"]))
        h_dat, _ = np.histogram(d[var][dsel], bins=edges)
        good = eps > 0
        n_corr = np.where(good, h_dat/np.where(good, eps, 1), 0)
        dn_corr = np.where(good, np.sqrt(np.maximum(h_dat,1))/np.where(good, eps, 1), 0)
        # densities
        y_d = n_corr/widths; dy_d = dn_corr/widths
        fw = np.diff(fine); fc = 0.5*(fine[:-1]+fine[1:])
        y_m = h_mod_f/fw
        sc = y_d[good].mean()/ (np.interp(ctr[good], fc, y_m).mean()) if good.any() else 1
        ax.plot(fc, y_m*sc, "r-", lw=2, label="Tuned model")
        ax.errorbar(ctr[good], y_d[good], yerr=dy_d[good], xerr=widths[good]/2,
                    fmt="ko", ms=5, capsize=2, label=r"Data/$\varepsilon$")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r"$dN^{corr}/d$" + xlabel.split()[0].strip("$") + " (shape)")
        ax.grid(alpha=0.3)
        if iax == 0:
            ax.legend(fontsize=9)
    fig.suptitle(r"Flux-folded cross-section shape vs $Q^2$ and $W$ "
                 r"(integrated over $|t|<3$)", fontsize=12)
    fig.tight_layout()
    fig.savefig(BASE / "jpsi_sigma_q2_w.png", dpi=130)
    plt.close(fig)
    print(f"Saved: jpsi_sigma_q2_w.png")

    with open(BASE / "jpsi_xsec_results.json", "w") as f:
        json.dump({"params_used": pars, "acceptance_mariana_bins": acc_rows,
                   "dsdt_bins": dsdt_out}, f, indent=2)
    print("Saved: jpsi_xsec_results.json")


if __name__ == "__main__":
    main()
