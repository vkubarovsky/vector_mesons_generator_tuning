#!/usr/bin/env python3
"""
jpsi_analyze.py — Acceptance-correct the 68 J/psi candidates and fit the model.

Data: pre-selected J/psi sample (counting, no mass fit).
Acceptance: from the MC caches produced by jpsi_pipeline.py.
  - per-event acceptance weight = mean accept over fastmc configs (in cache)
  - the two beam energies are combined so that ACCEPTED MC matches the
    data composition (52 @ 10.6 : 16 @ 10.2); the same per-energy scale
    is applied to generated MC -> mixture acceptance per bin.
Fit: reweight MC events with sigma_new/sigma_gen (dipole t-form model);
  free parameters from config (default alf2, mg2), chi2 over Q2, W, |t| bins.

Usage:  python3 jpsi_analyze.py [config.json]
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

BASE = Path(__file__).resolve().parent
CACHE_ROOT = Path.home() / "Downloads/volatile/clas12/vpk/fastmc/jpsi_tuning_cache"

MP  = 0.938272; MP2 = MP**2
MJP = 3.0969;   MJP2 = MJP**2
WTH2 = (MP + MJP)**2
MIN_GEN = 100   # min generated (unweighted) events per bin


def compute_tmin(q2, xb):
    w2 = MP2 + q2 * (1 - xb) / xb
    w = np.sqrt(np.maximum(w2, 0.01))
    ecm_i = (w2 + q2 + MP2) / (2 * w)
    pcm_i = np.sqrt(np.maximum(ecm_i**2 - MP2, 0))
    ecm_f = (w2 + MJP2 - MP2) / (2 * w)
    pcm_f = np.sqrt(np.maximum(ecm_f**2 - MJP2, 0))
    return (pcm_i - pcm_f)**2 - ((q2 + MJP2) / (2 * w))**2   # |t|_min > 0


def sigma_model(q2, xb, t_pos, tmin_pos, alf2, alf3, nuT, mg2, cR):
    """dsig/dt, dipole t-form. t_pos, tmin_pos are |t|, |t_min| (positive)."""
    w2 = MP2 + q2 * (1 - xb) / xb
    w = np.sqrt(np.maximum(w2, 0.01))
    thresh = np.maximum(1.0 - WTH2 / w2, 1e-30)**alf2
    sig_T = thresh * w**alf3 * (1.0 + q2 / MJP2)**(-nuT) \
            * 3.0 * (mg2 + tmin_pos)**3 / (mg2 + t_pos)**4
    return sig_T * (1.0 + cR * q2 / MJP2)


def read_data(cfg):
    """Read the combined J/psi candidate LUND file."""
    rows = {k: [] for k in ["ebeam", "q2", "w", "xb", "t", "t_tmin", "mee", "th_ft"]}
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
        th = np.degrees(np.arccos(cos_th))
        q2 = 2 * ebeam * eft[3] * (1 - cos_th)
        nu = ebeam - eft[3]
        xb = q2 / (2 * MP * nu)
        w2 = MP2 + 2 * MP * nu - q2
        w = np.sqrt(max(w2, 0))
        dE = pr[3] - MP
        t = abs(dE**2 - pr[0]**2 - pr[1]**2 - pr[2]**2)
        ee = ed + ep
        mee = np.sqrt(max(ee[3]**2 - ee[0]**2 - ee[1]**2 - ee[2]**2, 0))
        for k, v in [("ebeam", ebeam), ("q2", q2), ("w", w), ("xb", xb),
                     ("t", t), ("mee", mee), ("th_ft", th)]:
            rows[k].append(v)
        rows["t_tmin"].append(0.0)
    d = {k: np.asarray(v) for k, v in rows.items()}
    d["t_tmin"] = d["t"] - compute_tmin(d["q2"], d["xb"])
    return d


def main():
    cfg_path = BASE / "config.json" if len(sys.argv) < 2 else Path(sys.argv[1])
    with open(cfg_path) as f:
        cfg = json.load(f)
    cuts = cfg["cuts"]
    bins_def = cfg["bins"]
    mp = {k: v["value"] for k, v in cfg["cross_section_model"]["parameters"].items()}
    fp = cfg["fit_parameters"]

    # ─── Data ───────────────────────────────────────────────────
    d = read_data(cfg)
    d_cut = ((d["q2"] >= cuts["q2min"]) & (d["q2"] < cuts["q2max"]) &
             (d["w"] >= cuts["wmin"]) & (d["w"] < cuts["wmax"]) &
             (d["t"] < cuts["tmax"]) &
             (d["mee"] >= cuts["mee_min"]) & (d["mee"] < cuts["mee_max"]))
    print(f"Data: {len(d['q2'])} candidates, {d_cut.sum()} after cuts")
    comp = cfg["channel"]["data_composition"]

    # ─── MC caches ──────────────────────────────────────────────
    mcs = {}
    for eb in cfg["generation"]["energies"]:
        z = np.load(CACHE_ROOT / f"mc_e{eb}.npz")
        mc = {k: z[k] for k in z.files}
        cut = ((mc["q2"] >= cuts["q2min"]) & (mc["q2"] < cuts["q2max"]) &
               (mc["w"] >= cuts["wmin"]) & (mc["w"] < cuts["wmax"]) &
               (mc["t"] < cuts["tmax"]))
        mc = {k: v[cut] for k, v in mc.items()}
        # per-energy scale: accepted MC -> data count at this energy
        n_acc = mc["acc_weight"].sum()
        lam = comp[eb] / n_acc if n_acc > 0 else 0.0
        mcs[eb] = (mc, lam)
        print(f"MC e{eb}: {len(mc['q2'])} gen in cuts, "
              f"acc_weight sum={n_acc:.0f}, lambda={lam:.3e}")

    # ─── Binned acceptance + corrected yields ──────────────────
    plot_vars = ["q2", "w", "t", "t_tmin", "xb"]
    fit_vars = ["q2", "w", "t"]
    var_labels = {"q2": r"$Q^2$ [GeV$^2$]", "w": r"$W$ [GeV]",
                  "t": r"$|t|$ [GeV$^2$]", "t_tmin": r"$|t|-|t|_{min}$ [GeV$^2$]",
                  "xb": r"$x_B$"}

    hist = {}
    for var in plot_vars:
        edges = np.asarray(bins_def[var], dtype=float)
        h_gen = np.zeros(len(edges) - 1)
        h_acc = np.zeros(len(edges) - 1)
        h_gen_raw = np.zeros(len(edges) - 1)
        for eb, (mc, lam) in mcs.items():
            hg, _ = np.histogram(mc[var], bins=edges)
            ha, _ = np.histogram(mc[var], bins=edges, weights=mc["acc_weight"])
            h_gen += lam * hg
            h_acc += lam * ha
            h_gen_raw += hg
        with np.errstate(divide="ignore", invalid="ignore"):
            eps = np.where((h_gen > 0) & (h_gen_raw >= MIN_GEN), h_acc / h_gen, 0.0)
        h_dat, _ = np.histogram(d[var][d_cut], bins=edges)
        good = eps > 0
        n_corr = np.where(good, h_dat / np.where(good, eps, 1), 0.0)
        dn_corr = np.where(good, np.sqrt(np.maximum(h_dat, 1)) / np.where(good, eps, 1), 0.0)
        hist[var] = {"edges": edges, "eps": eps, "n_dat": h_dat,
                     "n_corr": n_corr, "dn_corr": dn_corr, "good": good}
        print(f"\n{var}: eps = {np.array2string(eps, precision=4)}")
        print(f"  N_dat  = {h_dat}")
        print(f"  N_corr = {np.array2string(n_corr, precision=0)}")

    # ─── Model fit via reweighting ──────────────────────────────
    gen_pars = dict(mp)
    free = fp["free"]
    bounds = [tuple(fp["bounds"][k]) for k in free]
    print(f"\nFitting: free={free}, fixed at gen values: "
          f"{[k for k in ['alf2','alf3','nuT','mg2','cR'] if k not in free]}")

    # precompute per-energy arrays and gen-model sigma
    pre = []
    for eb, (mc, lam) in mcs.items():
        tminp = compute_tmin(mc["q2"], mc["xb"])
        s_gen = sigma_model(mc["q2"], mc["xb"], mc["t"], tminp,
                            mp["alf2"], mp["alf3"], mp["nuT"], mp["mg2"], mp["cR"])
        # bin index per fit var (precomputed masks)
        masks = {}
        for var in fit_vars:
            edges = hist[var]["edges"]
            idx = np.digitize(mc[var], edges) - 1
            masks[var] = [(idx == i) for i in range(len(edges) - 1)]
        pre.append((mc, lam, tminp, s_gen, masks))

    def model_histograms(p):
        out = {var: np.zeros(len(hist[var]["edges"]) - 1) for var in fit_vars}
        for mc, lam, tminp, s_gen, masks in pre:
            s_new = sigma_model(mc["q2"], mc["xb"], mc["t"], tminp,
                                p["alf2"], p["alf3"], p["nuT"], p["mg2"], p["cR"])
            wgt = np.where(s_gen > 0, s_new / np.maximum(s_gen, 1e-300), 0.0)
            for var in fit_vars:
                for i, m in enumerate(masks[var]):
                    out[var][i] += lam * wgt[m].sum()
        return out

    def chi2_fn(x):
        p = dict(gen_pars)
        for k, v in zip(free, x):
            p[k] = v
        hm = model_histograms(p)
        # analytic normalization over all fit bins
        num = den = 0.0
        for var in fit_vars:
            g = hist[var]["good"]
            nc, dn = hist[var]["n_corr"][g], hist[var]["dn_corr"][g]
            m = hm[var][g]
            num += np.sum(nc * m / dn**2)
            den += np.sum(m**2 / dn**2)
        C = num / den if den > 0 else 0.0
        c2 = 0.0
        for var in fit_vars:
            g = hist[var]["good"]
            nc, dn = hist[var]["n_corr"][g], hist[var]["dn_corr"][g]
            c2 += np.sum((nc - C * hm[var][g])**2 / dn**2)
        return c2

    t0 = time.time()
    res = differential_evolution(chi2_fn, bounds, seed=42, tol=1e-8,
                                 maxiter=300, polish=True)
    print(f"Fit done ({time.time()-t0:.1f}s): chi2={res.fun:.2f}")
    best = dict(gen_pars)
    for k, v in zip(free, res.x):
        best[k] = v
    npts = sum(int(hist[v]["good"].sum()) for v in fit_vars)
    for k in ["alf2", "alf3", "nuT", "mg2", "cR"]:
        tag = "free" if k in free else "fixed"
        print(f"  {k:5s} = {best[k]:9.4f}  ({tag})")
    print(f"  ndf = {npts} pts - {len(free)} par = {npts - len(free)}")

    # ─── Plots: corrected data vs gen and fit models ───────────
    def model_hist_all(p):
        out = {}
        for var in plot_vars:
            edges = hist[var]["edges"]
            h = np.zeros(len(edges) - 1)
            for mc, lam, tminp, s_gen, _ in pre:
                s_new = sigma_model(mc["q2"], mc["xb"], mc["t"], tminp,
                                    p["alf2"], p["alf3"], p["nuT"], p["mg2"], p["cR"])
                wgt = np.where(s_gen > 0, s_new / np.maximum(s_gen, 1e-300), 0.0)
                hh, _ = np.histogram(mc[var], bins=edges, weights=lam * wgt)
                h += hh
            out[var] = h
        return out

    hm_gen = model_hist_all(gen_pars)
    hm_fit = model_hist_all(best)

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for iax, var in enumerate(plot_vars):
        ax = axes.flat[iax]
        h = hist[var]
        edges, good = h["edges"], h["good"]
        ctr = 0.5 * (edges[:-1] + edges[1:])
        wid = np.diff(edges)
        # normalize models to corrected data in good bins (per plot)
        for hm, color, label, style in [(hm_gen, "blue", "Gen model", "step"),
                                        (hm_fit, "red", "Fit model", "stepfilled")]:
            m = hm[var]
            sc = (h["n_corr"][good].sum() / m[good].sum()) if m[good].sum() > 0 else 1
            ax.stairs(m * sc, edges, color=color, lw=2,
                      linestyle="--" if color == "blue" else "-",
                      fill=(color == "red"), alpha=0.35 if color == "red" else 1,
                      label=label)
        ax.errorbar(ctr[good], h["n_corr"][good], yerr=h["dn_corr"][good],
                    xerr=wid[good]/2, fmt="ko", ms=5, capsize=2, label=r"Data/$\varepsilon$")
        ax.set_xlabel(var_labels[var]); ax.set_ylabel(r"$N(J/\psi)/\varepsilon$")
        ax.grid(alpha=0.3)
        if iax == 0:
            ax.legend(fontsize=9)
        if var in fit_vars:
            ax.set_title(f"{var} (fitted)", fontsize=10)
        else:
            ax.set_title(f"{var} (not fitted)", fontsize=10)
    # summary box
    ax = axes.flat[5]; ax.axis("off")
    lines = [r"J/psi tuning", "",
             f"{'Param':6s}{'Gen':>9s}{'Fit':>10s}"]
    for k in ["alf2", "alf3", "nuT", "mg2", "cR"]:
        tag = "free" if k in free else "fix"
        lines.append(f"{k:6s}{gen_pars[k]:9.3f}{best[k]:10.3f}  ({tag})")
    lines += ["", f"chi2 = {res.fun:.1f} / ndf {npts - len(free)}",
              f"N data (cuts) = {int(d_cut.sum())}"]
    ax.text(0.05, 0.95, "\n".join(lines), family="monospace", fontsize=11,
            va="top", transform=ax.transAxes,
            bbox=dict(boxstyle="round", fc="lightyellow"))
    fig.suptitle(r"J/$\psi$: acceptance-corrected counts — gen model (blue) vs fit (red)",
                 fontsize=13)
    fig.tight_layout()
    out_png = BASE / "jpsi_corrected.png"
    fig.savefig(out_png, dpi=130)
    print(f"Saved: {out_png}")

    # ─── Save results JSON ──────────────────────────────────────
    out = {
        "channel": "jpsi",
        "n_data": int(d_cut.sum()),
        "chi2": float(res.fun),
        "ndf": int(npts - len(free)),
        "free_params": free,
        "gen_params": gen_pars,
        "fit_params": best,
        "histograms": {
            var: {"edges": hist[var]["edges"].tolist(),
                  "eps": hist[var]["eps"].tolist(),
                  "n_dat": hist[var]["n_dat"].tolist(),
                  "n_corr": hist[var]["n_corr"].tolist(),
                  "dn_corr": hist[var]["dn_corr"].tolist(),
                  "model_gen": hm_gen[var].tolist(),
                  "model_fit": hm_fit[var].tolist()}
            for var in plot_vars
        },
    }
    out_json = BASE / "jpsi_results.json"
    with open(out_json, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved: {out_json}")


if __name__ == "__main__":
    main()
