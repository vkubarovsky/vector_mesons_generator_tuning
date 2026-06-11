#!/usr/bin/env python3
"""
fit_mkk_binned.py — Fit M(KK) in bins of Q², |t|, xB, W for each dataset.
Produces per-variable pages: 5 M(KK) fit panels + N(φ) vs variable.
Saves results to JSON for future generator tuning.
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path

# Constants
MK = 0.493677
THRESH = 2 * MK
MPHI = 1.0195   # fixed phi mass (PDG)
MPHI_REAL = 1.019412  # for tmin calculation
MPHI2 = MPHI_REAL**2
MP = 0.938272
MP2 = MP**2
WTH2 = (MP + MPHI_REAL)**2

BASE = Path("/Users/vpk/vector_mesons_generator_tuning/phi")
LUND_DIR = Path.home() / "Downloads/bhawani_phi_data"

DATASETS = ["rgafall18_inb", "rgafall18_outb", "rgasp18_inb", "rgasp18_outb", "rgasp19_inb"]

# M(KK) histogram settings: 2.5 MeV bins, 1019.5 at bin center
XMIN, XMAX, NBINS = 0.98325, 1.10075, 47

def compute_tmin(q2, xb):
    w2 = MP2 + q2 * (1 - xb) / xb
    w = np.sqrt(np.maximum(w2, 0.01))
    ecm_i = (w2 + q2 + MP2) / (2 * w)
    pcm_i = np.sqrt(np.maximum(ecm_i**2 - MP2, 0))
    ecm_f = (w2 + MPHI2 - MP2) / (2 * w)
    pcm_f = np.sqrt(np.maximum(ecm_f**2 - MPHI2, 0))
    return (ecm_i - ecm_f)**2 - (pcm_i - pcm_f)**2

# ─── Fit functions ───────────────────────────────────────────────
SIGMA_FIXED = 0.0047  # 4.7 MeV fixed resolution

def signal(m, A):
    return A * np.exp(-0.5 * ((m - MPHI) / SIGMA_FIXED)**2)

def background(m, B, n, b):
    dm = m - THRESH
    dm = np.where(dm > 0, dm, 0)
    return B * dm**n * np.exp(-b * m)

def total(m, A, B, n, b):
    return signal(m, A) + background(m, B, n, b)

# ─── Read LUND with full kinematics ─────────────────────────────
def read_lund_kinematics(lund_file, ebeam):
    """Read Bhawani 4-particle LUND, compute mkk, Q2, W, t, xB."""
    data = {"mkk": [], "q2": [], "w": [], "t": [], "xb": []}
    beam_e = np.array([0, 0, ebeam, ebeam])  # (px, py, pz, E)
    target = np.array([0, 0, 0, MP])

    with open(lund_file) as f:
        while True:
            header = f.readline()
            if not header:
                break
            parts = header.split()
            npart = int(parts[0])
            particles = {}
            for _ in range(npart):
                line = f.readline().split()
                pid = int(line[3])
                p4 = np.array([float(line[6]), float(line[7]),
                               float(line[8]), float(line[9])])
                if pid == 11:    particles["e"] = p4
                elif pid == 2212: particles["p"] = p4
                elif pid == 321:  particles["kp"] = p4
                elif pid == -321: particles["km"] = p4

            if len(particles) < 4:
                continue

            e, pro, kp, km = particles["e"], particles["p"], particles["kp"], particles["km"]

            # M(KK)
            pkk = kp + km
            mkk2 = pkk[3]**2 - pkk[0]**2 - pkk[1]**2 - pkk[2]**2
            mkk = np.sqrt(max(mkk2, 0))

            # q = beam - scattered electron
            q = beam_e - e
            Q2 = -(q[3]**2 - q[0]**2 - q[1]**2 - q[2]**2)
            nu = ebeam - e[3]

            # W
            w_vec = q + target
            W2 = w_vec[3]**2 - w_vec[0]**2 - w_vec[1]**2 - w_vec[2]**2
            W = np.sqrt(max(W2, 0))

            # t = (target - recoil)^2
            dt = target - pro
            t_val = dt[3]**2 - dt[0]**2 - dt[1]**2 - dt[2]**2
            t_abs = abs(t_val)

            # xB
            xB = Q2 / (2 * MP * nu) if nu > 0 else 0

            data["mkk"].append(mkk)
            data["q2"].append(Q2)
            data["w"].append(W)
            data["t"].append(t_abs)
            data["xb"].append(xB)

    result = {k: np.array(v) for k, v in data.items()}
    result["tmin"] = compute_tmin(result["q2"], result["xb"])
    result["t_tmin"] = result["t"] - result["tmin"]
    return result


def fit_mkk_bin(mkk_arr):
    """Fit M(KK) distribution, return fit results dict or None if too few events."""
    if len(mkk_arr) < 20:
        return None

    counts, edges = np.histogram(mkk_arr, bins=NBINS, range=(XMIN, XMAX))
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = edges[1] - edges[0]

    mask = counts > 0
    if mask.sum() < 10:
        return None

    x_fit = centers[mask]
    y_fit = counts[mask]
    y_err = np.sqrt(np.maximum(counts[mask], 1))

    try:
        p0 = [max(counts) * 0.5, 1e5, 0.5, 10]
        bounds_lo = [0, 0, 0.01, 0]
        bounds_hi = [1e5, 1e12, 5.0, 100]

        popt, pcov = curve_fit(total, x_fit, y_fit, p0=p0, sigma=y_err,
                               bounds=(bounds_lo, bounds_hi), maxfev=50000)
        perr = np.sqrt(np.diag(pcov))

        A, B, n, b = popt
        dA = perr[0]
        sig = SIGMA_FIXED

        N_phi = A * sig * np.sqrt(2 * np.pi) / bin_width
        dN_phi = (dA / A) * N_phi if A > 0 else 0

        y_model = total(x_fit, *popt)
        chi2 = np.sum(((y_fit - y_model) / y_err)**2)
        ndf = len(x_fit) - len(popt)

        return {
            "N_total": len(mkk_arr),
            "N_phi": float(N_phi), "dN_phi": float(dN_phi),
            "sigma": float(sig * 1000), "dsigma": 0.0,  # fixed
            "chi2": float(chi2), "ndf": int(ndf),
            "popt": popt.tolist(),
            "counts": counts, "centers": centers, "bin_width": bin_width,
        }
    except Exception as e:
        print(f"    Fit failed: {e}")
        return None


def plot_variable_page(dataset, var_name, var_label, var_unit, bin_edges, fit_results, out_dir):
    """Make page: 3×2 grid — 5 M(KK) fit panels + bottom-right N(φ) vs variable."""
    n_bins = len(bin_edges) - 1

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.subplots_adjust(hspace=0.35, wspace=0.30)

    x_smooth = np.linspace(XMIN, XMAX, 300)
    bin_centers_var = []
    n_phi_vals = []
    n_phi_errs = []

    # Map 5 bins to grid positions: (0,0),(0,1),(0,2),(1,0),(1,1)
    positions = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)]

    for i in range(n_bins):
        lo, hi = bin_edges[i], bin_edges[i + 1]
        bc = 0.5 * (lo + hi)
        bin_centers_var.append(bc)

        row, col = positions[i]
        ax = axes[row, col]
        r = fit_results[i]

        if r is not None:
            counts = r["counts"]
            centers = r["centers"]
            bw = r["bin_width"]
            popt = r["popt"]
            A, B, n, b = popt

            ax.bar(centers, counts, width=bw, color="lightskyblue",
                   edgecolor="none", label="Data")
            ax.errorbar(centers, counts, yerr=np.sqrt(np.maximum(counts, 1)),
                        fmt="ko", ms=2, capsize=0, lw=0.8, zorder=5)
            ax.fill_between(x_smooth, 0, background(x_smooth, B, n, b),
                            color="lightcoral", alpha=0.5)
            ax.plot(x_smooth, total(x_smooth, *popt), "r-", lw=1.5)

            n_phi_vals.append(r["N_phi"])
            n_phi_errs.append(r["dN_phi"])

            # Mini stat box
            stat = f"N(φ)={r['N_phi']:.0f}±{r['dN_phi']:.0f}\nσ={r['sigma']:.1f} MeV (fix)"
            ax.text(0.97, 0.95, stat, transform=ax.transAxes, fontsize=8,
                    va="top", ha="right", family="monospace",
                    bbox=dict(boxstyle="round", fc="lightyellow", ec="black", alpha=0.8))
        else:
            ax.text(0.5, 0.5, "Too few\nevents", transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, color="gray")
            n_phi_vals.append(0)
            n_phi_errs.append(0)

        ax.set_title(f"{lo:.1f} < {var_label} < {hi:.1f} {var_unit}", fontsize=10)
        ax.set_xlim(XMIN, XMAX)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("M(K⁺K⁻) [GeV]", fontsize=8)
        ax.set_ylabel("Events / 2.5 MeV", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.4, lw=0.5)

    # Bottom-right panel (1,2): N(phi) vs variable
    ax_bot = axes[1, 2]
    bc = np.array(bin_centers_var)
    nph = np.array(n_phi_vals)
    dnph = np.array(n_phi_errs)
    bw_var = np.array([bin_edges[i+1] - bin_edges[i] for i in range(n_bins)])

    valid = nph > 0
    if valid.any():
        ax_bot.errorbar(bc[valid], nph[valid], yerr=dnph[valid], xerr=bw_var[valid]/2,
                        fmt="ko", ms=5, capsize=3, lw=1.5)
    ax_bot.set_xlabel(f"{var_label} {var_unit}", fontsize=12)
    ax_bot.set_ylabel("N(φ)", fontsize=12)
    ax_bot.set_xlim(bin_edges[0], bin_edges[-1])
    ax_bot.set_ylim(bottom=0)
    ax_bot.grid(True, alpha=0.3)

    fig.suptitle(f"{dataset}: M(K⁺K⁻) fits in {var_label} bins", fontsize=14)

    out_path = out_dir / f"fit_mkk_vs_{var_name}.png"
    fig.savefig(str(out_path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(out_path)


# ─── Main ────────────────────────────────────────────────────────
VAR_INFO = {
    "q2":    ("Q²", "[GeV²]"),
    "t":     ("|t|", "[GeV²]"),
    "t_tmin": ("|t|-t_min", "[GeV²]"),
    "xb":    ("x_B", ""),
    "w":     ("W", "[GeV]"),
}

all_results = {}

for dataset in DATASETS:
    print(f"\n{'='*60}")
    print(f"  {dataset}")
    print(f"{'='*60}")

    cfg_path = BASE / "runs" / dataset / "config.json"
    with open(cfg_path) as f:
        cfg = json.load(f)

    ebeam = cfg["dataset"]["ebeam"]
    bins_def = cfg["mkk_bins"]
    out_dir = BASE / "runs" / dataset

    # Read data
    lund_file = LUND_DIR / cfg["dataset"]["data_file"]
    print(f"  Reading {lund_file}...")
    data = read_lund_kinematics(lund_file, ebeam)
    print(f"  Events: {len(data['mkk'])}")

    ds_results = {}

    for var_name, (var_label, var_unit) in VAR_INFO.items():
        bin_edges = bins_def[var_name]
        n_bins = len(bin_edges) - 1
        print(f"\n  {var_label} bins: {bin_edges}")

        fit_results = []
        bin_results = []

        for i in range(n_bins):
            lo, hi = bin_edges[i], bin_edges[i + 1]
            mask = (data[var_name] >= lo) & (data[var_name] < hi)
            mkk_bin = data["mkk"][mask]
            print(f"    [{lo:.1f}, {hi:.1f}): {len(mkk_bin)} events", end="")

            r = fit_mkk_bin(mkk_bin)
            fit_results.append(r)

            if r is not None:
                print(f"  → N(φ) = {r['N_phi']:.0f} ± {r['dN_phi']:.0f}")
                bin_results.append({
                    "bin_lo": lo, "bin_hi": hi,
                    "N_total": r["N_total"],
                    "N_phi": r["N_phi"], "dN_phi": r["dN_phi"],
                    "sigma_MeV": r["sigma"], "dsigma_MeV": r["dsigma"],
                    "chi2": r["chi2"], "ndf": r["ndf"],
                })
            else:
                print("  → skip (too few)")
                bin_results.append({
                    "bin_lo": lo, "bin_hi": hi,
                    "N_total": int(mask.sum()), "N_phi": 0, "dN_phi": 0,
                })

        # Plot
        plot_path = plot_variable_page(dataset, var_name, var_label, var_unit,
                                       bin_edges, fit_results, out_dir)
        print(f"  Saved: {plot_path}")

        ds_results[var_name] = {
            "bin_edges": bin_edges,
            "label": var_label,
            "unit": var_unit,
            "bins": bin_results,
        }

    # Save JSON
    json_path = out_dir / "mkk_binned_results.json"
    with open(json_path, "w") as f:
        json.dump(ds_results, f, indent=2)
    print(f"\n  Saved: {json_path}")

    all_results[dataset] = ds_results

print("\n" + "=" * 60)
print("  DONE — all datasets")
print("=" * 60)
