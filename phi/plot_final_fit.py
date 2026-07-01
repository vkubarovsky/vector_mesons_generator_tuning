#!/usr/bin/env python3
"""Fit-quality plots from the committed acceptance_results.json
(Roberts Eq.23b phi model: alf2=2, alf3=0.20 fixed; bt,nuT fit).
Data (acceptance-corrected n_corr) vs best-fit model, per dataset & variable."""
import json, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

DS = ["rgafall18_inb", "rgafall18_outb", "rgasp18_inb", "rgasp18_outb", "rgasp19_inb"]
VARS = [("w", "$W$ [GeV]"), ("q2", "$Q^2$ [GeV$^2$]"),
        ("t", "$|t|$ [GeV$^2$]"), ("xb", "$x_B$")]

fig, axes = plt.subplots(len(DS), len(VARS), figsize=(15, 15))
for i, d in enumerate(DS):
    r = json.load(open(f"runs/{d}/acceptance_results.json"))
    fp = r["fit_params"]
    for j, (v, xl) in enumerate(VARS):
        a = r["acceptance"][v]
        edges = np.array(a["bin_edges"]); c = 0.5 * (edges[:-1] + edges[1:])
        data = np.array(a["n_corr"]); derr = np.array(a["dn_corr"])
        model = np.array(a["h_fit_model"], dtype=float)
        # normalise both to unit area (shapes)
        data_n = data / data.sum(); derr_n = derr / data.sum()
        model_n = model / model.sum()
        ax = axes[i, j]
        ax.errorbar(c, data_n, yerr=derr_n, fmt="ko", ms=4, label="data (acc-corr)")
        ax.step(c, model_n, where="mid", color="tab:red", lw=1.8, label="fit model")
        ax.set_xlabel(xl, fontsize=9)
        if j == 0:
            ax.set_ylabel(d.replace("rga", "").replace("_", " "), fontsize=9)
        if i == 0 and j == 0:
            ax.legend(fontsize=8)
        chi2v = r["chi2_per_var"].get(v, {}).get("chi2")
        ttl = f"$\\chi^2$={chi2v:.1f}" if chi2v is not None else ""
        ax.set_title(ttl, fontsize=8)
        ax.grid(alpha=.3)
fig.suptitle("phi fit (Roberts Eq.23b: alf2=2, alf3=0.20 fixed; bt,nuT free) "
             "-- data vs best-fit model, 5 datasets x 4 variables", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.98])
fig.savefig("phi_final_fit_all.png", dpi=110)
print("wrote phi_final_fit_all.png")
