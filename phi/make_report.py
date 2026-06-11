#!/usr/bin/env python3
"""
make_report.py — Generate the LaTeX report for the phi tuning chain.
All generation parameters are read dynamically from the run configs and
acceptance_results.json (gen_params), so the report always matches the
actual iteration. Compiles the PDF if pdflatex is available.
"""

import json
import shutil
import subprocess
import numpy as np
from pathlib import Path

BASE = Path("/Users/vpk/vector_mesons_generator_tuning/phi")
TEX_DIR = BASE / "tex"
DATASETS = ["rgafall18_inb", "rgafall18_outb", "rgasp18_inb", "rgasp18_outb", "rgasp19_inb"]

DS_LABELS = {
    "rgafall18_inb":  "RGA Fall 2018, inbending",
    "rgafall18_outb": "RGA Fall 2018, outbending",
    "rgasp18_inb":    "RGA Spring 2018, inbending",
    "rgasp18_outb":   "RGA Spring 2018, outbending",
    "rgasp19_inb":    "RGA Spring 2019, inbending",
}
DS_SHORT = {
    "rgafall18_inb":  "Fall18 inb",
    "rgafall18_outb": "Fall18 outb",
    "rgasp18_inb":    "Sp18 inb",
    "rgasp18_outb":   "Sp18 outb",
    "rgasp19_inb":    "Sp19 inb",
}


def main():
    TEX_DIR.mkdir(parents=True, exist_ok=True)
    plots_dir = TEX_DIR / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    # ----- Load results -----
    all_acc, all_cfg = {}, {}
    for ds in DATASETS:
        p = BASE / "runs" / ds / "acceptance_results.json"
        if p.exists():
            with open(p) as f:
                all_acc[ds] = json.load(f)
        with open(BASE / "runs" / ds / "config.json") as f:
            all_cfg[ds] = json.load(f)

    # Generation parameters: taken from acceptance_results gen_params
    # (identical across datasets by construction)
    first = next(iter(all_acc.values()))
    gp = first["gen_params"]   # alf2, alf3, bt, nuT, cR

    # nev / Q2 ranges per dataset for the generator section
    gen_rows = []
    for ds in DATASETS:
        g = all_cfg[ds]["generation"]
        d = all_cfg[ds]["dataset"]
        gen_rows.append((DS_SHORT[ds], d["ebeam"], g["nev"],
                         g["gen_q2min"], g["gen_q2max"]))

    # ----- Copy plots -----
    plot_files = {}
    for ds in DATASETS:
        run_dir = BASE / "runs" / ds
        ds_plots = {}
        for png in sorted(run_dir.glob("*.png")):
            dest = plots_dir / f"{ds}_{png.name}"
            shutil.copy2(png, dest)
            ds_plots[png.stem] = f"plots/{ds}_{png.name}"
        plot_files[ds] = ds_plots

    # ----- Fitted parameter table data -----
    fit_vals = {k: [] for k in ["alf2", "alf3", "bt", "nuT"]}
    fit_rows = []
    for ds in DATASETS:
        if ds not in all_acc:
            continue
        r = all_acc[ds]
        fp = r["fit_params"]
        ndf = sum(v["npts"] for k, v in r["chi2_per_var"].items()
                  if k in ["q2", "t", "w"]) - 4
        fit_rows.append((DS_SHORT[ds], r["chi2_total"], ndf,
                         fp["alf2"], fp["alf3"], fp["bt"], fp["nuT"]))
        for k in fit_vals:
            fit_vals[k].append(fp[k])
    avg = {k: np.mean(v) for k, v in fit_vals.items()}
    std = {k: np.std(v) for k, v in fit_vals.items()}

    def subavg(dss):
        vals = {k: [] for k in fit_vals}
        for ds in dss:
            if ds in all_acc:
                fp = all_acc[ds]["fit_params"]
                for k in vals:
                    vals[k].append(fp[k])
        return ({k: np.mean(v) for k, v in vals.items()},
                {k: np.std(v) for k, v in vals.items()})

    inb_avg, inb_std = subavg([d for d in DATASETS if not d.endswith("outb")])
    outb_avg, outb_std = subavg([d for d in DATASETS if d.endswith("outb")])

    gen_line = (f"$\\alpha_2 = {gp['alf2']:.3f}$, $\\alpha_3 = {gp['alf3']:.3f}$, "
                f"$b_t = {gp['bt']:.3f}$, $\\nu_T = {gp['nuT']:.3f}$")

    # ================= LaTeX =================
    tex = r"""\documentclass[12pt,a4paper]{article}
\usepackage[margin=2.5cm]{geometry}
\usepackage{graphicx}
\usepackage{amsmath,amssymb}
\usepackage{booktabs}
\usepackage{caption}
\usepackage{hyperref}
\usepackage{float}

\title{Diffractive $\phi$ Meson Electroproduction at CLAS12:\\
Monte Carlo Tuning with Neural Network Fast Simulation}
\author{V.~Kubarovsky}
\date{\today}

\begin{document}
\maketitle

\begin{abstract}
We present results of Monte Carlo generator tuning for diffractive $\phi$ meson
electroproduction at CLAS12. The cross-section model parameters are determined
by fitting acceptance-corrected experimental yields in bins of $Q^2$, $|t|$, and $W$.
Five CLAS12 datasets (RGA Fall 2018 inbending/outbending, RGA Spring 2018
inbending/outbending, RGA Spring 2019 inbending) are analyzed using a neural
network-based fast Monte Carlo simulation for detector acceptance.
\end{abstract}

\tableofcontents
\newpage

%% ==================================================================
\section{Cross-Section Model}
\label{sec:model}

The differential cross section for exclusive $\phi$ meson electroproduction
$ep \to ep\phi$, $\phi \to K^+K^-$ is parametrized as:
\begin{equation}
\frac{d\sigma}{dt} = \frac{d\sigma_T}{dt} + \varepsilon\,\frac{d\sigma_L}{dt}
\end{equation}
where $\varepsilon$ is the virtual photon polarization parameter.

\subsection{Transverse Cross Section}
The transverse cross section is:
\begin{equation}
\frac{d\sigma_T}{dt} = \sigma_T(W,Q^2) \cdot b_t \cdot e^{b_t \cdot t}
\end{equation}
with
\begin{equation}
\sigma_T(W,Q^2) = \alpha_1 \left(1 - \frac{W_{\rm th}^2}{W^2}\right)^{\alpha_2}
W^{\alpha_3} \left(1 + \frac{Q^2}{m_\phi^2}\right)^{-\nu_T}
\end{equation}
where $W_{\rm th} = M_N + m_\phi = 1.96$~GeV is the production threshold,
$m_\phi = 1.019412$~GeV, and $M_N = 0.938272$~GeV.

\subsection{Longitudinal Cross Section}
The longitudinal cross section is related to the transverse via:
\begin{equation}
\frac{d\sigma_L}{dt} = c_R \cdot \frac{Q^2}{m_\phi^2} \cdot \frac{d\sigma_T}{dt}
\end{equation}

\subsection{Model Parameters}
The free parameters of the model are:
\begin{itemize}
\item $\alpha_1$ --- overall normalization [nb]
\item $\alpha_2$ --- threshold exponent: $(1 - W_{\rm th}^2/W^2)^{\alpha_2}$
\item $\alpha_3$ --- $W$ power law: $W^{\alpha_3}$
\item $b_t$ --- $t$-slope [GeV$^{-2}$]: $b_t \cdot e^{b_t \cdot t}$
\item $\nu_T$ --- $Q^2$ suppression power: $(1 + Q^2/m_\phi^2)^{-\nu_T}$
\item $c_R$ --- $\sigma_L/\sigma_T$ ratio coefficient (fixed to 1.0)
\end{itemize}

%% ==================================================================
\section{Monte Carlo Generator}
\label{sec:generator}

The Monte Carlo generator \texttt{diffrad\_gen\_exp} is based on the DIFFRAD code
by I.~Akushevich (1998) for diffractive vector meson electroproduction with QED
radiative corrections in the collinear approximation.

Events are generated in the kinematic space $(Q^2, y, t, \phi_h)$ using an
accept-reject method. The proposal distribution uses an exponential $t$-slope
for efficient sampling. The cross-section model parameters ($\alpha_1$, $\alpha_2$,
$\alpha_3$, $\nu_T$, $b_t$, $c_R$) are read from the input file, allowing
iteration without recompilation.

The generation parameters used for this iteration are:
\begin{center}
\begin{tabular}{lcl}
\toprule
Parameter & Value & Description \\
\midrule
"""
    tex += f"$\\alpha_2$ & ${gp['alf2']:.3f}$ & Threshold exponent \\\\\n"
    tex += f"$\\alpha_3$ & ${gp['alf3']:.3f}$ & $W$ power law \\\\\n"
    tex += f"$b_t$      & ${gp['bt']:.3f}$ & $t$-slope [GeV$^{{-2}}$] \\\\\n"
    tex += f"$\\nu_T$    & ${gp['nuT']:.3f}$ & $Q^2$ suppression power \\\\\n"
    tex += f"$c_R$      & ${gp['cR']:.3f}$ & $\\sigma_L/\\sigma_T$ ratio (fixed) \\\\\n"
    tex += r"""\bottomrule
\end{tabular}
\end{center}

The event statistics per dataset:
\begin{center}
\begin{tabular}{lccc}
\toprule
Dataset & $E_{\rm beam}$ [GeV] & $N_{\rm gen}$ & $Q^2$ range [GeV$^2$] \\
\midrule
"""
    for short, ebeam, nev, q2lo, q2hi in gen_rows:
        tex += f"{short} & {ebeam} & {nev//1000000}M & [{q2lo}, {q2hi}] \\\\\n"
    tex += r"""\bottomrule
\end{tabular}
\end{center}

The output is in LUND format with 7 particles per event:
beam $e^-$, target $p$, scattered $e^-$, $K^+$, $K^-$, recoil $p$, and
(for radiative events) a hard photon.

%% ==================================================================
\section{Neural Network Fast Monte Carlo}
\label{sec:fastmc}

Full GEANT4-based CLAS12 simulation and reconstruction is computationally
expensive. We use a neural network (NN) based fast Monte Carlo that provides:
\begin{itemize}
\item Acceptance probability: whether a particle with given $(p, \theta, \phi, v_z)$
      would be detected
\item Momentum smearing: realistic resolution effects
\end{itemize}

Separate NN models are trained for each particle species ($e^-$, $p$, $K^+$, $K^-$)
and each dataset configuration (beam energy, torus polarity). An event is
accepted if all four final-state particles ($e^-$, $p$, $K^+$, $K^-$) pass
their respective acceptance models.

Typical overall acceptance (all 4 particles detected) is 5--10\%, dominated
by the $K^-$ acceptance at backward angles.

"""

    for ds in DATASETS:
        lab = DS_LABELS[ds]
        if "step2_fastmc" in plot_files.get(ds, {}):
            tex += (f"\\begin{{figure}}[H]\n\\centering\n"
                    f"\\includegraphics[width=0.92\\textwidth]{{{plot_files[ds]['step2_fastmc']}}}\n"
                    f"\\caption{{Generated vs.\\ accepted distributions after fast MC for {lab}.}}\n"
                    f"\\end{{figure}}\n\n")

    tex += r"""
\clearpage
%% ==================================================================
\section{$M(K^+K^-)$ Signal Extraction}
\label{sec:mkk}

The $\phi$ meson yield $N(\phi)$ is extracted in each kinematic bin by fitting
the $K^+K^-$ invariant mass distribution $M(K^+K^-)$ with:
\begin{equation}
f(M) = A \cdot G(M;\, \mu_\phi, \sigma) \;+\; B\,(M - M_{\rm th})^n\, e^{-b\,M}
\end{equation}
where:
\begin{itemize}
\item $G(M;\, \mu_\phi, \sigma)$ is a Gaussian with fixed $\mu_\phi = 1019.5$~MeV
      and fixed $\sigma = 4.7$~MeV
\item $M_{\rm th} = 2 m_K = 987.354$~MeV is the $K^+K^-$ threshold
\item $A$, $B$, $n$, $b$ are free fit parameters
\end{itemize}

The signal yield is $N(\phi) = A \cdot \sigma \sqrt{2\pi} / \Delta M_{\rm bin}$.
Binning is performed in $Q^2$, $|t|$, $x_B$, and $W$ with 5 bins per variable.

\subsection{$M(K^+K^-)$ Fit Results}

"""

    var_labels = {"q2": "$Q^2$", "t": "$|t|$", "xb": "$x_B$", "w": "$W$"}
    for ds in DATASETS:
        lab = DS_LABELS[ds]
        tex += f"\\subsubsection*{{{lab}}}\n\n"
        for var in ["q2", "t", "xb", "w"]:
            key = f"fit_mkk_vs_{var}"
            if key in plot_files.get(ds, {}):
                tex += (f"\\begin{{figure}}[H]\n\\centering\n"
                        f"\\includegraphics[width=0.95\\textwidth]{{{plot_files[ds][key]}}}\n"
                        f"\\caption{{$M(K^+K^-)$ fits in bins of {var_labels[var]} for {lab}.}}\n"
                        f"\\end{{figure}}\n\n")
        tex += "\\clearpage\n\n"

    tex += r"""
%% ==================================================================
\section{Acceptance Correction}
\label{sec:acceptance}

The acceptance $\varepsilon_i$ in each kinematic bin $i$ is computed from the
Monte Carlo as:
\begin{equation}
\varepsilon_i = \frac{N_i^{\rm MC,\,accepted}}{N_i^{\rm MC,\,generated}}
\end{equation}
where $N_i^{\rm MC,\,generated}$ is the number of generated events in bin $i$
after kinematic cuts, and $N_i^{\rm MC,\,accepted}$ is the subset that passes
the NN fast MC acceptance.

The acceptance-corrected yield is:
\begin{equation}
N_i^{\rm corr} = \frac{N_i(\phi)}{\varepsilon_i}
\end{equation}
where $N_i(\phi)$ is the $\phi$ signal yield from the $M(K^+K^-)$ fit.

\subsection{Acceptance-Corrected Distributions with Model Fits}

"""

    for ds in DATASETS:
        lab = DS_LABELS[ds]
        if "acceptance_corrected" in plot_files.get(ds, {}):
            tex += (f"\\begin{{figure}}[H]\n\\centering\n"
                    f"\\includegraphics[width=0.95\\textwidth]{{{plot_files[ds]['acceptance_corrected']}}}\n"
                    f"\\caption{{Acceptance-corrected distributions for {lab}. "
                    f"Blue dashed histogram: generation model; red filled histogram: fitted model.}}\n"
                    f"\\end{{figure}}\n\n")

    tex += r"""
\clearpage
%% ==================================================================
\section{Model Fit Results}
\label{sec:fit}

The cross-section model parameters are determined by minimizing:
\begin{equation}
\chi^2 = \sum_{\rm var} \sum_{i} \frac{\left(N_i^{\rm corr} - N_i^{\rm model}\right)^2}
{\left(\delta N_i^{\rm corr}\right)^2}
\end{equation}
where the sum runs over the fitted kinematic variables ($Q^2$, $|t|$, $W$)
and their bins (5 bins each, 15 data points total).
The variables $x_B$ and $|t|-|t_{\rm min}|$ are not fitted since the model has
no independent parameters controlling them.

The model prediction $N_i^{\rm model}$ is obtained by reweighting the MC events:
\begin{equation}
N_i^{\rm model} = C \sum_{j \in {\rm bin}\,i} \frac{\sigma_{\rm new}(Q^2_j, x_{B,j}, t_j)}
{\sigma_{\rm gen}(Q^2_j, x_{B,j}, t_j)}
\end{equation}
where $\sigma_{\rm gen}$ is the cross section used in generation and $\sigma_{\rm new}$
uses the trial parameters. The normalization $C$ is determined analytically.

The optimization uses \texttt{scipy.optimize.differential\_evolution} with bounds:
$\alpha_2 \in [-2, 5]$, $\alpha_3 \in [-15, 15]$, $b_t \in [0.1, 8]$,
$\nu_T \in [0, 10]$. The parameter $c_R = 1$ is fixed.

\subsection{Results}

"""

    tex += (r"\begin{table}[H]" + "\n" + r"\centering" + "\n"
            + f"\\caption{{Fitted model parameters for each dataset. Generation parameters:\n{gen_line}.\n"
            + r"$\chi^2$ is for 15 data points (5 bins $\times$ 3 variables) minus 4 free parameters.}"
            + "\n" + r"\begin{tabular}{lccccc}" + "\n" + r"\toprule" + "\n"
            + r"Dataset & $\chi^2$/ndf & $\alpha_2$ & $\alpha_3$ & $b_t$ & $\nu_T$ \\" + "\n"
            + r"\midrule" + "\n")
    tex += (f"MC generation & --- & ${gp['alf2']:.3f}$ & ${gp['alf3']:.3f}$ & "
            f"${gp['bt']:.3f}$ & ${gp['nuT']:.3f}$ \\\\\n\\midrule\n")
    for short, chi2, ndf, a2, a3, bt, nt in fit_rows:
        tex += f"{short} & {chi2:.1f}/{ndf} & {a2:.3f} & {a3:.3f} & {bt:.3f} & {nt:.3f} \\\\\n"
    tex += r"\midrule" + "\n"
    tex += (f"Average & --- & ${avg['alf2']:.3f} \\pm {std['alf2']:.3f}$ & "
            f"${avg['alf3']:.3f} \\pm {std['alf3']:.3f}$ & "
            f"${avg['bt']:.3f} \\pm {std['bt']:.3f}$ & "
            f"${avg['nuT']:.3f} \\pm {std['nuT']:.3f}$ \\\\\n")
    tex += (f"Avg (inb) & --- & ${inb_avg['alf2']:.3f} \\pm {inb_std['alf2']:.3f}$ & "
            f"${inb_avg['alf3']:.3f} \\pm {inb_std['alf3']:.3f}$ & "
            f"${inb_avg['bt']:.3f} \\pm {inb_std['bt']:.3f}$ & "
            f"${inb_avg['nuT']:.3f} \\pm {inb_std['nuT']:.3f}$ \\\\\n")
    tex += (f"Avg (outb) & --- & ${outb_avg['alf2']:.3f} \\pm {outb_std['alf2']:.3f}$ & "
            f"${outb_avg['alf3']:.3f} \\pm {outb_std['alf3']:.3f}$ & "
            f"${outb_avg['bt']:.3f} \\pm {outb_std['bt']:.3f}$ & "
            f"${outb_avg['nuT']:.3f} \\pm {outb_std['nuT']:.3f}$ \\\\\n")
    tex += r"""\bottomrule
\end{tabular}
\end{table}

"""

    # chi2 per variable table
    tex += r"""\begin{table}[H]
\centering
\caption{$\chi^2$ per kinematic variable (5 bins each). Variables $x_B$ and $|t|-|t_{\rm min}|$
are not included in the fit but shown for comparison.}
\begin{tabular}{lccccc}
\toprule
Dataset & $Q^2$ & $|t|$ & $W$ & $x_B$ (not fit) & $|t|-|t_{\rm min}|$ (not fit) \\
\midrule
"""
    for ds in DATASETS:
        if ds in all_acc:
            cpv = all_acc[ds]["chi2_per_var"]
            cells = []
            for v in ["q2", "t", "w", "xb", "t_tmin"]:
                cells.append(f"{cpv[v]['chi2']:.1f}" if v in cpv else "---")
            tex += f"{DS_SHORT[ds]} & {' & '.join(cells)} \\\\\n"
    tex += r"""\bottomrule
\end{tabular}
\end{table}

\subsection{Pipeline Fit Plots}

"""
    for ds in DATASETS:
        lab = DS_LABELS[ds]
        if "step5_fit" in plot_files.get(ds, {}):
            tex += (f"\\begin{{figure}}[H]\n\\centering\n"
                    f"\\includegraphics[width=0.95\\textwidth]{{{plot_files[ds]['step5_fit']}}}\n"
                    f"\\caption{{Pipeline step 5: acceptance-corrected fit for {lab} "
                    f"(raw counts, not M(KK) signal extraction).}}\n"
                    f"\\end{{figure}}\n\n")

    tex += r"""
\clearpage
%% ==================================================================
\section{Summary}
\label{sec:summary}

The cross-section model for diffractive $\phi$ meson electroproduction has been
tuned using five CLAS12 datasets with the generation parameters
"""
    tex += f"{gen_line}.\n"
    tex += r"""
The acceptance-corrected $\phi$ yields, extracted via $M(K^+K^-)$ fits with fixed
Gaussian parameters ($\mu = 1019.5$~MeV, $\sigma = 4.7$~MeV), are fit using
reweighting of the Monte Carlo events.

The average fitted parameters across all 5 datasets are:
\begin{center}
\begin{tabular}{lc}
\toprule
Parameter & Average $\pm$ std \\
\midrule
"""
    tex += (f"$\\alpha_2$ & ${avg['alf2']:.3f} \\pm {std['alf2']:.3f}$ \\\\\n"
            f"$\\alpha_3$ & ${avg['alf3']:.3f} \\pm {std['alf3']:.3f}$ \\\\\n"
            f"$b_t$       & ${avg['bt']:.3f} \\pm {std['bt']:.3f}$ \\\\\n"
            f"$\\nu_T$    & ${avg['nuT']:.3f} \\pm {std['nuT']:.3f}$ \\\\\n")
    tex += r"""\bottomrule
\end{tabular}
\end{center}

These values can be used as input for the next iteration of the generator.

\appendix
\section{Generated Kinematics}

"""
    for ds in DATASETS:
        lab = DS_LABELS[ds]
        if "step1_generated" in plot_files.get(ds, {}):
            tex += (f"\\begin{{figure}}[H]\n\\centering\n"
                    f"\\includegraphics[width=0.92\\textwidth]{{{plot_files[ds]['step1_generated']}}}\n"
                    f"\\caption{{Generated kinematics for {lab}.}}\n"
                    f"\\end{{figure}}\n\n")

    tex += r"""
\section{MC vs Data Comparison}

"""
    for ds in DATASETS:
        lab = DS_LABELS[ds]
        if "step3_mc_vs_data" in plot_files.get(ds, {}):
            tex += (f"\\begin{{figure}}[H]\n\\centering\n"
                    f"\\includegraphics[width=0.92\\textwidth]{{{plot_files[ds]['step3_mc_vs_data']}}}\n"
                    f"\\caption{{MC vs.\\ data kinematic comparison for {lab}.}}\n"
                    f"\\end{{figure}}\n\n")

    tex += r"""
\section{Acceptance Maps}

"""
    for ds in DATASETS:
        lab = DS_LABELS[ds]
        if "step4_acceptance" in plot_files.get(ds, {}):
            tex += (f"\\begin{{figure}}[H]\n\\centering\n"
                    f"\\includegraphics[width=0.92\\textwidth]{{{plot_files[ds]['step4_acceptance']}}}\n"
                    f"\\caption{{Acceptance maps for {lab}.}}\n"
                    f"\\end{{figure}}\n\n")

    tex += "\n\\end{document}\n"

    tex_path = TEX_DIR / "phi_tuning_report.tex"
    with open(tex_path, "w") as f:
        f.write(tex)
    print(f"Written: {tex_path}")

    # Compile (twice for TOC)
    try:
        for _ in range(2):
            result = subprocess.run(
                ["pdflatex", "-interaction=nonstopmode", "phi_tuning_report.tex"],
                capture_output=True, text=True, timeout=120, cwd=str(TEX_DIR))
        if result.returncode == 0:
            print(f"PDF compiled: {TEX_DIR / 'phi_tuning_report.pdf'}")
        else:
            print(f"pdflatex exit {result.returncode}; check {TEX_DIR}/phi_tuning_report.log")
    except FileNotFoundError:
        print("pdflatex not found — .tex ready for manual compilation")


if __name__ == "__main__":
    main()
