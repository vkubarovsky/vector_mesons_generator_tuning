#!/usr/bin/env python3
"""
make_jpsi_report.py — LaTeX report for the J/psi tuning chain.
Reads config.json, jpsi_results.json, jpsi_xsec_results.json and the PNGs.
"""

import json
import shutil
import subprocess
from pathlib import Path

BASE = Path(__file__).resolve().parent
TEX_DIR = BASE / "tex"

# iteration history (gen -> fit); the last row is filled from jpsi_results.json
# Retune 2026-06-12: single-fastMC acceptance (F18in_45nA @10.6, S19in_50nA @10.2),
# generator diffrad_vm. Iter 1 reused the previous campaign's MC (gen at 4.105/3.106)
# with recomputed F18in_45nA-only weights.
ITER_HISTORY = [
    {"it": 1, "gen": (4.105, 3.106), "fit": (4.234, 4.065), "chi2": 20.4},
]


def main():
    TEX_DIR.mkdir(exist_ok=True)
    plots = TEX_DIR / "plots"
    plots.mkdir(exist_ok=True)

    with open(BASE / "config.json") as f:
        cfg = json.load(f)
    with open(BASE / "jpsi_results.json") as f:
        res = json.load(f)
    with open(BASE / "jpsi_xsec_results.json") as f:
        xs = json.load(f)

    gp = res["gen_params"]; fp = res["fit_params"]
    iters = list(ITER_HISTORY)
    iters.append({"it": len(iters) + 1,
                  "gen": (gp["alf2"], gp["mg2"]),
                  "fit": (fp["alf2"], fp["mg2"]),
                  "chi2": res["chi2"]})

    for png in ["jpsi_corrected.png", "step_gen_vs_acc.png",
                "jpsi_dsdt_mariana_bins.png", "jpsi_sigma_q2_w.png"]:
        if (BASE / png).exists():
            shutil.copy2(BASE / png, plots / png)

    nev106 = cfg["generation"]["energies"]["10.6"]["nev"]
    nev102 = cfg["generation"]["energies"]["10.2"]["nev"]
    nconf106 = len(cfg["generation"]["energies"]["10.6"]["fastmc_configs"])

    tex = r"""\documentclass[12pt,a4paper]{article}
\usepackage[margin=2.5cm]{geometry}
\usepackage{graphicx}
\usepackage{amsmath,amssymb}
\usepackage{booktabs}
\usepackage{hyperref}
\usepackage{float}

\title{Exclusive $J/\psi$ Electroproduction at CLAS12:\\
Monte Carlo Tuning with the Forward Tagger and NN Fast Simulation}
\author{V.~Kubarovsky}
\date{\today}

\begin{document}
\maketitle

\begin{abstract}
Tuning of the DIFFRAD Monte Carlo generator for exclusive $J/\psi$
electroproduction $ep \to e' p\, J/\psi$, $J/\psi \to e^+e^-$, with the
scattered electron in the CLAS12 Forward Tagger (quasi-photoproduction,
$Q^2 \lesssim 0.1$~GeV$^2$). The model uses a dipole $t$-dependence.
The combined CLAS12 RGA sample of 68 $J/\psi$ candidates
(52 at 10.6~GeV, 16 at 10.2~GeV) is acceptance-corrected with neural-network
fast simulation and fitted by MC reweighting.
\end{abstract}

\tableofcontents
\newpage

%% ==================================================================
\section{Cross-Section Model}

\begin{equation}
\frac{d\sigma_T}{dt} = \alpha_1 \left(1-\frac{W_{\rm th}^2}{W^2}\right)^{\alpha_2}
W^{\alpha_3} \left(1+\frac{Q^2}{m_\psi^2}\right)^{-\nu_T}
\cdot \frac{3\,(m_g^2 - t_{\rm min})^3}{(m_g^2-t)^4}
\end{equation}
\begin{equation}
\frac{d\sigma_L}{dt} = c_R\,\frac{Q^2}{m_\psi^2}\,\frac{d\sigma_T}{dt},
\qquad W_{\rm th} = M_N + m_\psi = 4.035~\text{GeV}.
\end{equation}
Here $t, t_{\rm min} < 0$ and the dipole factor is normalized:
$\int_{-\infty}^{t_{\rm min}} 3(m_g^2-t_{\rm min})^3/(m_g^2-t)^4\,dt = 1$,
so the $t$-integrated cross section is carried entirely by the $W$ and $Q^2$
factors. Near threshold $|t|_{\rm min}$ is large (2.7~GeV$^2$ at $W_{\rm th}$,
1.0~GeV$^2$ at $W=4.5$), which makes the correct $t_{\rm min}$ treatment
essential.

In Forward Tagger kinematics ($Q^2 \in [0.015, 0.1]$~GeV$^2$) the factors
$(1+Q^2/m_\psi^2)^{-\nu_T}$ and $c_R Q^2/m_\psi^2$ vary by less than 1\%:
the model has no usable $Q^2$ lever arm, and $\nu_T$, $c_R$, $\alpha_3$ are
kept fixed. The free parameters are $\alpha_2$ (threshold behaviour) and
$m_g^2$ (dipole $t$-slope).

%% ==================================================================
\section{Monte Carlo Generator}

\texttt{diffrad\_vm.f90} --- the combined all-meson production generator
(one-program policy; it supersedes the per-meson
\texttt{diffrad\_jpsi\_dipole.f90}, which is kept as frozen reference).
All model parameters ($\alpha_1,\alpha_2,\alpha_3,\nu_T,m_g^2,c_R$)
are read from the input file. The accept--reject proposal uses the same dipole
shape, so the sampling efficiency tracks $m_g^2$ automatically.

Two bugs were found and fixed during the original tuning campaign:
\begin{itemize}
\item \textbf{Dipole normalization sign}: the code evaluated
$(m_g^2 - |t_{\rm min}|)^3$ instead of $(m_g^2 + |t_{\rm min}|)^3$.
Since $|t|_{\rm min} > m_g^2$ below $W \approx 4.23$~GeV, the cross
section went negative and threshold events were silently rejected.
\item \textbf{32-bit overflow}: the attempt limit \texttt{maxtry = 100000*nev}
overflowed INTEGER*4 for $n_{\rm ev} \geq 2$M, terminating generation
instantly. Counters are now \texttt{integer*8} (also fixed in the $\phi$
generator).
\end{itemize}

"""
    tex += f"""Samples: {nev106//1000000}M events at 10.6 GeV, {nev102//1000000}M at 10.2 GeV.
Generation ranges: $Q^2 \\in [0.01, 0.12]$, $W \\in [4.05, 4.55]$,
$|t| < 6$~GeV$^2$, $\\theta_{{e'}} \\in [0^\\circ, 7.5^\\circ]$.

%% ==================================================================
\\section{{Data}}

Combined RGA $J/\\psi$ candidate sample (Mariana), 68 events
(52 at 10.6 GeV, 16 at 10.2 GeV), four detected particles:
$e'$ in the Forward Tagger ($\\theta = 2.7$--$6.0^\\circ$),
decay $e^+e^-$ and recoil proton in the Forward Detector.
$M(e^+e^-) \\in [2.93, 3.17]$~GeV --- the sample is pre-selected, the yield
is obtained by counting (no mass fit). 64 events pass the analysis cuts
($Q^2 \\in [{cfg['cuts']['q2min']}, {cfg['cuts']['q2max']}]$,
$W \\in [{cfg['cuts']['wmin']}, {cfg['cuts']['wmax']}]$, $|t| < {cfg['cuts']['tmax']}$).

%% ==================================================================
\\section{{Fast Monte Carlo and Acceptance}}

NN models (acceptance + smearing) trained per particle and per run
configuration: $e^-$ (FT, no vertex measurement, $v_z$ fixed at $-3$~cm),
decay $e^-$, $e^+$, $p$. A \\emph{{single}} fastMC configuration is used per
beam energy: \\textbf{{F18in\\_45nA}} for the 10.6 GeV MC and
\\textbf{{S19in\\_50nA}} for the 10.2 GeV MC. Rationale: the data LUND file
records the beam energy per event but not the torus polarity; by Mariana's
per-period yields in the $e'e^+e^-p$ topology the 10.6 GeV sample is
$\\approx 84\\%$ inbending ($N_{{J/\\psi}}$: F18in $16.0\\pm9.3$, S18in
$14.5\\pm5.5$, S18out $5.9\\pm5.0$, F18out consistent with zero), and the
10.2 GeV running (S19) was inbending only. A multi-configuration average
with equal weights would mis-weight periods that contributed few or no
events. The two beam energies are combined so that the \\emph{{accepted}}
MC matches the data composition 52:16.

\\begin{{figure}}[H]
\\centering
\\includegraphics[width=0.99\\textwidth]{{plots/step_gen_vs_acc.png}}
\\caption{{Generated vs fastMC-accepted distributions for both beam energies.}}
\\end{{figure}}

%% ==================================================================
\\section{{Tuning Iterations}}

Fit: acceptance-corrected counts in 3 $Q^2$ + 3 $W$ + 4 $|t|$ bins (10 points),
model predictions by event-by-event reweighting
$\\sigma_{{\\rm new}}/\\sigma_{{\\rm gen}}$, free parameters $\\alpha_2$ and $m_g^2$,
\\texttt{{differential\\_evolution}} minimization.

\\begin{{table}}[H]
\\centering
\\caption{{Retune iterations (single-fastMC acceptance, 2026-06-12).
Starting point: the previous multi-configuration tune
$\\alpha_2=4.105$, $m_g^2=3.106$; iteration 1 reuses that campaign's MC
with F18in\\_45nA-only acceptance weights, iteration 2 regenerates with
\\texttt{{diffrad\\_vm}} at the iteration-1 fit.}}
\\begin{{tabular}}{{ccccccc}}
\\toprule
Iter & $\\alpha_2^{{\\rm gen}}$ & $m_g^{{2,\\rm gen}}$ &
$\\alpha_2^{{\\rm fit}}$ & $m_g^{{2,\\rm fit}}$ & $\\chi^2$/ndf \\\\
\\midrule
"""
    for it in iters:
        tex += (f"{it['it']} & {it['gen'][0]:.3f} & {it['gen'][1]:.3f} & "
                f"{it['fit'][0]:.3f} & {it['fit'][1]:.3f} & {it['chi2']:.1f}/8 \\\\\n")
    tex += r"""\bottomrule
\end{tabular}
\end{table}

\begin{figure}[H]
\centering
\includegraphics[width=0.99\textwidth]{plots/jpsi_corrected.png}
\caption{Acceptance-corrected counts vs generation model (blue) and fitted
model (red). $Q^2$, $W$, $|t|$ are fitted; $|t|-|t|_{min}$ and $x_B$ are
cross-checks. The first $Q^2$ bin carries the dominant $\chi^2$ contribution;
the model $Q^2$ shape is flux-driven and not tunable in FT kinematics.}
\end{figure}

%% ==================================================================
\section{Cross-Section Shapes in the (W, Q$^2$) Analysis Bins}

Flux-folded $\gamma^* p \to J/\psi\,p$ shapes in the three correlated
$(W, Q^2)$ bins of the CLAS12 analysis (Table~16 of the analysis note).
Absolute normalization requires the integrated luminosity and is not
attempted here.

\begin{figure}[H]
\centering
\includegraphics[width=0.99\textwidth]{plots/jpsi_dsdt_mariana_bins.png}
\caption{$d\sigma/dt$ shape in the three $(W,Q^2)$ bins:
acceptance-corrected data (points) vs tuned model (red).}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=0.9\textwidth]{plots/jpsi_sigma_q2_w.png}
\caption{Cross-section shape vs $Q^2$ and $W$, integrated over $|t|<3$~GeV$^2$.}
\end{figure}

%% ==================================================================
\section{Acceptance Comparison with the CLAS12 Analysis}

"""
    tex += r"""\begin{table}[H]
\centering
\caption{Our NN-fastMC mixture acceptance vs the efficiencies of the CLAS12
analysis (Table~17): $\eta_1$ event-by-event, $\eta_{2a}$/$\eta_{2b}$
bin-averaged with two methods.}
\begin{tabular}{cccccccc}
\toprule
Bin & $Q^2$ [GeV$^2$] & $W$ [GeV] & $\varepsilon$ (this work) &
$\eta_1$ & $\eta_{2a}$ & $\eta_{2b}$ & $N_{\rm data}$ \\
\midrule
"""
    for r in xs["acceptance_mariana_bins"]:
        tex += (f"{r['bin']} & [{r['q2'][0]:.3f}, {r['q2'][1]:.3f}] & "
                f"[{r['w'][0]:.3f}, {r['w'][1]:.3f}] & {r['eps_ours']:.4f} & "
                f"{r['eta1']:.4f} & {r['eta2a']:.4f} & {r['eta2b']:.4f} & "
                f"{r['n_data']} \\\\\n")
    tex += r"""\bottomrule
\end{tabular}
\end{table}

%% ==================================================================
\section{Summary}

"""
    tex += (f"The tuned $J/\\psi$ generator parameters are "
            f"$\\alpha_2 = {fp['alf2']:.2f}$, $m_g^2 = {fp['mg2']:.2f}$~GeV$^2$ "
            f"($\\alpha_3 = {fp['alf3']}$, $\\nu_T = {fp['nuT']}$, "
            f"$c_R = {fp['cR']}$ fixed), $\\chi^2 = {res['chi2']:.1f}$ for "
            f"{res['ndf']}~ndf. ")
    tex += r"""The data prefer a much steeper threshold rise and a flatter
$t$-dependence than the starting (photoproduction-inspired) parametrization.
With 68 events the two-parameter description is adequate; the remaining
tension is confined to the lowest $Q^2$ bin, where the model shape is fixed
by the virtual photon flux.

\end{document}
"""

    tex_path = TEX_DIR / "jpsi_tuning_report.tex"
    with open(tex_path, "w") as f:
        f.write(tex)
    print(f"Written: {tex_path}")
    try:
        for _ in range(2):
            r = subprocess.run(["pdflatex", "-interaction=nonstopmode",
                                "jpsi_tuning_report.tex"],
                               capture_output=True, text=True, timeout=120,
                               cwd=str(TEX_DIR))
        if r.returncode == 0:
            print(f"PDF: {TEX_DIR/'jpsi_tuning_report.pdf'}")
        else:
            print(f"pdflatex exit {r.returncode}")
    except FileNotFoundError:
        print("pdflatex not found")


if __name__ == "__main__":
    main()
