#!/usr/bin/env python3
"""
Master script: full chain for all 5 datasets + LaTeX report.
Runs completely unattended.

1) Generate MC with new parameters (5M events × 5 datasets)
2) FastMC for each
3) Fit M(KK) in bins for all datasets
4) Acceptance correction + model fit
5) Write LaTeX report
"""

import os, sys, json, time, subprocess, shutil
import numpy as np
from pathlib import Path

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

BASE = Path("/Users/vpk/vector_mesons_generator_tuning/phi")
TEX_DIR = BASE / "tex"
TEX_DIR.mkdir(parents=True, exist_ok=True)

DATASETS = ["rgafall18_inb", "rgafall18_outb", "rgasp18_inb", "rgasp18_outb", "rgasp19_inb"]

LOG = open(BASE / "full_chain.log", "w")

def log(msg):
    print(msg, flush=True)
    LOG.write(msg + "\n")
    LOG.flush()

# =====================================================================
#  PHASE 1+2: Generate MC + FastMC via phi_pipeline.py
# =====================================================================
log("=" * 70)
log("  PHASE 1+2: Generate MC + FastMC for all datasets")
log("=" * 70)

for ds in DATASETS:
    log(f"\n>>> {ds} ...")
    t0 = time.time()
    cfg_path = BASE / "runs" / ds / "config.json"
    result = subprocess.run(
        [sys.executable, str(BASE / "phi_pipeline.py"), str(cfg_path)],
        capture_output=True, text=True, timeout=1800,
        cwd=str(BASE),
        env={**os.environ}
    )
    dt = time.time() - t0
    log(f"    exit={result.returncode}, time={dt:.0f}s")
    # Save stdout
    with open(BASE / "runs" / ds / "pipeline.log", "w") as f:
        f.write(result.stdout)
    if result.stderr:
        with open(BASE / "runs" / ds / "pipeline_err.log", "w") as f:
            f.write(result.stderr)
    if result.returncode != 0:
        log(f"    ERROR! stderr tail:\n{result.stderr[-500:]}")
        log("    Continuing to next dataset...")
    else:
        # Extract key info from stdout
        for line in result.stdout.split("\n"):
            if "PIPELINE COMPLETE" in line or "Result:" in line or "Chi2:" in line:
                log(f"    {line.strip()}")

# =====================================================================
#  PHASE 3: Fit M(KK) in bins
# =====================================================================
log("\n" + "=" * 70)
log("  PHASE 3: Fit M(KK) in bins")
log("=" * 70)

t0 = time.time()
result = subprocess.run(
    [sys.executable, str(BASE / "fit_mkk_binned.py")],
    capture_output=True, text=True, timeout=600,
    cwd=str(BASE),
    env={**os.environ}
)
dt = time.time() - t0
log(f"  exit={result.returncode}, time={dt:.0f}s")
with open(BASE / "mkk_fits.log", "w") as f:
    f.write(result.stdout)
if result.returncode != 0:
    log(f"  ERROR! {result.stderr[-500:]}")

# =====================================================================
#  PHASE 4: Acceptance correction + model fit
# =====================================================================
log("\n" + "=" * 70)
log("  PHASE 4: Acceptance correction + model fit")
log("=" * 70)

t0 = time.time()
result = subprocess.run(
    [sys.executable, str(BASE / "acceptance_correct.py")],
    capture_output=True, text=True, timeout=1800,
    cwd=str(BASE),
    env={**os.environ}
)
dt = time.time() - t0
log(f"  exit={result.returncode}, time={dt:.0f}s")
with open(BASE / "acceptance.log", "w") as f:
    f.write(result.stdout)
if result.returncode != 0:
    log(f"  ERROR! {result.stderr[-500:]}")
else:
    for line in result.stdout.split("\n"):
        if "SUMMARY" in line or "alf2" in line.lower() or "chi2" in line.lower() \
           or "---" in line or any(ds in line for ds in DATASETS):
            log(f"  {line.strip()}")

# =====================================================================
#  PHASE 5: Write LaTeX report (delegated to make_report.py)
# =====================================================================
log("\n" + "=" * 70)
log("  PHASE 5: Write LaTeX report")
log("=" * 70)

result = subprocess.run(
    [sys.executable, str(BASE / "make_report.py")],
    capture_output=True, text=True, timeout=600,
    cwd=str(BASE), env={**os.environ})
log(result.stdout.strip())
if result.returncode != 0:
    log(f"  REPORT ERROR! {result.stderr[-500:]}")

log("\n" + "=" * 70)
log("  FULL CHAIN COMPLETE")
log("=" * 70)
LOG.close()
