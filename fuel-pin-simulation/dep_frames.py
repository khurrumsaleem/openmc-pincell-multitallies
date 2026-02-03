import os
import matplotlib.pyplot as plt
import numpy as np
import openmc.deplete
from uncertainties import unumpy as unp

# =========================
# Read OpenMC results
# =========================
results = openmc.deplete.ResultsList.from_hdf5("depletion_results.h5")

time, k = results.get_eigenvalue()
openmc_keff = unp.uarray(k[:, 0], k[:, 1])

days = time / (24 * 60 * 60)

# =========================
# Reactor parameters
# =========================
hm_mass = 1.8898051683664847  # kgHM
power_MW = 0.06               # MW
burnup = (power_MW * days) / hm_mass

keff_nom = unp.nominal_values(openmc_keff)
keff_std = unp.std_devs(openmc_keff)

# =========================
# Axis limits 
# =========================
time_min, time_max = days.min(), days.max()
burn_min, burn_max = burnup.min(), burnup.max()

keff_min = (keff_nom - keff_std).min()
keff_max = (keff_nom + keff_std).max()

margin = 0.002  # 
keff_min -= margin
keff_max += margin

# =========================
# Frame folder
# =========================
frame_dir = "plot_frames"
os.makedirs(frame_dir, exist_ok=True)

# =========================
# Create frames ( STEP-BY-STEP)
# =========================
for i in range(1, len(days) + 1):

    # ---------- keff vs Time ----------
    plt.figure(figsize=(6, 4))
    plt.errorbar(
        days[:i],
        keff_nom[:i],
        yerr=keff_std[:i],
        fmt='o-',
        capsize=4
    )
    plt.xlim(time_min, time_max)
    plt.ylim(keff_min, keff_max)
    plt.xlabel("Time (days)")
    plt.ylabel("k-effective")
    plt.title(f"k-effective vs Time (Step {i})")
    plt.grid(True)
    plt.savefig(f"{frame_dir}/keff_time_{i:03d}.png", dpi=150)
    plt.close()

    # ---------- keff vs Burnup ----------
    plt.figure(figsize=(6, 4))
    plt.errorbar(
        burnup[:i],
        keff_nom[:i],
        yerr=keff_std[:i],
        fmt='o-',
        capsize=4
    )
    plt.xlim(burn_min, burn_max)
    plt.ylim(keff_min, keff_max)
    plt.xlabel("Burnup [MWd/kgHM]")
    plt.ylabel("k-effective")
    plt.title(f"k-effective vs Burnup (Step {i})")
    plt.grid(True)
    plt.savefig(f"{frame_dir}/keff_burnup_{i:03d}.png", dpi=150)
    plt.close()

print("✅ Frame-by-frame ")

