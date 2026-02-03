import matplotlib.pyplot as plt
import numpy as np
import openmc.deplete
#import serpentTools
#from tabulate import tabulate
from uncertainties import unumpy as unp

# Read OpenMC results ---------------------------------------------------------

# Open results file
results = openmc.deplete.ResultsList.from_hdf5("depletion_results.h5")

# Obtain K_eff as a function of time
time, k = results.get_eigenvalue()
openmc_keff = unp.uarray(k[:, 0], k[:, 1])
days = time/(24*60*60)
print(time)
# second → days
#Initial heavy metal mass [g] = 1889.8051683664846
#Initial heavy metal mass [kg] = 1.8898051683664847
hm_mass=1.8898051683664847
power_MW = 0.06 #MW
# Burnup [MWd/kgHM]
burnup = (power_MW * days) / hm_mass

# Nominal değerler ve belirsizlikler
keff_nom = unp.nominal_values(openmc_keff)
keff_std = unp.std_devs(openmc_keff)

plt.figure()
plt.errorbar(
    days,
    keff_nom,
    yerr=keff_std,
    fmt='o-',
    capsize=4
)

plt.xlabel("Time (days)")
plt.ylabel("k-effective")
plt.title("k-effective vs Time")
plt.grid(True)
plt.show()

plt.figure()
plt.errorbar(
    burnup,
    keff_nom,
    yerr=keff_std,
    fmt='o-',
    capsize=4
)

plt.xlabel("Burnup [MWd/kgHM]")
plt.ylabel("k-effective")
plt.title("k-effective vs Time")
plt.grid(True)
plt.show()
