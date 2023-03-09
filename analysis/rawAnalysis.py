import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

"""
Julia code to produce VMC energies

var_alpha = make_polaron(LinRange(0.01, 10, 1000), [0.1], [0.0]; ω=1.0, rtol = 1e-4, verbose = true, threads = true).F[:, 1]
CSV.write("var_alpha.csv",  Tables.table(var_alpha), writeheader=false)
"""


def perturbative(alpha, beta):
    return - alpha * (1 + (1/(4*beta)) + (9/(32*(beta**2))))


datadir = datadir = "../data/"

# Variational energies data
filename_VMC = "PolaronEnergyVMC.csv"
path_VMC = os.path.join(datadir, filename_VMC)
df_VMC = pd.read_csv(path_VMC)
#df_VMC = df_VMC.head(800)

# Raw data for PIMC T=0.1
alpha = [1, 2, 3, 4, 5, 6, 7, 8, 9]
energy =  [-0.511, -1.514, -2.703, -3.813, -5.046, -6.487, -8.01, -10.24, -12.24]
energy_err = [0.02, 0.097, 0.152, 0.05, 0.416, 0.0697, 0.118, 0.28, 0.45]

# Figure Parameters
#plt.rcParams.update({
#    "text.usetex": True,
#    "font.family": "sans-serif",
#    "font.sans-serif": "Dejavu Sans",
#    "font.size": 14
#})

plt.rcParams.update({
    "text.usetex": False,
    'font.sans-serif': 'CMU Sans Serif',
    'font.family': 'sans-serif',
    'text.color':'darkslategrey',
    'axes.labelcolor':'darkslategrey',
    "font.size": 16
})

fig, ax = plt.subplots(1, 1, figsize=(7, 4))


ax.tick_params(color='darkslategrey', labelcolor='darkslategrey')
for spine in ax.spines.values():
    spine.set_edgecolor('darkslategrey')

ax.errorbar(alpha, energy, yerr=energy_err, ls='none', marker='x', ms=1, color='chocolate')
ax.scatter(alpha, energy, marker='x', s=70, linewidth=2, color='chocolate', label='Path Integral Monte-Carlo')
ax.plot(df_VMC["alpha"], df_VMC["energy"], color="lightblue", zorder=0, label='Variational Method')
ax.plot(df_VMC["alpha"], perturbative(df_VMC["alpha"], 1/0.1), ls="dashed", color="lightsteelblue", zorder=0, label='Perturbative Theory')
ax.set_ylabel(r"Polaron Energy $E(\alpha)$ $T=0.1$")
ax.set_xlabel(r"Coupling Parameter $\alpha$")
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2))
ax.grid(color='darkslategrey', linestyle='dotted', linewidth=0.5)
plt.savefig('EnergyAlphaPlot.png', dpi=600, bbox_inches='tight', transparent=True)
plt.show()



# Variational energies data
filename_VMC = "PolaronEnergyVMC_Temp.csv"
path_VMC = os.path.join(datadir, filename_VMC)
df_VMC = pd.read_csv(path_VMC)
df_VMC = df_VMC.head(800)

T = np.arange(0.1, 0.8, 0.1)
PIMC_1 = [-0.531, -0.498, -0.337, -0.244, -0.18, -0.0484, 0.0128]
PIMC_1_err = [0.02, 0.019, 0.025, 0.027, 0.033, 0.03, 0.022]

fig, ax = plt.subplots(1, 1, figsize=(7, 4))

ax.tick_params(color='darkslategrey', labelcolor='darkslategrey')
for spine in ax.spines.values():
    spine.set_edgecolor('darkslategrey')

per = []
for temp in df_VMC["T"]:
    per.append(perturbative(0.1, 1/temp))

ax.errorbar(T, PIMC_1, yerr=PIMC_1_err, ls='none', marker='x', ms=1, color='chocolate')
ax.scatter(T, PIMC_1, marker='x', s=50, linewidth=2, color='chocolate', label='Path Integral Monte-Carlo')
ax.plot(df_VMC["T"], df_VMC["Energy"], color="lightblue", zorder=0, label='Variational Method')
#ax.plot(df_VMC["T"], per, ls="dashed", color="lightsteelblue", zorder=0)
ax.set_ylabel(r"Polaron Energy $E(T)$ $\alpha=1$")
ax.set_xlabel(r"Temperature $T$")
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2))
ax.grid(color='darkslategrey', linestyle='dotted', linewidth=0.5)
plt.savefig('EnergyTempPlot.png', dpi=600, bbox_inches='tight', transparent=True)
plt.show()

#%%

datadir = datadir = "../data/"

# Variational energies data
filename_pos = "Position_alpha4.0_T0.1.csv"
path_pos = os.path.join(datadir, filename_pos)
df_pos = pd.read_csv(path_pos)

plt.rcParams.update({
    "text.usetex": False,
    'font.sans-serif': 'CMU Sans Serif',
    'font.family': 'sans-serif',
    'text.color':'darkslategrey',
    'axes.labelcolor':'darkslategrey',
    "font.size": 14
})

fig, ax = plt.subplots(1, 1, figsize=(7, 3))


ax.tick_params(color='darkslategrey', labelcolor='darkslategrey')
for spine in ax.spines.values():
    spine.set_edgecolor('darkslategrey')

ax.hist(df_pos["product_id"], bins=100, density=True, color='#8FAADC')
ax.set_xlabel('Polaron Position')
plt.savefig('PosWavefunction.png', dpi=600, bbox_inches='tight', transparent=True)





