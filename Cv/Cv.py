import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator

rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = '16'
rcParams['font.sans-serif'] = 'Arial'
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

k_B = 8.617333262e-5 #eV/K

num_list = [108, 192, 256, 400, 500]
Cv, err = [], []

for index, n in enumerate(num_list):
    ke_square = np.loadtxt('./' + str(n) + '/ke.dat', skiprows=2, unpack=True)[1]
    ke_blocks = np.split(ke_square[1:], 100)
    t = np.loadtxt('./' + str(n) + '/temp.dat', skiprows=2, unpack=True)[1]
    avg_ke = np.mean(ke_blocks, axis=1)
    dK = np.var(ke_blocks, axis=1)
    cv = 3 * n * k_B / (2 - 4 * dK / (3 * k_B**2 * t**2 * n))
    Cv.append(np.mean(cv))
    err.append(2 * np.std(cv) / np.sqrt(len(cv)))

fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(num_list, Cv, marker='s', s=50, c='C0', edgecolors=(0.122, 0.467, 0.706, 0.6), linewidths=2, zorder=10)
ax.errorbar(num_list, Cv, marker='s', yerr=err, ls='none', ecolor='k', elinewidth=1, capsize=4, capthick=1, zorder=0)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('$N$')
ax.set_ylabel(r'$C_V$ (eV/K)')
plt.savefig('../Cv.pdf', bbox_inches='tight')