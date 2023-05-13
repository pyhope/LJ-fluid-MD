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

def structure_factor(q, r, g_r, rho):
    dr = r[1] - r[0]
    integrand = 4 * np.pi * r * rho * (g_r - 1) * np.sin(q * r) / q
    return 1 + np.sum(integrand) * dr

g2_exp = np.loadtxt('g2.dat', unpack=True)
Sq_exp = np.loadtxt('Sq.dat', unpack=True)

fig, ax = plt.subplots(figsize=(10, 8))
rho = 0.0207
q_values = np.linspace(0.01, 10, 500)
ax.set_xlabel('$q$ (Å$^{-1}$)')
ax.set_ylabel('$S(q)$')

clist = ['C1', 'C2', 'C4', 'C5', 'C0']
num_list = ['108', '192', '256', '400', '500']

ax.plot(Sq_exp[0], Sq_exp[1], c='k', label = 'Yarnell et al. (1973)', alpha = 0.5)
for index, i in enumerate(num_list):
    data = np.loadtxt('./' + i + '/rdf.dat', skiprows=4, unpack=True)
    r = data[1]
    g_r= data[2]
    s_q_values = [structure_factor(q, r, g_r, rho) for q in q_values]
    ax.plot(q_values, s_q_values, c = clist[index], label = '$N$ = ' + i, alpha = 0.8)

plt.legend(fancybox=False, edgecolor='black', fontsize=14)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.savefig('../Sq.pdf', bbox_inches='tight')

plt.cla()
ax.plot(g2_exp[0], g2_exp[1], c='k', label = 'Yarnell et al. (1973)', alpha = 0.5)
for index, i in enumerate(num_list):
    data = np.loadtxt('./' + i + '/rdf.dat', skiprows=4, unpack=True)
    r = data[1]
    g_r = data[2]
    ax.plot(r, g_r, c = clist[index], label = '$N$ = ' + i, alpha = 0.8)
ax.set_xlabel('$r$ (Å)')
ax.set_ylabel('$g_2(r)$')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

plt.legend(fancybox=False, edgecolor='black', fontsize=14)
plt.savefig('../g2.pdf', bbox_inches='tight')