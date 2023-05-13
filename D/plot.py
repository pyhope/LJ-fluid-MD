import numpy as np
from matplotlib import pyplot as plt
import pickle as pkl
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

def fit(x, y, deg):
    z = np.polyfit(x, y, deg=deg, full=True)
    pn = np.poly1d(z[0])
    R_square = 1 - z[1][0] / sum((y - np.mean(y))** 2)
    return pn, R_square

def vacftoPS(vacf, w, t):
    dt = t[1] - t[0]
    integrand = vacf * np.cos(w * t)
    return np.sum(integrand * dt)

with open('data.pkl', 'rb') as file:
    D = pkl.load(file)

with open('vacf.pkl', 'rb') as file:
    vacflist = pkl.load(file)

num_list = ['108', '192', '256', '400', '500']
clist = ['C1', 'C2', 'C4', 'C5', 'C0']
l = [17.34045, 23.1206, 28.90075]
length = np.array([(l[i]*l[j]*l[k])**(1/3) for i, j, k in [(0,0,0),(1,1,0),(1,1,1),(2,2,1),(2,2,2)]])
k_B = 1.380649e-23 #J/K

# for i in num_list:
#     p = [np.mean(D[i][j]) for j in ['msd', 'vacf', 'temp']]
#     print(i, p[0], p[1], p[2])

for i in num_list:
    p = [np.mean(D[i][j]) for j in ['msd', 'vacf', 'temp']]
    p = [np.std(D[i][j]) for j in ['msd', 'vacf', 'temp']]
T = np.mean([np.mean(D[i]['temp']) for i in num_list])

fig, ax = plt.subplots(figsize=(8, 8))

for sys in ['msd', 'vacf']:
    Diffusion = np.array([np.mean(D[i][sys]) for i in num_list]) * 1e9
    err = np.array([2 * np.std(D[i][sys]) / np.sqrt(len(D[i][sys])) for i in num_list]) * 1e9
    ax.scatter(1/length, Diffusion, marker='s', s=50, c='C0', edgecolors=(0.122, 0.467, 0.706, 0.6), linewidths=2, zorder=10)
    ax.errorbar(1/length, Diffusion, marker='s', yerr=err, ls='none', ecolor='k', elinewidth=1, capsize=4, capthick=1, zorder=0)
    linear, Rs = fit(1/length, Diffusion, deg=1)
    x = np.linspace(0, 1 / np.min(length), 1000)
    ax.plot(x, linear(x), ls='--', c = 'gray')
    ax.set_xlim(0, 0.06)
    ax.set_ylim(2.02, 2.72)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlabel('$1/L$ ($\mathrm{\AA^{-1}}$)')
    ax.set_ylabel(r'$D$ ($10^{-9}$ m$^2$/s)')
    plt.savefig('../D_' + sys + '.pdf', bbox_inches='tight')
    plt.cla()
    eta = - k_B * T * 2.837 / (6 * np.pi * linear[1] * 1e-19) * 1e6
    print(sys, eta, 'uPa·s', linear[0], '*1e-9 m2/s')

# msd  183.9 uPa·s 2.702 *1e-9 m2/s
# vacf 192.3 uPa·s 2.699 *1e-9 m2/s

w = np.linspace(0, 0.04, 10000)

for index, i in enumerate(num_list):
    vacf = vacflist[index][:5000]
    t = np.arange(len(vacf))
    ps = np.array([vacftoPS(vacf, j, t) for j in w]) / (300 * Diffusion[index])
    ax.plot(w, ps, c = clist[index], label = '$N$ = ' + i, alpha = 0.8)
ax.plot(w, 1e-4 / (w**2 + 1e-4), c = 'k', ls = '--', label = 'Brownian fluid', alpha = 0.5)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.set_xlabel('$\omega$ (fs$^{-1}$)')
ax.set_ylabel('$f(\omega$)')
ax.legend(fancybox=False, edgecolor='black', fontsize=14)
plt.savefig('../ps.pdf', bbox_inches='tight')
plt.show()