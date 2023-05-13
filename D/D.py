import numpy as np
import pickle as pkl
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import LogLocator

rcParams['font.family'] = 'sans-serif'
rcParams['font.size'] = '16'
rcParams['font.sans-serif'] = 'Arial'
rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

def vacftoD(vacf, t):
    dt = t[1] - t[0]
    integrand = vacf
    return np.sum(integrand * dt)

def fit(x, y, deg):
    z = np.polyfit(x, y, deg=deg, full=True)
    pn = np.poly1d(z[0])
    R_square = 1 - z[1][0] / sum((y - np.mean(y))** 2)
    return pn, R_square

fig, ax = plt.subplots(figsize=(10, 8))

q_values = np.linspace(0.01, 10, 500)
ax.set_xlabel('$t$ (fs)')
ax.set_ylabel('$\Delta R^2$ (Å$^{2}$)')

clist = ['C1', 'C2', 'C4', 'C5', 'C0']
num_list = ['108', '192', '256', '400', '500']
D = dict()
for index, i in enumerate(num_list):
    msd = np.zeros(10001)
    t = np.arange(len(msd))
    D[i] = dict()
    D[i]['msd'], D[i]['vacf'], D[i]['temp'] = [], [], []
    for j in range(1, 1001):
        msd_tmp = np.loadtxt('./' + i + '/nve/msd/msd.' + str(j) + '.dat', skiprows=2, unpack=True)[1]
        temp = np.loadtxt('./' + i + '/nve/temp/temp.' + str(j) + '.dat', skiprows=2, unpack=True)[1]
        pn, Rs = fit(t[2000:], msd_tmp[2000:], 1)
        D[i]['msd'].append(pn[1]/6 * 1e-5)
        D[i]['temp'].append(np.mean(temp))
        msd += msd_tmp
    msd /= 1000
    ax.loglog(t, msd, c = clist[index], label = '$N$ = ' + i, alpha = 0.8)

ax.legend(fancybox=False, edgecolor='black', fontsize=14)
ll = LogLocator(base=10.0,subs=([0.1 * i for i in range(1,10)]),numticks=12)
ax.xaxis.set_minor_locator(ll)
ax.yaxis.set_minor_locator(ll)
plt.savefig('../msd.pdf', bbox_inches='tight')

vacflist = []
plt.cla()
for index, i in enumerate(num_list):
    vacf = np.zeros(10001)
    t = np.arange(len(vacf))
    for j in range(1, 1001):
        vacf_tmp = np.loadtxt('./' + i + '/nve/vacf/vacf.' + str(j) + '.dat', skiprows=2, unpack=True)[1]
        D[i]['vacf'].append(vacftoD(vacf_tmp[:5000], t[:5000]) / 3 * 1e-11)
        vacf += vacf_tmp
    vacf /= 1000
    vacflist.append(vacf)
    ax.plot(t[:5000], vacf[:5000], c = clist[index], label = '$N$ = ' + i, alpha = 0.8)

ax.legend(fancybox=False, edgecolor='black', fontsize=14)
ax.set_xlabel('$t$ (fs)')
ax.set_ylabel('VACF ($10^{-6}$ Å$^{2}$fs$^{-2}$)')
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
plt.savefig('../vacf.pdf', bbox_inches='tight')

with open('data.pkl', 'wb') as file:
    pkl.dump(D, file)
with open('vacf.pkl', 'wb') as file:
    pkl.dump(vacflist, file)