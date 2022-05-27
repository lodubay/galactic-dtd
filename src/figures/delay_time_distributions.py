"""
Plot the Type Ia supernova delay time distributions (DTDs)
"""

import sys
import os
sys.path.append(os.path.abspath('../../migration'))
from src.simulations.dtd import PowerLaw, BrokenPowerLaw, Exponential, Bimodal
import matplotlib.pyplot as plt

plaw = {
    'name': 'powerlaw',
    'func': PowerLaw(slope=-1.1),
    'label': r'Power-Law ($\alpha=-1.1$)',
    'color': 'k',
    'line': '-',
}
plaw_steep = {
    'name': 'powerlaw_steep',
    'func': PowerLaw(slope=-1.4),
    'label': r'Power-Law ($\alpha=-1.4$)',
    'color': '#aa3377',
    'line': '--',
}
plaw_broken = {
    'name': 'powerlaw_broken',
    'func': BrokenPowerLaw(),
    'label': 'Broken Power-Law',
    'color': '#ee6677',
    'line': '-.',
}
exp = {
    'name': 'exponential',
    'func': Exponential(timescale=1.5),
    'label': r'Exponential ($\tau=1.5$ Gyr)',
    'color': '#66ccee',
    'line': '-',
}
exp_long = {
    'name': 'exponential_long',
    'func': Exponential(timescale=3),
    'label': r'Exponential ($\tau=3$ Gyr)',
    'color': '#4477aa',
    'line': '--',
}
bimodal = {
    'name': 'bimodal',
    'func': Bimodal(),
    'label': 'Bimodal',
    'color': '#228833',
    'line': '-.',
}

dtds = [bimodal, plaw_steep, plaw, plaw_broken, exp, exp_long]

fig, ax = plt.subplots(figsize=(10, 10), tight_layout=True)
time = [0.001*i for i in range(40, 13201)]
for dtd in dtds:
    func = dtd['func']
    ax.plot(time, [func(t) for t in time], label=dtd['label'], c=dtd['color'], ls=dtd['line'])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim((1e-12, 3e-8))
ax.set_xlabel('Time [Gyr]')
ax.set_ylabel(r'Normalized SN Ia Rate [$\rm{M}_\odot^{-1}$ yr$^{-1}$]')
ax.legend(frameon=False, loc='upper right')
plt.savefig('delay_time_distributions.pdf')
plt.show()
