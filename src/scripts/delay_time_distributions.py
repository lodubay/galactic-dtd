"""
Plot the Type Ia supernova delay time distributions (DTDs) as a function of time.
"""

import sys
import os
import pickle
sys.path.append(os.path.abspath(
    '/mnt/c/Users/dubay.11/Repos/galactic-dtd-migration/'))
import matplotlib.pyplot as plt
import paths

PARENT_DIR = paths.data / 'migration' / 'diffusion' / 'insideout'

plaw = {
    'name': 'powerlaw',
    'label': r'Power-Law ($\alpha=-1.1$)',
    'color': 'k',
    'line': '-',
}
plaw_steep = {
    'name': 'powerlaw_steep',
    'label': r'Power-Law ($\alpha=-1.4$)',
    'color': '#aa3377',
    'line': '--',
}
plaw_broken = {
    'name': 'powerlaw_broken',
    'label': 'Broken Power-Law',
    'color': '#ee6677',
    'line': '-.',
}
exp = {
    'name': 'exponential',
    'label': r'Exponential ($\tau=1.5$ Gyr)',
    'color': '#66ccee',
    'line': '-',
}
exp_long = {
    'name': 'exponential_long',
    'label': r'Exponential ($\tau=3$ Gyr)',
    'color': '#4477aa',
    'line': '--',
}
bimodal = {
    'name': 'bimodal',
    'label': 'Bimodal',
    'color': '#228833',
    'line': '-.',
}

dtds = [bimodal, plaw_steep, plaw, plaw_broken, exp, exp_long]


def main():
    fig, ax = plt.subplots(figsize=(3.25, 3.25), tight_layout=True)
    time = [0.001*i for i in range(40, 13200)]
    for dtd in dtds:
        func = import_dtd(dtd['name'])
        ax.plot(time, [func(t) for t in time], label=dtd['label'],
                c=dtd['color'], ls=dtd['line'], lw=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim((1e-12, 3e-8))
    ax.set_xlabel('Time [Gyr]')
    ax.set_ylabel(r'Normalized SN Ia Rate [$\rm{M}_\odot^{-1}$ yr$^{-1}$]')
    ax.legend(frameon=False, loc='upper right', fontsize=7)
    fig.savefig(paths.figures / 'delay_time_distributions.pdf')
    plt.close()


def import_dtd(name, zone=0):
    output = PARENT_DIR / ('%s.vice' % name)
    attributes_path = output / ('zone%s.vice' % zone) / 'attributes'
    with open(attributes_path / 'RIa.obj', 'rb') as f:
        RIa = pickle.load(f)
    return RIa


if __name__ == '__main__':
    main()
