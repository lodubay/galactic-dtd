"""
Plot the Type Ia supernova delay time distributions (DTDs) as a function of time.
"""

import sys
import matplotlib.pyplot as plt
import paths
sys.path.append(str(paths.root))
from migration.src.simulations import dtds
from colormaps import paultol
plt.rcParams['axes.prop_cycle'] = plt.cycler('color', paultol.bright.colors)

class styles:
    plaw = {
        'name': 'powerlaw',
        'func': dtds.powerlaw(),
        'label': r'Power-Law ($\alpha=-1.1$)',
        'color': 'k',
        'line': '-',
    }
    plaw_steep = {
        'name': 'powerlaw_steep',
        'func': dtds.powerlaw(slope=-1.4),
        'label': r'Power-Law ($\alpha=-1.4$)',
        'color': paultol.bright.colors[5],
        'line': '--',
    }
    plaw_broken = {
        'name': 'powerlaw_broken',
        'func': dtds.powerlaw_broken(),
        'label': 'Broken Power-Law',
        'color': paultol.bright.colors[1],
        'line': '-.',
    }
    exp = {
        'name': 'exponential',
        'func': dtds.exponential(),
        'label': r'Exponential ($\tau=1.5$ Gyr)',
        'color': paultol.bright.colors[4],
        'line': '--',
    }
    exp_long = {
        'name': 'exponential_long',
        'func': dtds.exponential(timescale=3),
        'label': r'Exponential ($\tau=3$ Gyr)',
        'color': paultol.bright.colors[0],
        'line': '-',
    }
    bimodal = {
        'name': 'bimodal',
        'func': dtds.bimodal(),
        'label': 'Bimodal',
        'color': paultol.bright.colors[2],
        'line': '-.',
    }

distributions = [styles.bimodal, styles.plaw_steep, styles.plaw,
                 styles.plaw_broken, styles.exp, styles.exp_long]

def main():
    fig, ax = plt.subplots(figsize=(3.25, 3.25), tight_layout=True)
    time = [0.001*i for i in range(40, 13200)]
    for dtd in distributions:
        func = dtd['func']
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

if __name__ == '__main__':
    main()
