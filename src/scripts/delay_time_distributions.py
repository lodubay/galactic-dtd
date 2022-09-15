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

DELAY = 0.04 # minimum Ia delay timie in Gyr

class styles:
    plaw = {
        'func': dtds.powerlaw(slope=-1.1, tmin=DELAY),
        'label': r'Power-Law ($\alpha=-1.1$)',
        'color': paultol.vibrant.colors[4],
        'line': '-',
    }
    plaw_steep = {
        'func': dtds.powerlaw(slope=-1.4, tmin=DELAY),
        'label': r'Power-Law ($\alpha=-1.4$)',
        'color': paultol.bright.colors[5],
        'line': '--',
    }
    plateau = {
        'func': dtds.plateau(width=0.3, slope=-1.1, tmin=DELAY),
        'label': r'Power-Law with 300 Myr plateau',
        'color': paultol.vibrant.colors[0],
        'line': '-.',
    }
    plateau_long = {
        'func': dtds.plateau(width=1., slope=-1.1, tmin=DELAY),
        'label': r'Plateau ($W=1$ Gyr)',
        'color': paultol.bright.colors[3],
        'line': '-.',
    }
    exp = {
        'func': dtds.exponential(timescale=1.5, tmin=DELAY),
        'label': r'Exponential ($\tau=1.5$ Gyr)',
        'color': paultol.bright.colors[4],
        'line': '--',
    }
    exp_long = {
        'func': dtds.exponential(timescale=3, tmin=DELAY),
        'label': r'Exponential ($\tau=3$ Gyr)',
        'color': paultol.vibrant.colors[1],
        'line': '--',
    }
    prompt = {
        'func': dtds.prompt(peak=0.05, stdev=0.015, timescale=3, tmin=DELAY),
        'label': r'Exponential with prompt component',
        'color': paultol.vibrant.colors[2],
        'line': ':',
    }

# distributions = [styles.prompt, styles.plaw_steep, styles.plaw, styles.plateau, 
#                  styles.exp, styles.exp_long, styles.plateau_long]
distributions = [styles.plaw, styles.plateau, styles.exp_long, styles.prompt]

def main():
    fig, ax = plt.subplots(figsize=(3.25, 3.25), tight_layout=True)
    time = [0.001*i for i in range(40, 13200)]
    for dtd in distributions:
        func = dtd['func']
        ax.plot([t * 1e9 for t in time], [func(t) for t in time], 
                label=dtd['label'], c=dtd['color'], ls=dtd['line'], lw=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim((1e-12, 5e-8))
    ax.set_xlabel('Time After Starburst [yr]')
    ax.set_ylabel(r'Normalized SN Ia Rate [$\rm{M}_\odot^{-1}$ yr$^{-1}$]')
    ax.legend(frameon=False, loc='upper right', fontsize=8, handlelength=1.25)
    fig.savefig(paths.figures / 'delay_time_distributions.png')
    plt.close()

if __name__ == '__main__':
    main()
