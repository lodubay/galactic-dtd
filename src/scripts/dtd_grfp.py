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

DELAY = 0.04 # minimum Ia delay time in Gyr
DTD_LIST = [
    dtds.powerlaw(slope=-1.1, tmin=DELAY),
    dtds.greggio05_single(tmin=DELAY),
    dtds.greggio05_approximate.from_defaults('wide'),
    dtds.prompt(peak=0.05, stdev=0.015, timescale=5, tmin=DELAY)
]
LABELS = [
    "Cosmic Ia rate (MG17)",
    "Single degenerate (G05)",
    "Double degenerate (G05)",
    "Two-population (MDP06)"
]
# COLORS = [
    # paultol.muted.colors[0],
    # paultol.muted.colors[1],
    # paultol.muted.colors[3],
    # paultol.muted.colors[6]
# ]
LINE_STYLES = [
    '-',
    '--',
    ':',
    '-.'
]

def main():
    fig, ax = plt.subplots(figsize=(3.25, 2.5), tight_layout=True)
    time = [0.001*i for i in range(40, 13200)]
    for i in range(len(DTD_LIST)):
        func = DTD_LIST[i]
        ax.plot([t * 1e9 for t in time], [func(t) for t in time], 
                label=LABELS[i], ls=LINE_STYLES[i], lw=1.5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim((5e-12, 3e-8))
    ax.set_xlabel('Years after starburst')
    ax.set_ylabel(r'Normalized Type Ia supernova rate')
    ax.legend(frameon=False, loc='upper right', fontsize=8, handlelength=1.75)
    fig.savefig(paths.figures / 'dtd_grfp.png')
    plt.close()

if __name__ == '__main__':
    main()
