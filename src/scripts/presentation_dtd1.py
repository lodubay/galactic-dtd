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

# Custom presentation plot settings
plt.style.use('presentation.mplstyle')

DELAY = 0.04 # minimum Ia delay time in Gyr
DTD_LIST = [
    dtds.powerlaw(slope=-1.1, tmin=DELAY),
    # dtds.prompt(peak=0.05, stdev=0.015, timescale=5, tmin=DELAY),
    dtds.greggio05_single(tmin=DELAY),
    dtds.greggio05_approximate.from_defaults('wide'),
]
LABELS = [
    "Cosmic Ia rate\n(Maoz & Graur 2017)",
    # "Two-population model\n(Mannucci+ 2006)",
    "Analytic single degenerate\n(Greggio 2005)",
    "Analytic double degenerate\n(Greggio 2005)",
]

def main():
    # Set up plot
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    plt.subplots_adjust(left=0.13, bottom=0.15, right=0.6, top=0.95)
    
    # for i in range(len(DTD_LIST)):
    #     func = DTD_LIST[i]
    #     ax.plot([t * 1e9 for t in time], [func(t) for t in time], 
    #             label=LABELS[i], ls=LINE_STYLES[i], lw=2)
    # Observational DTDs
    plaw, = plot_dtd(dtds.powerlaw(slope=-1.1, tmin=DELAY), 
                     ax, label='Power law', ls='-')
    # prompt, = plot_dtd(dtds.prompt(peak=0.05, stdev=0.015, timescale=5, tmin=DELAY),
    #                    ax, label='Two population', ls='--')
    # Theoretical DTDs
    gsd, = plot_dtd(dtds.greggio05_single(tmin=DELAY),
                    ax, label='SD', ls='--')
    gdd, = plot_dtd(dtds.greggio05_approximate.from_defaults('wide'),
                    ax, label='DD', ls=':')
        
    # Configure axes
    ax.tick_params(axis='x', pad=6)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim((4e7, 1.32e10))
    ax.set_ylim((5e-12, 1e-8))
    ax.set_xlabel('Years after star formation')
    ax.set_ylabel(r'Normalized Type Ia SN rate')
    
    # Legend
    # leg1 = fig.legend([plaw, prompt], 
    #                   ["Cosmic Ia rate\n(Maoz & Graur 2017)", 
    #                    "Two-population model\n(Mannucci+ 2006)"],
    #                   title="Observational DTDs", frameon=False, 
    #                   loc='upper left', bbox_to_anchor=(0.6, 0.95))
    leg1 = fig.legend([plaw], ["Cosmic Ia rate\n(Maoz & Graur 2017)"],
                      title="Simple Functional Forms", frameon=False, 
                      loc='upper left', bbox_to_anchor=(0.6, 0.95))
    leg1._legend_box.align = 'left'
    leg2 = fig.legend([gsd, gdd], ["Single degenerate", "Double degenerate"],
                      title="Analytic models\n(Greggio 2005)", frameon=False, 
                      loc='upper left', bbox_to_anchor=(0.6, 0.6),)
    leg2._legend_box.align = 'left'
    
    fig.savefig(paths.figures / 'presentation_dtd1.png')
    plt.close()
    

def plot_dtd(func, ax, label='', ls=None, c=None):
    time = [0.001*i for i in range(40, 13200)]
    out = ax.plot([t * 1e9 for t in time], [func(t) for t in time], 
                  label=label, ls=ls, lw=2, c=c)
    return out
    

if __name__ == '__main__':
    main()
