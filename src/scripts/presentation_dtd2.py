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

def main():
    # Set up plot
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    plt.subplots_adjust(left=0.13, bottom=0.15, right=0.6, top=0.95)
    
    plaw, = plot_dtd(dtds.powerlaw(slope=-1.1, tmin=DELAY), 
                     ax, label=r'Power law ($t^{-1.1}$)', 
                     ls='-', lw=4, c=paultol.bright.colors[0])
    plateau, = plot_dtd(dtds.plateau(width=0.3, slope=-1.1, tmin=DELAY),
                        ax, label='Plateau (300 Myr)', 
                        ls='--', lw=4, c=paultol.bright.colors[1])
    gdd, = plot_dtd(dtds.greggio05_approximate.from_defaults('close'),
                    ax, label='Double degenerate\n(Greggio 2005)', 
                    ls='-.', lw=2, c=paultol.bright.colors[2])
        
    # Configure axes
    ax.tick_params(axis='x', pad=6)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim((4e7, 1.32e10))
    ax.set_ylim((5e-12, 2e-8))
    ax.set_xlabel('Years after star formation')
    ax.set_ylabel(r'Normalized Type Ia SN rate')
    
    # Legend
    # leg1 = fig.legend([plaw, plateau], 
    #                   [r"Power law ($t^{-1.1}$)",
    #                    "Plateau (300 Myr)"],
    #                   title="Simple Functional Forms", frameon=False, 
    #                   loc='upper left', bbox_to_anchor=(0.6, 0.95))
    # leg1._legend_box.align = 'left'
    # leg2 = fig.legend([gsd, gdd], ["Single degenerate", "Double degenerate"],
    #                   title="Analytic models\n(Greggio 2005)", frameon=False, 
    #                   loc='upper left', bbox_to_anchor=(0.6, 0.6),)
    # leg2._legend_box.align = 'left'
    fig.legend(loc='upper left', bbox_to_anchor=(0.6, 0.95), frameon=False)
    
    fig.savefig(paths.figures / 'presentation_dtd2.png')
    plt.close()
    

def plot_dtd(func, ax, label='', ls=None, c=None, lw=2):
    time = [0.001*i for i in range(40, 13200)]
    out = ax.plot([t * 1e9 for t in time], [func(t) for t in time], 
                  label=label, ls=ls, lw=lw, c=c)
    return out
    

if __name__ == '__main__':
    main()
