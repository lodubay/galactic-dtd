"""
Compare [O/Fe]-[Fe/H] plots for the Solar annulus for VICE outputs with
different star formation histories and delay time distributions
"""

from tqdm import tqdm
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
# from matplotlib.colors import Normalize
# from matplotlib.cm import ScalarMappable
import vice
from multizone_stars import MultizoneStars
from scatter_plot_grid import setup_colorbar
from apogee_tools import import_apogee, gen_kde
from _globals import ZONE_WIDTH, TWO_COLUMN_WIDTH, MAX_SF_RADIUS
import paths

FEH_LIM = (-1.3, 0.7)
OFE_LIM = (-0.15, 0.55)
GALR_LIM = (7, 9)
ABSZ_LIM = (0, 2)

SFH_LIST = ['insideout', 'lateburst', 'earlyburst', 'twoinfall']
SFH_LABELS = ['Inside-out', 'Late-burst', 'Early-burst', 'Two-infall']
DTD_LIST = ['prompt', 
            'powerlaw_slope11', 
            'exponential_timescale15', 
            'plateau_width10', 
            'triple']
DTD_LABELS = ['Prompt', 
              'Power-law\n($\\alpha=-1.1$)', 
              'Exponential\n($\\tau=1.5$ Gyr)',
              'Plateau\n($W=1$ Gyr)',
              'Triple-system']

CMAP_NAME = 'winter'

def main():
    # Set up figure
    plt.style.use(paths.styles / 'paper.mplstyle')
    width = TWO_COLUMN_WIDTH
    fig, axs = plt.subplots(5, 4, sharex=True, sharey=True,
                            figsize=(width, 5/4*width))
    plt.subplots_adjust(top=0.95, wspace=0., hspace=0.)
    cbar = setup_colorbar(fig, cmap=CMAP_NAME, vmin=0, vmax=MAX_SF_RADIUS, 
                          width=0.04, pad=0.02, 
                          label=r'Birth $R_{\rm{gal}}$ [kpc]')
    cbar.ax.yaxis.set_major_locator(MultipleLocator(2))
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
    apogee_data = import_apogee()
    
    with tqdm(total=len(SFH_LIST) * len(DTD_LIST)) as t:
        for i, dtd in enumerate(DTD_LIST):
            for j, sfh in enumerate(SFH_LIST):
                output_name = '/'.join(['gaussian', sfh, dtd, 'diskmodel'])
                # Import multioutput stars data
                mzs = MultizoneStars.from_output(output_name)
                mzs.model_uncertainty(inplace=True)
                mzs.region(GALR_LIM, ABSZ_LIM, inplace=True)
                # Plot sample of star particle abundances
                mzs.scatter_plot(axs[i,j], '[fe/h]', '[o/fe]', color='galr_origin',
                                 cmap=CMAP_NAME, norm=cbar.norm, markersize=0.1)
                # Plot abundance tracks
                zone = int(0.5 * (GALR_LIM[0] + GALR_LIM[1]) / ZONE_WIDTH)
                zone_path = str(mzs.fullpath / ('zone%d' % zone))
                hist = vice.history(zone_path)
                axs[i,j].plot(hist['[fe/h]'], hist['[o/fe]'], c='k', ls='-', 
                        linewidth=0.5)
                # Plot APOGEE contours
                xx, yy, logz = gen_kde(apogee_data, bandwidth=0.02,
                                       galr_lim=GALR_LIM, absz_lim=ABSZ_LIM)
                # scale the linear density to the max value
                scaled_density = np.exp(logz) / np.max(np.exp(logz))
                # contour levels at 1, and 2 sigma
                levels = np.exp(-0.5 * np.array([2, 1])**2)
                axs[i,j].contour(xx, yy, scaled_density, levels, colors='r',
                                 linewidths=0.5, linestyles=['--', '-'])
                t.update()
    
    # Set x-axis ticks
    axs[0,0].xaxis.set_major_locator(MultipleLocator(0.5))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.1))
    # Set y-axis ticks
    axs[0,0].yaxis.set_major_locator(MultipleLocator(0.2))
    axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    # Set axis limits
    axs[0,0].set_xlim(FEH_LIM)
    axs[0,0].set_ylim(OFE_LIM)
    # Set axis labels
    for ax in axs[-1]:
        ax.set_xlabel('[Fe/H]')
    for ax in axs[:,0]:
        ax.set_ylabel('[O/Fe]')
    # Axis titles
    for j, ax in enumerate(axs[0]):
        ax.set_title(SFH_LABELS[j])
    # Axis row labels
    for i, ax in enumerate(axs[:,0]):
        ax.text(0.9, 0.9, DTD_LABELS[i], va='top', ha='right', 
                transform=axs[i,j].transAxes)
    
    plt.savefig(paths.figures / 'ofe_feh_summary.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
