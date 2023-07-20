"""
Compare [O/Fe]-[Fe/H] plots for the Solar annulus for VICE outputs with
different star formation histories.
"""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import vice
from multizone_stars import MultizoneStars
from scatter_plot_grid import setup_colorbar
from _globals import ZONE_WIDTH, ONE_COLUMN_WIDTH
import paths
from ofe_feh_dtd import apogee_contours
from apogee_tools import import_apogee

FEH_LIM = (-1.3, 0.6)
OFE_LIM = (-0.1, 0.55)
GALR_LIM = (7, 9)
ABSZ_LIM = (0, 2)

SFH_LABELS = ['Inside-out', 'Late-burst', 'Early-burst', 'Two-infall']
DTD_MODEL = 'exponential_timescale15'

CMAP_NAME = 'winter'

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True,
                            figsize=(ONE_COLUMN_WIDTH, 0.9*ONE_COLUMN_WIDTH))
    plt.subplots_adjust(top=0.98, right=0.93, wspace=0., hspace=0.)
    cbar = setup_colorbar(fig, cmap=CMAP_NAME, vmin=0, vmax=15.5, width=0.04,
                          label=r'Birth $R_{\rm{gal}}$ [kpc]', pad=0.02)
    cbar.ax.yaxis.set_major_locator(MultipleLocator(2))
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
    apogee_data = import_apogee()
    
    for label, ax in zip(SFH_LABELS, axs.flatten()):
        sfh = ''.join(label.split('-')).lower()
        output_name = '/'.join(['gaussian', sfh, DTD_MODEL, 'diskmodel'])
        # Import multioutput stars data
        mzs = MultizoneStars.from_output(output_name)
        mzs.model_uncertainty(inplace=True)
        mzs.region(GALR_LIM, ABSZ_LIM, inplace=True)
        # Plot sample of star particle abundances
        mzs.scatter_plot(ax, '[fe/h]', '[o/fe]', color='galr_origin',
                         cmap=CMAP_NAME, norm=cbar.norm, markersize=0.1)
        apogee_contours(ax, apogee_data, GALR_LIM, ABSZ_LIM)
        # Plot abundance tracks
        zone = int(0.5 * (GALR_LIM[0] + GALR_LIM[1]) / ZONE_WIDTH)
        # Import post-processed output for the given annulus
        zone_path = str(mzs.fullpath / ('zone%d' % zone))
        hist = vice.history(zone_path)
        ax.plot(hist['[fe/h]'], hist['[o/fe]'], c='k', ls='-', 
                linewidth=0.5)
        # Label axis
        ax.text(0.9, 0.9, label, va='top', ha='right', transform=ax.transAxes)
    
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
    
    plt.savefig(paths.figures / 'ofe_feh_sfh.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
