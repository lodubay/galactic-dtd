"""
Plot a single panel of [O/Fe] vs [Fe/H] for a galactic region.
"""

# import sys
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import vice
from multizone_stars import MultizoneStars
from scatter_plot_grid import setup_axes, setup_colorbar
from _globals import ABSZ_BINS, ZONE_WIDTH, ONE_COLUMN_WIDTH
import paths

FEH_LIM = (-1.3, 0.7)
OFE_LIM = (-0.15, 0.55)
GALR_LIM = (7, 9)
ABSZ_LIM = (0, 2)

def main(output_name, cmap='winter', uncertainties=True, tracks=True,
         data_dir=paths.data/'migration'):
    # Import multioutput stars data
    mzs = MultizoneStars.from_output(output_name, data_dir=data_dir)
    # Model observational uncertainties
    if uncertainties:
        mzs.model_uncertainty(inplace=True)
    fig, ax = plt.subplots(figsize=(ONE_COLUMN_WIDTH, 0.8*ONE_COLUMN_WIDTH))
    cbar = setup_colorbar(fig, cmap=cmap, vmin=0, vmax=15.5, width=0.04,
                          label=r'Birth $R_{\rm{Gal}}$ [kpc]', pad=0.02)
    cbar.ax.yaxis.set_major_locator(MultipleLocator(2))
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
    subset = mzs.region(GALR_LIM, ABSZ_LIM)
    subset.scatter_plot(ax, '[fe/h]', '[o/fe]', color='galr_origin',
                        cmap=cmap, norm=cbar.norm, markersize=0.5)
    if tracks:
        zone = int(0.5 * (GALR_LIM[0] + GALR_LIM[1]) / ZONE_WIDTH)
        # Import post-processed output for the given annulus
        zone_path = str(mzs.fullpath / ('zone%d' % zone))
        zone_path = zone_path.replace('diffusion', 'post-process')
        zone_path = zone_path.replace('gaussian', 'post-process')
        hist = vice.history(zone_path)
        ax.plot(hist['[fe/h]'], hist['[o/fe]'], c='k', ls='-', 
                linewidth=0.5)
    
    # Set x-axis ticks
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    # Set y-axis ticks
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    # Set labels
    ax.set_xlabel('[Fe/H]')
    ax.set_ylabel('[O/Fe]')
    ax.set_title(r'$%.1f < R_{\rm Gal} < %.1f$ kpc, $%.1f < |z| < %.1f$ kpc' % (
        GALR_LIM[0], GALR_LIM[1], ABSZ_LIM[0], ABSZ_LIM[1]))
    # Set axis limits
    ax.set_xlim(FEH_LIM)
    ax.set_ylim(OFE_LIM)
    
    fname = output_name.replace('/', '_')
    if uncertainties: fname += '_errors'
    fname += '.png'
    plt.savefig(paths.debug / 'ofe_feh_single' / fname, dpi=300)
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='ofe_feh.py',
        description='Generate a grid of [O/Fe] vs [Fe/H] scatterplots ' + \
            'from a VICE multizone run.'
    )
    parser.add_argument('output_name', metavar='NAME',
                        help='Name of VICE multizone output')
    parser.add_argument('-u', '--uncertainties', action='store_true',
                        help='Model APOGEE uncertainties in VICE output')
    parser.add_argument('-t', '--tracks', action='store_true',
                        help='Plot ISM tracks in addition to stellar abundances')
    parser.add_argument('-c', '--cmap', metavar='COLORMAP', type=str,
                        default='winter',
                        help='Name of colormap for color-coding VICE ' + \
                             'output (default: winter)')
    args = parser.parse_args()
    main(**vars(args))
