"""
Plot a grid of [O/Fe] vs [Fe/H] at varying Galactic radii and z-heights.
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
from _globals import ABSZ_BINS, ZONE_WIDTH
import paths

FEH_LIM = (-1.3, 0.7)
OFE_LIM = (-0.15, 0.55)
GALR_BINS = [3, 5, 7, 9, 11, 13]

def main(output_name, cmap='winter', uncertainties=True, tracks=True,
         data_dir=paths.data/'migration'):
    # Import multioutput stars data
    mzs = MultizoneStars.from_output(output_name, data_dir=data_dir)
    # Model observational uncertainties
    if uncertainties:
        mzs.model_uncertainty(inplace=True)
    fig, axs = setup_axes(xlim=FEH_LIM, ylim=OFE_LIM, xlabel='[Fe/H]', 
                          ylabel='[O/Fe]')
    cbar = setup_colorbar(fig, cmap=cmap, vmin=0, vmax=15.5, 
                          label=r'Birth $R_{\rm{Gal}}$ [kpc]', pad=0.)
    cbar.ax.yaxis.set_major_locator(MultipleLocator(2))
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
        
    for i, row in enumerate(axs):
        absz_lim = (ABSZ_BINS[-(i+2)], ABSZ_BINS[-(i+1)])
        for j, ax in enumerate(row):
            galr_lim = (GALR_BINS[j], GALR_BINS[j+1])
            subset = mzs.region(galr_lim, absz_lim)
            subset.scatter_plot(ax, '[fe/h]', '[o/fe]', color='galr_origin',
                                cmap=cmap, norm=cbar.norm)
            if tracks:
                zone = int(0.5 * (galr_lim[0] + galr_lim[1]) / ZONE_WIDTH)
                # Import post-processed output for the given annulus
                zone_path = str(mzs.fullpath / ('zone%d' % zone))
                zone_path = zone_path.replace('diffusion', 'post-process')
                zone_path = zone_path.replace('gaussian', 'post-process')
                hist = vice.history(zone_path)
                ax.plot(hist['[fe/h]'], hist['[o/fe]'], c='k', ls='-', 
                        linewidth=0.5)
    
    # Set x-axis ticks
    axs[0,0].xaxis.set_major_locator(MultipleLocator(0.5))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.1))
    # Set y-axis ticks
    axs[0,0].yaxis.set_major_locator(MultipleLocator(0.2))
    axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    
    fname = output_name.replace('/', '_')
    if uncertainties: fname += '_errors'
    fname += '.png'
    plt.savefig(paths.debug / 'ofe_feh_grid' / fname, dpi=300)
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
