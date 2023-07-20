"""
Plot a grid of [O/Fe] vs [Fe/H] at varying Galactic radii and z-heights.
"""

# import sys
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
# from matplotlib.colors import Normalize
# from matplotlib.cm import ScalarMappable
import numpy as np
import vice
from multizone_stars import MultizoneStars
from apogee_tools import import_apogee, gen_kde
from scatter_plot_grid import setup_axes, setup_colorbar
# from utils import kde2D
from _globals import ABSZ_BINS, ZONE_WIDTH, MAX_SF_RADIUS
import paths

FEH_LIM = (-1.3, 0.7)
OFE_LIM = (-0.15, 0.55)
GALR_BINS = [3, 5, 7, 9, 11, 13]

def main(output_name, cmap='winter', uncertainties=True, tracks=True, 
         contours=True):
    plt.style.use(paths.styles / 'paper.mplstyle')
    # Import multioutput stars data
    mzs = MultizoneStars.from_output(output_name)
    # Import APOGEE data
    if contours:
        apogee_data = import_apogee()
    # Model observational uncertainties
    if uncertainties:
        mzs.model_uncertainty(inplace=True)
    fig, axs = setup_axes(xlim=FEH_LIM, ylim=OFE_LIM, xlabel='[Fe/H]', 
                          ylabel='[O/Fe]', row_label_pos=(0.33, 0.88),
                          title=output_name)
    cbar = setup_colorbar(fig, cmap=cmap, vmin=0, vmax=MAX_SF_RADIUS, 
                          label=r'Birth $R_{\rm{Gal}}$ [kpc]')
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
                zone_path = str(mzs.fullpath / ('zone%d' % zone))
                hist = vice.history(zone_path)
                ax.plot(hist['[fe/h]'], hist['[o/fe]'], c='k', ls='-', 
                        linewidth=0.5)
            if contours:
                xx, yy, logz = gen_kde(apogee_data, bandwidth=0.02,
                                       galr_lim=galr_lim, absz_lim=absz_lim)
                # scale the linear density to the max value
                scaled_density = np.exp(logz) / np.max(np.exp(logz))
                # contour levels at 1, 2, and 3 sigma
                levels = np.exp(-0.5 * np.array([3, 2, 1])**2)
                ax.contour(xx, yy, scaled_density, levels, colors='k',
                           linewidths=0.5, linestyles=[':', '--', '-'])
    
    # Set x-axis ticks
    axs[0,0].xaxis.set_major_locator(MultipleLocator(0.5))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.1))
    # Set y-axis ticks
    axs[0,0].yaxis.set_major_locator(MultipleLocator(0.2))
    axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    
    # Save
    fname = output_name.replace('diskmodel', 'ofe_feh_grid.png')
    fullpath = paths.figures / 'supplementary' / fname
    if not fullpath.parents[0].exists():
        fullpath.parents[0].mkdir(parents=True)
    plt.savefig(fullpath, dpi=300)
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='ofe_feh_grid.py',
        description='Generate a grid of [O/Fe] vs [Fe/H] scatterplots ' + \
            'from a VICE multizone run.'
    )
    parser.add_argument('output_name', metavar='NAME',
                        help='Name of VICE multizone output')
    parser.add_argument('-u', '--uncertainties', action='store_true',
                        help='Model APOGEE uncertainties in VICE output')
    parser.add_argument('-t', '--tracks', action='store_true',
                        help='Plot ISM tracks in addition to stellar abundances')
    parser.add_argument('-c', '--contours', action='store_true',
                        help='Plot contour lines from APOGEE data')
    parser.add_argument('--cmap', metavar='COLORMAP', type=str,
                        default='winter',
                        help='Name of colormap for color-coding VICE ' + \
                             'output (default: winter)')
    args = parser.parse_args()
    main(**vars(args))
