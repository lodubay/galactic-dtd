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
import numpy as np
import vice
from multizone_stars import MultizoneStars
from apogee_tools import import_apogee, apogee_region
from scatter_plot_grid import setup_axes, setup_colorbar
from utils import kde2D
from _globals import ABSZ_BINS, ZONE_WIDTH
import paths

FEH_LIM = (-1.3, 0.7)
OFE_LIM = (-0.15, 0.55)
GALR_BINS = [3, 5, 7, 9, 11, 13]

def main(output_name, cmap='winter', uncertainties=True, tracks=True, 
         contours=True, data_dir=paths.data/'migration'):
    # Import multioutput stars data
    mzs = MultizoneStars.from_output(output_name, data_dir=data_dir)
    # Import APOGEE data
    if contours:
        apogee_data = import_apogee()
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
    
    fname = output_name.replace('/', '_')
    if uncertainties: fname += '_errors'
    fname += '.png'
    plt.savefig(paths.debug / 'ofe_feh_grid' / fname, dpi=300)
    plt.close()
    

def plot_contours(ax, data, bandwidth=0.02, absz_lim=(0, 5), galr_lim=(0, 20),
                  overwrite=False, **kwargs):
    """
    Plot APOGEE density contours in the specified region.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object on which to draw contours.
    data : pandas.DataFrame
        APOGEE data containing columns 'FE_H' and 'O_FE'.
    bandwidth : float
        Kernel density estimate bandwidth. A larger number will produce
        smoother contour lines. The default is 0.02.
    absz_lim : tuple
        Limits on absolute Galactic z-height in kpc. The default is (0, 5).
    galr_lim : tuple
        Limits on Galactocentric radius in kpc. The default is (0, 20).
    overwrite : bool
        If True, force re-generate the 2D KDE and save the output.
    **kwargs passed to matplotlib.axes.Axes.contour().
    
    Returns
    -------
    contours: matplotlib.contour.QuadContourSet
    """
    return contours


def gen_kde(data, bandwidth=0.02, absz_lim=(0, 5), galr_lim=(0, 20), 
            overwrite=False):
    """
    Generate kernel density estimate (KDE) of APOGEE data, or import previously
    saved KDE if it already exists.
    
    Parameters
    ----------
    data : pandas.DataFrame
        APOGEE data containing columns 'FE_H' and 'O_FE'.
    bandwidth : float
        Kernel density estimate bandwidth. A larger number will produce
        smoother contour lines. The default is 0.02.
    absz_lim : tuple
        Limits on absolute Galactic z-height in kpc. The default is (0, 5).
    galr_lim : tuple
        Limits on Galactocentric radius in kpc. The default is (0, 20).
    overwrite : bool
        If True, force re-generate the 2D KDE and save the output.
    
    Returns
    -------
    xx, yy, logz: tuple of numpy.array
        Outputs of kde2D()
    """    
    # Path to save 2D KDE for faster plot times
    path = kde_path(galr_lim, absz_lim, savedir=paths.data/'APOGEE/kde/ofe_feh/')
    if path.exists() and not overwrite:
        xx, yy, logz = read_kde(path)
    else:
        # Limit to specified Galactic region
        subset = apogee_region(data, galr_lim, absz_lim)
        subset = subset.copy().dropna(axis=0, subset=['FE_H', 'O_FE'])
        xx, yy, logz = kde2D(subset['FE_H'], subset['O_FE'], bandwidth)
        save_kde(xx, yy, logz, path)
    return xx, yy, logz


def read_kde(path):
    """
    Read a text file generated by save_kde()
    """
    arr2d = np.genfromtxt(path)
    nrows = int(arr2d.shape[0]/3)
    xx = arr2d[:nrows]
    yy = arr2d[nrows:2*nrows]
    logz = arr2d[2*nrows:]
    return xx, yy, logz


def save_kde(xx, yy, logz, path):
    """
    Generate a text file containing the KDE of the given region along with its
    corresponding x and y coordinates.
    """
    with open(path, 'w') as f:
        for arr in [xx, yy, logz]:
            f.write('#\n')
            np.savetxt(f, arr)


def kde_path(galr_lim, absz_lim, savedir=paths.data/'APOGEE/kde'):
    """
    Generate file name for the KDE of the given region.
    """
    filename = 'r%s-%s_z%s-%s.dat' % (galr_lim[0], galr_lim[1],
                                      absz_lim[0], absz_lim[1])
    return Path(savedir) / filename


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
    parser.add_argument('-c', '--contours', action='store_true',
                        help='Plot contour lines from APOGEE data')
    parser.add_argument('--cmap', metavar='COLORMAP', type=str,
                        default='winter',
                        help='Name of colormap for color-coding VICE ' + \
                             'output (default: winter)')
    args = parser.parse_args()
    main(**vars(args))
