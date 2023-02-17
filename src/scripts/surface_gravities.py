import numpy as np
import matplotlib.pyplot as plt
from utils import import_apogee, select_giants, get_bin_centers, \
    discrete_colormap, setup_discrete_colorbar, get_color_list
from ofe_feh_apogee import apogee_region
from _globals import GALR_BINS, ABSZ_BINS
import paths

BIN_WIDTH = 0.1

def main(cmap_name='plasma_r'):
    data = import_apogee()
    fig, axs = setup_axes()
    cmap, norm = discrete_colormap(cmap_name, GALR_BINS)
    cax = setup_discrete_colorbar(fig, cmap, norm, 
                                  label=r'Galactocentric radius [kpc]')
    colors = get_color_list(cmap, GALR_BINS)
    plot_apogee_logg(data, axs, colors=colors, label=None)
    fig.savefig(paths.figures / 'surface_gravities.png', dpi=300)
    

def plot_apogee_logg(data, axs, colors=[], label='APOGEE', 
                     logg_bin_width=BIN_WIDTH, logg_lim=(-0.5, 4.),
                     absz_bins=ABSZ_BINS, galr_bins=GALR_BINS):
    """
    Plot surface gravity distributions from APOGEE data.
    
    Parameters
    ----------
    data : pandas.DataFrame
        Combined data from APOGEE and astroNN.
    axs : list of matplotlib.axes.Axes
        Axes on which to plot the age distributions, the length of which must
        correspond to len(absz_bins)-1 (3 by default); usually a single
        column from a larger array of axes.
    colors : list, optional
        List of colors corresponding to galactic radius. If len(colors) !=
        len(galr_bins)-1, the default matplotlib color scheme will be used.
        The default is [].
    label : str, optional
        Axis column label. The default is 'astroNN'.
    age_bin_width : float, optional
        Width of age bins in Gyr. The default is 1.
    max_age : float, optional
        Maximum age in Gyr to include. The default is 14.
    absz_bins : list, optional
        Bin edges of galactic z-height in kpc. The default is [0, 0.5, 1, 2].
    galr_bins : list, optional
        Bin edges of galactic radius in kpc. The default is
        [3, 5, 7, 9, 11, 13, 15].
    """
    if len(colors) != len(galr_bins) - 1:
        colors = [None for galr in galr_bins[:-1]]
    if len(axs) == len(absz_bins) - 1:
        for i, ax in enumerate(axs.flatten()):
            absz_lim = absz_bins[-(i+2):len(absz_bins)-i]
            for j in range(len(galr_bins)-1):
                galr_lim = galr_bins[j:j+2]
                subset = apogee_region(data, galr_lim, absz_lim)
                logg_bins = np.arange(logg_lim[0], 
                                      logg_lim[1] + logg_bin_width, 
                                      logg_bin_width)
                logg_hist, _ = np.histogram(subset['LOGG'], 
                                            bins=logg_bins, density=True)
                ax.plot(get_bin_centers(logg_bins), logg_hist, 
                        color=colors[j], linewidth=1)
        axs[0].set_title(label)
    else:
        raise ValueError('Mismatch between axes and z-height bins.')
        
def setup_axes():
    fig, axs = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(3.25, 6))
    axs[-1].set_xlabel(r'$\log(g)$')
    axs[0].set_ylabel('PDF')
    return fig, axs

if __name__ == '__main__':
    main()
