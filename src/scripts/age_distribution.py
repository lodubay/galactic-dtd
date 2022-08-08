"""
This script plots the distribution of stellar ages as predicted by VICE.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import vice
import paths
from utils import multioutput_to_pandas, filter_multioutput_stars, import_astroNN
from _globals import DT, GALR_BINS, ZONE_WIDTH, ABSZ_BINS, END_TIME
from mdf_9panel import setup_colorbar, discrete_colormap, box_smooth
from ofe_feh_apogee import apogee_region

global AGE_LIM
global SMOOTH_WIDTH

AGE_LIM = (0, 14)
SMOOTH_WIDTH = 1
BIN_WIDTH = 1

def main(evolution, RIa, cmap_name='plasma_r'):
    output = '%s/%s' % (evolution, RIa)
    stars_pp = multioutput_to_pandas(paths.data / 'migration' / 'post-process' / output)
    stars_diff = multioutput_to_pandas(paths.data / 'migration' / 'diffusion' / output)
    astroNN_data = import_astroNN()
    # diff_subset = filter_multioutput_stars(stars_diff,
    #                                         (0, 15), (0, 5),
    #                                         ZONE_WIDTH, min_mass=0)
    # ages, n_stars = age_distribution(diff_subset)
    # ages, n_stars = age_distribution(stars_diff)
    # plt.plot(ages, n_stars)
    # plt.xlabel('Age [Gyr]')
    # plt.ylabel('dN/dAge')
    # plt.savefig(paths.figures / 'adf.png', dpi=300)
    # plt.show()
    
    fig, axs = setup_axes(xlim=AGE_LIM)
    cmap, norm = discrete_colormap(cmap_name, GALR_BINS)
    cax = setup_colorbar(fig, cmap, norm, label=r'Galactocentric radius [kpc]')
    # Define color scheme
    rmin, rmax = GALR_BINS[0], GALR_BINS[-2]
    colors = cmap([(r-rmin)/(rmax-rmin) for r in GALR_BINS[:-1]])
    for i in range(len(ABSZ_BINS)-1):
        absz_lim = ABSZ_BINS[-(i+2):len(ABSZ_BINS)-i]
        axs[i,0].set_ylabel(r'$|z| = %s - %s$' % tuple(absz_lim))
        for j in range(len(GALR_BINS)-1):
            galr_lim = GALR_BINS[j:j+2]
            # Plot VICE post-process in left panels
            # TODO bin in 1 Gyr age bins instead
            # TODO bin all 10+ Gyr stars in one bin
            pp_subset = filter_multioutput_stars(stars_pp,
                                                 galr_lim, absz_lim,
                                                 ZONE_WIDTH, min_mass=0)
            dndt, bins = age_distribution(pp_subset)
            adf_smooth = box_smooth(dndt, bins, SMOOTH_WIDTH)
            axs[i,0].plot(bins[:-1], adf_smooth, color=colors[j], linewidth=1)
            
            # Plot VICE diffusion in center panels
            diff_subset = filter_multioutput_stars(stars_diff,
                                                   galr_lim, absz_lim,
                                                   ZONE_WIDTH, min_mass=0)
            dndt, bins = age_distribution(diff_subset)
            adf_smooth = box_smooth(dndt, bins, SMOOTH_WIDTH)
            axs[i,1].plot(bins[:-1], adf_smooth, color=colors[j], linewidth=1)

            # Plot APOGEE in right panels
            apogee_subset = apogee_region(astroNN_data, galr_lim, absz_lim)
            apogee_adf, _ = np.histogram(apogee_subset['ASTRONN_AGE'], bins=bins,
                                         density=True)
            apogee_smooth = box_smooth(apogee_adf, bins, SMOOTH_WIDTH)
            axs[i,2].plot(bins[:-1], apogee_smooth, color=colors[j], linewidth=1)
            
    for ax in axs[-1]:
        ax.set_xlabel('Age [Gyr]')
    axs[0,0].set_title('Post-process')
    axs[0,1].set_title('Diffusion')
    axs[0,2].set_title('astroNN')
    plt.savefig(paths.figures / ('adf.png'), dpi=300)
    plt.close()


def setup_axes(xlim=AGE_LIM):
    fig, axs = plt.subplots(3, 3, figsize=(4.5, 4),
                            sharex=True, sharey=True)
    fig.subplots_adjust(left=0.07, top=0.93, right=0.97, bottom=0.1,
                        wspace=0.07, hspace=0.)
    # axs[0,0].set_xlim((xlim[0]-0.09, xlim[1]+0.09))
    axs[0,0].set_xlim(xlim)
    axs[0,0].xaxis.set_major_locator(MultipleLocator(5))
    axs[0,0].xaxis.set_minor_locator(MultipleLocator(1))
    # axs[0,0].yaxis.set_major_locator(MultipleLocator(1))
    # axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.2))
    # Refine axes
    for ax in axs.flatten():
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # ax.spines['bottom'].set_bounds(xlim[0], xlim[1])
        ax.yaxis.set_ticks_position('none')
        ax.yaxis.set_ticklabels([])
        ax.patch.set_alpha(0)
        ax.tick_params(top=False, which='both')
    # Turn on y-axis spine and ticks for left-hand panels
    # for ax in axs[:,0]:
    #     ax.spines['left'].set_visible(True)
    #     ax.yaxis.set_ticks_position('left')
    #     ax.spines['left'].set_bounds(0, 2.6)
    return fig, axs


def age_distribution(stars, end_time=END_TIME, dt=DT, **kwargs):
    """
    Calculate the distribution of ages in a VICE multizone output.
    
    Parameters
    ----------
    stars : DataFrame
        Pandas conversion of a VICE multizone output
    end_time : float, optional
        Simulation end time in Gyr. The default is 13.2
    dt : float, optional
        Simulation time step in Gyr. The default is 0.01
    kwargs : dict, optional
        Keyword arguments passed to mean_stellar_mass
        
    Returns
    -------
    dndt : 1xn numpy.ndarray
        Fraction of stars in each age bin
    bins : 1x(n+1) numpy.ndarray
        Age bins derived from simulation timesteps
    """
    # Create dummy entries to count at least 0 mass at every age
    bins = np.arange(0, end_time + dt, dt)
    temp_df = pd.DataFrame({'age': bins[:-1], 'mass': np.zeros(bins[:-1].shape)})
    stars = pd.concat([stars, temp_df])
    stars['age'] = np.round(stars['age'], decimals=2)
    # Sum stellar mass in each timestep
    mass_total = stars.groupby(['age']).sum()['mass']
    ages = np.array(mass_total.index)
    mass_total = np.array(mass_total)
    # Calculate remaining stellar mass today
    mass_remaining = mass_total * (1 - np.array(
        [vice.cumulative_return_fraction(age) for age in ages]))
    # Average mass of a star of that particular age
    mass_average = np.array([mean_stellar_mass(age, **kwargs) for age in ages])
    # Number of stars in each age bin
    nstars = np.around(mass_remaining / mass_average)
    # Fraction of stars in each age bin
    dndt = nstars / (dt * nstars.sum())
    return dndt, bins


def mean_stellar_mass(age, imf=vice.imf.kroupa, mlr=vice.mlr.larson1974,
                      m_lower=0.08, m_upper=100, dm=0.01):
    """
    Calculate the mean mass of a stellar population of a given age.

    Parameters
    ----------
    age : float
        Stellar age in Gyr
    imf : <function>, optional
        Initial mass function which takes mass in solar masses as an argument.
        The default is vice.imf.kroupa
    mlr : <function>, optional
        Mass-lifetime relation which takes age in Gyr as an argument. The
        default is vice.mlr.larson1974
    m_lower : float, optional
        Lower mass limit on IMF in solar masses. The default is 0.08
    m_upper : float, optional
        Upper mass limit on IMF in solar masses. The default is 100
    dm : float, optional
        IMF integration step in solar masses. The default is 0.01

    Returns
    -------
    float
        Mean mass of stars with lifetime greater than or equal to the given age
        weighted by the IMF
    """
    m_max = min((mlr(age, which='age'), m_upper))
    masses = np.arange(m_lower, m_max + dm, dm)
    dndm = np.array([imf(m) for m in masses])
    weighted_mean = np.average(masses, weights=dndm)
    return weighted_mean


# def f_survive(age, mlr='larson1974', imf='kroupa', m_lower=0.08, m_upper=100,
#               dm=0.01):
#     """
#     Calculate the surviving mass fraction of a stellar population with a given
#     age.

#     Parameters
#     ----------
#     age : float
#         Age of the stellar population in Gyr
#     mlr : str
#         Mass-lifetime relation (MLR) to use. The default is 'larson1974'.
#     imf : str
#         Which IMF to use. Options are 'kroupa' or 'salpeter'.
#     m_lower : float
#         Lower limit of stellar mass. The default is 0.08 solar masses.
#     m_upper : float
#         Upper limit of stellar mass. The default is 100 solar masses.
#     dm : float
#         Integration step size in solar masses. The default is 0.01.

#     Returns
#     -------
#     float
#         Stellar population surviving mass fraction.
#     """
#     mlr_select = {
#         'larson1974': vice.mlr.larson1974,
#         'mm1989': vice.mm1989,
#         'pm1993': vice.mlr.pm1993,
#         'ka1997': vice.mlr.ka1997,
#         'hpt2000': vice.mlr.hpt2000,
#         'vincenzo2016': vice.mlr.vincenzo2016,
#         'powerlaw': vice.mlr.powerlaw
#     }
#     if mlr in mlr_select.keys():
#         # Mass of a star with a lifetime equal to age
#         mass = mlr_select[mlr](age, which='age')
#         m_arr = np.array(m_lower, min((m_upper, mass)), dm)
#         normal_imf = NormalIMF(which=imf, m_lower=m_lower, m_upper=m_upper,
#                                dm=dm)
#         f_survive = sum([normal_imf(m) for m in m_arr])
#         return f_survive
#     else:
#         raise ValueError('MLR must be in acceptable list.')


if __name__ == '__main__':
    main('insideout_johnson21', 'powerlaw_slope11_delay040')
