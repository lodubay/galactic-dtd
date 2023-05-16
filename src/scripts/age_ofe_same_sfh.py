"""
Compare plots of [O/Fe] vs age in a single Galactic region for various DTDs.
"""

from utils import import_apogee, apogee_region, multioutput_to_pandas, \
    filter_multioutput_stars, model_uncertainty
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from age_ofe import plot_vice_medians, plot_astroNN_medians
from scatter_plot_grid import setup_colorbar, plot_vice_sample
from _globals import ZONE_WIDTH, ONE_COLUMN_WIDTH
import paths

SFH = 'insideout'
DTD_LIST = ['powerlaw_slope11', 'exponential_timescale15', 
            'plateau_width1000_slope11', 'triple_delay040']
LABEL_LIST = [r'Power law ($\alpha=-1.1$)', r'Exponential ($\tau=1.5$ Gyr)',
              r'Plateau ($W=1.0$ Gyr)', 'Triple system evolution']
AGE_LABELS = {'ASTRONN_AGE': 'Mackereth et al. 2019',
              'LATENT_AGE': 'Leung et al. 2023'}
AGE_LIM = (0.3, 20)
OFE_LIM = (-0.15, 0.55)

def main(cmap='winter', galr_lim=(7, 9), absz_lim=(0, 0.5), age_col='LATENT_AGE'):
    apogee_data = import_apogee()
    apogee_subset = apogee_region(apogee_data, galr_lim, absz_lim)
    
    fig, axs = plt.subplots(2, 2, figsize=(3.25, 3.),
                            sharex=True, sharey=True)
    plt.subplots_adjust(right=0.94, left=0.1, bottom=0.1, top=0.98,
                        wspace=0, hspace=0)
    cbar = setup_colorbar(fig, cmap=cmap, vmin=0, vmax=15.5, 
                          label=r'Birth $R_{\rm{Gal}}$ [kpc]', pad=0.02, width=0.04)
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
    for ax, dtd, label in zip(axs.flatten(), DTD_LIST, LABEL_LIST):
        ax.text(0.07, 0.88, label, transform=ax.transAxes, fontsize=7)
        output_name = '/'.join(['diffusion', SFH, dtd])
        vice_stars = multioutput_to_pandas(output_name)
        # [O/Fe] uncertainty 
        ofe_err = apogee_data['O_FE_ERR'].median()
        vice_stars['[o/fe]'] = model_uncertainty(vice_stars['[o/fe]'], ofe_err)
        # Age uncertainty
        age_err = apogee_data['LOG_LATENT_AGE_ERR'].median()
        vice_stars['age'] = model_uncertainty(vice_stars['age'], age_err,
                                              how='logarithmic')
        vice_subset = filter_multioutput_stars(vice_stars, galr_lim, absz_lim,
                                                zone_width=ZONE_WIDTH)
        plot_vice_sample(ax, vice_subset, 'age', '[o/fe]', 
                         cmap=cmap, norm=cbar.norm)
        plot_astroNN_medians(ax, apogee_subset, age_col=age_col, 
                              label=AGE_LABELS[age_col], 
                              plot_low_count_bins=False)
        plot_vice_medians(ax, vice_subset, label='VICE',
                          plot_low_mass_bins=False)
    
    # Set x-axis scale and ticks
    axs[0,0].set_xlim(AGE_LIM)
    axs[0,0].set_xscale('log')
    axs[0,0].xaxis.set_major_formatter(FormatStrFormatter('%d'))
    
    # Set y-axis ticks
    axs[0,0].set_ylim(OFE_LIM)
    axs[0,0].yaxis.set_major_locator(MultipleLocator(0.2))
    axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    
    # Axis labels
    for ax in axs[-1]:
        ax.set_xlabel('Age [Gyr]')
    for i, ax in enumerate(axs[:,0]):
        ax.set_ylabel('[O/Fe]', labelpad=2)
        
    # Legend
    axs[0,1].legend(loc='upper left', frameon=False, 
                    bbox_to_anchor=(0.02, 0.89), handlelength=0.7)
    
    fig.savefig(paths.figures / 'age_ofe_same_sfh.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
