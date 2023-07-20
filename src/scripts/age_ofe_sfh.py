"""
Compare [O/Fe]-Age plots for the Solar annulus for VICE outputs with
different star formation histories.
"""

from utils import multioutput_to_pandas, filter_multioutput_stars, model_uncertainty
from apogee_tools import import_apogee, apogee_region
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from age_ofe import plot_vice_medians, plot_astroNN_medians
from scatter_plot_grid import setup_colorbar, plot_vice_sample
from _globals import ZONE_WIDTH, ONE_COLUMN_WIDTH, MAX_SF_RADIUS
import paths

SFH_LIST = ['insideout', 'lateburst', 'conroy22_JW20yields', 'twoinfall']
DTD = 'powerlaw_slope11'
LABEL_LIST = ['Inside-Out', 'Late-Burst', 'Early-Burst', 'Two-Infall']
AGE_COL = 'LATENT_AGE'
AGE_LABEL = 'Leung et al. (2023)'
AGE_LIM = (0.3, 20)
OFE_LIM = (-0.15, 0.55)
CMAP = 'winter'
GALR_LIM = (7, 9)
ABSZ_LIM = (0, 0.5)

def main():
    apogee_data = import_apogee()
    apogee_subset = apogee_region(apogee_data, GALR_LIM, ABSZ_LIM)
    
    fig, axs = plt.subplots(2, 2, figsize=(ONE_COLUMN_WIDTH, 3.),
                            sharex=True, sharey=True)
    plt.subplots_adjust(right=0.94, left=0.1, bottom=0.1, top=0.98,
                        wspace=0, hspace=0)
    cbar = setup_colorbar(fig, cmap=CMAP, vmin=0, vmax=MAX_SF_RADIUS, 
                          label=r'Birth $R_{\rm{Gal}}$ [kpc]', pad=0.02, width=0.04)
    cbar.ax.yaxis.set_major_locator(MultipleLocator(2))
    cbar.ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    
    for ax, sfh, label in zip(axs.flatten(), SFH_LIST, LABEL_LIST):
        ax.text(0.07, 0.88, label, transform=ax.transAxes, fontsize=7)
        output_name = '/'.join(['diffusion', sfh, DTD])
        vice_stars = multioutput_to_pandas(output_name)
        # [O/Fe] uncertainty 
        ofe_err = apogee_data['O_FE_ERR'].median()
        vice_stars['[o/fe]'] = model_uncertainty(vice_stars['[o/fe]'], ofe_err)
        # Age uncertainty
        age_err = apogee_data['LOG_LATENT_AGE_ERR'].median()
        vice_stars['age'] = model_uncertainty(vice_stars['age'], age_err,
                                              how='logarithmic')
        vice_subset = filter_multioutput_stars(vice_stars, GALR_LIM, ABSZ_LIM,
                                                zone_width=ZONE_WIDTH)
        plot_vice_sample(ax, vice_subset, 'age', '[o/fe]', 
                         cmap=CMAP, norm=cbar.norm)
        plot_astroNN_medians(ax, apogee_subset, age_col=AGE_COL, 
                             label=AGE_LABEL, 
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
    
    fig.savefig(paths.figures / 'age_ofe_same_dtd.pdf', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
