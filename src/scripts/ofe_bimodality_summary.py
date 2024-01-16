"""
This script plots a comparison of the [O/Fe] bimodality for several DTDs.
It takes the [O/Fe] distribution within a narrow slice of [Fe/H]
for each VICE multi-zone output, and runs a peak-finding algorithm to determine
whether the output exhibits alpha-element bimodality.
"""
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from multizone_stars import MultizoneStars
import paths
from utils import get_bin_centers, highlight_panels
from apogee_tools import import_apogee, apogee_region, apogee_mdf
from colormaps import paultol
from _globals import TWO_COLUMN_WIDTH

OFE_LIM = (-0.15, 0.55)
PROMINENCE = 0.1
FEH_BINS = [(-0.6, -0.4), (-0.4, -0.2)]
COLORS = [paultol.highcontrast.colors[0], paultol.highcontrast.colors[2]]
LINESTYLES = ['--', '-']
GALR_LIM = (7, 9)
ABSZ_LIM = (0, 2)
SMOOTH_WIDTH = 0.05

SFH_LIST = ['insideout', 'lateburst', 'earlyburst', 'twoinfall']
SFH_LABELS = ['Inside-out', 'Late-burst', 'Early-burst', 'Two-infall']
DTD_LIST = ['prompt', 
            'powerlaw_slope11', 
            'exponential_timescale15', 
            'plateau_width10', 
            'triple']
DTD_LABELS = ['Two-population', 
              'Power law\n($\\alpha=-1.1$)', 
              'Exponential\n($\\tau=1.5$ Gyr)',
              'Plateau\n($W=1.0$ Gyr)', 
              'Triple-system']

def main():
    plt.style.use(paths.styles / 'paper.mplstyle')
    figwidth = TWO_COLUMN_WIDTH
    title_size = plt.rcParams['figure.titlesize']
    # create subfigures for top and bottom rows
    # fig = plt.figure(layout='constrained', figsize=(figwidth, figwidth*0.4))
    # subfigs = fig.subfigures(2, 1, hspace=0.1)
    fig, axs = plt.subplots(2, 5, figsize=(figwidth, figwidth*0.45), 
                            sharex=True, sharey='row')
    fig.subplots_adjust(left=0.04, top=0.83, right=0.98, bottom=0.12,
                        wspace=0.07, hspace=0.52)
    # Remove spines and y-axis labels
    for ax in axs.flatten():
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('none')
        ax.yaxis.set_ticklabels([])
        ax.patch.set_alpha(0)
        ax.tick_params(top=False, which='both')
        # Set bottom ticks pointing out
        ax.tick_params(axis='x', which='both', direction='out')
    
    apogee_data = import_apogee()
    apogee_subset = apogee_region(apogee_data, galr_lim=GALR_LIM, 
                                  absz_lim=ABSZ_LIM)
    
    print('Plotting DTDs...')
    for i, RIa in enumerate(tqdm(DTD_LIST)):
        axs[0,i].set_title(DTD_LABELS[i], va='top', pad=18)
        # Import VICE multi-zone output data
        output_name = '/'.join(['gaussian', 'lateburst', RIa, 'diskmodel'])
        plot_bimodality(axs[0,i], output_name, apogee_subset, uncertainties=True)
    axs[0,0].set_ylabel('Normalized PDF')
    fig.text(0.51, 0.98, 'Late-burst SFH',
             ha='center', va='top', size=title_size)
    # axs[0,0].set_ylabel('Late-burst SFH')
    
    print('Plotting SFHs...')
    for j, evolution in enumerate(tqdm(SFH_LIST)):
        axs[1,j].set_title(SFH_LABELS[j], pad=-8)
        output_name = '/'.join(['gaussian', evolution, 'exponential_timescale15', 'diskmodel'])
        plot_bimodality(axs[1,j], output_name, apogee_subset, uncertainties=True)
    axs[1,0].set_ylabel('Normalized PDF')
    fig.text(0.42, 0.5, 'Exponential DTD ($\\tau=1.5$ Gyr)',
             ha='center', va='top', size=title_size)
    # axs[1,0].set_ylabel('Exponential DTD\n($\\tau=1.5$ Gyr)')
        
    # Plot APOGEE
    print('Plotting APOGEE...')
    for i, feh_bin in enumerate(FEH_BINS):
        subset_slice = apogee_subset[(apogee_subset['FE_H'] >= feh_bin[0]) &
                                     (apogee_subset['FE_H'] < feh_bin[1])]
        mdf, bin_edges = apogee_mdf(subset_slice, col='O_FE', 
                                    smoothing=SMOOTH_WIDTH, 
                                    bins=np.arange(-0.15, 0.56, 0.01))
        mdf /= mdf.max()
        bin_centers = get_bin_centers(bin_edges)
        axs[1,-1].plot(bin_centers, mdf, ls=LINESTYLES[i], label=feh_bin, c=COLORS[i])
    axs[1,-1].set_title('APOGEE', pad=8, size=title_size)
    
    for ax in axs[-1]:
        ax.set_xlabel('[O/Fe]')
    for ax in axs[:,0]:
        ax.xaxis.set_major_locator(MultipleLocator(0.2))
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    axs[0,0].set_xlim(OFE_LIM)
    axs[0,0].set_ylim((0, 1.1))
    axs[1,0].set_ylim((0, 1.1))
    axs[0,0].legend(loc='upper left', frameon=False, bbox_to_anchor=(0.65, 1.),
                    title='[Fe/H] bin')
    
    # Add gray box underneath APOGEE subpanel for visual separation
    highlight_panels(fig, axs, (1,-1))
    
    # Save
    plt.savefig(paths.figures / 'ofe_bimodality_summary.pdf', dpi=300)
    plt.close()
    print('Done!')
    

def plot_bimodality(ax, output_name, apogee_data, smoothing=SMOOTH_WIDTH, 
                    uncertainties=True):
    mzs = MultizoneStars.from_output(output_name)
    if uncertainties:
        mzs.model_uncertainty(inplace=True, apogee_data=apogee_data)
    subset = mzs.region(galr_lim=GALR_LIM, absz_lim=ABSZ_LIM)
    subset.resample_zheight(20000, apogee_data, inplace=True)
    for i, feh_bin in enumerate(FEH_BINS):
        subset_slice = subset.filter({'[fe/h]': feh_bin})
        mdf, bin_edges = subset_slice.mdf('[o/fe]', smoothing=smoothing,
                                          bins=np.arange(-0.15, 0.56, 0.01))
        mdf /= mdf.max()
        bin_centers = get_bin_centers(bin_edges)
        ax.plot(bin_centers, mdf, ls=LINESTYLES[i], label=feh_bin, c=COLORS[i])


if __name__ == '__main__':
    main()
    