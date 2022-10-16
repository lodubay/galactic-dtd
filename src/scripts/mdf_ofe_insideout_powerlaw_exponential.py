"""
This script plots multiple [O/Fe] distribution functions (MDFs) for 
multizone runs which assume the same star formathion history (SFH).
"""

from ofe_distribution import plot_multiple_comparison

MIGRATION = 'diffusion'
SFH_NAME = 'insideout'
DTD_LIST = ['powerlaw_slope14',
            'exponential_timescale30',]
LABELS = ['Power-Law DTD\n($\\alpha=-1.4$)',
          'Exponential DTD\n($\\tau=3$ Gyr)',]

def main(verbose=False):
    outputs = ['%s/%s/%s' % (MIGRATION, SFH_NAME, dtd) for dtd in DTD_LIST]
    plot_multiple_comparison(outputs, LABELS, 
                             fname='mdf_ofe_insideout_powerlaw_exponential.png', 
                             verbose=verbose, double_line_titles=True,
                             figure_width=3.25, panel_aspect_ratio=1.6,
                             cbar_width=0.75)

if __name__ == '__main__':
    main(verbose=True)
