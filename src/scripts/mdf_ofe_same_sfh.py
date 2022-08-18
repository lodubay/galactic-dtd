"""
This script plots multiple [O/Fe] distribution functions (MDFs) for 
multizone runs which assume the same star formathion history (SFH).
"""

from ofe_distribution import plot_multiple_comparison

MIGRATION = 'diffusion'
SFH_NAME = 'insideout'
DTD_LIST = ['powerlaw_slope11',
            'powerlaw_slope14',
            'plateau_width300_slope11',
            'plateau_width1000_slope11']
LABELS = ['Power-Law\n($\\alpha=-1.1$)',
          'Power-Law\n($\\alpha=-1.4$)',
          'Broken Power-Law\n($W=300$ Myr)',
          'Broken Power-Law\n($W=1$ Gyr)']

def main(verbose=False):
    outputs = ['%s/%s/%s' % (MIGRATION, SFH_NAME, dtd) for dtd in DTD_LIST]
    plot_multiple_comparison(outputs, LABELS, fname='mdf_ofe_same_sfh.pdf', 
                             verbose=verbose, double_line_titles=True)

if __name__ == '__main__':
    main(verbose=True)
