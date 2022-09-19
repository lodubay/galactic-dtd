"""
This script plots multiple [O/Fe] distribution functions (MDFs) for 
multizone runs which assume the same star formathion history (SFH).
"""

from ofe_distribution import plot_multiple_comparison

MIGRATION = 'diffusion'
SFH_NAME = 'twoinfall'
PRIMARY_DTD_LIST = ['powerlaw_slope11',
                    'plateau_width300_slope11',
                    'exponential_timescale15',
                    'prompt_peak050_stdev015_timescale30']
PRIMARY_LABELS = ['Power-Law\n($\\alpha=-1.1$)',
                  'Broken Power-Law\n($W=300$ Myr)',
                  'Exponential\n($\\tau=1.5$ Gyr)',
                  'Prompt ($\\mu=50$ Myr)\n+ Exponential ($\\tau=3$ Gyr)']
SECONDARY_DTD_LIST = ['powerlaw_slope14',
                      'plateau_width1000_slope11',
                      'exponential_timescale30']
SECONDARY_LABELS = ['Power-Law\n($\\alpha=-1.4$)',
                    'Broken Power-Law\n($W=1$ Gyr)',
                    'Exponential\n($\\tau=3$ Gyr)']

def main(verbose=False):
    primary_outputs = ['%s/%s/%s' % (MIGRATION, SFH_NAME, dtd) for dtd in \
                       PRIMARY_DTD_LIST]
    plot_multiple_comparison(primary_outputs, PRIMARY_LABELS, 
                             fname='mdf_ofe_twoinfall_1.pdf', 
                             verbose=verbose, double_line_titles=True)
    secondary_outputs = ['%s/%s/%s' % (MIGRATION, SFH_NAME, dtd) for dtd in \
                         SECONDARY_DTD_LIST]
    plot_multiple_comparison(secondary_outputs, SECONDARY_LABELS, 
                             fname='mdf_ofe_twoinfall_2.pdf', 
                             verbose=verbose, double_line_titles=True)

if __name__ == '__main__':
    main(verbose=True)
