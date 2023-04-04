"""
This script plots multiple [O/Fe] distribution functions (MDFs) for 
multizone runs which assume the same star formathion history (SFH).
"""

from ofe_distribution import plot_multiple_comparison

MIGRATION = 'diffusion'
SFH_NAME = 'insideout'
PRIMARY_DTD_LIST = [
    'exponential_timescale15',
    'powerlaw_slope11',
    'plateau_width300_slope11',
    # 'exponential_timescale15',
    # 'triple_delay040',
]
PRIMARY_LABELS = [
    'Exponential DTD\n($\\tau=1.5$ Gyr)',
    'Power law DTD\n($\\alpha=-1.1$)',
    'Plateau DTD\n($W=300$ Myr)',
    # 'Exponential\n($\\tau=1.5$ Gyr)',
    # 'Triple-System Evolution\n(Rajamuthukumar+ 2022)',
]

def main(verbose=False):
    outputs = ['%s/%s/%s' % (MIGRATION, SFH_NAME, dtd) for dtd in PRIMARY_DTD_LIST]
    plot_multiple_comparison(outputs, PRIMARY_LABELS, 
                             fname='mdf_ofe_insideout_alt.png', 
                             verbose=verbose, double_line_titles=True)

if __name__ == '__main__':
    main(verbose=True)
