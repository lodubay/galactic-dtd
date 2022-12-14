"""
Plot age vs [O/Fe] from all multizone runs
"""

from age_ofe_compare import main

SFH_LIST = [
    'insideout', 
    'lateburst', 
    'conroy22', 
    'twoinfall'
]
DTD_LIST = [
    'powerlaw_slope11', 
    'powerlaw_slope14', 
    'exponential_timescale15',
    'exponential_timescale30', 
    'plateau_width300_slope11', 
    'plateau_width1000_slope11', 
    'prompt_peak050_stdev015_timescale30'
]

for evolution in SFH_LIST:
    for RIa in DTD_LIST:
        output_name = '/'.join(['diffusion', evolution, RIa])
        print('\n%s' % output_name)
        main('/'.join(['diffusion', evolution, RIa]))