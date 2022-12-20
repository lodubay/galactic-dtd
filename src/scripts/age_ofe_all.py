"""
Plot age vs [O/Fe] from all multizone runs
"""

from age_ofe import main

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
MIGRATION = 'diffusion'

for evolution in SFH_LIST:
    for RIa in DTD_LIST:
        output_name = '/'.join([MIGRATION, evolution, RIa])
        print('\n%s.vice' % output_name)
        main(evolution, RIa, migration=MIGRATION, log=True, verbose=True)
