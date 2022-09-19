"""
Plot MDFs from all multizone runs
"""

from ofe_distribution import main

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
        main(evolution, RIa)