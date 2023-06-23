"""
This script plots [O/Fe] MDFs for all VICE multizone runs.
"""

from tqdm import tqdm
from ofe_distribution import main

MIGR_LIST = ['diffusion', 'gaussian']
SFH_LIST = [
    'insideout', 
    'lateburst', 
    'earlyburst', 
    'twoinfall'
]
DTD_LIST = [
    'powerlaw_slope11', 
    'powerlaw_slope14', 
    'exponential_timescale15',
    'exponential_timescale30', 
    'plateau_width300_slope11', 
    'plateau_width1000_slope11', 
    'prompt_peak050_stdev015_timescale30',
    'triple_delay040'
]

with tqdm(total=len(MIGR_LIST) * len(SFH_LIST) * len(DTD_LIST)) as t:
    for migration in MIGR_LIST:
        for evolution in SFH_LIST:
            for RIa in DTD_LIST:
                # Import VICE multi-zone output data
                output_name = '/'.join([migration, evolution, RIa])
                main(output_name, uncertainties=True)
                t.update()
