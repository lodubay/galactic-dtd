"""
Plot age vs [O/Fe] from all multizone runs and generate a summary table of
RMS median-difference scores.
"""

from tqdm import tqdm
import pandas as pd
from age_ofe import plot_age_ofe
from utils import import_apogee, multioutput_to_pandas
import paths

AGE_SOURCES = ['M19', 'L23']
SFH_LIST = [
    'insideout', 
    'lateburst', 
    'conroy22_JW20yields', 
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
MIGRATION = 'diffusion'
DATA_DIR = '../data/migration'

def main():
    # Import APOGEE and astroNN data
    apogee_data = import_apogee(verbose=True)
    
    summary_table = pd.DataFrame([], 
        index=pd.MultiIndex.from_product([DTD_LIST, SFH_LIST], 
                                         names=['DTD', 'SFH']),
        columns=AGE_SOURCES,
    )
    
    with tqdm(total=len(SFH_LIST) * len(DTD_LIST) * len(AGE_SOURCES)) as t:
        for evolution in SFH_LIST:
            for RIa in DTD_LIST:
                # Import VICE multi-zone output data
                output_name = '/'.join([MIGRATION, evolution, RIa])
                vice_stars = multioutput_to_pandas(output_name, DATA_DIR, 
                                                   verbose=True)
                scores = []
                for age_source in AGE_SOURCES:
                    print(evolution, RIa, age_source)
                    fname = '%s_%s_%s.png' % (evolution, RIa, age_source)
                    score = plot_age_ofe(vice_stars, apogee_data, fname=fname, 
                                         ages=age_source, log=True, score=True,
                                         uncertainties=True, verbose=True)
                    scores.append(score)
                    t.update()
                summary_table.loc[RIa, evolution] = scores
    
    print('\n-----------------------------------\n')
    print(summary_table)
    summary_table.to_csv(paths.debug / 'age_ofe/scores.csv')


if __name__ == '__main__':
    main()
