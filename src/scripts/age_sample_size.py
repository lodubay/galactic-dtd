"""
This script prints the size of the Leung et al. (2023) stellar age sample 
to the file src/tex/output/age_sample_size.txt
"""

import pandas as pd
from apogee_tools import import_apogee
import paths

data = import_apogee()
ages = data[data['LATENT_AGE'].notna()]

with open(paths.output / 'age_sample_size.txt', 'w') as f:
    f.write('\\num{%s}' % str(ages.shape[0]))
