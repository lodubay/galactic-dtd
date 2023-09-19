"""
This script prints the size of the APOGEE sample to  the file
src/tex/output/sample_size.txt
"""

import pandas as pd
from apogee_tools import import_apogee
import paths

data = import_apogee()

with open(paths.output / 'sample_size.txt', 'w') as f:
    f.write('\\num{%s}' % str(data.shape[0]))
