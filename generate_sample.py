"""
This script generates the APOGEE sample with Leung et al. (2023) ages from scratch.
The resulting CSV file will be located at src/data/APOGEE/sample.csv
"""

# Append scripts directory to path to access multizone directory
import sys
import os
sys.path.append(os.path.abspath('./src/scripts'))
from apogee_tools import import_apogee

import_apogee(overwrite=True, verbose=True)
