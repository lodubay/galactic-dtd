"""
This script extracts all multizone outputs from the multizone.tar.gz archive
to src/data/multizone/.
"""

import tarfile
import paths

with tarfile.open(paths.data / 'multizone.tar.gz') as f:
    f.extractall(paths.data / 'multizone')
