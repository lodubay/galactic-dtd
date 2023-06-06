"""
File description
"""

from pathlib import Path
import numpy as np
import pandas as pd
import vice
import paths
from _globals import ZONE_WIDTH, END_TIME

SOLAR_Z_TOTAL = 0.014

def main():
    # test
    mzs = MultizoneStars('powerlaw_slope11', data_dir=paths.data)
    print(mzs.stars)


class MultizoneStars:
    def __init__(self, name, data_dir=paths.data/'migration', 
                 zone_width=ZONE_WIDTH):
        self._name = str(name)
        self._fullpath = Path(data_dir) / (self.name + '.vice')
        # Import star tracer data
        # stars = pd.DataFrame(vice.stars(str(self.fullpath)).todict())
        stars = pd.read_csv(self.fullpath / 'tracers.out', 
                            comment='#', sep='\t', 
                            names=['formation_time', 'zone_origin', 
                                   'zone_final', 'mass', 'z(fe)', 'z(o)'],
                            index_col=False
        )
        # Convert radial zone indices to Galactic radii in kpc
        stars['galr_origin'] = zone_width * stars['zone_origin']
        stars['galr_final'] = zone_width * stars['zone_final']
        # Convert Z to bracket notation
        np.seterr(divide='ignore')
        stars['[fe/h]'] = np.log10(stars['z(fe)'] / vice.solar_z['fe'])
        stars['[o/h]'] = np.log10(stars['z(o)'] / vice.solar_z['o'])
        stars['[o/fe]'] = stars['[o/h]'] - stars['[fe/h]']
        # Calculate total Z and [M/H]
        stars['z'] = SOLAR_Z_TOTAL * (stars['z(fe)'] + stars['z(o)']) \
                                     / (vice.solar_z['fe'] + vice.solar_z['o'])
        stars['[m/h]'] = np.log10(stars['z'] / SOLAR_Z_TOTAL)
        stars['age'] = END_TIME - stars['formation_time']
        # Import star particle analogue (z-height) data
        analogdata = pd.read_csv(self.fullpath / 'analogdata.out',
                                 comment='#', sep='\t',
                                 names=['zone_origin', 'time_origin', 
                                        'analog_id', 'zfinal']
        )
        # Limit analogdata to same max time as stars data
        tmax = stars['formation_time'].max()
        analogdata = analogdata[analogdata['time_origin'] <= tmax]
        assert(stars.shape[0] == analogdata.shape[0])
        # Combine
        stars[['analog_id', 'zfinal']] = analogdata[['analog_id', 'zfinal']]
        self._stars = stars
        
    @property
    def name(self):
        return self._name
    
    @property
    def fullpath(self):
        return self._fullpath
    
    @property
    def stars(self):
        return self._stars
        

if __name__ == '__main__':
    main()
