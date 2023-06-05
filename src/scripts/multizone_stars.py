"""
File description
"""

from pathlib import Path
import numpy as np
import pandas as pd
import vice
import paths
from _globals import ZONE_WIDTH

class MultizoneStars:
    def __init__(self, name, data_dir=paths.data/'migration', 
                 zone_width=ZONE_WIDTH):
        self._name = str(name)
        self._fullpath = Path(data_dir) / (self.name + '.vice')
        data = pd.read_csv(self.fullpath / 'tracers.out', 
                           comment='#', sep='\t', 
                           names=['formation_time', 'zone_origin', 
                                  'zone_final', 'mass', 'z_fe', 'z_o']
        )
        data['galr_origin'] = zone_width * data['zone_origin']
        data['galr_final'] = zone_width * data['zone_final']
        
        analogdata = pd.read_csv(self.fullpath / 'analogdata.out',
                                 comment='#', sep='\t',
                                 names=['zone_origin', 'time_origin', 
                                        'analog_id', 'zfinal']
        )
        # Limit analogdata to same max time as stars data
        tmax = data['formation_time'].max()
        analogdata = analogdata[analogdata['time_origin'] <= tmax]
        assert(data.shape[0] == analogdata.shape[0])
        
        self._data = data
        
    @property
    def name(self):
        return self._name
    
    @property
    def fullpath(self):
        return self._fullpath
    
    @property
    def data(self):
        return self._data
        
