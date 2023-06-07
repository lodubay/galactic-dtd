"""
File description
"""

from pathlib import Path
import numpy as np
import pandas as pd
import vice
import paths
from _globals import ZONE_WIDTH

SOLAR_Z_TOTAL = 0.014

def main():
    # test
    from utils import multioutput_to_pandas, filter_multioutput_stars
    output_name = 'diffusion/insideout/powerlaw_slope11'
    mzs = MultizoneStars(output_name)
    print(mzs())
    old_stars = multioutput_to_pandas(output_name)
    print(old_stars)
    galr_lim = (3, 5)
    absz_lim = (0, 0.5)
    print(mzs.region(galr_lim, absz_lim))
    print(filter_multioutput_stars(old_stars, galr_lim, absz_lim))


class MultizoneStars:
    def __init__(self, name, data_dir=paths.data/'migration', 
                 zone_width=ZONE_WIDTH):
        self._name = str(name)
        self._fullpath = Path(data_dir) / (self.name + '.vice')
        self._zone_width = zone_width
        # Import star tracer data
        # stars = pd.DataFrame(vice.stars(str(self.fullpath)).todict())
        stars = self.import_tracers()
        self._end_time = stars['formation_time'].max()
        # Import star particle analogue (z-height) data
        analogdata = pd.read_csv(self.fullpath / 'analogdata.out',
                                 comment='#', sep='\t',
                                 names=['zone_origin', 'time_origin', 
                                        'analog_id', 'zfinal']
        )
        # Limit analogdata to same max time as stars data
        analogdata = analogdata[analogdata['time_origin'] <= self.end_time]
        assert(stars.shape[0] == analogdata.shape[0])
        # Combine
        stars[['analog_id', 'zfinal']] = analogdata[['analog_id', 'zfinal']]
        self._stars = stars
        
    def __call__(self, cols=[]):
        """
        Return the ``stars`` dataframe or a subset of the dataframe.
        
        Parameters
        ----------
        cols : str or list of strings, optional
            If an empty list, return the entire DataFrame. If a string, return
            that column of the DataFrame. If a list of strings, return a subset
            of those columns in the DataFrame. The default is [].
            
        Returns
        -------
        pandas.DataFrame or pandas.Series
            Star particle data or subset of that data.
        """
        if cols == []:
            return self.stars
        else:
            # Error handling
            if isinstance(cols, str):
                if cols not in self.stars.columns:
                    raise ValueError('Parameter "cols" must be an element ' + \
                                     'of stars.columns.')
            elif isinstance(cols, list):
                if all([isinstance(c, str) for c in cols]):
                    if not all([c in self.stars.columns for c in cols]):
                        raise ValueError('Each element of "cols" must be ' + \
                                         'an element of stars.columns.')
                else:
                    raise TypeError('Each element of "cols" must be a string.')
            else:
                raise TypeError('Parameter "cols" must be a string or list ' +\
                                'of strings.')
            return self.stars[cols]
        
    def region(self, galr_lim=(0, 20), absz_lim=(0, 5), min_mass=1.0, 
               origin=False):
        """
        Slice DataFrame of stars within a given Galactic region of radius and
        z-height.

        Parameters
        ----------
        stars : pandas DataFrame
            VICE multizone star data.
        galr_lim : tuple
            Minimum and maximum Galactic radius in kpc. The default is (0, 20).
        absz_lim : tuple
            Minimum and maximum of the absolute value of z-height in kpc. The
            default is (0, 5).
        min_mass : float, optional
            Minimum mass of stellar particle. The default is 1.
        origin : bool, optional
            If True, filter by star's original radius instead of final radius. 
            The default is False.

        Returns
        -------
        pandas DataFrame
            Re-indexed DataFrame of stellar parameters
        """
        galr_min, galr_max = galr_lim
        absz_min, absz_max = absz_lim
        if origin:
            galr_col = 'galr_origin'
        else:
            galr_col = 'galr_final'
        # Select subset
        subset = self.stars[(self.stars[galr_col] >= galr_min) &
                            (self.stars[galr_col] < galr_max) &
                            (self.stars['zfinal'].abs() >= absz_min) &
                            (self.stars['zfinal'].abs() < absz_max) &
                            (self.stars['mass'] >= min_mass)]
        subset.reset_index(inplace=True)
        return subset
        
    @property
    def name(self):
        """
        str
            Multizone output name, excluding ``.vice`` extension or parent dir.
        """
        return self._name
    
    @property
    def fullpath(self):
        """
        pathlib.Path
            Full path to multizone output.
        """
        return self._fullpath
    
    @property
    def stars(self):
        """
        pandas.DataFrame
            Complete star particle data.
        """
        return self._stars
    
    @property
    def zone_width(self):
        """
        float
            Width of each zone in kpc.
        """
        return self._zone_width
    
    @property
    def end_time(self):
        """
        float
            Multizone simulation end time in Gyr.
        """
        return self._end_time
    
    def import_tracers(self, solar_z_total=SOLAR_Z_TOTAL):
        """
        Import star particle data from ``tracers.out``.
        
        Parameters
        ----------
        solar_z_total : float, optional
            Solar metallicity. The default is 0.014.
        
        Returns
        -------
        pandas.DataFrame
            Star particle tracer data.
        """        
        stars = pd.read_csv(self.fullpath / 'tracers.out', 
                            comment='#', sep='\t', 
                            names=['formation_time', 'zone_origin', 
                                   'zone_final', 'mass', 'z(fe)', 'z(o)'],
                            index_col=False
        )
        # Convert radial zone indices to Galactic radii in kpc
        stars['galr_origin'] = self.zone_width * stars['zone_origin']
        stars['galr_final'] = self.zone_width * stars['zone_final']
        # Convert Z to bracket notation
        np.seterr(divide='ignore')
        stars['[fe/h]'] = np.log10(stars['z(fe)'] / vice.solar_z['fe'])
        stars['[o/h]'] = np.log10(stars['z(o)'] / vice.solar_z['o'])
        stars['[o/fe]'] = stars['[o/h]'] - stars['[fe/h]']
        # Calculate total Z and [M/H]
        stars['z'] = solar_z_total * (stars['z(fe)'] + stars['z(o)']) \
                                     / (vice.solar_z['fe'] + vice.solar_z['o'])
        stars['[m/h]'] = np.log10(stars['z'] / solar_z_total)
        stars['age'] = stars['formation_time'].max() - stars['formation_time']
        return stars
        

if __name__ == '__main__':
    main()
