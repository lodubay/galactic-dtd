"""
File description
"""

import math as m
from numbers import Number
from pathlib import Path
import numpy as np
import pandas as pd
import vice
import paths
from _globals import ZONE_WIDTH

SOLAR_Z_TOTAL = 0.014

def main():
    # test
    output_name = 'gaussian/conroy22_JW20yields/exponential_timescale15'
    mzs = MultizoneStars.from_output(output_name)
    print(mzs.galr_lim)
    print(mzs.absz_lim)
    region = mzs.region((7, 9), (0, 0.5))
    print(region)
    print(region.galr_lim)
    print(region.absz_lim)
    print(mzs('age'))
    noisy = mzs.model_uncertainty()
    print(noisy('age'))


class MultizoneStars:
    """
    Star particle data from VICE multizone outputs.
    
    Note
    ----
    In all but rare circumstances, new instances should be generated using
    the MultizoneStars.from_output() class method.
    """
    def __init__(self, stars, name='', fullpath='', zone_width=ZONE_WIDTH,
                 galr_lim=None, absz_lim=None, noisy=False):
        self.stars = stars
        self.name = name
        self.fullpath = fullpath
        self.zone_width = zone_width
        # Automatically calculate bounds of galactic region if none are given
        if galr_lim is None:
            galr_lim = (m.floor(10 * self.stars['galr_final'].min()) / 10.,
                        m.ceil(10 * self.stars['galr_final'].max()) / 10.)
        self.galr_lim = galr_lim
        if absz_lim is None:
            absz_lim = (m.floor(10 * self.stars['zfinal'].abs().min()) / 10.,
                        m.ceil(10 * self.stars['zfinal'].abs().max()) / 10.)
        self.absz_lim = absz_lim
        self.noisy = noisy
        
    @classmethod
    def from_output(cls, name, data_dir=paths.data/'migration', 
                    zone_width=ZONE_WIDTH, verbose=False):
        """
        Generate an instance of MultizoneStars from a VICE multizone output.
        
        Parameters
        ----------
        name : str
            Name of VICE output, excluding ``.vice`` extension.
        data_dir : str or pathlib.Path, optional
            Parent directory of VICE output. The default is 
            ``../data/migration``.
        zone_width : float, optional
            Width of simulation zones in kpc. The default is 0.1.
        verbose : bool, optional
            If True, print some updates to the terminal. The default is False.
        
        Returns
        -------
        MultizoneStars instance
        """
        fullpath = Path(data_dir) / (name + '.vice')
        # Import star tracer data
        if verbose: 
            print('Importing VICE multizone data from', str(fullpath))
        stars = cls.import_tracers(fullpath, zone_width=zone_width)
        # stars = pd.DataFrame(vice.stars(str(self.fullpath)).todict())
        # Import star particle analogue (z-height) data
        analogdata = pd.read_csv(fullpath / 'analogdata.out',
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
        return cls(stars, name=name, fullpath=fullpath, zone_width=zone_width)
        
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
               origin=False, inplace=False):
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
        inplace : bool, optional
            If True, update the current class instance. If False, returns a 
            new class instance with the limited subset. The default is False.

        Returns
        -------
        MultizoneStars instance or None
        """
        galr_min, galr_max = galr_lim
        absz_min, absz_max = absz_lim
        if origin:
            galr_col = 'galr_origin'
        else:
            galr_col = 'galr_final'
        # Select subset
        subset = self.stars[(self.stars[galr_col]       >= galr_min) &
                            (self.stars[galr_col]       <  galr_max) &
                            (self.stars['zfinal'].abs() >= absz_min) &
                            (self.stars['zfinal'].abs() <  absz_max) &
                            (self.stars['mass']         >= min_mass)]
        subset.reset_index(inplace=True)
        if inplace:
            self.stars = subset
        else:
            return MultizoneStars(subset, name=self.name, 
                                  fullpath=self.fullpath, 
                                  zone_width=self.zone_width, 
                                  galr_lim=galr_lim, absz_lim=absz_lim, 
                                  noisy=self.noisy)
    
    def model_uncertainty(self, apogee_data=None, age_source='L23',
                          inplace=False):
        """
        Forward-model observational uncertainties from median data errors.
        Star particle data are modified in-place, so only run this once!
        
        Parameters
        ----------
        apogee_data : pandas.DataFrame or NoneType, optional
            Full APOGEE data. If None, will be imported from ``sample.csv``.
            The default is None.
        age_source : str, optional
            Source for age uncertainty. Options are 'F19' (Feuillet+ 2019),
            'M19' (Mackereth+ 2019, e.g. astroNN), or 'L23' (Leung+ 2023, e.g.
            latent ages). The default is 'L23'.
        inplace : bool, optional
            If True, update the current class instance. If False, returns a 
            new class instance with the noisy outputs. The default is False.
            
        Returns
        -------
        MultizoneStars instance or None
        """
        noisy_stars = self.stars.copy()
        if apogee_data is None:
            apogee_data = pd.read_csv(paths.data/'APOGEE/sample.csv')
        rng = np.random.default_rng()
        # age uncertainty depends on the age source
        if age_source == 'F19':
            log_age_err = 0.15 # Feuillet et al. (2016), ApJ, 817, 40
            log_age_noise = rng.normal(scale=log_age_err, 
                                       size=noisy_stars.shape[0])
            noisy_stars['age'] *= 10 ** log_age_noise
        elif age_source == 'M19':
            frac_age_err = (apogee_data['ASTRONN_AGE_ERR'] / 
                            apogee_data['ASTRONN_AGE']).median()
            frac_age_noise = rng.normal(scale=frac_age_err, 
                                        size=noisy_stars.shape[0])
            noisy_stars['age'] *= (1 + frac_age_noise)
        elif age_source == 'L23': # L23
            log_age_err = apogee_data['LOG_LATENT_AGE_ERR'].median()
            log_age_noise = rng.normal(scale=log_age_err, 
                                       size=noisy_stars.shape[0])
            noisy_stars['age'] *= 10 ** log_age_noise
        else:
            raise ValueError('Invalid value for "age_source".')
        # [Fe/H] uncertainty
        feh_med_err = apogee_data['FE_H_ERR'].median()
        feh_noise = rng.normal(scale=feh_med_err, size=noisy_stars.shape[0])
        noisy_stars['[fe/h]'] += feh_noise
        # [O/Fe] uncertainty 
        ofe_med_err = apogee_data['O_FE_ERR'].median()
        ofe_noise = rng.normal(scale=ofe_med_err, size=noisy_stars.shape[0])
        noisy_stars['[o/fe]'] += ofe_noise
        if inplace:
            self.stars = noisy_stars
        else:
            return MultizoneStars(noisy_stars, name=self.name, 
                                  fullpath=self.fullpath, 
                                  zone_width=self.zone_width, 
                                  galr_lim=self.galr_lim, 
                                  absz_lim=self.absz_lim, 
                                  noisy=True)
    
    @staticmethod
    def import_tracers(fullpath, zone_width=ZONE_WIDTH, 
                       solar_z_total=SOLAR_Z_TOTAL):
        """
        Import star particle data from ``tracers.out``.
        
        Parameters
        ----------
        fullpath : str or pathlib.Path
            Full path to VICE output directory.
        zone_width : float, optional
            Width of simulation zones in kpc. The default is 0.1.
        solar_z_total : float, optional
            Solar metallicity. The default is 0.014.
        
        Returns
        -------
        pandas.DataFrame
            Star particle tracer data.
        """        
        stars = pd.read_csv(Path(fullpath) / 'tracers.out', 
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
        stars['z'] = solar_z_total * (stars['z(fe)'] + stars['z(o)']) \
                                     / (vice.solar_z['fe'] + vice.solar_z['o'])
        stars['[m/h]'] = np.log10(stars['z'] / solar_z_total)
        stars['age'] = stars['formation_time'].max() - stars['formation_time']
        return stars
    
    def __str__(self):
        return self.stars.__str__()
        
    @property
    def name(self):
        """
        str
            Multizone output name, excluding ``.vice`` extension or parent dir.
        """
        return self._name
    
    @name.setter
    def name(self, value):
        if isinstance(value, str):
            self._name = value
        else:
            raise TypeError('Attribute "name" must be a string. Got:', 
                            type(value))
    
    @property
    def fullpath(self):
        """
        pathlib.Path
            Full path to multizone output.
        """
        return self._fullpath
    
    @fullpath.setter
    def fullpath(self, value):
        if isinstance(value, (str, Path)):
            if '.vice' in str(value) and Path(value).is_dir():
                self._fullpath = Path(value)
            else:
                raise ValueError('Value is not a valid VICE output directory.')
        else:
            raise TypeError('Attribute "fullpath" must be a string or Path. Got:',
                            type(value))
    
    @property
    def stars(self):
        """
        pandas.DataFrame
            Complete star particle data.
        """
        return self._stars
    
    @stars.setter
    def stars(self, value):
        if isinstance(value, pd.DataFrame):
            self._stars = value
        else:
            raise TypeError('Attribute "stars" must be a DataFrame. Got:',
                            type(value))
    
    @property
    def zone_width(self):
        """
        float
            Width of each zone in kpc.
        """
        return self._zone_width
    
    @zone_width.setter
    def zone_width(self, value):
        if isinstance(value, float):
            self._zone_width = value
        else:
            raise TypeError('Attribute "zone_width" must be a float. Got:',
                            type(value))
            
    @property
    def galr_lim(self):
        """
        tuple
            Minimum and maximum bounds on the Galactic radius in kpc.
        """
        return self._galr_lim
    
    @galr_lim.setter
    def galr_lim(self, value):
        if isinstance(value, tuple):
            if len(value) == 2:
                if all([isinstance(x, Number) for x in value]):
                    self._galr_lim = value
                else:
                    raise TypeError('Each item in "galr_lim" must be a number.')
            else:
                raise ValueError('Attribute "galr_lim" must have length 2.')
        else:
            raise TypeError('Attribute "galr_lim" must be a tuple. Got:',
                            type(value))
            
    @property
    def absz_lim(self):
        """
        tuple
            Minimum and maximum bounds on the absolute z-height in kpc.
        """
        return self._absz_lim
    
    @absz_lim.setter
    def absz_lim(self, value):
        if isinstance(value, tuple):
            if len(value) == 2:
                if all([isinstance(x, Number) for x in value]):
                    self._absz_lim = value
                else:
                    raise TypeError('Each item in "absz_lim" must be a number.')
            else:
                raise ValueError('Attribute "absz_lim" must have length 2.')
        else:
            raise TypeError('Attribute "absz_lim" must be a tuple. Got:',
                            type(value))
    
    @property
    def noisy(self):
        """
        bool
            If True, indicates that forward-modeled observational uncertainty
            has been applied to certain columns of the VICE output.
        """
        return self._noisy
    
    @noisy.setter
    def noisy(self, value):
        if isinstance(value, bool):
            self._noisy = value
        else:
            raise TypeError('Attribute "noisy" must be a Boolean. Got:',
                            type(value))
        

if __name__ == '__main__':
    main()
