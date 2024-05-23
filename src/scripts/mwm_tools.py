"""
This file contains functions for importing and handling Milky Way Mapper and
related datasets. If run as a script, it will re-generate the main sample file
(src/data/MWM/sample.csv).
"""

from pathlib import Path
import numpy as np
import pandas as pd
import paths
from utils import fits_to_pandas, quad_add

# Data file names
ALLSTAR_FNAME = 'astraAllStarASPCAP-0.5.0.fits.gz'

def main():
    apogee_data = import_mwm(overwrite=True, verbose=True)
    print(apogee_data)


def import_mwm(name='sample.csv', parent_dir=paths.data/'MWM', 
               verbose=False, overwrite=False):
    """
    Import combined and cut sample of MWM data, generating first if needed.
    
    Parameters
    ----------
    name : str, optional
        Name of CSV file containing sample data. The default is 'sample.csv'.
    parent_dir : str or pathlib.Path, optional
        The parent directory containing APOGEE data files. The default is
        '../data/MWM/'.
    verbose : bool, optional
        Whether to print verbose output to terminal. The default is False.
    overwrite : bool, optional
        If True, re-generates the sample file even if it already exists.
        The default is False.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing combined and cut sample of APOGEE data.
    """
    sample_file_path = parent_dir / name
    if sample_file_path.exists() and not overwrite:
        if verbose:
            print('Reading MWM sample from', sample_file_path)
        df = pd.read_csv(sample_file_path)
    else:
        if verbose:
            print('Sample file at', sample_file_path, 'not found.\n' + \
                  'Importing MWM catalog and generating sample...')
        df = gen_mwm_sample(parent_dir=parent_dir, verbose=verbose)
        df.to_csv(sample_file_path, index=False)
        if verbose:
            print('Done.')
    return df


def gen_mwm_sample(parent_dir=paths.data/'MWM', verbose=False):
    """
    Make selection cuts to MWM sample.
    
    Parameters
    ----------
    parent_dir : str or pathlib.Path, optional
        The parent directory containing APOGEE data files. The default is
        '../data/MWM/'.
    verbose : bool, optional
        Whether to print verbose output to terminal. The default is False.
    
    Returns
    -------
    sample : pandas.DataFrame
    """
    # Make parent dir if needed
    if not Path(parent_dir).is_dir():
        parent_dir.mkdir(parents=True)
    mwm_catalog_path = parent_dir / ALLSTAR_FNAME
    if verbose: print('Importing allStar file...')
    mwm_catalog = fits_to_pandas(mwm_catalog_path, hdu=3)
    if verbose: print('Implementing quality cuts...')
    sample = mwm_quality_cuts(mwm_catalog)
    # Calculate [O/Fe] ratio and errors
    sample['O_FE'] = sample['o_h'] - sample['fe_h']
    sample['O_FE_ERR'] = quad_add(sample['e_o_h'], sample['e_fe_h'])
    # add duplicate columns with consistent names from APOGEE
    sample['FE_H'] = sample['fe_h'].copy()
    sample['FE_H_ERR'] = sample['e_fe_h'].copy()
    # Calculate galactocentric coordinates based on galactic l, b and Gaia dist
    galr, galphi, galz = sky_to_galactocentric(
        sample['ra'], sample['dec'], sample['r_med_photogeo']/1000
    )
    sample['GALR'] = galr # kpc
    sample['GALPHI'] = galphi # deg
    sample['GALZ'] = galz # kpc
    # Drop unneeded columns
    return sample.copy()


def mwm_quality_cuts(df):
    """
    Make quality cuts on the APOGEE catalog.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Full APOGEE catalog
    
    Returns
    -------
    pandas.DataFrame
    """
    # Limit to main red star sample
    # df = df[df['EXTRATARG'] == 0]
    # Weed out bad flags
    # fatal_flags = (2**23) # STAR_BAD
    # df = df[df['ASPCAPFLAG'] & fatal_flags == 0]
    # Cut low-S/N targets
    df = df[df['snr'] > 80]
    # Limit to giants
    df = df[(df['logg'] > 1) & (df['logg'] < 3.8) & 
            (df['teff'] > 3500) & (df['teff'] < 5500)]
    # Replace NaN stand-in values with NaN
    # df.replace(99.999, np.nan, inplace=True)
    # Limit to stars with measurements of both [Fe/H] and [O/Fe]
    df.dropna(subset=['fe_h', 'o_fe'], inplace=True)
    df.reset_index(inplace=True, drop=True)
    return df
    

def sky_to_galactocentric(ra, dec, distance):
    r"""
    Use astropy's SkyCoord to convert Galactic (l, b, distance) coordinates
    to galactocentric (r, phi, z) coordinates.

    Parameters
    ----------
    ra : array-like
        Right ascension in degrees
    dec : array-like
        Declination in degrees
    distance : array-like
        Distance (from Sun) in kpc

    Returns
    -------
    galr : numpy array
        Galactocentric radius in kpc
    galphi : numpy array
        Galactocentric phi-coordinates in degrees
    galz : numpy arraay
        Galactocentric z-height in kpc
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord, Galactocentric
    ra = np.array(ra)
    dec = np.array(dec)
    d = np.array(distance)
    if ra.shape == dec.shape == d.shape:
        if not isinstance(ra, u.quantity.Quantity):
            ra *= u.deg
        if not isinstance(dec, u.quantity.Quantity):
            dec *= u.deg
        if not isinstance(d, u.quantity.Quantity):
            d *= u.kpc
        galactic = SkyCoord(ra=ra, dec=dec, distance=d, frame='icrs')
        galactocentric = galactic.transform_to(frame=Galactocentric())
        galactocentric.representation_type = 'cylindrical'
        galr = galactocentric.rho.to(u.kpc).value
        galphi = galactocentric.phi.to(u.deg).value
        galz = galactocentric.z.to(u.kpc).value
        return galr, galphi, galz
    else:
        raise ValueError('Arrays must be of same length.')


if __name__ == '__main__':
    main()
