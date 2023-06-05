"""
This file contains functions for importing and handling APOGEE and related
datasets. If run as a script, it will re-generate the main sample file
(src/data/APOGEE/sample.csv).
"""

from pathlib import Path
import numpy as np
import pandas as pd
import paths
from utils import galactic_to_galactocentric, fits_to_pandas

def main():
    import_apogee(overwrite=True, verbose=True)


def apogee_region(data, galr_lim=(0, 20), absz_lim=(0, 5)):
    """
    Slice APOGEE data within a given Galactic region of radius and z-height.

    Parameters
    ----------
    data : pandas DataFrame
        APOGEE data
    galr_lim : tuple
        Minimum and maximum Galactic radius in kpc
    absz_lim : tuple
        Minimum and maximum of the absolute value of z-height in kpc

    Returns
    -------
    pandas DataFrame
        Re-indexed DataFrame of stellar parameters
    """
    galr_min, galr_max = galr_lim
    absz_min, absz_max = absz_lim
    # Select subset
    subset = data[(data['GALR'] >= galr_min) &
                  (data['GALR'] < galr_max) &
                  (data['GALZ'].abs() >= absz_min) &
                  (data['GALZ'].abs() < absz_max)].copy()
    subset.reset_index(inplace=True, drop=True)
    return subset


def import_apogee(name='sample.csv', parent_dir=paths.data/'APOGEE', 
                  verbose=False, overwrite=False):
    """
    Import combined and cut sample of APOGEE data, generating first if needed.
    
    Parameters
    ----------
    name : str, optional
        Name of CSV file containing sample data. The default is 'sample.csv'.
    parent_dir : str or pathlib.Path, optional
        The parent directory containing APOGEE data files. The default is
        '../data/APOGEE/'.
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
            print('Reading APOGEE sample from', sample_file_path)
        df = pd.read_csv(sample_file_path)
    else:
        if verbose:
            print('Sample file at', sample_file_path, 'not found.\n' + \
                  'Importing APOGEE catalog and generating sample...')
        df = gen_apogee_sample(parent_dir=parent_dir, verbose=verbose)
        df.to_csv(sample_file_path, index=False)
        if verbose:
            print('Done.')
    return df


def gen_apogee_sample(parent_dir=paths.data/'APOGEE', verbose=False):
    """
    Make selection cuts to APOGEE sample and combine with age catalogs.
    
    Parameters
    ----------
    parent_dir : str or pathlib.Path, optional
        The parent directory containing APOGEE data files. The default is
        '../data/APOGEE/'.
    verbose : bool, optional
        Whether to print verbose output to terminal. The default is False.
    
    Returns
    -------
    sample : pandas.DataFrame
    """
    if verbose: print('Importing allStar file...')
    apogee_catalog = fits_to_pandas(parent_dir / 
                                    'allStarLite-dr17-synspec.fits', hdu=1)
    # Add ages from row-matched datasets BEFORE any cuts
    # Add ages from astroNN (Leung & Bovy 2019)
    if verbose: print('Joining with astroNN age catalog...')
    astroNN_catalog = fits_to_pandas(parent_dir / 'apogee_astroNN-DR17.fits')
    full_catalog = join_astroNN_ages(apogee_catalog, astroNN_catalog)
    # Add ages from Leung et al. (2023)
    if verbose: print('Joining with latent age catalog...')
    leung23_catalog = pd.read_csv(parent_dir / 'nn_latent_age_dr17.csv')
    full_catalog = join_latent_ages(full_catalog, leung23_catalog)
    if verbose: print('Implementing quality cuts...')
    sample = apogee_quality_cuts(full_catalog)
    sample = apogee_galactocentric_coords(sample)
    sample = drop_apogee_columns(sample)
    return sample


def join_astroNN_ages(apogee_df, astroNN_df):
    """
    Join the recommended age from astroNN to the row-matched APOGEE dataset.
    
    Parameters
    ----------
    apogee_df : pandas.DataFrame
        Full APOGEE dataset without cuts
    astroNN_df : pandas.DataFrame
        astroNN dataset
    
    Returns
    -------
    joined_df : pandas.DataFrame
        APOGEE dataset with astroNN ages
    """
    cols = ['age_lowess_correct', 'age_total_error']
    astroNN_ages = astroNN_df[cols].copy()
    astroNN_ages.columns = ['ASTRONN_AGE', 'ASTRONN_AGE_ERR']
    joined = apogee_df.join(astroNN_ages)
    return joined


def join_latent_ages(apogee_df, leung23_df):
    """
    Join ages from Leung et al. (2023) to the row-matched APOGEE dataset.
    
    Parameters
    ----------
    apogee_df : pandas.DataFrame
        Full APOGEE dataset without cuts
    leung23_df : pandas.DataFrame
        Dataset from Leung et al. (2023)
    
    Returns
    -------
    joined_df : pandas.DataFrame
        APOGEE dataset with astroNN ages
    """
    cols = ['LogAge', 'LogAge_Error', 'Age', 'Age_Error']
    latent_ages = leung23_df[cols].copy()
    latent_ages.columns = ['LOG_LATENT_AGE', 'LOG_LATENT_AGE_ERR', 
                           'LATENT_AGE', 'LATENT_AGE_ERR']
    joined = apogee_df.join(latent_ages)
    return joined


def apogee_quality_cuts(df):
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
    df = df[df['EXTRATARG'] == 0]
    # Weed out bad flags
    fatal_flags = (2**23) # STAR_BAD
    df = df[df['ASPCAPFLAG'] & fatal_flags == 0]
    # Cut low-S/N targets
    df = df[df['SNREV'] > 80]
    # Limit to giants
    df = df[(df['LOGG'] > 1) & (df['LOGG'] < 3.8) & 
            (df['TEFF'] > 3500) & (df['TEFF'] < 5500)]
    # Replace NaN stand-in values with NaN
    df.replace(99.999, np.nan, inplace=True)
    df.reset_index(inplace=True, drop=True)
    return df


def apogee_galactocentric_coords(df):
    """
    Add columns to the APOGEE dataset with galactocentric coordinates
    
    Parameters
    ----------
    df : pandas.DataFrame
        APOGEE dataset containing galactic longitude, latitutde, and Gaia
        photogeometric distances
    
    Returns
    -------
    df : pandas.DataFrame
        Same as input with three new columns:
          - 'GALR': galactocentric radius in kpc
          - 'GALPHI': galactocentric azimuth in degrees
          - 'GALZ': height above the Galactic midplane in kpc
    """
    # Calculate galactocentric coordinates based on galactic l, b and Gaia dist
    galr, galphi, galz = galactic_to_galactocentric(
        df['GLON'], df['GLAT'], df['GAIAEDR3_R_MED_PHOTOGEO']/1000
    )
    df['GALR'] = galr # kpc
    df['GALPHI'] = galphi # deg
    df['GALZ'] = galz # kpc
    return df


def drop_apogee_columns(df):
    """
    Drop unneeded columns from the APOGEE dataset
    """
    cols = ['APOGEE_ID', 'RA', 'DEC', 'GALR', 'GALPHI', 'GALZ', 'SNREV',
            'TEFF', 'TEFF_ERR', 'LOGG', 'LOGG_ERR', 'FE_H', 'FE_H_ERR',
            'O_FE', 'O_FE_ERR', 'ASTRONN_AGE', 'ASTRONN_AGE_ERR', 
            'LATENT_AGE', 'LATENT_AGE_ERR', 'LOG_LATENT_AGE', 
            'LOG_LATENT_AGE_ERR']
    return df[cols].copy()


if __name__ == '__main__':
    main()
