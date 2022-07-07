"""
Functions used by many plotting scripts.
"""

from pathlib import Path
import numpy as np
from numpy.random import default_rng
import pandas as pd
from astropy.table import Table
import vice
import paths

# =============================================================================
# DATA IMPORT AND UTILITY FUNCTIONS
# =============================================================================

def import_allStar(verbose=False):
    """
    Import APOGEE allStar data and, if necessary, export for faster re-use.
    """
    clean_df_path = paths.data / 'APOGEE' / 'allStarLite-dr17-synspec-clean.csv'
    try:
        if verbose:
            print('Importing from %s' % clean_df_path)
        df = pd.read_csv(clean_df_path, dtype={'MEMBER': 'object'})
    except FileNotFoundError:
        if verbose:
            print('Clean allStar file not found, generating from source')
        df = clean_allStar()
        df.to_csv(clean_df_path, index=False)
    return df

def clean_allStar(name='allStarLite-dr17-synspec.fits'):
    """
    Import APOGEE AllStar file and convert it to a pandas DataFrame.

    Parameters
    ----------
    name : star, optional
        The name of the AllStar data file.

    Returns
    -------
    df : pandas DataFrame
    """
    df = fits_to_pandas(paths.data / 'APOGEE' / name, hdu=1)
    # Limit to main red star sample
    df = df[df['EXTRATARG'] == 0]
    df.drop(columns=['EXTRATARG'], inplace=True)
    # Weed out bad flags
    fatal_flags = (2**23) # STAR_BAD
    df = df[df['ASPCAPFLAG'] & fatal_flags == 0]
    # Replace NaN stand-in values with NaN
    df.replace(99.999, np.nan, inplace=True)
    # Replace NaN in certain columns with empty string
    df.replace({'ASPCAPFLAGS': np.nan, 'MEMBER': np.nan}, '', inplace=True)
    # Calculate galactocentric coordinates based on galactic l, b and Gaia dist
    galr, galphi, galz = galactic_to_galactocentric(
        df['GLON'], df['GLAT'], df['GAIAEDR3_R_MED_PHOTOGEO']/1000
    )
    df['GALR'] = galr
    df['GALPHI'] = galphi
    df['GALZ'] = galz
    df.reset_index(inplace=True, drop=True)
    return df

def fits_to_pandas(path, **kwargs):
    """
    Import a table in the form of a FITS file and convert it to a pandas
    DataFrame.

    Parameters
    ----------
    path : Path or str
        Path to fits file
    Other keyword arguments are passed to astropy.table.Table

    Returns
    -------
    df : pandas DataFrame
    """
    # Read FITS file into astropy table
    table = Table.read(path, format='fits', **kwargs)
    # Filter out multidimensional columns
    cols = [name for name in table.colnames if len(table[name].shape) <= 1]
    # Convert byte-strings to ordinary strings and convert to pandas
    df = decode(table[cols].to_pandas())
    return df

def decode(df):
    """
    Decode DataFrame with byte strings into ordinary strings.

    Parameters
    ----------
    df : pandas DataFrame
    """
    str_df = df.select_dtypes([object])
    str_df = str_df.stack().str.decode('utf-8').unstack()
    for col in str_df:
        df[col] = str_df[col]
    return df

def galactic_to_galactocentric(l, b, distance):
    r"""
    Use astropy's SkyCoord to convert Galactic (l, b, distance) coordinates
    to galactocentric (r, phi, z) coordinates.

    Parameters
    ----------
    l : array-like
        Galactic longitude in degrees
    b : array-like
        Galactic latitude in degrees
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
    from astropy.coordinates import SkyCoord, Galactic, Galactocentric
    l = np.array(l)
    b = np.array(b)
    d = np.array(distance)
    if l.shape == b.shape == d.shape:
        if not isinstance(l, u.quantity.Quantity):
            l *= u.deg
        if not isinstance(b, u.quantity.Quantity):
            b *= u.deg
        if not isinstance(d, u.quantity.Quantity):
            d *= u.kpc
        galactic = SkyCoord(l=l, b=b, distance=d, frame=Galactic())
        galactocentric = galactic.transform_to(frame=Galactocentric())
        galactocentric.representation_type = 'cylindrical'
        galr = galactocentric.rho.to(u.kpc).value
        galphi = galactocentric.phi.to(u.deg).value
        galz = galactocentric.z.to(u.kpc).value
        return galr, galphi, galz
    else:
        raise ValueError('Arrays must be of same length.')

def multioutput_to_pandas(output_name, data_dir='../data/migration_outputs'):
    """
    Convert VICE multizone stars output to pandas DataFrame (slow).

    Parameters
    ----------
    output_name : str
        Path to the .vice directory containing the migration simulation output
    data_dir : str, optional
        Path to the parent directory of all migration outputs. The default is
        '../data/migration_outputs'.

    Returns
    -------
    pandas DataFrame
        Parameters of simulated stellar populations including galactic z-height
    """
    full_path = Path(data_dir) / output_name
    output = vice.output(str(full_path))
    stars = pd.DataFrame(dict(output.stars))
    analogdata = pd.read_csv('%s_analogdata.out' % full_path, sep='\t')
    # Limit analogdata to same max time as stars data
    tmax = max(output.stars['formation_time'])
    analogdata = analogdata[analogdata['time_origin'] <= tmax]
    # Combine relevant data
    stars[['analog_id', 'zfinal']] = analogdata[['analog_id', 'zfinal']]
    stars.dropna(how='any', inplace=True)
    return stars

def filter_multioutput_stars(stars, galr_lim=(0, 20), absz_lim=(0, 5),
                             zone_width=0.1, min_mass=1.0):
    """
    Slice DataFrame of stars within a given Galactic region of radius and
    z-height.

    Parameters
    ----------
    stars : pandas DataFrame
        Output from stars_dataframe()
    galr_lim : tuple
        Minimum and maximum Galactic radius in kpc
    absz_lim : tuple
        Minimum and maximum of the absolute value of z-height in kpc
    zone_width : float
        Width of each simulation zone in kpc

    Returns
    -------
    pandas DataFrame
        Re-indexed DataFrame of stellar parameters
    """
    galr_min, galr_max = galr_lim
    absz_min, absz_max = absz_lim
    # Convert Galactic radius in kpc to zone #
    zone_min = galr_min / zone_width
    zone_max = galr_max / zone_width
    # Select subset
    subset = stars[(stars['zone_final'] >= zone_min) &
                   (stars['zone_final'] < zone_max) &
                   (stars['zfinal'].abs() >= absz_min) &
                   (stars['zfinal'].abs() < absz_max) &
                   (stars['mass'] >= min_mass)]
    subset.reset_index(inplace=True)
    return subset

def sample_dataframe(df, n):
    """
    Randomly sample n unique rows from a pandas DataFrame.

    Parameters
    ----------
    df : pandas DataFrame
    n : int
        Number of random samples to draw

    Returns
    -------
    pandas DataFrame
        Re-indexed DataFrame of n sampled rows
    """
    if isinstance(df, pd.DataFrame):
        # Initialize default numpy random number generator
        rng = default_rng()
        # Randomly sample without replacement
        rand_indices = rng.choice(df.index, size=n, replace=False)
        sample = df.loc[rand_indices]
        return sample.reset_index()
    else:
        raise TypeError('Expected pandas DataFrame.')
