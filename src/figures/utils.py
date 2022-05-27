"""
Functions used by many plotting scripts.
"""

from pathlib import Path
from numpy.random import default_rng
import pandas as pd
import vice

def migration_stars(output_name, data_dir='.'):
    """
    Convert VICE milkyway stars output to pandas DataFrame.
    Includes z-height from analogdata.
    Takes a few seconds to run.
    """
    output = vice.output(Path(data_dir) / output_name)
    stars = pd.DataFrame(dict(output.stars))
    analogdata = pd.read_csv('%s_analogdata.out' % output_name, sep='\t')
    # Limit analogdata to same max time as stars data
    tmax = max(output.stars['formation_time'])
    analogdata = analogdata[analogdata['time_origin'] <= tmax]
    # Combine relevant data
    stars[['analog_id', 'zfinal']] = analogdata[['analog_id', 'zfinal']]
    # Remove massless particles
    stars = stars[stars['mass'] > 0]
    return stars

def galactic_region(stars, galr_lim=(0, 20), absz_lim=(0, 5), zone_width=0.1):
    """
    Slice DataFrame of stars within a given Galactic region of radius and z-height.

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
                   (stars['zfinal'].abs() < absz_max)]
    subset.reset_index(inplace=True)
    return subset

def sample_stars(stars, n):
    """
    Randomly sample n stars from VICE output.

    Parameters
    ----------
    stars : pandas DataFrame
        Output of stars_dataframe()
    n : Number of random samples

    Returns
    -------
    pandas DataFrame
        Re-indexed DataFrame of stellar parameters
    """
    # Initialize default numpy random number generator
    rng = default_rng()
    # Randomly sample without replacement
    rand_indices = rng.choice(stars.index, size=n, replace=False)
    sample = stars.loc[rand_indices]
    return sample.reset_index()