"""
Functions used by many plotting scripts.
"""

from pathlib import Path
import numpy as np
from numpy.random import default_rng
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, BoundaryNorm, LogNorm
from matplotlib.cm import ScalarMappable
from astropy.table import Table
import vice
import paths
from _globals import ZONE_WIDTH

# =============================================================================
# APOGEE / ASTRONN DATA IMPORT AND UTILITY FUNCTIONS
# =============================================================================

allStar_file_name = 'allStarLite-dr17-synspec.fits'

def apogee_region(data, galr_lim=(0, 20), absz_lim=(0, 5)):
    """
    Slice APOGEE data within a given Galactic region of radius and z-height.

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
    # Select subset
    subset = data[(data['GALR'] >= galr_min) &
                  (data['GALR'] < galr_max) &
                  (data['GALZ'].abs() >= absz_min) &
                  (data['GALZ'].abs() < absz_max)]
    subset.reset_index(inplace=True)
    return subset

def select_giants(data, logg_col='LOGG', teff_col='TEFF'):
    """
    Select stars in the RGB and red clump by cuts in log(g) and Teff.
    
    Parameters
    ----------
    data : pandas.DataFrame
        APOGEE allStar data
    logg_col : str, optional
        Name of column containing surface gravities. The default is 'LOGG'.
    teff_col : str, optional
        Name of column containing effective temperatures. The default is 'TEFF'.
    
    Returns
    -------
    giants : pandas.DataFrame
        APOGEE allStar data with cuts
    """
    giants = data[(data[logg_col] > 1) & (data[logg_col] < 3.8) & 
                  (data[teff_col] > 3500) & (data[teff_col] < 5500) & 
    # exclude the lower-left corner which may contain main sequence stars
                  ~((data[logg_col] > 3) & (data[teff_col] < 4000))]
    giants.reset_index(inplace=True)
    return giants

def import_astroNN(verbose=False, allStar_name=allStar_file_name, 
                   only_giants=True):
    """
    Import the joint data table of APOGEE allStar and astroNN.
    
    If no joint table exists, import each data set separately and join them, 
    then export the result as a CSV file.
    
    Parameters
    ----------
    verbose : bool, optional
        If True, print verbose progress messages. The default is False.
    allStar_name : str, optional
        Name of the APOGEE allStar file. The default is 
        'allStarLite-dr17-synspec.fits'
    only_giants : bool
        If True, makes cuts in log(g) and Teff to return only the giant sample.
        If False, returns all main sample stars. The default is True.
        
    Returns
    -------
    pandas.DataFrame
    """
    joined_path = paths.data / 'APOGEE' / 'allStarLite-astroNN-dr17-clean.csv'
    try:
        if verbose:
            print('Importing APOGEE+astroNN data from %s' % joined_path)
        joined = pd.read_csv(joined_path, dtype={'MEMBER': 'object'})
    except FileNotFoundError:
        if verbose:
            print('No complete CSV found, importing and joining from source')
        # Important! Data can only be cleaned *after* joining the two sets
        allStar = fits_to_pandas(paths.data / 'APOGEE' / allStar_name, hdu=1)
        joined = join_astroNN(allStar)
        joined = clean_allStar(joined, only_giants=only_giants)
        joined.to_csv(joined_path, index=False)
    return joined

def join_astroNN(allStar, name='apogee_astroNN-DR17.fits'):
    """
    Import the astroNN dataset and join it with the APOGEE allStar dataset.
    
    This relies on allStar and astroNN having exactly the same number of rows
    and in the same order, which should be the case if they haven't been
    changed.
    
    Parameters
    ----------
    allStar : pandas.DataFrame
        The APOGEE allStar dataset
    name : str, optional
        Name of the astroNN fits file
    
    Returns
    -------
    pandas.DataFrame
    """
    astroNN = fits_to_pandas(paths.data / 'APOGEE' / name)
    astroNN = rename_astroNN_columns(astroNN)
    # Add a few more abundance combinations
    subtract_abundances(astroNN, 'O_H', 'FE_H', prefix='ASTRONN_')
    subtract_abundances(astroNN, 'MG_H', 'FE_H', prefix='ASTRONN_')
    subtract_abundances(astroNN, 'TI_H', 'FE_H', prefix='ASTRONN_')
    subtract_abundances(astroNN, 'C_H', 'N_H', prefix='ASTRONN_')
    joined = allStar.join(astroNN.drop('APOGEE_ID', axis=1))
    return joined
    
def subtract_abundances(data, abund1, abund2, prefix=''):
    """
    Subtract one abundance measurement from another to calculate the ratio
    of the numerators.
    
    Parameters
    ----------
    data : pandas.DataFrame
        DataFrame containing abundance data
    abund1 : str
        Name of abundance from which abund2 is subtracted, in the format 'X_Z'
        (e.g., 'O_FE').
    abund2 : str
        Name of abundance to subtract from abund1, in the format 'Y_Z' 
        (e.g. 'FE_H').
    prefix : str, optional
        Prefix to abundance column names. The default is ''.
        
    Returns
    -------
    pandas.DataFrame
        Abundance data with new column of combined abundances
    """
    abund1_name = prefix + abund1
    abund2_name = prefix + abund2
    result_name = prefix + abund1.split('_')[0] + '_' + abund2.split('_')[0]
    data[result_name] = data[abund1_name] - data[abund2_name]
    abund1_err = abund1_name + '_ERR'
    abund2_err = abund2_name + '_ERR'
    result_err = result_name + '_ERR'
    data[result_err] = quad_add(data[abund1_err], data[abund2_err])
    return data

def quad_add(arr1, arr2):
    """
    Add two input arrays in quadrature.
    """
    return np.sqrt(arr1**2 + arr2**2)

def rename_astroNN_columns(astroNN):
    """
    Select relevant astroNN columns and rename them to avoid conflicts
    
    Parameters
    ----------
    astroNN : pandas.DataFrame
        DataFrame generated from astroNN fits file
    
    Returns
    -------
    pandas.DataFrame
    """
    cols = ['APOGEE_ID'] + astroNN.columns[5:49].tolist()
    cols += ['weighted_dist', 'weighted_dist_error', 'age_lowess_correct', 
             'age_total_error']
    cols += ['galr', 'galphi', 'galz']
    astroNN = astroNN[cols].copy()
    cols_new = ['APOGEE_ID'] + ['ASTRONN_'+name for name in cols[1:45]]
    cols_new += ['ASTRONN_DIST', 'ASTRONN_DIST_ERR', 'ASTRONN_AGE', 
                 'ASTRONN_AGE_ERR']
    cols_new += ['ASTRONN_'+name.upper() for name in cols[49:]]
    astroNN.columns = cols_new
    return astroNN

def import_allStar(verbose=False, name=allStar_file_name, only_giants=True):
    """
    Import APOGEE allStar data and, if necessary, export for faster re-use.
    
    Parameters
    ----------
    only_giants : bool
        If True, makes cuts in log(g) and Teff to return only the giant sample.
        If False, returns all main sample stars. The default is True.
    """
    clean_df_path = paths.data / 'APOGEE' / 'allStarLite-dr17-synspec-clean.csv'
    try:
        if verbose:
            print('Importing from %s' % clean_df_path)
        df = pd.read_csv(clean_df_path, dtype={'MEMBER': 'object'})
    except FileNotFoundError:
        if verbose:
            print('Clean allStar file not found, generating from source')
        raw = fits_to_pandas(paths.data / 'APOGEE' / name, hdu=1)
        df = clean_allStar(raw.copy(), only_giants=only_giants)
        df.to_csv(clean_df_path, index=False)
    return df

def clean_allStar(df, snr=80, only_giants=True):
    """
    Import APOGEE AllStar file and convert it to a pandas DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        The unmodified, freshly converted DataFrame of the allStar file
    snr : float
        The minimum signal-to-noise threshold
    only_giants : bool
        If True, makes cuts in log(g) and Teff to return only the giant sample.
        If False, returns all main sample stars. The default is True.

    Returns
    -------
    df : pandas DataFrame
    """
    # df = fits_to_pandas(paths.data / 'APOGEE' / name, hdu=1)
    # Limit to main red star sample
    df = df[df['EXTRATARG'] == 0]
    df.drop(columns=['EXTRATARG'], inplace=True)
    # Weed out bad flags
    fatal_flags = (2**23) # STAR_BAD
    df = df[df['ASPCAPFLAG'] & fatal_flags == 0]
    # Cut low-S/N targets
    df = df[df['SNREV'] > snr]
    # Limit to giants
    if only_giants: df = select_giants(df)
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
        
# =============================================================================
# VICE MULTIZONE INPUT AND UTILITY FUNCTIONS
# =============================================================================

def multioutput_to_pandas(output_name, data_dir=paths.data/'migration'):
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
    # Convert zone to radius
    stars['galr_origin'] = stars['zone_origin'] * ZONE_WIDTH
    stars.dropna(how='any', inplace=True)
    return stars

def filter_multioutput_stars(stars, galr_lim=(0, 20), absz_lim=(0, 5),
                             zone_width=ZONE_WIDTH, min_mass=1.0, 
                             zone_origin=False):
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
    min_mass : float, optional
        Minimum mass of stellar particle
    zone_origin : bool, optional
        If True, filter by star's original zone instead of final zone

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

def sample_dataframe(df, n, weights=None):
    """
    Randomly sample n unique rows from a pandas DataFrame.

    Parameters
    ----------
    df : pandas DataFrame
    n : int
        Number of random samples to draw
    weights : array
        Probability weights of the given DataFrame

    Returns
    -------
    pandas DataFrame
        Re-indexed DataFrame of n sampled rows
    """
    if isinstance(df, pd.DataFrame):
        # Initialize default numpy random number generator
        rng = default_rng()
        # Randomly sample without replacement
        rand_indices = rng.choice(df.index, size=n, replace=False, p=weights)
        sample = df.loc[rand_indices]
        return sample.reset_index()
    else:
        raise TypeError('Expected pandas DataFrame.')

# =============================================================================
# FILE NAMING AND STRING FORMATTING
# =============================================================================

_MIGRATION_MODELS = ['diffusion', 'post-process']
_EVOLUTION_MODELS = ['insideout', 'lateburst']
_EFFICIENCY_MODELS = ['johnson21', 'conroy22']
_RIA_MODELS = ['powerlaw', 'exponential', 'plateau', 'prompt']
_RIA_SETTING_DEFAULTS = {
    'powerlaw': -1.1,
    'exponential': 1.5,
    'plateau': 0.2,
    'prompt': 0.05}

def multizone_output_path(migration='diffusion', evolution='insideout',
                          efficiency='johnson21', RIa='powerlaw',
                          RIa_setting=None, minimum_delay=0.04,
                          parent=paths.data/'migration'):
    r"""
    Generate the name of a VICE multizone output based on its simulation
    parameters.

    Parameters
    ----------
    migration : str [default: 'diffusion']
        A keyword denoting the stellar migration model.
        Allowed values:

        - 'diffusion'
        - 'post-process'

    evolution : str [default: 'insideout']
        A keyword denoting the star formation history.
        Allowed values:

        - 'insideout'
        - 'lateburst'

    efficiency : str [default: 'johnson21']
        A keyword denoting the star formation efficiency prescription.
        Allowed values:

        - 'johnson21'
        - 'conroy22'

    RIa : str [default: 'powerlaw']
        A keyword denoting the general form of the delay-time distribution.
        Allowed values:

        - 'powerlaw'
        - 'exponential'
        - 'plateau'
        - 'prompt'

    RIa_setting : float [default: None]
        The primary setting for the given delay-time distribution. If None,
        the default value for the given DTD is assumed. The primary setting
        and default for each DTD is:

        - 'powerlaw': the power-law slope [default: -1.1]
        - 'exponential': the exponential timescale in Gyr [default: 1.5]
        - 'plateau': the plateau width in Gyr [default: 0.2]
        - 'prompt': the peak of the prompt component in Gyr [default: 0.05]

    minimum_delay : float [default: 0.04]
        The minimum delay time in Gyr before the first SNe Ia explode.

    Returns
    -------
    pathlib.Path
        Path to VICE multizone output directory
    """
    if not isinstance(migration, str):
        raise TypeError('Parameter "migration" must be a string. Got: %s'
                        % type(migration))
    if migration not in _MIGRATION_MODELS:
        raise ValueError('Parameter "migration" must be one of %s. Got: %s'
                         % (_MIGRATION_MODELS, migration))
    if not isinstance(evolution, str):
        raise TypeError('Parameter "evolution" must be a string. Got: %s'
                        % type(evolution))
    if evolution not in _EVOLUTION_MODELS:
        raise ValueError('Parameter "evolution" must be one of %s. Got: %s'
                         % (_EVOLUTION_MODELS, evolution))
    if not isinstance(RIa, str):
        raise TypeError('Parameter "RIa" must be a string. Got: %s' % type(RIa))
    if RIa not in _RIA_MODELS:
        raise ValueError('Parameter "RIa" must be one of %s. Got: %s'
                         % (_RIA_MODELS, RIa))
    if RIa_setting is None:
        RIa_setting = _RIA_SETTING_DEFAULTS[RIa]
    elif not isinstance(RIa_setting, float):
        raise TypeError('Parameter "RIa_setting" must be a float. Got: %s'
                        % type(RIa_setting))
    if not isinstance(minimum_delay, float):
        raise TypeError('Parameter "minimum_delay" must be a float. Got: %s'
                        % type(minimum_delay))
    if RIa in ['powerlaw', 'exponential']:
        RIa_name = f'{RIa}{int(abs(RIa_setting)*10):02d}'
    else:
        RIa_name = f'{RIa}{int(RIa_setting*1000):03d}'
    RIa_name += f'_delay{int(minimum_delay*1000):03d}.vice'
    evol_name = '_'.join((evolution, efficiency))
    path = parent / migration / evol_name / RIa_name
    return path

def bracket_string_to_allstar_column(bracket):
    """
    Convert a bracket notation string (e.g. '[Fe/H]') to its corresponding
    allStar column name (e.g. FE_H).

    Parameters
    ----------
    bracket : str
        Bracket notation for chemical abundance (e.g. '[Fe/H]')

    Returns
    -------
    col : str
        allStar column name (e.g. 'FE_H')
    """
    col = bracket.upper()
    col = col.replace('[', '').replace(']', '')
    col = col.replace('/', '_')
    return col

def format_bracket_string(bracket):
    """
    Capitalize the proper characters in a bracket notation string.

    Parameters
    ----------
    bracket : str
        Bracket notation for chemical abundance (e.g. '[Fe/H]')

    Returns
    -------
    str
    """
    bracket = bracket.lower()
    for i in range(len(bracket)):
        if bracket[i-1] in ['[', '/']:
            bracket[i] = bracket[i].upper()
    return bracket

# =============================================================================
# ONE-ZONE MODELS
# =============================================================================

def run_singlezone(name, simtime, overwrite=False, **kwargs):
    """
    Run a VICE one-zone simulation, or, if it has already been run, import
    it instead.

    Parameters
    ----------
    name : str
        Path to the one-zone output directory, not necessarily including the
        '.vice' extension
    simtime : array-like
        Array of simulation times in Gyr
    overwrite : bool, optional
        If True, re-run the model regardless of whether an output already exists
    kwargs : dict, optional
        Keyword arguments passed to vice.singlezone

    Returns
    -------
    vice.singlezone object
    """
    if overwrite:
        sz = vice.singlezone(name=name, **kwargs)
        sz.run(simtime, overwrite=overwrite)
    else:
        try:
            sz = vice.singlezone.from_output(name)
        except OSError:
            sz = vice.singlezone(name=name, **kwargs)
            sz.run(simtime, overwrite=overwrite)
    return sz

# =============================================================================
# VICE FUNCTION WRAPPERS
# =============================================================================

class NormalIMF:
    """
    A normalized initial mass function (IMF).
    """
    def __init__(self, which='kroupa', m_lower=0.08, m_upper=100, dm=0.01):
        """
        Initialize the IMF.

        Parameters
        ----------
        which : string, optional
            Which version of the IMF to use. Must be one of 'salpeter' or
            'kroupa'. The default is 'kroupa'
        m_lower : float, optional
            Lower limit of integration in solar masses. The default is 0.08.
        m_upper : TYPE, optional
            Upper limit of integration in solar masses. The default is 100.
        dm : TYPE, optional
            Integration step in solar masses. The default is 0.01.
        """
        select = {
            'salpeter': vice.imf.salpeter,
            'kroupa': vice.imf.kroupa
        }
        if which in select.keys():
            self.dm = dm
            self.masses = np.arange(m_lower, m_upper + dm, dm)
            self._imf = select[which]
            self.norm = 1
            self.norm = 1 / self.integrate()
        else:
            raise ValueError('IMF must be either "kroupa" or "salpeter".')

    def __call__(self, mass):
        """
        Calculate the normalized IMF at a given stellar mass.

        Parameters
        ----------
        mass : float
            Stellar mass in solar masses.

        Returns
        -------
        float
            The normalized value of the IMF at that stellar mass.
        """
        return self.norm * self._imf(mass)

    def integrate(self):
        """
        float
            The integral of the IMF
        """
        integral = sum([self.__call__(m) * self.dm for m in self.masses])
        return integral

    def weighted_mean(self):
        """
        Calculate the average stellar mass of the IMF.
        """
        weights = np.array([self.__call__(m) for m in self.masses])
        weighted_mean = np.average(self.masses, weights=weights)
        return weighted_mean
    
# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def scatter_hist(ax, x, y, xlim=None, ylim=None, log_norm=True, cmap='gray',
                 cmin=10, vmin=None, vmax=None, nbins=50, color='k',
                 rasterized=True):
    """
    Generate a scatter plot and overlayed 2D histogram for dense data.

    Parameters
    ----------
    ax : matplotlib.axis.Axes
        Axes object on which to plot the data.
    x : array-like
        Horizontal coordinates of the data points.
    y : array-like
        Vertical coordinates of the data points.
    xlim : float, optional
        Bounds for x-axis. The default is None.
    ylim : float, optional
        Bounds for y-axis. The default is None.
    log_norm : bool, optional
        Shade the 2D histogram on a logarithmic scale. The default is True.
    cmap : str, optional
        Colormap for 2D histogram. The default is'gray'.
    cmin : int, optional
        Minimum counts per bin; any number below this will show individual points.
        The default is 10.
    vmin : float or None, optional
        Value to map to minimum of histogram normalization. The default is None.
    vmax : float or None, optional
        Value to map to maximum of histogram normalization. The default is None.
    nbins : int or tuple of ints, optional
        Number of histogram bins. If a tuple, presumed to be (xbins, ybins).
        The default is 50.
    color : str, optional
        Color of individual points. The default is 'k'.
    rasterized : bool, optional [default: True]
        Whether to rasterize the scattered points
    """
    # Set automatic plot bounds
    if not xlim:
        xlim = (np.min(x), np.max(x))
    if not ylim:
        ylim = (np.min(y), np.max(y))
    # Set bin edges
    if type(nbins) == 'tuple':
        xbins, ybins = nbins
    else:
        xbins = ybins = nbins
    xbins = np.linspace(xlim[0], xlim[1], num=xbins, endpoint=True)
    ybins = np.linspace(ylim[0], ylim[1], num=ybins, endpoint=True)
    # Histogram normalization
    if log_norm:
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)
    # Plot
    ax.scatter(x, y, c=color, s=0.5, rasterized=rasterized, edgecolor='none')
    return ax.hist2d(x, y, bins=[xbins, ybins], cmap=cmap, norm=norm, cmin=cmin)


def get_bin_centers(bin_edges):
    """
    Calculate the centers of bins defined by the given bin edges.
    
    Parameters
    ----------
    bin_edges : array-like of length N
        Edges of bins, including the left-most and right-most bounds.
     
    Returns
    -------
    bin_centers : numpy.ndarray of length N-1
        Centers of bins
    """
    bin_edges = np.array(bin_edges, dtype=float)
    if len(bin_edges) > 1:
        return 0.5 * (bin_edges[:-1] + bin_edges[1:])
    else:
        raise ValueError('The length of bin_edges must be at least 2.')


def get_color_list(cmap, bins):
    """
    Split a discrete colormap into a list of colors based on bin edges.
    
    Parameters
    ----------
    cmap : matplotlib colormap
    bins : array-like
        Bin edges, including left- and right-most edges
    
    Returns
    -------
    list
        List of colors of length len(bins) - 1
    """
    rmin, rmax = bins[0], bins[-2]
    colors = cmap([(r-rmin)/(rmax-rmin) for r in bins[:-1]])
    return colors


def discrete_colormap(cmap_name, bounds):
    """
    Convert a continuous colormap into a discrete one.
    
    Parameters
    ----------
    cmap_name : str
        Name of matplotlib colormap
    bounds : array-like
        Bounds of discrete colormap
    
    Returns
    -------
    cmap : matplotlib colormap
    norm : colormap normalization
    """
    cmap = plt.get_cmap(cmap_name)
    norm = BoundaryNorm(bounds, cmap.N)
    return cmap, norm


def setup_discrete_colorbar(fig, cmap, norm, label='', width=0.6):
    """
    Adds a colorbar for a discrete colormap to the bottom of the figure.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
    cmap : matplotlib colormap
    norm : matplotlib colormap normalization
    label : str, optional
        Colorbar label. The default is ''.
    width : float, optional
        Width of colorbar as fraction of whole figure. The default is 0.6.
    
    Returns
    -------
    cax : matplotlib.axes.Axes
        Colorbar axes
    """
    fig.subplots_adjust(bottom=0.2)
    cax = plt.axes([0.5 - (width / 2), 0.09, width, 0.02])
    # Add colorbar
    cbar = fig.colorbar(ScalarMappable(norm, cmap), cax,
                        orientation='horizontal')
    cbar.set_label(label)
    return cax


# =============================================================================
# MATHEMATICAL / DATA UTILITY FUNCTIONS
# =============================================================================

def weighted_quantile(df, val, weight, quantile=0.5):
    """
    Calculate the quantile of a pandas column weighted by another column.
    
    Parameters
    ----------
    df : pandas.DataFrame
    val : str
        Name of values column.
    weight : str
        Name of weights column.
    quantile : float, optional
        The quantile to calculate. Must be in [0,1]. The default is 0.5.
    
    Returns
    -------
    wq : float
        The weighted quantile of the dataframe column.
    """
    if quantile >= 0 and quantile <= 1:
        if df.shape[0] == 0:
            return np.nan
        else:
            df_sorted = df.sort_values(val)
            cumsum = df_sorted[weight].cumsum()
            cutoff = df_sorted[weight].sum() * quantile
            wq = df_sorted[cumsum >= cutoff][val].iloc[0]
            return wq
    else:
        raise ValueError("Quantile must be in range [0,1].")
