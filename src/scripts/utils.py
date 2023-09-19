"""
Functions used by many plotting scripts.
"""

# import math as m
from pathlib import Path
import numpy as np
from numpy.random import default_rng
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, BoundaryNorm, LogNorm
# from matplotlib.ticker import MultipleLocator
from matplotlib.cm import ScalarMappable
from astropy.table import Table
from astropy.io import fits
import vice
import paths
from _globals import ZONE_WIDTH, RANDOM_SEED

# =============================================================================
# DATA IMPORT
# =============================================================================        

def feuillet2019_data(filename):
    r"""
    Obtain the Feuillet et al. (2019) [1]_ data.

    Parameters
    ----------
    filename : ``str``
        The relative path to the file containing the data for a given region.

    Returns
    -------
    age : ``list``
        The mean ages of stars in Gyr in bins of abundance, assuming a gaussian
        distribution in log-age.
    abundance : ``list``
        The abundances at which the mean ages are measured. Same length as
        ``age``.
    age_disp : ``list``
        The standard deviation of the age in Gyr distribution in each bin of
        abundance, assuming a gaussian distribution in log-age. Same length as
        ``age``.
    abundance_disp : ``list``
        The width of the bin in abundance, centered on each element of the
        ``abundance`` array.

    .. [1] Feuillet et al. (2019), MNRAS, 489, 1724
    """
    raw = fits.open(filename)
    abundance = len(raw[1].data) * [0.]
    abundance_disp = len(raw[1].data) * [0.]
    age = len(raw[1].data) * [0.]
    age_disp = [len(raw[1].data) * [0.], len(raw[1].data) * [0.]]
    for i in range(len(raw[1].data)):
        if raw[1].data["nstars"][i] > 15:
            abundance[i] = (raw[1].data["bin_ab"][i] +
                raw[1].data["bin_ab_max"][i]) / 2.
            abundance_disp[i] = (raw[1].data["bin_ab_max"][i] -
                raw[1].data["bin_ab"][i]) / 2.
            age[i] = 10**(raw[1].data["mean_age"][i] - 9) # converts yr to Gyr
            age_disp[0][i] = age[i] - 10**(raw[1].data["mean_age"][i] -
                raw[1].data["age_disp"][i] - 9)
            age_disp[1][i] = 10**(raw[1].data["mean_age"][i] +
                raw[1].data["age_disp"][i] - 9) - age[i]
        else:
            abundance[i] = abundance_disp[i] = float("nan")
            age[i] = age_disp[0][i] = age_disp[1][i] = float("nan")
    return [age, abundance, age_disp, abundance_disp]


# =============================================================================
# DATA UTILITY FUNCTIONS
# =============================================================================

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


def quad_add(arr1, arr2):
    """
    Add two input arrays in quadrature.
    """
    return np.sqrt(arr1**2 + arr2**2)


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
        

def group_by_bins(df, bin_col, bins=10):
    """
    Bin a DataFrame column and group data by those bins.
    
    Parameters
    ----------
    df : pandas.DataFrame
    bin_col : str
        Column name by which to bin the data.
    bins : int or array-like
        If an int is provided, the number of bins between the min and max
        data values. If an array, the bin edges to use. The default is 10.
    
    Returns
    -------
    grouped : pandas.DataFrameGroupBy
        A groupby object that contains information about the groups.
    """
    df = df.copy()
    # Handle different types for "bins" parameter
    if isinstance(bins, int):
        bin_edges = np.linspace(df[bin_col].min(), df[bin_col].max(), bins)
    elif isinstance(bins, (list, np.ndarray, pd.Series)):
        bin_edges = np.array(bins)
    else:
        raise ValueError('Parameter "bins" must be int or array-like.')
    # segment and sort data into bins
    bin_centers = get_bin_centers(bin_edges)
    df.insert(len(df.columns), 
              'bin', 
              pd.cut(df[bin_col], bin_edges, labels=bin_centers))
    # group data by bins
    return df.groupby('bin')


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


def sample_dataframe(df, n, weights=None, reset=True, seed=RANDOM_SEED):
    """
    Randomly sample n unique rows from a pandas DataFrame.

    Parameters
    ----------
    df : pandas DataFrame
    n : int
        Number of random samples to draw
    weights : array, optional
        Probability weights of the given DataFrame
    reset : bool, optional
        If True, reset sample DataFrame index

    Returns
    -------
    pandas DataFrame
        Re-indexed DataFrame of n sampled rows
    """
    if isinstance(df, pd.DataFrame):
        # Number of samples can't exceed length of DataFrame
        n = min(n, df.shape[0])
        # Initialize default numpy random number generator
        rng = default_rng(seed)
        # Randomly sample without replacement
        rand_indices = rng.choice(df.index, size=n, replace=False, p=weights)
        sample = df.loc[rand_indices]
        if reset:
            sample.reset_index(inplace=True, drop=True)
        return sample
    else:
        raise TypeError('Expected pandas DataFrame.')


def median_standard_error(x, B=1000, seed=RANDOM_SEED):
    """
    Use bootstrapping to calculate the standard error of the median.
    
    Parameters
    ----------
    x : array-like
        Data array.
    B : int, optional
        Number of bootstrap samples. The default is 1000.
    
    Returns
    -------
    float
        Standard error of the median.
    """
    rng = np.random.default_rng(seed)
    # Randomly sample input array *with* replacement, all at once
    samples = rng.choice(x, size=len(x) * B, replace=True).reshape((B, len(x)))
    medians = np.median(samples, axis=1)
    # The standard error is the standard deviation of the medians
    return np.std(medians)


def error_fit(df, col, deg, err_col='', bins=30, range=None):
    """
    Fit a polynomial to the error in a parameter as a function of that parameter.
    
    Parameters
    ----------
    df : pandas.DataFrame
    col : str
        Name of parameter column
    deg : int
        Degree of polynomial to fit
    err_col : str, optional
        Name of parameter error column. If '', assumed to be col+'_ERR'. The
        default is None.
    bins : int
        Number of bins to divide data by the parameter. The default is 30.
    range : tuple or NoneType, optional
        Parameter range to include data. If None, include the entire data range.
        The default is None.
    
    Returns
    -------
    p : numpy.ndarray
        Polynomial coefficients, highest power first.
    """
    df = df.copy()
    if err_col in ('', None):
        err_col = col + '_ERR'
    if range is not None:
        df = df[(df[col] >= range[0]) & (df[col] < range[1])]
    
    grouped = group_by_bins(df, col, bins=bins)
    # index labels are bin centers
    medians = grouped[err_col].median()
    # standard error of the medians
    median_errs = grouped[err_col].apply(median_standard_error)
    
    # model fit
    p = np.polyfit(medians.index, medians, deg, w=1/median_errs)
    return p


def model_uncertainty(x, err, how='linear', seed=RANDOM_SEED):
    """
    Apply Gaussian uncertainty to the given data array.
    
    Parameters
    ----------
    x : array-like
        Input (clean) data.
    err : float
        Standard deviation of the Gaussian.
    how : str, optional
        How the uncertainty should be applied to the data. Options are 'linear',
        'logarithmic' or 'log', and 'fractional' or 'frac'. The default is 
        'linear'.
    
    Returns
    -------
    y : array-like
        Noisy data, with same dimensions as x.
    """
    rng = np.random.default_rng(seed)
    noise = rng.normal(loc=0, scale=err, size=x.shape[0])
    if how.lower() == 'linear':
        y = x + noise
    elif how.lower() in ['logarithmic', 'log']:
        y = x * 10 ** noise
    elif how.lower() in ['fractional', 'frac']:
        y = x * (1 + noise)
    else:
        raise ValueError('Parameter "how" must be one of ("linear", ' + 
                         '"logarithmic", "log", "fractional", "frac")')
    return y    

        
# =============================================================================
# VICE MULTIZONE INPUT AND UTILITY FUNCTIONS
# =============================================================================

def multioutput_to_pandas(output_name, data_dir=paths.simulation_outputs, 
                          verbose=False, zone_width=ZONE_WIDTH):
    """
    Convert VICE multizone stars output to pandas DataFrame (slow).

    Parameters
    ----------
    output_name : str
        Path to the .vice directory containing the migration simulation output
    data_dir : str, optional
        Path to the parent directory of all migration outputs. The default is
        '../data/migration_outputs'.
    verbose : bool, optional
        If True, print verbose output to terminal.

    Returns
    -------
    pandas DataFrame
        Parameters of simulated stellar populations including galactic z-height
    """
    full_path = Path(data_dir) / output_name
    if verbose: 
        print('Importing VICE multizone data from %s.vice' % full_path)
    stars = pd.DataFrame(vice.stars(str(full_path)).todict())
    analogdata = pd.read_csv('%s_analogdata.out' % full_path, sep='\t')
    # Limit analogdata to same max time as stars data
    tmax = stars['formation_time'].max()
    analogdata = analogdata[analogdata['time_origin'] <= tmax]
    # Combine relevant data
    stars[['analog_id', 'zfinal']] = analogdata[['analog_id', 'zfinal']]
    # Convert zone to radius
    stars['galr_origin'] = stars['zone_origin'] * zone_width
    stars['galr_final'] = stars['zone_final'] * zone_width
    stars.dropna(how='any', inplace=True)
    return stars


def filter_multioutput_stars(stars, galr_lim=(0, 20), absz_lim=(0, 5),
                             min_mass=1.0, origin=False):
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
        If True, filter by star's original zone instead of final zone. The
        default is False.

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
    subset = stars[(stars[galr_col] >= galr_min) &
                   (stars[galr_col] < galr_max) &
                   (stars['zfinal'].abs() >= absz_min) &
                   (stars['zfinal'].abs() < absz_max) &
                   (stars['mass'] >= min_mass)]
    subset.reset_index(inplace=True)
    return subset


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
                          parent=paths.simulation_outputs):
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
# PLOTTING FUNCTIONS
# =============================================================================


def axes_grid(rows, cols, width=8, xlim=None, ylim=None):
    """
    Set up a blank grid of axes plus a colorbar axis.

    Parameters
    ----------
    rows : int
        Number of rows of axes
    cols : int
        Number of columns of axes
    width : float, optional
        Width of the figure in inches. The default is 8 in.
    xlim : tuple or None, optional
        Limits of x-axis for all axes
    ylim : tuple or None, optional
        Limits of y-axis for all axes

    Returns
    -------
    fig : matplotlib figure
    axs : list of axes
    cax : axis object for colorbar
    """
    fig, axs = plt.subplots(rows, cols, figsize=(width, (width/cols)*rows),
                            sharex=True, sharey=True)
    # Configure plot dimensions
    plt.subplots_adjust(right=0.98, left=0.05, bottom=0.09, top=0.94,
                        wspace=0.05, hspace=0.05)
    # Configure axis limits and ticks (will be applied to all axes)
    axs[0,0].set_xlim(xlim)
    axs[0,0].set_ylim(ylim)
    # axs[0,0].xaxis.set_major_locator(MultipleLocator(0.5))
    # axs[0,0].xaxis.set_minor_locator(MultipleLocator(0.1))
    # axs[0,0].yaxis.set_major_locator(MultipleLocator(0.2))
    # axs[0,0].yaxis.set_minor_locator(MultipleLocator(0.05))
    return fig, axs


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


def highlight_panels(fig, axs, idx, color='#cccccc'):
    """
    Add a colored box behind subplots to make them stand out.
    
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Main plot figure.
    axs : list of matplotlib.axes.Axes
        All subplots of the figure.
    idx : tuple or list of tuples
        Index of the subplot(s) to highlight.
    color : str, optional
        Color to highlight the panel. The default is '#cccccc'.
    """
    # Get spacing between subplots (assuming all are identical)
    # Note: bbox coordinates converted from display to figure
    bbox00 = axs[0,0].get_window_extent().transformed(fig.transFigure.inverted())
    bbox01 = axs[0,1].get_window_extent().transformed(fig.transFigure.inverted())
    bbox10 = axs[1,0].get_window_extent().transformed(fig.transFigure.inverted())
    pad_h = bbox01.x0 - bbox00.x0 - bbox00.width
    pad_v = bbox00.y0 - bbox10.y0 - bbox10.height
    if not isinstance(idx, list):
        idx = [idx]
    for i in idx:
        bbox = axs[i].get_tightbbox().transformed(fig.transFigure.inverted())
        fig.patches.extend([plt.Rectangle((bbox.x0 - pad_h/2, bbox.y0 - pad_v/4),
                                          bbox.x1 - bbox.x0 + pad_h, # width
                                          bbox.y1 - bbox.y0 + pad_v/2, # height
                                          fill=True, color=color, zorder=-1,
                                          transform=fig.transFigure, figure=fig)])
               

# =============================================================================
# DISTRIBUTION STATISTICS
# =============================================================================

def cross_entropy(pk, qk):
    """
    Calculate the cross entropy between two distributions.
    
    The cross entropy is defined as CE = -sum(pk * log(qk)). This function will 
    normalize pk and qk to 1 if needed.
    
    Parameters
    ----------
    pk : numpy.ndarray
        The discrete probability distribution.
    qk : numpy.ndarray
        The probability distribution against which to compute the cross entropy.
        
    Returns
    -------
    CE : float
        The cross entropy of the input distributions
    """
    if pk.shape != qk.shape:
        raise ValueError('Arrays logp and logq must have the same shape.')
    # Normalize distributions
    pk /= np.sum(pk)
    qk /= np.sum(qk)
    # Mask array 0s with smallest non-zero value
    qk[qk == 0] = np.min(qk[qk > 0])
    return -np.sum(pk * np.log(qk))


def kl_divergence(pk, qk, dx):
    r"""
    Calculate the Kullback-Leibler (KL) divergence between two distributions.
    
    For a continuous random variable, KL divergence is defined to be
    $D_{\rm{KL}}(P\parallel Q) = \int_{-\infty}^{\infty} p(x)\log(p(x)/q(x))dx$
    
    Parameters
    ----------
    pk : numpy.ndarray
        Probability density of the observed (true) distribution.
    qk : numpy.ndarray
        Probability density of the model distribution.
    dx : float
        Integration step of the observed variable.
    
    Returns
    -------
    kl : float
        The KL divergence between the two distributions, which is 0 if they
        are identical and positive otherwise.
    """
    # mask zeroes with smallest non-zero value
    pk_nz = np.where(pk != 0, pk, np.min(pk[pk > 0]))
    qk_nz = np.where(qk != 0, qk, np.min(qk[qk > 0]))
    return np.sum(np.where(pk != 0, pk * np.log(pk_nz / qk_nz) * dx, 0))


def kl_div_2D(x, y):
    """
    Compute the Kullback-Leibler divergence between two multivariate samples.
    
    Parameters
    ----------
    x : 2D array (n,d)
        Samples from distribution P, which typically represents the true
        distribution.
    y : 2D array (m,d)
        Samples from distribution Q, which typically represents the approximate
        distribution.
        
    Returns
    -------
    out : float
        The estimated Kullback-Leibler divergence D(P||Q).
        
    References
    ----------
    Pérez-Cruz, F. Kullback-Leibler divergence estimation of
        continuous distributions IEEE International Symposium on Information
        Theory, 2008.
    Source: https://mail.python.org/pipermail/scipy-user/2011-May/029521.html
    """
    from scipy.spatial import cKDTree as KDTree

    # Check the dimensions are consistent
    x = np.atleast_2d(x)
    y = np.atleast_2d(y)

    n,d = x.shape
    m,dy = y.shape

    assert(d == dy)

    # Build a KD tree representation of the samples and find the nearest neighbour
    # of each point in x.
    xtree = KDTree(x)
    ytree = KDTree(y)

    # Get the first two nearest neighbours for x, since the closest one is the
    # sample itself.
    r = xtree.query(x, k=2, eps=.01, p=2)[0][:,1]
    s = ytree.query(x, k=1, eps=.01, p=2)[0]

    return np.log(s/r).sum() * d / n + np.log(m / (n - 1.))


def kde2D(x, y, bandwidth, xbins=100j, ybins=100j, **kwargs):
    """Build 2D kernel density estimate (KDE).

    Parameters
    ----------
    x : array-like
    y : array-like
    bandwidth : float
    xbins : complex, optional [default: 100j]
    ybins : complex, optional [default: 100j]

    Other keyword arguments are passed to sklearn.neighbors.KernelDensity

    Returns
    -------
    xx : MxN numpy array
        Density grid x-coordinates (M=xbins, N=ybins)
    yy : MxN numpy array
        Density grid y-coordinates
    logz : MxN numpy array
        Grid of log-likelihood density estimates
    """
    from sklearn.neighbors import KernelDensity
    # Error handling for xbins and ybins
    if type(xbins) == np.ndarray and type(ybins) == np.ndarray:
        if xbins.shape == ybins.shape:
            if len(xbins.shape) == 2 and len(ybins.shape) == 2:
                xx = xbins
                yy = ybins
            else:
                raise ValueError('Input xbins and ybins must have dimension 2.')
        else:
            raise ValueError('Got xbins and ybins of different shape.')
    elif type(xbins) == complex and type(ybins) == complex:
        # create grid of sample locations (default: 100x100)
        xx, yy = np.mgrid[x.min():x.max():xbins,
                          y.min():y.max():ybins]
    else:
        raise TypeError('Input xbins and ybins must have type complex ' + \
                        '(e.g. 100j) or numpy.ndarray.')

    xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
    xy_train  = np.vstack([y, x]).T

    kde_skl = KernelDensity(kernel='gaussian', bandwidth=bandwidth, **kwargs)
    kde_skl.fit(xy_train)

    # score_samples() returns the log-likelihood of the samples
    logz = kde_skl.score_samples(xy_sample)
    return xx, yy, np.reshape(logz, xx.shape)


def box_smooth(hist, bins, width):
    """
    Box-car smoothing function for a pre-generated histogram.

    Parameters
    ----------
    bins : array-like
        Bins dividing the histogram, including the end. Length must be 1 more
        than the length of hist, and bins must be evenly spaced.
    hist : array-like
        Histogram of data
    width : float
        Width of the box-car smoothing function in data units
    """
    bin_width = bins[1] - bins[0]
    box_width = int(width / bin_width)
    box = np.ones(box_width) / box_width
    hist_smooth = np.convolve(hist, box, mode='same')
    return hist_smooth


def gaussian_smooth(hist, bins, width):
    """
    Box-car smoothing function for a pre-generated histogram.

    Parameters
    ----------
    bins : array-like
        Bins dividing the histogram, including the end. Length must be 1 more
        than the length of hist, and bins must be evenly spaced.
    hist : array-like
        Histogram of data
    width : float
        Standard deviation of the Gaussian in data units
    """
    from scipy.stats import norm
    bin_width = bins[1] - bins[0]
    sigma = int(width / bin_width)
    gaussian = norm.pdf(np.arange(-5*sigma, 5*sigma), loc=0, scale=sigma)
    hist_smooth = np.convolve(hist, gaussian, mode='same')
    return hist_smooth

# =============================================================================
# SCIENCE FUNCTIONS
# =============================================================================

def mean_stellar_mass(age, imf=vice.imf.kroupa, mlr=vice.mlr.larson1974,
                      m_lower=0.08, m_upper=100, dm=0.01):
    """
    Calculate the mean mass of a stellar population of a given age.

    Parameters
    ----------
    age : float
        Stellar age in Gyr
    imf : <function>, optional
        Initial mass function which takes mass in solar masses as an argument.
        The default is vice.imf.kroupa
    mlr : <function>, optional
        Mass-lifetime relation which takes age in Gyr as an argument. The
        default is vice.mlr.larson1974
    m_lower : float, optional
        Lower mass limit on IMF in solar masses. The default is 0.08
    m_upper : float, optional
        Upper mass limit on IMF in solar masses. The default is 100
    dm : float, optional
        IMF integration step in solar masses. The default is 0.01

    Returns
    -------
    float
        Mean mass of stars with lifetime greater than or equal to the given age
        weighted by the IMF
    """
    m_max = min((mlr(age, which='age'), m_upper))
    masses = np.arange(m_lower, m_max + dm, dm)
    dndm = np.array([imf(m) for m in masses])
    weighted_mean = np.average(masses, weights=dndm)
    return weighted_mean
