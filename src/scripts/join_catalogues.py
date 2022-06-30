"""
This script imports relevant columns from APOGEE value-added catalogues for
astroNN, APOKASC-2, StarHorse, and BACCHUS, joins them to the main APOGEE DR17
catalog, and exports result to 'all_data.csv'. A subset with cuts
based on APOGEE parameters is also exported to 'cut_data.csv'.
"""
from astropy.table import Table
import numpy as np
import pandas as pd
from pathlib import Path

data_dir = '../src/data/APOGEE'
data_path = Path(data_dir)

aspcap_file = 'allStarLite-dr17-synspec.fits'
apokasc_file = 'APOKASC_cat_v6.7.2.fits'
starhorse_file = 'APOGEE_DR17_EDR3_STARHORSE_v2.fits'
astroNN_file = 'apogee_astroNN-DR17.fits'
bacchus_file = 'dr17_nc_abund_v1_0.fits'
apok2_file = 'new_apo_k2_updated.csv'

# Name for output CSV
all_data_file = 'dr17_all_data.csv'
cut_file = 'dr17_cut_data.csv'

def decode(df):
    """
    Decode DataFrame with byte strings into ordinary strings.
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


print('Importing ASPCAP catalog...')
# Import APOGEE allStar data
aspcap_table = Table.read(data_path / aspcap_file, format='fits', hdu=1)
# Combine select abundances
aspcap_table['C_N'] = aspcap_table['C_FE'] - aspcap_table['N_FE']
aspcap_table['C_N_ERR'] = quad_add(aspcap_table['C_FE_ERR'], aspcap_table['N_FE_ERR'])
# Separate paramflags into individual columns
for i in range(len(aspcap_table['PARAMFLAG'][0])):
    aspcap_table['PARAMFLAG' + str(i)] = aspcap_table['PARAMFLAG'][:,i]
# Filter out multidimensional columns
cols = [name for name in aspcap_table.colnames if len(aspcap_table[name].shape) <= 1]
aspcap_df = decode(aspcap_table[cols].to_pandas())
# Relace NaN
aspcap_df.replace(99.999, np.nan, inplace=True)
# Replace '' with 'none' in columns of type 'object'
aspcap_df.replace('', 'none', inplace=True)
# Weed out bad flags
fatal_flags = (2**23) # STAR_BAD
aspcap_df = aspcap_df[aspcap_df["ASPCAPFLAG"] & fatal_flags == 0]
aspcap_df.reset_index(inplace=True)

print('Importing astroNN catalog...')
# Import astroNN data
astroNN_table = Table.read(data_path / astroNN_file, format='fits')
# Select and rename astroNN columns
cols = ['APOGEE_ID'] + astroNN_table.colnames[5:49]
cols += ['weighted_dist', 'weighted_dist_error', 'age_lowess_correct', 'age_total_error']
cols += ['galr', 'galphi', 'galz']
cols_new = ['APOGEE_ID'] + ['ASTRONN_'+name for name in cols[1:45]]
cols_new += ['ASTRONN_DIST', 'ASTRONN_DIST_ERR', 'ASTRONN_AGE', 'ASTRONN_AGE_ERR']
cols_new += ['ASTRONN_'+name.upper() for name in cols[49:]]
# Convert to pandas DataFrame
astroNN_df = decode(astroNN_table[cols].to_pandas())
astroNN_df.columns = cols_new
# Combine select abundances
astroNN_df['ASTRONN_O_FE'] = astroNN_df['ASTRONN_O_H'] - astroNN_df['ASTRONN_FE_H']
astroNN_df['ASTRONN_O_FE_ERR'] = quad_add(astroNN_df['ASTRONN_O_H_ERR'], astroNN_df['ASTRONN_FE_H_ERR'])
astroNN_df['ASTRONN_MG_FE'] = astroNN_df['ASTRONN_MG_H'] - astroNN_df['ASTRONN_FE_H']
astroNN_df['ASTRONN_MG_FE_ERR'] = quad_add(astroNN_df['ASTRONN_MG_H_ERR'], astroNN_df['ASTRONN_FE_H_ERR'])
astroNN_df['ASTRONN_C_N'] = astroNN_df['ASTRONN_C_H'] - astroNN_df['ASTRONN_N_H']
astroNN_df['ASTRONN_C_N_ERR'] = quad_add(astroNN_df['ASTRONN_C_H_ERR'], astroNN_df['ASTRONN_N_H_ERR'])
astroNN_df['ASTRONN_TI_FE'] = astroNN_df['ASTRONN_TI_H'] - astroNN_df['ASTRONN_FE_H']
astroNN_df['ASTRONN_TI_FE_ERR'] = quad_add(astroNN_df['ASTRONN_TI_H_ERR'], astroNN_df['ASTRONN_FE_H_ERR'])

print('\tJoining...')
# Join APOGEE and astroNN tables row by row
cat = aspcap_df.join(astroNN_df.drop('APOGEE_ID', axis=1))
# Drop duplicate APOGEE_IDs flagged by EXTRATARG bitmask
duplicate_bitmask = 16
cat = cat[(cat['EXTRATARG'] & duplicate_bitmask) == 0]
# Drop completely identical entries (takes a while)
cat = cat.drop_duplicates()
# cat = cat.sort_values(['APOGEE_ID', 'SNREV'])\
#          .drop_duplicates(subset='APOGEE_ID', keep='last')\
#          .sort_index()
# Drop calibration fields
cat = cat[cat['FIELD'].str.contains('calibration')==False]

print('Importing StarHorse catalog...')
# Import StarHorse data (duplicates are already removed)
starhorse_table = Table.read(data_path / starhorse_file, format='fits')
# Convert percentiles to minus and plus errors
cols = starhorse_table.colnames[8:-2]
median_cols = [name for name in cols if '50' in name]
cols_new = []
for col in median_cols:
    param = col[:-2]
    mean_col = 'STARHORSE_' + param.upper()
    starhorse_table[mean_col] = starhorse_table[col]
    merr_col = 'STARHORSE_' + param.upper() + '_MERR'
    starhorse_table[merr_col] = starhorse_table[col] - starhorse_table[param+'16']
    perr_col = 'STARHORSE_' + param.upper() + '_PERR'
    starhorse_table[perr_col] = starhorse_table[param+'84'] - starhorse_table[col]
    cols_new += [mean_col, merr_col, perr_col]
cols_new = ['APOGEE_ID'] + cols_new
# Convert to DataFrame
starhorse_df = decode(starhorse_table[cols_new].to_pandas()).set_index('APOGEE_ID')
print('\tJoining...')
# Join StarHorse to ASPCAP and astroNN
cat = cat.join(starhorse_df, on='APOGEE_ID')

print('Importing BACCHUS catalog...')
# Import BACCHUS carbon isotope ratios
bacchus_table = Table.read(data_path / bacchus_file, format='fits', hdu=1)
cols = ['APOGEE_ID', 'FIELD', 'SNR', 'C12C13', 'C12C13_ERR_MEAS', 'C12C13_ERR_EMP']
bacchus_df = decode(bacchus_table[cols].to_pandas())
# Drop duplicate APOGEE_IDs, keeping highest SNR
bacchus_df = bacchus_df.sort_values(['APOGEE_ID', 'SNR'])\
                       .drop_duplicates(subset='APOGEE_ID', keep='last')\
                       .sort_index()
bacchus_df = bacchus_df.set_index(['APOGEE_ID', 'FIELD']).drop('SNR', axis=1)
print('\tJoining...')
# Join BACCHUS to other catalogues
cat = cat.join(bacchus_df, on=['APOGEE_ID', 'FIELD'])

print('Importing APOKASC2 catalog...')
# Import APOKASC2 data
apokasc_table = Table.read(data_path / apokasc_file, format='fits')
# Select columns to keep
cols = ['KEPLER_INT', '2MASS_ID'] + \
       [name for name in apokasc_table.colnames if 'APOKASC2_' in name]
apokasc_df = decode(apokasc_table[cols].to_pandas())
# Rename index column to match
apokasc_df.rename(columns={'2MASS_ID': 'APOGEE_ID'}, inplace=True)
apokasc_df.set_index('APOGEE_ID', inplace=True)
# Replace NaN values
apokasc_df.replace([np.inf, -np.inf, -9999., -9999.99, -999., -999.99],
                   np.nan, inplace=True)
apokasc_df.replace('-9999', 'none', inplace=True)
print('\tJoining...')
# Join APOKASC2 to other catalogues
cat = cat.join(apokasc_df, on='APOGEE_ID')

print('Importing APO-K2 catalog...')
# Import APO-K2 catalogue
cols = ['APOGEE_ID', 'SNREV',
        'new_mass', 'new_mass_err', 'new_radius', 'new_radius_err']
apok2_df = pd.read_csv(data_path / apok2_file, usecols=cols)
apok2_df.rename(columns={'new_mass': 'APOK2_MASS',
                         'new_mass_err': 'APOK2_MASS_ERR',
                         'new_radius': 'APOK2_RADIUS',
                         'new_radius_err': 'APOK2_RADIUS_ERR'},
                inplace=True)
# Drop duplicate APOGEE_IDs, keeping highest SNREV
apok2_df = apok2_df.sort_values(['APOGEE_ID', 'SNREV'])\
                   .drop_duplicates(subset='APOGEE_ID', keep='last')\
                   .sort_index()
apok2_df = apok2_df.drop('SNREV', axis=1).set_index('APOGEE_ID')
print('\tJoining...')
# Join APO-K2 to other catalogues
cat = cat.join(apok2_df, on='APOGEE_ID')

# After all joins, replace NaN in 'object' columns with 'none' for consistency
cat[cat.select_dtypes(include=object).columns] = cat.select_dtypes(include=object).replace(np.nan, 'none')

print('Exporting CSVs...')
# Export CSV
cat.to_csv(data_path / all_data_file, index=False)

# Cut in log(g) and T_eff to remove discrepant metallicities
cut_data = cat[(cat['LOGG'] < 4) & (cat['TEFF'] < 6000) &
               ~((cat['LOGG'] > 3) & (cat['TEFF'] < 4000))]
cut_data.to_csv(data_path / cut_file, index=False)

print('Done!')
