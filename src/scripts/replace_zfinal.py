# -*- coding: utf-8 -*-
"""
This script replaces the old zfinal values in VICE multizone outputs with
the new calculation (Gaussian migration scheme only).
"""

from pathlib import Path
import paths
import pandas as pd
import numpy as np
from utils import multioutput_to_pandas
import random
from tqdm import tqdm

def main():
    p = paths.data / 'migration' / 'gaussian'
    for filepath in list(p.glob('**/*_analogdata.out')):
        print(filepath)
        output_name = str(filepath).split('/')[-3].split('_analogdata')[0]
        stars = multioutput_to_pandas(output_name)
        analogdata = pd.read_csv(filepath, sep='\t')
        zfinal_arr = np.arange(-5, 5.01, 0.01)
        samples = []
        for i in tqdm(range(stars.shape[0])):
            hz = sech_scale(data['age'].iloc[i], data['rfinal_sample'].iloc[i])
            s = random.choices(zfinal_arr, cum_weights=sech_cdf(zfinal_arr, hz), k=1)
            samples.append(s)
        analogdata['zfinal'] = np.array(samples)
        
        
def sech_scale(age, rfinal):
    return 0.25 * np.exp((age - 5)/7.0) * np.exp((rfinal-8)/6.0)

def sech_dist(z, scale):
    return 1/(4 * scale) * np.cosh(z / (2 * scale))**-2

def sech_cdf(z, scale):
    return 1 / (1 + np.exp(-z / scale))
        
    
if __name__ == '__main__':
    main()