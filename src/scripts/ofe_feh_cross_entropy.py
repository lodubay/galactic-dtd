"""
Calculate the cross entropy in the [O/Fe]-[Fe/H] plane between VICE stellar
populations from multizone runs and APOGEE.
"""

import numpy as np
import math as m

def main(verbose=False):
    pass


def cross_entropy(pk, qk, base=m.e):
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
    base : float, optional
        Logarithmic base to compute the cross entropy. The default is e.
        
    Returns
    -------
    CE : float
        The cross entropy of the input distributions
    """
    if pk.shape != qk.shape:
        raise ValueError('Arrays pk and qk must have the same shape.')
    pk /= np.sum(pk)
    qk /= np.sum(qk)
    return -np.sum(pk * np.log(qk)) / np.log(base)


if __name__ == '__main__':
    main()
