"""
This file defines the Gaussian SN Ia DTDs from Strolger et al. (2005).
"""

from .utils import gaussian
from ..._globals import END_TIME

class strolger05(gaussian):
    """
    The SN Ia DTDs from Strolger et al. (2005).
    
    Inherits functionality from utils.gaussian()
    """
    def __init__(self, case='narrow', tmin=0.04, tmax=END_TIME):
        """
        Initialize the Gaussian DTD.
        
        Parameters
        ----------
        case : str, optional
            Which case from Strolger (2004) to use, either `wide` or `narrow`.
            The default is `narrow`.
        tmin, tmax : float, optional
            Dummy parameters to make this class work with multi-zone models.
        """
        if not type(case) == str:
            raise TypeError('Parameter `case` must be a string. Got: %s' 
                            % str(type(case)))
        case = case.lower()
        if case == 'narrow':
            center = 3.4 # Gyr
            stdev = 0.2 * center
        elif case == 'wide':
            center = 3.0 # Gyr
            stdev = 0.5 * center
        else:
            raise ValueError('Parameter `case` must be either `wide` or `narrow`.')
        super().__init__(center=center, stdev=stdev, normalize=True)
        self._name = 'gaussian_%s' % case
        
    @property
    def name(self):
        return self._name
    