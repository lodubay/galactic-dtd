"""
This file contains classes and functions related to star formation efficiency.
"""

import math as m
from utils import Gaussian

class Conroy2022:
    def __init__(self, t1=2.5, t2=3.7, slope=3, tau_star_init=50, tau_star_final=2.36):
        self.t1 = t1
        self.t2 = t2
        self.slope = slope
        self.tau_star_init = tau_star_init
        self.tau_star_final = tau_star_final
        
    def __call__(self, time):
        if time < self.t1:
            return self.tau_star_init
        elif time >= self.t1 and time <= self.t2:
            return self.tau_star_init / ((1 + self.slope*(time-self.t1))**2)
        else:
            return self.tau_star_final


class DoubleBurst:
    def __init__(self, tau_star_init=50, tau_star_final=2.36, tau_star_min=0.5, t1=2.8, slope1=5, 
                 t2=11.5, width2=1, amp2=0.8):
        logistic_amp = tau_star_init / tau_star_final
        self.logistic = Logistic(tau_star_init=logistic_amp, tau_star_final=1, slope=slope1, center=t1)
        self.burst = EfficiencyBurst(plateau=tau_star_final, center=t2, width=width2, amp=amp2)
    
    def __call__(self, time):
        return self.logistic(time) * self.burst(time)
        
    
class EfficiencyBurst:
    def __init__(self, plateau=2, center=11.5, width=1, amp=0.75):
        self.burst = Gaussian(center=center, stdev=width, coeff=amp, normalize=False)
        self.plateau = plateau
        
    def __call__(self, time):
        return self.plateau * (1 - self.burst(time))
        


class Logistic:
    def __init__(self, tau_star_init=2, tau_star_final=1, slope=2, center=1):
        self.plateau = tau_star_init
        self.drop = tau_star_init - tau_star_final
        self.center = center
        self.slope = slope
    
    def __call__(self, time):
        return self.plateau - self.drop / (1 + m.exp(-self.slope * (time - self.center)))
