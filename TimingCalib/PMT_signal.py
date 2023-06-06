# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 12:49:46 2019

@author: karlen
"""
import numpy as np
from scipy import stats

class PMT_signal:
    """A signal from a PMT"""

    def __init__(self, pmt, n_photons, delay):
        self.pmt = pmt
        self.amplitude = None
        self.time = None
        self.time_sig = None
        
        mpmt = pmt.controller

        n_pe = n_photons * pmt.prop_true['qe']
        setattr(self, 'amplitude', n_pe*pmt.prop_true['gain'])

        sig_timing = 4./np.sqrt(n_pe)   # simple model for now
        arrival = delay + pmt.prop_true['delay'] + \
            stats.norm.rvs(0, pmt.prop_true['delay_jitter']) + \
            stats.norm.rvs(0, sig_timing) - \
            mpmt.prop_true['clock_offset']
        setattr(self, 'time', arrival)
        setattr(self, 'time_sig', sig_timing)
