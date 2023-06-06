# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 13:32:35 2019

@author: karlen
"""
from scipy import stats

class LED_pulse:
    """A pulse of photons emitted at the end of a fiber connected to an LED"""

    def __init__(self, led, voltage):
        self.led = led
        self.intensity = None
        self.delay = None
        self.pulse_width = None

        # make pulse based on the LED properties:
        # intensity:
        mean = led.prop_true['intensity']*voltage
        sigma = mean*led.prop_true['intensity_rel_sig']
        val = mean
        if sigma > 0.:
            val = stats.norm.rvs(mean, sigma)
        setattr(self, 'intensity', val)

        # delay:
        mean = led.prop_true['delay']
        sigma = led.prop_true['delay_jitter']
        val = mean
        if sigma > 0.:
            val = stats.norm.rvs(mean, sigma)
        setattr(self, 'delay', val)
        
        # pulse width:
        mean = led.prop_true['pulse_width']
        sigma = led.prop_true['pulse_width_jitter']
        val = mean
        if sigma > 0.:
            val = stats.norm.rvs(mean, sigma)
        setattr(self, 'pulse_width', val)