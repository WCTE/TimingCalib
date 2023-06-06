# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 15:23:50 2019

@author: karlen
"""

from TimingCalib.Device import Device
from TimingCalib.PMT import PMT
from TimingCalib.LED import LED

class MPMT(Device):
    """A multi-PMT module"""

    # class properties:
    prop_mean = {} # dictionary of mean properties of different kinds of mPMTs
    prop_scale = {} # scale of variations of properties used to create objects
    prop_var = {} # distribution of variations

    # A dictionary of pmt kinds and placements in the mPMT:
    pmts_design = {}

    # A dictionary of led kinds and placements in the mPMT:
    leds_design = {}

    # default properties of MPMTs
    # all properties are defined by primitives, so shallow dictionary copy works
    # if new properties are needed, be sure to add to def_prop_mean 
    # and def_prop_scale
    def_prop_mean = {'clock_offset':0.,  # ns
                     'adc_cf':2., # 2 ADC per mV
                     }
    def_prop_scale = {'clock_offset':100.,
                      'adc_cf':0.1
                      }
    def_prop_var = {'clock_offset':'uniform'}

    def_pmts = []
    # pattern of PMTs:
    for i in range(-1, 2):
        for j in range(-1, 2):
            # if i*j == 0:  # 5 PMTs per mPMT (cross pattern)
            if i*j != 10:  # 9 PMTs per mPMT
                def_pmts.append({'kind':'P3',
                                 'loc':[100.*i, 100.*j, 200.,],
                                 'loc_sig':[1.0, 1.0, 1.0],
                                 'rot_axes':'XZ',
                                 'rot_angles':[0., 0.],
                                 'rot_angles_sig':[0.01, 100.]})

    def_leds = []
    def_leds.append({'kind':'L1',
                     'loc':[50., 0., 200.,],
                     'loc_sig':[1.0, 1.0, 1.0],
                     'rot_axes':'XZ',
                     'rot_angles':[0., 0.],
                     'rot_angles_sig':[0.01, 100.]})
    def_leds.append({'kind':'L2',
                     'loc':[0., 50., 200.,],
                     'loc_sig':[1.0, 1.0, 1.0],
                     'rot_axes':'XZ',
                     'rot_angles':[0., 0.],
                     'rot_angles_sig':[0.01, 100.]})

    # Standard MPMT:
    m1_prop_mean = def_prop_mean.copy()
    m1_prop_mean['adc_cf'] = 2.2  # 2.2 ADC channels per mV
    m1_prop_scale = def_prop_scale.copy()
    m1_prop_var = def_prop_var.copy()

    prop_mean['M1'] = m1_prop_mean
    prop_scale['M1'] = m1_prop_scale
    prop_var['M1'] = m1_prop_var

    pmts_design['M1'] = def_pmts
    leds_design['M1'] = def_leds

    def __init__(self, controller, kind, place_design, place_true):
        super().__init__()
        self.controller = controller
        self.kind = kind

        self.set_properties(MPMT, kind)
        self.set_placement(place_design, place_true)

        # create the set of PMTs
        self.pmts = self.place_devices(PMT, self.pmts_design, kind)

        # create the set of LEDs
        self.leds = self.place_devices(LED, self.leds_design, kind)
