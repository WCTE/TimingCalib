# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 21:16:14 2019

@author: karlen
"""

from TimingCalib.Device import Device

class LED(Device):
    """ A LED-fibre light source """

    # class properties:
    prop_mean = {} # dictionary of mean properties of different kinds of LEDs
    prop_scale = {} # scale of variations of properties used to create objects
    prop_var = {} # distribution for variations

    # default properties of LEDs
    # all properties are defined by primitives, so shallow dictionary copy works
    # if new properties are needed, be sure to add to def_prop_mean
    # and def_prop_scale
    def_prop_mean = {'cone_angle':1.0,      # 1 radian cone
                     'cone_ring_width':0.0, # if zero, then full cone, otherwise specifies ring width (radians)
                     'intensity':1.E8,    # 0.1 billion photons per ns per volt
                     'intensity_rel_sig':0.01, # relative standard deviation
                     'delay':4.,          # 1 ns delay for leading edge
                     'delay_jitter':0.05,    # 0.05 ns jitter
                     'rise_time':1.0,     # 1 ns rise time
                     'pulse_width':2.0,   # 2 ns mean pulse width (flat top)
                     'pulse_width_jitter':0.05 # standard deviation of above
                     }
    def_prop_scale = {'cone_angle':0.01,
                      'cone_ring_width':0.0,
                      'intensity':2.E7,
                      'intensity_rel_sig':0.,
                      'delay':1.0,
                      'delay_jitter':0.01,
                      'rise_time':0.1,
                      'pulse_width':0.1,
                      'pulse_width_jitter':0.01
                      }
    def_prop_var = {}

    # Standard LED:
    l1_prop_mean = def_prop_mean.copy()
    l1_prop_scale = def_prop_scale.copy()
    l1_prop_var = def_prop_var.copy()

    prop_mean['L1'] = l1_prop_mean
    prop_scale['L1'] = l1_prop_scale
    prop_var['L1'] = l1_prop_var

    # Weak LED:
    l2_prop_mean = def_prop_mean.copy()
    l2_prop_mean['intensity'] = 1.E7  # 0.01 billion photons per ns per volt
    l2_prop_scale = def_prop_scale.copy()
    l2_prop_scale['intensity'] = 2.E6
    l2_prop_var = def_prop_var.copy()

    prop_mean['L2'] = l2_prop_mean
    prop_scale['L2'] = l2_prop_scale
    prop_var['L2'] = l2_prop_var

    def __init__(self, controller, kind, place_design, place_true):
        super().__init__()
        self.kind = kind
        self.controller = controller

        self.set_properties(LED, kind)
        self.set_placement(place_design, place_true)
