# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 15:15:51 2019

@author: karlen
"""

from TimingCalib.Device import Device

class PMT(Device):
    """ A photomultiplier tube """

    # class properties:
    prop_mean = {} # dictionary of properties of different kinds of PMTs
    prop_scale = {} # scales of variations of properties used to create objects
    prop_var = {} # distribution of variations (Normal by default, "uniform")

    # default properties of PMTs
    # all properties are defined by primitives, so shallow dictionary copy works
    # if new properties are needed, be sure to add to def_prop_mean
    # and def_prop_scale
    def_prop_mean = {'fov':1.0,  # radian field of view (cone half angle)
                     'size':75, # mm diameter
                     'qe':0.5,   # quantum efficiency
                     'gain':2.,  # mV per photoelectron
                     'delay':4.0, # ns mean delay of signal wrt pulse
                     'delay_jitter':0.05 # standard deviation of variation of delay (jitter)
                     }
    def_prop_scale = {'fov':0.001,
                      'size':0.1,
                      'qe':0.01,
                      'gain':0.1,
                      'delay':0.5,
                      'delay_jitter':0.01
                      }
    def_prop_var = {}

    # Standard 3inch PMT:
    p3_prop_mean = def_prop_mean.copy()
    p3_prop_mean['gain'] = 2.1  # 2.1 mV per photoelectron
    p3_prop_scale = def_prop_scale.copy()
    p3_prop_var = def_prop_var.copy()

    prop_mean['P3'] = p3_prop_mean
    prop_scale['P3'] = p3_prop_scale
    prop_var['P3'] = p3_prop_var

    # Modified 3inch PMT:
    p31_prop_mean = def_prop_mean.copy()
    p31_prop_mean['gain'] = 2.5  # 2.5 mV per photoelectron
    p31_prop_scale = def_prop_scale.copy()
    p31_prop_var = def_prop_var.copy()

    prop_mean['P31'] = p31_prop_mean
    prop_scale['P31'] = p31_prop_scale
    prop_var['P31'] = p31_prop_var

    def __init__(self, controller, kind, place_design, place_true):
        super().__init__()
        self.controller = controller
        self.kind = kind

        self.set_properties(PMT, kind)
        self.set_placement(place_design, place_true)
