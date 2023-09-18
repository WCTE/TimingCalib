# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 15:23:50 2019

@author: karlen
"""

from TimingCalib.Device import Device
from TimingCalib.PMT import PMT
from TimingCalib.LED import LED
import numpy as np
from scipy.spatial.transform import Rotation as R


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

    # default properties of (rectangular) MPMTs
    # all properties are defined by primitives, so shallow dictionary copy works
    # if new properties are needed, be sure to add to def_prop_mean 
    # and def_prop_scale
    def_prop_mean = {'clock_offset':0.,  # ns
                     'adc_cf':2., # 2 ADC per mV
                     'size': 500.,  # mm diameter
                     }
    def_prop_scale = {'clock_offset':100.,
                      'adc_cf':0.1,
                      'size': 0.5
                      }
    def_prop_var = {'clock_offset':'uniform'}

    def_pmts = []
    # rectangular pattern of PMTs:
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

    # Standard rectangular MPMT:
    m1_prop_mean = def_prop_mean.copy()
    m1_prop_mean['adc_cf'] = 2.2  # 2.2 ADC channels per mV
    m1_prop_scale = def_prop_scale.copy()
    m1_prop_var = def_prop_var.copy()

    prop_mean['M1'] = m1_prop_mean
    prop_scale['M1'] = m1_prop_scale
    prop_var['M1'] = m1_prop_var

    pmts_design['M1'] = def_pmts
    leds_design['M1'] = def_leds

    dome_pmts = []
    # dome pattern of PMTs:
    number_by_row = [1, 6, 12] # number of PMTs per row
    angle_by_row = [0., -0.297, -0.593] # radians
    dz_by_row = [0., -14.242, -55.724] # mm wrt PMT 0
    distance_by_row = [0., 96.355, 190.594] # mm distance to PMT centres
    #transverse_radius_by_row = [np.sqrt(distance_by_row[i]**2 - dz_by_row[i]**2) for i in range(len(number_by_row))]
    transverse_radius_by_row = []
    for i in range(len(number_by_row)):
        val = np.sqrt(distance_by_row[i]**2 - dz_by_row[i]**2)
        transverse_radius_by_row.append(val)

    # baseplate top surface definition
    dz_to_pmt0 = 246.8 # mm from top surface to PMT0
    long_edge = 264. # mm length of long edge
    long_edge_separation = 528.0 # mm separation of the long edges
    halfs = [long_edge_separation/2., long_edge/2.]
    signs = [1.,-1.]
    base_xy_points = []
    for i in range(8):
        x = halfs[[0, 1, 1, 0, 0, 1, 1, 0][i]] * signs[[0, 0, 1, 1, 1, 1, 0, 0][i]]
        y = halfs[[1, 0, 0, 1, 1, 0, 0, 1][i]] * signs[[1, 1, 1, 1, 0, 0, 0, 0][i]]
        base_xy_points.append([x,y,0.])
    # feedthrough hole definition (to show orientation clearly)
    ft_xy = [195.26, -43.29] # mm centre of feedthrough hole
    ft_diameter = 43. # mm as seen from outside
    feedthough_xy_points = []
    nft = 20
    for i in range(nft):
        theta = 2.*np.pi * i/nft
        x = ft_xy[0] + ft_diameter / 2. * np.cos(theta)
        y = ft_xy[1] + ft_diameter / 2. * np.sin(theta)
        feedthough_xy_points.append([x, y, 0.])

    for i_row,number in enumerate(number_by_row):
        if i_row == 0:
            dome_pmts.append({'kind': 'P3',
                             'loc': [0., 0., dz_to_pmt0],
                             'loc_sig': [1.0, 1.0, 1.0],
                             'rot_axes': 'XZ',
                             'rot_angles': [0., 0.],
                             'rot_angles_sig': [0.01, 100.]})
        else:
            for i_pmt in range(number):
                # start with a PMT located on the mpmt y axis
                loc = [0., transverse_radius_by_row[i_row], dz_by_row[i_row]+dz_to_pmt0]
                # now rotate it around the mpmt z axis
                phi_angle = 2.*np.pi*i_pmt/number
                rot_phi = R.from_euler('Z', phi_angle)
                rot_loc = rot_phi.apply(loc)
                # rotations of the normal defined by 2 extrinsic rotations
                rot_angles = [angle_by_row[i_row],phi_angle]
                dome_pmts.append({'kind': 'P3',
                                  'loc': rot_loc,
                                  'loc_sig': [1.0, 1.0, 1.0],
                                  'rot_axes': 'xz',
                                  'rot_angles': rot_angles,
                                  'rot_angles_sig': [0.01, 0.01]})


    dome_leds = []
    # The dome LED holes are located with respect to the outer top flat surface of the matrix
    matrix_z = 115.85 # mm in zm coordinate of outer top flat surface of the matrix
    # The LED diffuser location is the end of the LED diffuser holder
    diffuser_holder_length = 66.7 # mm distance from matrix surface to end of diffuser holder
    # dome pattern of LED hole locations on surface of matrix (z is wrt outer top flat surface of the matrix)
    led_number_by_row = [3, 3, 6] # number of LED holes per row
    led_angle_by_row = [-0.17, -0.388, -0.707] # radians
    led_dz_by_row = [68.709, 52.644, 8.504] # mm wrt outer top flat surface of the matrix
    led_xm_by_row = [-22.645, -101.328, -44.963]  # mm xm coordinate for first LED hole in the row (numbering azimuthally)
    led_ym_by_row = [39.221, 0., 167.804]  # mm ym coordinate for first LED hole in the row (numbering azimuthally)
    led_transverse_radius_by_row = []
    for i in range(len(number_by_row)):
        val = np.sqrt(led_xm_by_row[i]**2 + led_ym_by_row[i]**2)
        led_transverse_radius_by_row.append(val)

    for i_row,number in enumerate(led_number_by_row):
        if i_row == 0:
            kind = 'L1'
        else:
            kind = 'L2'

        # azimuthal angle of first LED hole in row
        phi_0 = np.arctan2(led_ym_by_row[i_row], led_xm_by_row[i_row])
        for i_led in range(number):
            # start with the first LED hole located at positive azimuth angle wrt ym axis
            # start with a vertically oriented diffuser holder
            loc = [0.,0.,diffuser_holder_length]
            # rotate it about the y-axis
            rot_holder = R.from_euler('Y', led_angle_by_row[i_row])
            rot_loc = rot_holder.apply(loc)
            # translate to the mPMT coordinates on xm axis (had it been located on the xm axis)
            trans_loc = [rot_loc[0] + led_transverse_radius_by_row[i_row], rot_loc[1], rot_loc[2] + matrix_z + led_dz_by_row[i_row]]
            # now rotate it about the mpmt z axis

            phi_angle = 2.*np.pi*i_led/number + phi_0
            rot_phi = R.from_euler('Z', phi_angle)
            rot_trans = rot_phi.apply(trans_loc)
            # rotations of the normal defined by 2 extrinsic rotations
            rot_angles = [led_angle_by_row[i_row], phi_angle]
            dome_leds.append({'kind': kind,
                              'loc': rot_trans,
                              'loc_sig': [1.0, 1.0, 1.0],
                              'rot_axes': 'xz',
                              'rot_angles': rot_angles,
                              'rot_angles_sig': [0.01, 0.01]})

    # Standard dome MPMT:
    m2_prop_mean = def_prop_mean.copy()
    m2_prop_mean['adc_cf'] = 2.2  # 2.2 ADC channels per mV
    m2_prop_scale = def_prop_scale.copy()
    m2_prop_var = def_prop_var.copy()

    prop_mean['M2'] = m2_prop_mean
    prop_scale['M2'] = m2_prop_scale
    prop_var['M2'] = m2_prop_var

    pmts_design['M2'] = dome_pmts
    leds_design['M2'] = dome_leds

    def get_xy_points(self, place_info, feature='base'):
        """Return set of points that shows features on x-y plane (z=0)
        To show feedthrough, set feature='feedthrough'
        """
        device_place = getattr(self, 'place_' + place_info, None)
        location = device_place['loc']

        rot = R.from_euler(device_place['rot_axes'], device_place['rot_angles'])
        xy_points = self.base_xy_points
        if feature == 'feedthrough':
            xy_points = self.feedthough_xy_points
        points = []
        for point in xy_points:
            rotated_point = rot.apply(point)
            points.append(list(np.add(location, rotated_point)))

        return points

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
