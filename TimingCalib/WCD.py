# -*- coding: utf-8 -*-
"""
WCD: A water Cherenkov detector consisting of MPMTs and calibration devices

properties:
    - mpmts: an array of MPMTs

methods:
    - pulse_LED: produce one LED pulse from one LED in one MPMT and
    transmit that light to all PMTs in all MPMTs

The WCD has a cartesian coordinate system with its origin near the centre of the
array of the MPMTS and "z" axis pointing up. Orientations of all PMTs are
calculated according to this coordinate system by using device placement information.

@author: karlen
"""
import numpy as np

from TimingCalib.Device import Device
from TimingCalib.MPMT import MPMT
from TimingCalib.LED import LED
from TimingCalib.LED_pulse import LED_pulse
from TimingCalib.PMT_signal import PMT_signal
from TimingCalib.WCD_event import WCD_event

from scipy.spatial.transform import Rotation as R

class WCD(Device):
    """A water cherenkov detector"""

    # class properties:
    light_velocity = 2.99792E2 # speed of light in vacuum (mm/ns)

    prop_mean = {} # dictionary of mean properties of different kinds of WCDs
    prop_scale = {} # scale of variations of properties used to create objects
    prop_var = {} # distribution of variations

    # A dictionary of mpmt kinds and placements in the WCD:
    mpmts_design = {}

    # A dictionary of calib_source kinds and placements in the WCD:
    calibs_design = {}

    # default properties of WCDs
    # all properties are defined by primitives, so shallow dictionary copy works
    # if new properties are needed, be sure to add to def_prop_mean
    # and def_prop_scale
    def_prop_mean = {'clock_offset': 0., # ns offset to reference clock
                     'refraction_index':1.4,
                     'absorption_length':80.E3 # 80 m
                     }
    def_prop_scale = {'clock_offset':0.,
                      'refraction_index':0.,
                      'absorption_length':1.E3
                      }
    def_prop_var = {}
    
    def_calibs = []
    # one calibration beacon located somewhere
    def_calibs.append({'kind':'L1',
                       'loc':[1000., 1500., 2000.,],
                       'loc_sig':[0.0, 0.0, 0.0],
                       'rot_axes':'XZ',
                       'rot_angles':[0., 0.],
                       'rot_angles_sig':[0.0, 0.]})

    def_mpmts = []
    # two pairs of mPMTs facing each other:
    for i in range(2):
        for j in range(2):
            sign_x = 2*i - 1
            x_mpmt = sign_x*4000.
            z_mpmt = (2*j-1)*1000.
            def_mpmts.append({'kind':'M1',
                              'loc':[x_mpmt, 0., z_mpmt],
                              'loc_sig':[1.0, 1.0, 1.0],
                              'rot_axes':'XYZ',
                              'rot_angles':[np.pi/2., -1.*sign_x*np.pi/2., 0.],
                              'rot_angles_sig':[0.01, 0.01, 0.01]})

    def3_mpmts = def_mpmts.copy()
    # 3 pairs of mPMTs
    def3_mpmts.append({'kind':'M1',
                       'loc':[0., 0., -4000.],
                       'loc_sig':[1.0, 1.0, 1.0],
                       'rot_axes':'XYZ',
                       'rot_angles':[0., 0., 0.],
                       'rot_angles_sig':[0.01, 0.01, 0.01]})
    def3_mpmts.append({'kind':'M1',
                       'loc':[0., 0., 4000.],
                       'loc_sig':[1.0, 1.0, 1.0],
                       'rot_axes':'XYZ',
                       'rot_angles':[np.pi, 0., 0.],
                       'rot_angles_sig':[0.01, 0.01, 0.01]})
    def3_mpmts.append({'kind':'M1',
                       'loc':[0., -4000., 0.],
                       'loc_sig':[1.0, 1.0, 1.0],
                       'rot_axes':'XYZ',
                       'rot_angles':[-np.pi/2., 0., 0.],
                       'rot_angles_sig':[0.01, 0.01, 0.01]})
    def3_mpmts.append({'kind':'M1',
                       'loc':[0., 4000., 0.],
                       'loc_sig':[1.0, 1.0, 1.0],
                       'rot_axes':'XYZ',
                       'rot_angles':[np.pi/2., 0., 0.],
                       'rot_angles_sig':[0.01, 0.01, 0.01]})

    # Standard WCD:
    w1_prop_mean = def_prop_mean.copy()
    w1_prop_scale = def_prop_scale.copy()
    w1_prop_var = def_prop_var.copy()

    prop_mean['W1'] = w1_prop_mean
    prop_scale['W1'] = w1_prop_scale
    prop_var['W1'] = w1_prop_var

    mpmts_design['W1'] = def_mpmts
    calibs_design['W1'] = def_calibs

    # 6 mPMT Standard WCD:
    w3_prop_mean = def_prop_mean.copy()
    w3_prop_scale = def_prop_scale.copy()
    w3_prop_var = def_prop_var.copy()

    prop_mean['W3'] = w3_prop_mean
    prop_scale['W3'] = w3_prop_scale
    prop_var['W3'] = w3_prop_var

    mpmts_design['W3'] = def3_mpmts
    calibs_design['W3'] = def_calibs

    # Standard WCD in air:
    w2_prop_mean = def_prop_mean.copy()
    w2_prop_mean['refraction_index'] = 1.0
    w2_prop_mean['absorption_length'] = 80.E6
    w2_prop_scale = def_prop_scale.copy()
    w2_prop_var = def_prop_var.copy()

    prop_mean['W2'] = w2_prop_mean
    prop_scale['W2'] = w2_prop_scale
    prop_var['W2'] = w2_prop_var

    mpmts_design['W2'] = def_mpmts
    calibs_design['W2'] = def_calibs

    # WCTE
    ######

    wcte_prop_mean = def_prop_mean.copy()
    wcte_prop_scale = def_prop_scale.copy()
    wcte_prop_var = def_prop_var.copy()

    prop_mean['WCTE'] = w1_prop_mean
    prop_scale['WCTE'] = w1_prop_scale
    prop_var['WCTE'] = w1_prop_var

    wcte_height = 3450.  # mm separation of mPMT opposite baseplates top-bottom
    wcte_diameter = 3654.  # mm separation of mPMT baseplates on opposite wall locations
    wall_vertical_pitch = 580.  # mm separation of mPMT centres for rows
    tb_pitch = 580.  # mm separation of mPMT centres on top and bottom (x and y the same)

    loc_sig = [1.0, 1.0, 1.0]  # mm positioning accuracy
    rot_angles_sig = [0.01, 0.01]  # rad rotational angle positioning accuracy

    # WCTE x-axis aligned with beam, z-axis vertical origin is centre
    # three separate groups of mPMTs: bottom, wall, top

    # Start with mPMT at centre of bottom, then order like the PMTs in an mPMT
    # Next do rows starting at phi = 0

    bottom_mpmts = []
    #################

    loc_centre = [0., -wcte_height / 2., 0.]
    offsets = [[0., 0., 0.]]

    offs = [-tb_pitch, 0, tb_pitch]
    for i in range(8):
        offset = [offs[[1, 2, 2, 2, 1, 0, 0, 0][i]], 0., offs[[2, 2, 1, 0, 0, 0, 1, 2][i]]]
        offsets.append(offset)

    offs = [-2. * tb_pitch, -tb_pitch, 0, tb_pitch, 2. * tb_pitch]
    for i in range(12):
        offset = [offs[[2, 3, 4, 4, 4, 3, 2, 1, 0, 0, 0, 1][i]], 0., offs[[4, 4, 3, 2, 1, 0, 0, 0, 1, 2, 3, 4][i]]]
        offsets.append(offset)

    for offset in offsets:
        location = np.add(loc_centre, offset)
        bottom_mpmts.append({
            'kind': 'M2',
            'loc': location,
            'loc_sig': loc_sig,
            'rot_axes': 'xy',
            'rot_angles': [-np.pi/2., 0.],
            'rot_angles_sig': rot_angles_sig
        })

    wall_mpmts = []
    ###############

    n_col = 18
    for i_row in range(-2, 3):
        loc = [0., i_row * wall_vertical_pitch, wcte_diameter / 2.]
        for j_col in range(n_col):
            phi_angle = 2. * np.pi * j_col / n_col
            rot_phi = R.from_euler('Y', phi_angle)
            rot_loc = rot_phi.apply(loc)
            # rotations of the normal defined by 2 extrinsic rotations
            rot_angles = [np.pi, phi_angle]
            wall_mpmts.append({'kind': 'M2',
                               'loc': rot_loc,
                               'loc_sig': loc_sig,
                               'rot_axes': 'xy',
                               'rot_angles': rot_angles,
                               'rot_angles_sig': rot_angles_sig})

    top_mpmts = []
    #################

    loc_centre = [0., wcte_height / 2., 0.]

    for offset in offsets:
        location = np.add(loc_centre, offset)
        top_mpmts.append({
            'kind': 'M2',
            'loc': location,
            'loc_sig': loc_sig,
            'rot_axes': 'xy',
            'rot_angles': [np.pi / 2., 0.],
            'rot_angles_sig': rot_angles_sig
        })

    mpmts_design['WCTE'] = bottom_mpmts + wall_mpmts + top_mpmts


    def __init__(self, controller, kind, place_design, place_true):
        super().__init__()
        self.controller = controller
        self.kind = kind

        self.set_properties(WCD, kind)
        self.set_placement(place_design, place_true)

        # create the set of mPMTs
        self.mpmts = self.place_devices(MPMT, self.mpmts_design, kind)

        # create some calibration sources (LED type only for now)
        self.calibs = self.place_devices(LED, self.calibs_design, kind)

    def get_pmt_sep_residual(self):
        """ Calculate the mean and standard deviations of the
        led-pmt distance residuals (true-design)
        Do not include pmts within same mpmt as the LED"""
        sum1 = 0. ; sum2 = 0. ; count = 0.
        for l in range(len(self.mpmts)):
            l_mpmt = self.mpmts[l]
            led = l_mpmt.leds[0]
            led_true, direction = led.get_z_orientation('true')
            led_design, direction = led.get_z_orientation('design')
            for j in range(len(self.mpmts)):
                if j>l:
                    j_mpmt = self.mpmts[j]
                    for k in range(len(j_mpmt.pmts)):
                        pmt = j_mpmt.pmts[k]
                        pmt_true, direction = pmt.get_z_orientation('true')
                        pmt_design, direction = pmt.get_z_orientation('design')
                        light_true = np.subtract(pmt_true, led_true)
                        light_design = np.subtract(pmt_design, led_design)
                        dist_true = np.linalg.norm(light_true)
                        dist_design = np.linalg.norm(light_design)
                        resid = dist_true - dist_design
                        sum1 += resid
                        sum2 += resid**2
                        count += 1.
        mean = sum1/count
        std = np.sqrt(sum2/count-mean*mean)
        return mean, std                   

    def pulse_LED(self, led, volts):
        """Pulse an LED
        and propagate its light in the WCD to all PMTs, assigning a signal
        to each one, and returning these data as WCD_event"""

        wcd_event = WCD_event(self, 'LED', led)

        vc = self.light_velocity/self.prop_true['refraction_index']
        mpmt_led = led.controller
        led_pulse = LED_pulse(led, volts)

        # find vector to true location of the fiber endpoint
        # and unit vector that represents true z-axis of LED fiber

        led_loc, led_z_axis = led.get_z_orientation('true')

        for j_mpmt in range(len(self.mpmts)):
            mpmt = self.mpmts[j_mpmt]
            for k_pmt in range(len(mpmt.pmts)):
                pmt = mpmt.pmts[k_pmt]
                pmt.signal = None

                # find vector to true location of the PMT centre
                # and unit vector that represents true z-axis of PMT

                pmt_loc, pmt_z_axis = pmt.get_z_orientation('true')

                # find distance between the two, the LED emission angle
                # and the PMT arrival angle

                light_vec = np.subtract(pmt_loc, led_loc)
                dist = np.linalg.norm(light_vec)
                cos_emission = np.dot(light_vec, led_z_axis)/dist
                cos_arrival = -1. * np.dot(light_vec, pmt_z_axis)/dist

                if cos_emission > np.cos(led.prop_true['cone_angle']) and \
                    cos_arrival > np.cos(pmt.prop_true['fov']):
                    # simple model for light intensity vs cos_emission_angle
                    falloff = 1.
                    # beacon is uniform intensity
                    if led.prop_true['cone_angle'] < np.pi:
                        falloff = cos_emission**2
                    n_photons = (pmt.prop_true['size']/dist)**2 /16. * \
                                falloff * led_pulse.intensity
                    tof = dist/vc
                    delay = mpmt_led.prop_true['clock_offset'] + \
                            led_pulse.delay + tof
                    pmt_signal = PMT_signal(pmt, n_photons, delay)
                    wcd_event.add_PMT_signal(pmt_signal, j_mpmt, k_pmt)

        return wcd_event                    