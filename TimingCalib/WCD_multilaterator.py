# -*- coding: utf-8 -*-
"""


@author: karlen
"""

import numpy as np
from scipy.optimize import least_squares

class WCD_multilaterator():
    """A water cherenkov detector multilaterator
    The position of a pulsed light source is estimated by the
    arrival times in mPMTs"""

    def __init__(self, wcd):
        self.wcd = wcd
        self.beacon_event = None
        self.vc = wcd.light_velocity/wcd.prop_design['refraction_index']

    def multilaterate(self, i_beacon):
        """Pulse the beacon and record arrival times at all PMTs.
        From this data, perform the multilateration and return
        position estimates. The beacon firing time is presumed
        unknown as is treated as a nuiscance parameter."""

        wcd = self.wcd
        beacon_volts = 16.
        beacon = wcd.calibs[i_beacon]
        self.beacon_event = wcd.pulse_LED(beacon, beacon_volts)
        
        sum_sigma_t = 0.
        count_sigma_t = 0
        n_mpmt = len(self.wcd.mpmts)
        for j in range(n_mpmt):
            mpmt = self.wcd.mpmts[j]
            m_pmt = len(mpmt.pmts)
            for k in range(m_pmt):
                pmt = self.wcd.mpmts[j].pmts[k]
                pmt_signal = self.beacon_event.get_PMT_signal(j, k)
                if pmt_signal is not None:
                    pmt_jitter = pmt.prop_true['delay_jitter']
                    sigma_t = np.sqrt(pmt_signal.time_sig**2 +\
                                      pmt_jitter**2) 
                    sum_sigma_t += sigma_t
                    count_sigma_t += 1
        ave_sigma_t = -1.
        if count_sigma_t > 0:
            ave_sigma_t = sum_sigma_t/count_sigma_t
        
        # start with guess in centre of tank
        pars0 = np.array([0., 0., 0., 0.])
#        x_scale = np.array([100.,100.,100.,0.5])
        bounds = ([-4000.,-4000.,-4000.,-1000.], [4000.,4000.,4000.,1000.])
#        result = least_squares(self.rho, pars0, self.jac, bounds,'trf',\
#                               x_scale = x_scale)
        result = least_squares(self.rho, pars0, self.jac, bounds,'trf')
        
        jacob = result.jac
        jacob_t = np.transpose(jacob)
        hess = np.matmul(jacob_t,jacob)
        hess_inv = np.linalg.inv(hess)
        
        pulls = self.rho(result.x)
        chi_2 = 0.
        for pull in pulls:
            chi_2 += pull**2
        n_dof = len(pulls)-4

        return result.x, hess_inv, ave_sigma_t, chi_2, n_dof

        
    def rho(self,pars):
        """ the loss function as defined in scipy.optimize.least_squares"""
        loc = pars[0:3] # presumed location of beacon
        eps = pars[3] # presumed time offset of beacon pulse wrt master
        
        n_mpmt = len(self.wcd.mpmts)
        pulls = []
        for j in range(n_mpmt):
            mpmt = self.wcd.mpmts[j]
            m_pmt = len(mpmt.pmts)
            clock_offset = mpmt.prop_est['clock_offset']
            for k in range(m_pmt):
                pmt = self.wcd.mpmts[j].pmts[k]
                pmt_signal = self.beacon_event.get_PMT_signal(j, k)
                if pmt_signal is not None:
                    response_delay = pmt.prop_est['delay']
                    time = pmt_signal.time + clock_offset - response_delay

                    pmt_loc, pmt_z_axis = pmt.get_z_orientation('design')
                    light_vec = np.subtract(pmt_loc, loc)
                    dist = np.linalg.norm(light_vec)
                    
                    dev = time - eps - dist/self.vc
                    # total std of timing signal:
                    # from pmt jitter and pe stats
                    pmt_jitter = pmt.prop_design['delay_jitter']
                    sigma_t = np.sqrt(pmt_signal.time_sig**2 +\
                                      pmt_jitter**2)                        
                        
                    pulls.append(dev/sigma_t)

        return np.array(pulls)
    
    def jac(self,pars):
        """the Jacobian as defined by scipy.optimize.least_squares"""
        loc = pars[0:3] # presumed location of beacon
        
        n_mpmt = len(self.wcd.mpmts)
        rows = []
        for j in range(n_mpmt):
            mpmt = self.wcd.mpmts[j]
            m_pmt = len(mpmt.pmts)
            for k in range(m_pmt):
                pmt = self.wcd.mpmts[j].pmts[k]
                pmt_signal = self.beacon_event.get_PMT_signal(j, k)
                if pmt_signal is not None:

                    pmt_loc, pmt_z_axis = pmt.get_z_orientation('design')
                    light_vec = np.subtract(pmt_loc, loc)
                    dist = np.linalg.norm(light_vec)

                    # total std of timing signal:
                    # from pmt jitter and pe stats
                    pmt_jitter = pmt.prop_design['delay_jitter']
                    sigma_t = np.sqrt(pmt_signal.time_sig**2 +\
                                      pmt_jitter**2)

                    row = []
                    for i in range(3):
                        row.append(-1.*(loc[i]-pmt_loc[i])/dist/self.vc/sigma_t)
                    row.append(-1./sigma_t)
                    rows.append(row)
                    
        return np.array(rows)
