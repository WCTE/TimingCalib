# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 06:05:07 2020

@author: karlen
"""

import numpy as np

from TimingCalib.LED import LED
from TimingCalib.PMT import PMT

class WCD_calibrator():
    """A water cherenkov detector timing calibrator"""

    def __init__(self, wcd):
        self.wcd = wcd

    def calibrate(self):
        """Pulse LED i_led on all mPMTs and record arrival times at all PMTs.
        From this data, perform calibration, and save calibrations as
        estimates for the device constants: prop_est"""
        
        # force delta_0 = 0. - using variable intended to be read only!
        # calibration does not adjust this constant, 
        # and it is assumed to be zero in the set of linear equations
        self.wcd.mpmts[0]._prop_true['clock_offset'] = 0.

        wcd = self.wcd
        led_volts = 8.
        i_led = 0
        sum_sigma_t = 0.
        count_sigma_t = 0.
        vc = wcd.light_velocity/wcd.prop_design['refraction_index']
        n_mpmt = len(wcd.mpmts)
        m_pmt = len(wcd.mpmts[0].pmts)
        h = []; t = []; t_sig = []; d = []
        for l in range(n_mpmt):
            hl = []; tl = []; tl_sig = []; dl = []
            led = wcd.mpmts[l].leds[i_led]
            led_jitter = led.prop_design['delay_jitter']
            led_loc, led_z_axis = led.get_z_orientation('design')
            wcd_event = wcd.pulse_LED(led, led_volts)
            for j in range(n_mpmt):
                for k in range(m_pmt):
                    pmt = wcd.mpmts[j].pmts[k]

                    pmt_signal = wcd_event.get_PMT_signal(j, k)
                    if pmt_signal is None:
                        hl.append(0.)
                        tl.append(np.inf)
                        tl_sig.append(np.inf)
                        dl.append(np.inf)
                    else:
                        pmt_loc, pmt_z_axis = pmt.get_z_orientation('design')
                        light_vec = np.subtract(pmt_loc, led_loc)
                        dist = np.linalg.norm(light_vec)

                        hl.append(1.)
                        tl.append(pmt_signal.time)
                        # total std of timing signal:
                        # from led jitter, pmt jitter, and pe stats
                        pmt_jitter = pmt.prop_design['delay_jitter']
                        sigma_t = np.sqrt(pmt_signal.time_sig**2 +\
                                          led_jitter**2 + pmt_jitter**2)
                        tl_sig.append(sigma_t)
                        sum_sigma_t += sigma_t
                        count_sigma_t += 1.
                        dl.append(dist)

            h.append(hl); t.append(tl); t_sig.append(tl_sig); d.append(dl)

        # form linear equations to be solved: calculate coefficients for
        # each parameter - parameters are e, d, a
        aa = []; bb = []
        for l in range(n_mpmt):
            e_sum = [0.]*n_mpmt
            d_sum = [0.]*n_mpmt
            a_sum = [0.]*n_mpmt*m_pmt
            bb_sum = 0.
            for j in range(n_mpmt):
                for k in range(m_pmt):
                    p = j*m_pmt + k
                    if h[l][p] > 0:
                        wt = np.power(t_sig[l][p], -2)
                        e_sum[l] += wt
                        d_sum[l] += wt
                        d_sum[j] -= wt
                        a_sum[p] += wt
                        bb_sum += (t[l][p] - d[l][p]/vc)*wt
            led = wcd.mpmts[l].leds[i_led]
            e_mean = LED.prop_mean[led.kind]['delay']
            e_sig = LED.prop_scale[led.kind]['delay']
            we = np.power(e_sig, -2)
            e_sum[l] -= we
            aa_row = e_sum + d_sum[1:] + a_sum
            aa.append(aa_row)
            bb_sum -= e_mean*we
            bb.append(bb_sum)
        # remove one equation, force delta_0 = 0
        for l in range(1, n_mpmt):
            e_sum = [0.]*n_mpmt
            d_sum = [0.]*n_mpmt
            a_sum = [0.]*n_mpmt*m_pmt
            bb_sum = 0.
            for j in range(n_mpmt):
                for k in range(m_pmt):
                    p = j*m_pmt + k
                    if h[l][p] > 0 and j != l:
                        wt = np.power(t_sig[l][p], -2)
                        e_sum[l] += wt
                        d_sum[l] += wt
                        d_sum[j] -= wt
                        a_sum[p] += wt
                        bb_sum += (t[l][p] - d[l][p]/vc)*wt
            for i in range(n_mpmt):
                for k in range(m_pmt):
                    p = l*m_pmt + k
                    if h[i][p] > 0 and i != l:
                        wt = np.power(t_sig[i][p], -2)
                        e_sum[i] -= wt
                        d_sum[i] -= wt
                        d_sum[l] += wt
                        a_sum[p] -= wt
                        bb_sum -= (t[i][p] - d[i][p]/vc)*wt

            aa_row = e_sum + d_sum[1:] + a_sum
            aa.append(aa_row)
            bb.append(bb_sum)

        for j in range(n_mpmt):
            for k in range(m_pmt):
                p = j*m_pmt + k
                e_sum = [0.]*n_mpmt
                d_sum = [0.]*n_mpmt
                a_sum = [0.]*n_mpmt*m_pmt
                bb_sum = 0.
                for i in range(n_mpmt):
                    if h[i][p] > 0:
                        wt = np.power(t_sig[i][p], -2)
                        e_sum[i] += wt
                        d_sum[i] += wt
                        d_sum[j] -= wt
                        a_sum[p] += wt
                        bb_sum += (t[i][p] - d[i][p]/vc)*wt

                pmt = wcd.mpmts[j].pmts[k]
                a_mean = PMT.prop_mean[pmt.kind]['delay']
                a_sig = PMT.prop_scale[pmt.kind]['delay']
                wa = np.power(a_sig, -2)
                a_sum[p] -= wa
                aa_row = e_sum + d_sum[1:] + a_sum
                aa.append(aa_row)
                bb_sum -= a_mean*wa
                bb.append(bb_sum)

        a = np.array(aa)
        b = np.array(bb)
        pars = np.linalg.solve(a, b)

        # Store pars in the parameter estimates
        for i in range(n_mpmt):
            wcd.mpmts[i].leds[0].prop_est['delay'] = pars[i]

        wcd.mpmts[0].prop_est['clock_offset'] = 0.# by definition
        for i in range(1, n_mpmt):
            wcd.mpmts[i].prop_est['clock_offset'] = pars[n_mpmt+i-1]

        for j in range(n_mpmt):
            for k in range(m_pmt):
                p = j*m_pmt + k
                wcd.mpmts[j].pmts[k].prop_est['delay'] = pars[n_mpmt*2-1+p]

        # calculate chi^2:
        chi2_t = 0.
        n_term_t = 0
        for i in range(n_mpmt):
            for j in range(n_mpmt):
                for k in range(m_pmt):
                    p = j*m_pmt + k
                    if h[i][p] > 0:
                        wt = np.power(t_sig[i][p], -2)
                        dev = t[i][p] - d[i][p]/vc  \
                            - wcd.mpmts[i].prop_est['clock_offset'] \
                            + wcd.mpmts[j].prop_est['clock_offset'] \
                            - wcd.mpmts[i].leds[0].prop_est['delay'] \
                            - wcd.mpmts[j].pmts[k].prop_est['delay']
                        n_term_t += 1
                        chi2_t += wt*dev**2

        chi2_e = 0.
        n_term_e = 0
        for l in range(n_mpmt):
            led = wcd.mpmts[l].leds[i_led]
            e_mean = LED.prop_mean[led.kind]['delay']
            e_sig = LED.prop_scale[led.kind]['delay']
            dev = wcd.mpmts[l].leds[0].prop_est['delay'] - e_mean
            we = np.power(e_sig, -2)
            n_term_e += 1
            chi2_e += we*dev**2

        chi2_a = 0.
        n_term_a = 0
        for j in range(n_mpmt):
            for k in range(m_pmt):
                pmt = wcd.mpmts[j].pmts[k]
                a_mean = PMT.prop_mean[pmt.kind]['delay']
                a_sig = PMT.prop_scale[pmt.kind]['delay']
                dev = wcd.mpmts[j].pmts[k].prop_est['delay'] - a_mean
                wa = np.power(a_sig, -2)
                n_term_a += 1
                chi2_a += wa*dev**2

        chi2 = chi2_t + chi2_e + chi2_a
        n_dof = n_term_t + n_term_e + n_term_a - len(b)
        mean_sigma_t = sum_sigma_t/count_sigma_t

        return chi2, n_dof, pars, mean_sigma_t, n_term_t
    