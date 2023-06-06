# -*- coding: utf-8 -*-
"""

@author: karlen
"""
import numpy as np

from TimingCalib.WCD import WCD
from TimingCalib.WCD_calibrator import WCD_calibrator
from TimingCalib.LED import LED
from TimingCalib.PMT import PMT
from TimingCalib.MPMT import MPMT

my_wcd = WCD(None, 'W3', {}, {})

#for i in range(10):
#    pmt_signal = my_wcd.pulse_LED()
#
#   print(pmt_signal.time)

# force delta_0 = 0. - using a variable intended to be read only!
my_wcd.mpmts[0]._prop_true['clock_offset'] = 0.

n_mpmt = len(my_wcd.mpmts)
m_pmt = len(my_wcd.mpmts[0].pmts)

my_wcd_calibrator = WCD_calibrator(my_wcd)
chi2, n_dof, pars, mean_sigma_t, n_term_t  = my_wcd_calibrator.calibrate()

print('chi2, n_dof, n_pars, mean_sigma_t, n_term_t')
print(chi2, n_dof, len(pars), mean_sigma_t, n_term_t)

print('MPMT delays')
mpmt = my_wcd.mpmts[0]
print('prior mean = {0:.2f}  prior half-range = {1:.2f}' \
      .format(MPMT.prop_mean[mpmt.kind]['clock_offset'],
              MPMT.prop_scale[mpmt.kind]['clock_offset']))
diff = []
for i in range(1,n_mpmt):
    diff.append(pars[n_mpmt+i-1] - my_wcd.mpmts[i].prop_true['clock_offset'])
print('Residuals mean = {0:.2f}  Residual sig = {1:.2f}' \
      .format(np.mean(diff), np.std(diff)))

print('LED delays')
led = my_wcd.mpmts[0].leds[0]
print('prior mean = {0:.2f}  prior sig = {1:.2f}' \
      .format(LED.prop_mean[led.kind]['delay'], LED.prop_scale[led.kind]['delay']))
diff = []
for i in range(n_mpmt):
    diff.append(pars[i] - my_wcd.mpmts[i].leds[0].prop_true['delay'])
print('Residuals mean = {0:.2f}  Residual sig = {1:.2f}' \
      .format(np.mean(diff), np.std(diff)))
    
print('PMT delays')
pmt = my_wcd.mpmts[0].pmts[0]
print('prior mean = {0:.2f}  prior sig = {1:.2f}' \
      .format(PMT.prop_mean[pmt.kind]['delay'], PMT.prop_scale[pmt.kind]['delay']))
diff = []
for j in range(n_mpmt):
    for k in range(m_pmt):
        p = j*m_pmt + k
        diff.append(pars[n_mpmt*2-1+p] - my_wcd.mpmts[j].pmts[k].prop_true['delay'])
print('Residuals mean = {0:.2f}  Residual sig = {1:.2f}' \
      .format(np.mean(diff), np.std(diff)))