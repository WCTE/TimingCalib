# -*- coding: utf-8 -*-
"""
Systematic study of various WCDs and the resulting errors in the
timing calibration constants.

WCD made of mPMTs on the 6 surfaces of a cube (cube diameter 8m)

WCD design parameters to be varied:
    - number of PMTs in each mPMT 1 - 5 - 9
    - number of mPMTs on each surface 1 - 5 - 9
    - positional variation of PMTs/LEDs in mPMTs: 0. - 1. - 10.
    - rotational variation of PMTs/LEDs in mPMTs: 0. - 0.01 - 0.10
    - positional variation of mPMTs in WCD: 0. - 5. - 20.
    - rotational variation of mPMTs in WCD: 0. - 0.001 - 0.02
    - variation of LED emission delays: 0.01 - 0.5 - 2.
    - variation of PMT response delays: 0.01 - 0.5 - 2.
    - LED delay jitter: 0.01 - 0.05 - 0.5

Quantities to follow:
     - standard deviation of led-pmt distance residuals (true-design)
     - standard deviation of led-pmt tof residuals (meas-design)
     - standard deviation of the 3 types of timing constant residuals
     - mean deviation of the 3 types of timing constant residuals
     - standard deviations of beacon position measurement

@author: karlen
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from TimingCalib.WCD_calibrator import WCD_calibrator
from TimingCalib.WCD_multilaterator import WCD_multilaterator
from TimingCalib.WCD import WCD
from TimingCalib.MPMT import MPMT
from TimingCalib.PMT import PMT
from TimingCalib.LED import LED

def plot_WCD_layout(wcd, led, wcd_event):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    n_mpmt = len(wcd.mpmts)
    for j_mpmt in range(n_mpmt):
        m_pmt = len(wcd.mpmts[j_mpmt].pmts)
        for k_pmt in range(m_pmt):
            pmt = wcd.mpmts[j_mpmt].pmts[k_pmt]
            [x,y,z], direction = pmt.get_z_orientation('true')
            color = 'r'
            if wcd_event is not None:
                if wcd_event.get_PMT_signal(j_mpmt, k_pmt) is not None:
                    color = 'g'
            ax.scatter(x, y, z, c=color, marker='o', s=1)
    if led is not None:
        [x,y,z], direction = led.get_z_orientation('true')
        ax.scatter(x, y, z, c='b', marker='x', s=5)
    ax.set_xlim3d(-4000, 4000); ax.set_ylim3d(-4000, 4000); ax.set_zlim3d(-4000, 4000);
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    plt.show()

def define_PMT(pmt_kind, delay_sig):
    PMT.prop_mean[pmt_kind] = PMT.def_prop_mean.copy()
    PMT.prop_scale[pmt_kind] = PMT.def_prop_scale.copy()
    PMT.prop_var[pmt_kind] = PMT.def_prop_var.copy()

    PMT.prop_scale[pmt_kind]['delay'] = delay_sig
#    PMT.prop_mean[pmt_kind]['fov'] = np.pi/2.

def define_LED(led_kind, delay_sig):
    LED.prop_mean[led_kind] = LED.def_prop_mean.copy()
    LED.prop_scale[led_kind] = LED.def_prop_scale.copy()
    LED.prop_var[led_kind] = LED.def_prop_var.copy()

    LED.prop_scale[led_kind]['delay'] = delay_sig

def define_MPMT(mpmt_kind, mpmt_conf, pmt_kind, led_kind, loc_sig, rot_angles_sig):
    MPMT.prop_mean[mpmt_kind] = MPMT.def_prop_mean.copy()
    MPMT.prop_scale[mpmt_kind] = MPMT.def_prop_scale.copy()
    MPMT.prop_var[mpmt_kind] = MPMT.def_prop_var.copy()
    
    pmts = []
    for i in range(-1, 2):
        for j in range(-1, 2):
            include = False
            if mpmt_conf == '1':
                include = (i==0 and j==0)
            elif mpmt_conf == '5+':
                include = (i*j == 0)
            elif mpmt_conf == '9':
                include = True
            if include:
                pmts.append({'kind':pmt_kind,
                             'loc':[100.*i, 100.*j, 200.,],
                             'loc_sig':loc_sig,
                             'rot_axes':'XZ',
                             'rot_angles':[0., 0.],
                             'rot_angles_sig':rot_angles_sig})
                
    leds = []
    leds.append({'kind':led_kind,
                     'loc':[50., 0., 200.,],
                     'loc_sig':loc_sig,
                     'rot_axes':'XZ',
                     'rot_angles':[0., 0.],
                     'rot_angles_sig':rot_angles_sig})
    MPMT.pmts_design[mpmt_kind] = pmts
    MPMT.leds_design[mpmt_kind] = leds

def define_cubic_WCD(wcd_kind, n_mpmt_per_side, mpmt_kind, loc_sig, \
                     rot_angles_sig, calib_loc):
    WCD.prop_mean[wcd_kind] = WCD.def_prop_mean.copy()
    WCD.prop_scale[wcd_kind] = WCD.def_prop_scale.copy()
    WCD.prop_var[wcd_kind] = WCD.def_prop_var.copy()

    # make arrangement of up to 9 mpmts on the 3 faces of a cube
    mpmts = []
    rotations = [ [0., np.pi/2., 0.], \
                  [0., -np.pi/2., 0.], \
                  [-np.pi/2., 0., 0.], \
                  [np.pi/2., 0., 0.], \
                  [0., 0., 0.], \
                  [0., np.pi, 0.] ]
    face = -1
    for plane in range(3):
        for sign in range(-1,2,2):
            face += 1
            for j in range(-3,4,1):
                for k in range(-3,4,1):
                    include = False
                    if n_mpmt_per_side == '1':
                        include = (j==0 and k==0)
                    elif n_mpmt_per_side == '5+':
                        include = (j*k == 0 and j*j<2 and k*k<2)
                    elif n_mpmt_per_side == '9+':
                        include = (j*k == 0 and j*j<5 and k*k<5)
                    elif n_mpmt_per_side == '9x':
                        include = (j*j == k*k and j*j<5 and k*k<5)
                    elif n_mpmt_per_side == '9':
                        include = (j*j<2 and k*k<2)
                    elif n_mpmt_per_side == '25':
                        include = (j*j<5 and k*k<5)
                    elif n_mpmt_per_side == '49':
                        include = True
                    if include:
                        j_pnt = (plane+1)%3
                        k_pnt = (plane+2)%3
                        location = [0.]*3
                        location[plane] = sign*4000.
                        location[j_pnt] = 1000.*j
                        location[k_pnt] = 1000.*k
                        mpmts.append({  \
    'kind':mpmt_kind, \
    'loc':location, \
    'loc_sig':loc_sig, \
    'rot_axes':'XYZ', \
    'rot_angles':rotations[face], \
    'rot_angles_sig':rot_angles_sig})

    WCD.mpmts_design[wcd_kind] = mpmts
    
    calibs = []
    # one calibration beacon located somewhere
    define_LED('beacon_LED', 0.)
    LED.prop_mean['beacon_LED']['cone_angle'] = np.pi
    LED.prop_scale['beacon_LED']['cone_angle'] = 0.
    LED.prop_mean['beacon_LED']['delay_jitter'] = 0.
    LED.prop_scale['beacon_LED']['delay_jitter'] = 0.
    
    calibs.append({'kind':'beacon_LED',
                   'loc':calib_loc,
                   'loc_sig':[0.0, 0.0, 0.0],
                   'rot_axes':'XY',
                   'rot_angles':[0., 0.],
                   'rot_angles_sig':[0., 0.]})
    WCD.calibs_design[wcd_kind] = calibs

def calibrate(wcd, print_results):
    """ Run calibration and optionally print results"""
    my_wcd = wcd
    my_wcd_calibrator = WCD_calibrator(my_wcd)
    chi2, n_dof, pars, mean_sigma_t, n_term_t = my_wcd_calibrator.calibrate()

    if print_results:
        n_mpmt = len(my_wcd.mpmts)
        m_pmt = len(my_wcd.mpmts[0].pmts)
        frac_signal = 1. * n_term_t / (n_mpmt * n_mpmt * m_pmt)
        print('chi2 = {0:.1f}, n_dof = {1:d}, n_pars = {2:d}, '.\
              format(chi2, n_dof, len(pars)) + \
              'mean_s_t = {0:.3f}, frac_signal = {1:.3f}, n_t = {2:d}'.\
              format(mean_sigma_t, frac_signal, n_term_t))
        
        print('MPMT delays:')
        mpmt = my_wcd.mpmts[0]
        print('prior mean = {0:.2f}  prior half-range = {1:.2f}' \
              .format(MPMT.prop_mean[mpmt.kind]['clock_offset'], \
                      MPMT.prop_scale[mpmt.kind]['clock_offset']))
        diff = []
        for i in range(1,n_mpmt):
            diff.append(pars[n_mpmt+i-1] - my_wcd.mpmts[i].prop_true['clock_offset'])
        print('Residuals mean = {0:.2f}  Residual sig = {1:.2f}' \
              .format(np.mean(diff), np.std(diff)))
        
        print('LED delays:')
        led = my_wcd.mpmts[0].leds[0]
        print('prior mean = {0:.2f}  prior sig = {1:.2f}' \
              .format(LED.prop_mean[led.kind]['delay'], LED.prop_scale[led.kind]['delay']))
        diff = []
        for i in range(n_mpmt):
            diff.append(pars[i] - my_wcd.mpmts[i].leds[0].prop_true['delay'])
        print('Residuals mean = {0:.2f}  Residual sig = {1:.2f}' \
              .format(np.mean(diff), np.std(diff)))
            
        print('PMT delays:')
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
            
    return chi2, n_dof, pars, mean_sigma_t, n_term_t

def multilaterate(wcd, print_results):
    my_wcd = wcd
    my_wcd_ml = WCD_multilaterator(my_wcd)
    result, cov, ave_sigma_t, chi_2, n_dof = my_wcd_ml.multilaterate(0)
    if print_results:
        calib_loc = my_wcd.calibs[0].place_true['loc']
        print('Beacon location estimate:')
        print('Average t_sigma = {0:0.3f}, chi2 = {1:0.1f}, n_dof = {2:d}'\
              .format(ave_sigma_t, chi_2, n_dof))
        for i in range(3):
            print('x{0:1d}: true={1:.1f}   est={2:.1f} +/- {3:.1f}'\
                  .format(i,calib_loc[i],result[i],np.sqrt(cov[i][i])))
        print('correlations:')
        for i in range(3):        
            buff = ''
            for j in range(3):
                    buff += '{0:0.3f} '.format(cov[i][j]/np.sqrt(cov[i][i]*cov[j][j]))
            print(buff)
        print('eps  true={0:.3f}   est={1:.3f} +/- {2:.3f}'\
                  .format(my_wcd.calibs[0].prop_true['delay'],result[3],np.sqrt(cov[3,3])))
        
    return result, cov, ave_sigma_t, chi_2, n_dof

def repeat_multilaterate(wcd, n_rep):
    my_wcd = wcd
    calib_loc = my_wcd.calibs[0].place_true['loc']
    sum = [0.]*4
    sum2 = [0.]*4
    for i in range(n_rep):
        result, cov, ave_sigma_t, chi_2, n_dof = multilaterate(my_wcd, False)
        for j in range(4):
            sum[j] += result[j]
            sum2[j] += result[j]**2
    print('{0:d} repeated multilaterations:'.format(n_rep))
    # may want to change locations?
    for j in range(3):
        mean = sum[j]/n_rep
        std = np.sqrt(sum2[j]/n_rep - mean**2)
        bias = mean - calib_loc[j]
        err_mean = std/np.sqrt(1.*n_rep)
        print('Ave x[{0:d}] = {1:0.1f}, std = {2:0.1f}, bias = {3:0.1f} +/- {4:0.1f}'\
              .format(j, mean, std, bias, err_mean))
            
def shift_plane(wcd, dist):
    #Systematic shift:
    #move x for all mpmts on first face towards centre, using "read-only" variables!
    my_wcd = wcd
    print('Systematic shift of -x plane applied!')
    n_mpmt = len(my_wcd.mpmts)
    for j in range(int(n_mpmt/6)):
        loc = my_wcd.mpmts[j].place_true['loc']
        loc[0] += dist
        my_wcd.mpmts[j]._place_true['loc'] = loc

def draw_wcd(wcd, led):
    wcd_event = wcd.pulse_LED(led, 8.)
    plot_WCD_layout(wcd,led, wcd_event)    

def case_study():
# CASE STUDY combinations

    photo_MPMT={'perfect':{'loc':0.,'rot':0.}, \
                'good':{'loc':1.,'rot':0.01}, \
                'poor':{'loc':10.,'rot':0.1}}
        
    photo_WCD={'perfect':{'loc':0.,'rot':0.}, \
                'good':{'loc':5.,'rot':0.001}, \
                'poor':{'loc':20.,'rot':0.02}}
    
    timing_LED={'perfect':0.01, \
                'good':0.5, \
                'poor':2.0}
    
    timing_PMT={'perfect':0.01, \
                'good':0.5, \
                'poor':2.0}
        
    # different WCD configurations
    
#    wcd_configs = [['1','1'],['1','5+'],['5+','1'],['5+','5+'],\
#                  ['5+','9'],['9','5+'],['9','9'],['25','5+']]

#    wcd_configs = [['1','1'],['5+','5+'],\
#                   ['9','5+'],['9','9'],['25','5+']]

    wcd_configs = [['9','9'],['25','5+'],['49','5+']]
    wcd_configs = [['49','9']]

    # positions of beacon
    #calib_locs = [[500.,-1200.,2000.], 
    #              [3500.,-1200.,2000.], 
    #              [-3500.,-1200.,2000.]]
    calib_locs = [[500.,-600.,1000.], 
                  [3000.,-600.,1000.], 
                  [-3000.,-600.,1000.]]

    
    calib_rots = [[0.,0.], [0., np.pi/2.], [0., -np.pi/2.]]
    
    # angles of beacon
    
    #qualities = ['good','poor']
    qualities = ['good']
    
    for wcd_config in wcd_configs:
        for photo in qualities:
            #for timing in qualities:
            if 1==1:
                timing = photo
                #for x_shift in [False,True]:
                for x_shift in [True]:
                    print('wcd_config:',wcd_config,\
                          ' photo:',photo,\
                          ' timing:',timing,\
                          ' x-shift:',x_shift)
                    define_PMT('pmt_a', timing_PMT[timing])
                    define_LED('led_a', timing_LED[timing])
                    d_loc = [photo_MPMT[photo]['loc']]*3
                    d_rot = photo_MPMT[photo]['rot']
                    mpmt_config = wcd_config[1]
                    define_MPMT('mpmt_a',mpmt_config,'pmt_a','led_a',\
                                d_loc,[d_rot, 100.])
                    n_mpmt_per_face = wcd_config[0]
                    d_loc = [photo_WCD[photo]['loc']]*3
                    d_rot = [photo_WCD[photo]['rot']]*3
                    define_cubic_WCD('wcd_a',n_mpmt_per_face,'mpmt_a',\
                     d_loc,d_rot,[0.,0.,0.])
    
                    my_wcd = WCD(None, 'wcd_a', {}, {})
    
                    if x_shift:
                        shift_plane(my_wcd, 100.)
    
#                    mean, std = my_wcd.get_pmt_sep_residual()
#                    print('LED-PMT distance residual mean = {0:.2f} mm std = {1:.2f} mm'\
#                          .format(mean,std))
    
                    #draw_wcd(my_wcd, my_wcd.mpmts[0].leds[0])
    
                    chi2, n_dof, pars, mean_sigma_t, n_term_t = calibrate(my_wcd, True)
                    
                    for calib_loc in calib_locs:
                        for calib_rot in calib_rots:
                            my_wcd.calibs[0]._place_true['loc'] = calib_loc
                            # non-isotropic beacon:
                            my_wcd.calibs[0]._prop_true['cone_angle'] = np.pi/2.
                            my_wcd.calibs[0]._place_true['rot_angles'] = calib_rot
                            print('beacon angles =',calib_rot)
                        
                            result, cov, ave_sigma_t, chi_2, n_dof = multilaterate(my_wcd, True)
    
                            repeat_multilaterate(my_wcd, 100)



case_study()

# single case study:

n_mpmt_per_face = '49'
mpmt_config = '5+'

# position of beacon
calib_loc = [3000.,-600.,1000.]

define_PMT('pmt_a', 0.5)
define_LED('led_a', 0.5)
define_MPMT('mpmt_a',mpmt_config,'pmt_a','led_a',[1.,1.,1.],[0.01, 100.])
define_cubic_WCD('wcd_a',n_mpmt_per_face,'mpmt_a',\
                 [1.,1.,1.],[0.01,0.01,0.01],calib_loc)

my_wcd = WCD(None, 'wcd_a', {}, {})

if 1==2:
    shift_plane(my_wcd, 100.)

#mean, std = my_wcd.get_pmt_sep_residual()
#print('LED-PMT distance residual mean = {0:.2f} mm std = {1:.2f} mm'\
#      .format(mean,std))
    
my_wcd.calibs[0]._prop_true['cone_angle'] = np.pi/2.
my_wcd.calibs[0]._place_true['rot_angles'] = [0.,np.pi/2.]

#draw_wcd(my_wcd, my_wcd.calibs[0])

#chi2, n_dof, pars, mean_sigma_t, n_term_t = calibrate(my_wcd, True)
#result, cov, ave_sigma_t, chi_2, n_dof = multilaterate(my_wcd, True)

#repeat_multilaterate(my_wcd, 100)
    