# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 05:35:00 2020

@author: karlen
"""

class WCD_event():
    """A water cherenkov detector event consisting of a list of PMT_signals"""

    def __init__(self, wcd, event_type, device):
        self.wcd = wcd
        self.event_type = event_type
        self.device = device
        n_mpmt = len(wcd.mpmts)
        self.pmt_signals = [None]*n_mpmt
        for i_mpmt in range(n_mpmt):
            m_pmt = len(wcd.mpmts[i_mpmt].pmts)
            self.pmt_signals[i_mpmt] = [None]*m_pmt

    def add_PMT_signal(self, pmt_signal, index_mpmt, index_pmt):
        """add a PMT signal to list"""
        self.pmt_signals[index_mpmt][index_pmt] = pmt_signal

    def get_PMT_signal(self, index_mpmt, index_pmt):
        """get a PMT signal"""
        return self.pmt_signals[index_mpmt][index_pmt]
