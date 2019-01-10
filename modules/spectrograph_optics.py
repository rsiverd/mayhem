#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Routines related to spectrograph optical properties and calculations.
#
# Rob Siverd
# Created:       2018-12-29
# Last modified: 2019-01-10
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.1.2"

## Python version-agnostic module reloading:
try:
    reload                              # Python 2.7
except NameError:
    try:
        from importlib import reload    # Python 3.4+
    except ImportError:
        from imp import reload          # Python 3.0 - 3.3

## Modules:
import os
import sys
import time
import numpy as np
#from numpy.lib.recfunctions import append_fields
#import datetime as dt
#from dateutil import parser as dtp
#import scipy.linalg as sla
#import scipy.signal as ssig
#import scipy.ndimage as ndi
#import scipy.optimize as opti
#import scipy.interpolate as stp
#import scipy.spatial.distance as ssd
#from functools import partial
#from collections import OrderedDict
#import multiprocessing as mp
#np.set_printoptions(suppress=True, linewidth=160)
#import pandas as pd
#import statsmodels.api as sm
#import statsmodels.formula.api as smf
#from statsmodels.regression.quantile_regression import QuantReg
#import theil_sen as ts
#import window_filter as wf
#import itertools as itt

##--------------------------------------------------------------------------##
## Calculating index of refraction for several glass types:
class Glass(object):

    def __init__(self):
        self._bcoeffs = {
            'SiO2':np.array([0.67071081e0, 0.433322857e0, 0.877379057e0]),
            'LLF1':np.array([1.21640125e0, 1.33664540e-1, 8.83399468e-1]),
            'PBM2':np.array([1.39446503e0, 1.59230985e-1, 2.45470216e-1]),
             'LF5':np.array([1.28035628e0, 1.6350597e-1,  8.93930112e-1]), }
        self._ccoeffs = {
            'SiO2':np.array([0.00449192312e0, 0.0132812976e0, 95.8899878e0]),
            'LLF1':np.array([8.57807248e-3,   4.20143003e-2,   1.07593060e+2]),
            'PBM2':np.array([1.10571872e-2,   5.07194882e-2,   3.14440142e1]),
             'LF5':np.array([9.29854416e-3,   4.49135769e-2,   1.10493685e2]), }
        return

    def _unknown_glass(self, glasstype):
        if not glasstype in self._bcoeffs.keys():
            sys.stderr.write("Unknown glass type: %s\n" % glasstype)
            return True
        else:
            return False

    # Squared index of refraction for specified wavelengths:
    def refraction_index_squared(self, wlen_um, glasstype):
        if self._unknown_glass(glasstype):
            raise
        lam_um_sq = wlen_um**2
        n_squared = np.ones_like(wlen_um, dtype='float')
        for bb,cc in zip(self._bcoeffs[glasstype], self._ccoeffs[glasstype]):
            n_squared += (lam_um_sq * bb) / (lam_um_sq - cc)
            #n_squared += (wlen_um**2 * bb) / (wlen_um**2 - cc)
        return n_squared

    #def glass_nn_vs_lambda(self, wlen_um, glasstype):
    def refraction_index(self, wlen_um, glasstype):
        if self._unknown_glass(glasstype):
            raise
        return np.sqrt(self.refraction_index_squared(wlen_um, glasstype))
    
    #def deprecated():
    #    bvals = bcoeffs[glasstype]
    #    cvals = ccoeffs[glasstype]
    #    n2 = np.ones_like(wlen_um, dtype='float')
    #    for bb,cc in zip(bvals, cvals):
    #        n2 += (wlen_um**2 * bb) / (wlen_um**2 - cc)
    #    return np.sqrt(n2)
    #    #lamsq = wlen_um[:, None]**2
    #    #tmpnn = np.sum((lamsq * bvals) / (lamsq - cvals), axis=1)
    #    #return np.sqrt(tmpnn + 1.0)

    def glass_dn_dlambda_easy(self, wlen_um, glasstype, epsfrac=1e-5):
        if self._unknown_glass(glasstype):
            raise

        wlen_lower = wlen_um * (1.0 - epsfrac)
        wlen_upper = wlen_um * (1.0 + epsfrac)

        nn_lower = self.refraction_index(wlen_lower, glasstype)
        nn_upper = self.refraction_index(wlen_upper, glasstype)

        return (nn_upper - nn_lower) / (wlen_upper - wlen_lower)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Prism deflection:
def prism_deflection_n(incid_r, apex_r, n):
    ptemp = np.sqrt(n**2 - np.sin(incid_r)**2) * np.sin(apex_r) \
                - np.cos(apex_r) * np.sin(incid_r)
    return incid_r - apex_r + np.arcsin(ptemp)

def wiki_prism_deflection_n(i, A, n):
    return i - A + np.arcsin(n * np.sin(A - np.arcsin(np.sin(i) / n)))

##--------------------------------------------------------------------------##
## Notes on notation, relations, identities, etc.:

# theta_B  --> blaze angle
# lambda_B --> blaze wavelength
# theta    --> facet illumination angle
# alpha    --> angle of incidence
# beta     --> angle of diffraction
# gamma    --> angle of incidence relative to facet normal SPECIFICALLY
#                   in the plane that is parallel to grooves

# Related by:
# alpha = theta_B + theta
# beta  = theta_B - theta

# Note that:
# * theta is SMALL in Littrow configuration, alpha ~= beta

## Grating methods implemented:
class GratingTools(object):

    #groove_spacing_units = {'mm':1.0, 'um':1e3, 'nm':1e6}

    def __init__(self, gratio=None, lines_per_mm=None):
        # Placeholders:
        self._blaze_angle_rad = 0.0
        self._g_spacing_units = {'mm':1.0, 'um':1e3, 'nm':1e6}
        self._groove_spacing = {x:0.0 for x in self._g_spacing_units.keys()}

        # Handle useful inputs:
        if gratio:
            self._blaze_angle_rad = self._calc_blaze_angle(gratio)
        if lines_per_mm:
            for kk,mult in self._g_spacing_units.items():
                self._groove_spacing[kk] = mult / lines_per_mm

        return

    # --------------------------------------
    # Blaze angle and wavelength:

    @staticmethod
    def _calc_blaze_angle(ratio):
        """For R4 grating, use ratio = 4.  Blaze angle returns in RADIANS."""
        return np.arctan(ratio)

    @staticmethod
    def _calc_blaze_wavelength(spec_orders, groove_spacing,
            theta_B_rad, facet_rad, gamma_rad):
        """
        Compute blaze wavelength for an array of order numbers.
        NOTE: wavelength computed in same units as groove spacing.
        """
        ord_vec = np.atleast_1d(spec_orders)
        if not np.all(ord_vec == np.int_(ord_vec)):
            sys.stderr.write("Error: non-integer spec_orders provided!\n")
            return np.nan * spec_orders

        ang_tmp = np.sin(theta_B_rad) * np.cos(facet_rad) * np.cos(gamma_rad)
        return 2.0 * groove_spacing * ang_tmp / np.float_(spec_orders)

    def get_blaze_wavelengths(self, spec_orders, facet_rad=0.0, gamma_rad=0.0,
            units='nm'):
        return self._calc_blaze_wavelength(spec_orders,
                self._groove_spacing[units], self._blaze_angle_rad, 
                facet_rad=facet_rad, gamma_rad=gamma_rad)

    @staticmethod
    def _calc_order_FSR(wlen_cen, spec_orders):
        return wlen_cen / spec_orders
    
    @staticmethod
    def _calc_order_width(wlen_cen, groove_spacing, 
            theta_B_rad, facet_rad, gamma_rad):
        """
        Compute angular span of the order (free spectral range) with central
        wavelength wlen_cen. 
        
        NOTES: 
        * groove_spacing and wlen_cen need same units for correct result!
        * diff_angle_rad should correspond to central wavelength (wlen_cen)
        * wlen_cen supports numpy arrays
        """
        diff_angle_rad = theta_B_rad - facet_rad
        cos_diff_gamma = np.cos(diff_angle_rad) * np.cos(gamma_rad)
        return wlen_cen / (groove_spacing * cos_diff_gamma)

    def get_order_params(self, spec_orders, facet_rad=0.0, gamma_rad=0.0,
            units='nm'):
        """
        Computes central (blaze) wavelength and angular extent of specified
        orders using chosen geometry/angles.

        Outputs:
        central_wlen -- corresponding central wavelengths in requested units
        order_FSR    -- order free spectral ranges (units of central_wlen)
        order_angwid -- order angular width (RADIANS)
        """
        use_spacing = self._groove_spacing[units]
        central_wlen = self._calc_blaze_wavelength(spec_orders, use_spacing,
                self._blaze_angle_rad, facet_rad, gamma_rad)
        order_FSR = self._calc_order_FSR(central_wlen, spec_orders)
        order_angwid = self._calc_order_width(central_wlen, use_spacing,
                self._blaze_angle_rad, facet_rad, gamma_rad)
        return (central_wlen, order_FSR, order_angwid)


##--------------------------------------------------------------------------##




##--------------------------------------------------------------------------##




######################################################################
# CHANGELOG (spectrograph_optics.py):
#---------------------------------------------------------------------
#
#  2018-12-29:
#     -- Increased __version__ to 0.1.0.
#     -- First created spectrograph_optics.py.
#
