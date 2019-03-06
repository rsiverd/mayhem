#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# First crack at a wavelength solution system.
#
# Rob Siverd
# Created:       2019-02-12
# Last modified: 2019-03-06
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.2.5"

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
#import scipy.linalg as sla
import scipy.signal as ssig
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
import window_filter as wf
import itertools as itt

from math import ceil
import robust_stats as rs

## NRES extraction tools:
import nres_extraction
reload(nres_extraction)
nrex = nres_extraction


##--------------------------------------------------------------------------##
## Calculate conservative (greedy) wavelength range from spectroscopic order
## numbers and central wavelengths. 
## NOTES:
## * The spread in wavelength per order is (in a single free spectral range): 
##          Δλ_FSR = λ / m, where m is spectroscopic order
## * Half the FSR occurs on each side of the central wavelength.
class ApproxFSR(object):

    def __init__(self):
        return

    @staticmethod
    def greedy_wl_limits(lam_cen, spec_ord, nFSR=1.5):
        lam_full_spread = lam_cen / spec_ord
        lam_difference  = nFSR * 0.5 * lam_full_spread
        return (lam_cen - lam_difference, lam_cen + lam_difference)


##--------------------------------------------------------------------------##
## Find lines in input spectrum for wavelength verification:
class LineFinder(object):

    def __init__(self, boxwid=3):
        self._boxwid = boxwid
        #self._kernel = self.make_pix_kernel(boxwid)
        return

    # Clean/smooth spectrum for line-finding:
    def _make_smooth_flux(self, xpix, flux, boxwid=3, pctile=False):
        smooth_flux = nrex.boxcar_smooth(flux, boxwid)
        if pctile:
            running_25th = wf.window_filter(xpix, xpix, smooth_flux,
                    xwidth=10*boxwid, pctile=25)
            very_clean = nrex.boxcar_smooth(flux / running_25th, boxwid)
        else:
            #masked_flux = self._make_masked_flux(smooth_flux, 9.0)
            #running_25th = nrex.boxcar_smooth(masked_flux, 
            running_mean = nrex.boxcar_smooth(smooth_flux, 
                    self._get_upper_oddval(100*boxwid))
            very_clean = nrex.boxcar_smooth(flux - running_mean, boxwid)
        return very_clean
        #smooth_flux = nrex.boxcar_smooth(flux / running_25th, boxwid)
        #return smooth_flux

    # Helper routines:
    @staticmethod
    def _get_upper_oddval(val):
        upper = ceil(val)
        return upper if ((upper % 2) == 1.) else upper+1
 
    # Masked flux where outliers are replaced by average:
    @staticmethod
    def _make_masked_flux(data, sig_thresh):
        msk_data = np.copy(data)
        med, sig = rs.calc_ls_med_IQR(data)
        junk_pts = ((np.abs(data - med) / sig) <= sig_thresh)
        msk_data[junk_pts] = med
        return msk_data

    # Creates small kernels for 1-D (line) centroids:
    def _make_pix_kernel(self, npix):
        return np.arange(self._get_upper_oddval(npix)) - (npix // 2)

    # Compute centroid offset from peak index based on neighbor flux:
    def _calc_offsets(self, flux, centers, ksize):
        offset = np.int_(self._make_pix_kernel(ksize))
        ntotal = np.array([np.sum(offset*flux[offset+x]) for x in centers])
        ftotal = np.array([np.sum(flux[offset+x]) for x in centers])
        return ntotal / ftotal

    # Driver routine to find and return line positions for wavelength checks:
    def extract_lines_xpix(self, xpix, flux, boxwid=3, shallow=0.01, 
            border=20, pctile=False):
        smooth_flux = self._make_smooth_flux(xpix, flux, 
                boxwid=boxwid, pctile=pctile)
        maxima_idx = ssig.argrelmax(smooth_flux, order=3)[0]
        peak_fluxes = np.array([smooth_flux[x] for x in maxima_idx])
        highest_flx = np.max(peak_fluxes)
        #sys.stderr.write("Peaks found: %6d\n" % len(peak_fluxes))
        #sys.stderr.write("highest_flx: %10.3f\n" % highest_flx)
        flx_cut = shallow * highest_flx
        #sys.stderr.write("flx_cut: %10.5f\n" % flx_cut)
        keepers = (peak_fluxes >= flx_cut).nonzero()[0]
        #sys.stderr.write("Uncut peaks: %6d\n" % keepers.size)
        peak_center_idx = maxima_idx[keepers]
        x_lo, x_hi = border, xpix.max() - border
        safe_edge = (x_lo < peak_center_idx) & (peak_center_idx < x_hi)
        keep_center_idx = peak_center_idx[safe_edge]
        keep_center_pix = xpix[keep_center_idx]
        ctr_offsets = \
                self._calc_offsets(smooth_flux, keep_center_idx, boxwid)
        return (keep_center_pix + ctr_offsets)




##--------------------------------------------------------------------------##
# dummy_x = f0_thar_data[50]['xpix']
# dummy_f = f0_thar_data[50]['spec']
# niceflx = make_smooth_flux(dummy_x, dummy_f)
# linepos, lineidx = find_maxima_xpix(dummy_x, dummy_f)


# clf()
# plot(niceflx)
# [axvline(x) for x in linepos]




##--------------------------------------------------------------------------##




######################################################################
# CHANGELOG (wl_solve_test.py):
#---------------------------------------------------------------------
#
#  2019-02-12:
#     -- Increased __version__ to 0.1.0.
#     -- First created wl_solve_test.py.
#
