#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# First crack at a wavelength solution system.
#
# Rob Siverd
# Created:       2019-02-12
# Last modified: 2019-02-12
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.1.0"

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
## Segment matcher for spectral lines:
class WLSegMatch(object):

    def __init__(self, vlevel=0):
        self._vlevel = vlevel
        self._MAXSEG = 1e7

        # Source information:
        self._nsource  = [0, 0]
        self._nsegment = [0, 0]
        self._raw_data = [None, None]
        self._seg_data = [None, None]

        # Bin varieties:
        self._nb_sizes = [30, 35, 40, 45, 50, 55, 60]
        return

    # --------------------------------------------------------
    # Reset functions:
    def _reset_results(self):
        self.all_diffs = None
        self.par_votes = None
        self.seg_votes = None
        self.seg_tally = None
        self.obj_votes = None
        self.obj_tally = None
        self.obj_pairs = None
        return

    # --------------------------------------------------------
    # Set line positions/lists:
    def set_measured_peaks(self, xcoords):
        self._reset_results()
        return self._set_cat_data(0, xcoords)

    def set_reference_lines(self, xcoords):
        self._reset_results()
        return self._set_cat_data(1, xcoords)

    def _set_cat_data(self, which, xcoords):
        self._nsource[which]  = xcoords.size
        self._raw_data[which] = xcoords.astype('float32')
        self._seg_data[which] = self._list_segments(xcoords)
        self._nsegment[which] = self._seg_data[which]['seg'].shape[0]
        return

    def _get_ndata(self):
        return (self._nsegment[0] * self._nsegment[1])

    # --------------------------------------------------------
    # Set line positions/lists:
    def _list_segments(self, xcoords):
        pairs = itt.combinations(range(len(xcoords)), 2)
        ii,jj = zip(*pairs)
        x_sep = np.log10(xcoords[jj,] - xcoords[ii,])
        params = np.column_stack((xcoords[ii,], x_sep))
        indices = np.column_stack((ii, jj)).astype('uint16')
        return {'idx':indices, 'seg':params}

    # Check for existing all-vs-all differences:
    def _have_all_diffs(self):
        if not isinstance(self.all_diffs, np.ndarray):
            return False
        else:
            return True

    # --------------------------------------------------------
    # Create empty vote arrays:
    def _empty_seg_obj_vote_arrays(self):
        seg_votes = np.zeros(self._nsegment).astype('uint16')
        obj_votes = np.zeros(self._nsource, dtype='uint16')
        return (seg_votes, obj_votes)

    # --------------------------------------------------------
    @staticmethod
    def _vote_weighted_avg_params(votes, params):
        tvotes = np.array(votes)
        tparam = np.array(params).T
        return (np.sum(tvotes * tparam, axis=1) / np.sum(tvotes))

    def _print_dither_hist_results(self):
        sys.stderr.write("%7s %7s %10s %10s %10s\n" \
                % ("votes", "margin", "deltaL", "deltaT", "deltaM"))
        for pars,votes,margin in self.dh_results:
            pstring = " ".join(["%10.5f"%x for x in pars])
            sys.stderr.write("%7d %7d %s\n" % (votes, margin, pstring))
        return

    def dither_hist_best_fit(self, ranges, n_bins, n_divs, mode='weighted'):
        # Compute all diffs if necessary:
        if not self._have_all_diffs():
            self.all_diffs = self._calc_all_diffs()

        # Initialize result and vote counters:
        t_results = []
        cumu_seg_votes, cumu_obj_votes = self._empty_seg_obj_vote_arrays()

        # Bin widths, tolerances, and bin center adjustments:
        htols = self._tols_from_ranges_bins(ranges, n_bins)
        bfracs = [np.linspace(0.0, 1.0, nd, endpoint=False) for nd in n_divs]
        nudges = [x*bf for (x, bf) in zip(htols, bfracs)]
        nudge_sets = [np.array(x) for x in itt.product(*nudges)]

        # Analyze the requested combinations:
        for adjustment in nudge_sets:
            adj_diffs = self.all_diffs.reshape(-1, 3) + adjustment[None, :]
            adj_diffs[(adj_diffs > 180.0)] -= 360.0
            dh_pars, mvotes, margin = \
                    self._params_from_hbins(adj_diffs, ranges, n_bins)
            dh_pars = np.array(dh_pars) - adjustment
            t_results.append((dh_pars, mvotes, margin))

            tmp_seg_votes, tmp_obj_votes = \
                    self._make_seg_obj_vtest(dh_pars, htols)
            cumu_seg_votes += tmp_seg_votes
            cumu_obj_votes += tmp_obj_votes

        # Record 
        self.dh_seg_votes = cumu_seg_votes
        self.dh_obj_votes = cumu_obj_votes
        self.dh_results = sorted(t_results, key=lambda x:x[1], reverse=True)

        # Pick or calculate best parameters:
        p_list, v_list, _ = zip(*self.dh_results)
        if (mode == 'weighted'):
            best_params = self._vote_weighted_avg_params(v_list, p_list)
            self.seg_votes = self.dh_seg_votes
            self.obj_votes = self.dh_obj_votes
        elif (mode == 'best'):
            best_params = p_list[0]     # pick best
            self.seg_votes, self.obj_votes = \
                    self._make_seg_obj_vtest(best_params, htols)
        return best_params

##--------------------------------------------------------------------------##




######################################################################
# CHANGELOG (wl_solve_test.py):
#---------------------------------------------------------------------
#
#  2019-02-12:
#     -- Increased __version__ to 0.1.0.
#     -- First created wl_solve_test.py.
#
