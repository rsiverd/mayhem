#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# First crack at a wavelength solution system.
#
# Rob Siverd
# Created:       2019-02-12
# Last modified: 2019-02-14
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.2.0"

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

    # Calculate parameter range from nbins and tolerance:
    @staticmethod
    def bintol_range(nbins, tolerance):
        return (-(0.5 * nbins) * tolerance, +(0.5 * nbins) * tolerance)

    # Calculate nbins from range and tolerance:
    @staticmethod
    def rangetol_bins(ranges, tolerance):
        return (-(0.5 * nbins) * tolerance, +(0.5 * nbins) * tolerance)

    # Get bin widths and param tols from param range/binning:
    @staticmethod
    def _tols_from_ranges_bins(ranges, n_bins):
        widths = np.array([(hi-lo)/nn for (lo,hi),nn in zip(ranges, n_bins)])
        tols   = 0.5 * widths
        return tols

    # Select matching segments given transform parameters and tolerances:
    def _segments_from_diffs(self, pars, tols):
        if not self._have_all_diffs():
            sys.stderr.write("Differences not yet calculated!\n")
            self.all_diffs = self._calc_all_diffs()
            sys.stderr.write("Diffs calculated on-demand.\n")
            #return None
        #return np.all(np.abs(self.all_diffs - pars) < tols, axis=-1)
        sys.stderr.write("self.all_diffs.shape: %s\n" 
                % str(self.all_diffs.shape))
        #diff_minus_pars = self.all_diffs.reshape(-1, 1) - pars[None, :]
        #ttols = tols[None, :]
        #return np.all(np.abs(diff_minus_pars) < ttols, axis=-1).nonzero()
        return np.all(np.abs(self.all_diffs - pars) < tols, axis=-1).nonzero()

    # --------------------------------------------------------- #
    #               Mid-level Workhorse Routines:               #
    # --------------------------------------------------------- #

    # Look for matches in scale+rotation+mag (numpy method):
    def _count_matches_srm(self, cdata1, cdata2, params, tols):
        segs1 = cdata1['seg']
        segs2 = cdata2['seg']
        hits = self._match_segments_srm(segs1, segs2, params, tols)
        return np.sum(hits)

    # Create empty vote arrays:
    def _empty_seg_obj_vote_arrays(self):
        seg_votes = np.zeros(self._nsegment).astype('uint16')
        obj_votes = np.zeros(self._nsource, dtype='uint16')
        return (seg_votes, obj_votes)

    # Identify corresponding stars:
    #def _make_seg_obj_votes(self, cdata1, cdata2, params, tols):
    def _make_seg_obj_votes(self, params, tols):
        """Voting results reported using RAW catalog indexes."""
        segs1, sidx1 = self._seg_data[0]['seg'], self._seg_data[0]['idx']
        segs2, sidx2 = self._seg_data[1]['seg'], self._seg_data[1]['idx']
        hits = self._match_segments_srm(segs1, segs2, params, tols)
        #sys.stderr.write("hits.shape: %s\n" % str(hits.shape))
        #sys.stderr.write("hits: %s\n" % str(hits))
        sc1, sc2 = hits.nonzero()   # segment indexes in original catalogs
        #sys.stderr.write("sc1: %s\n" % str(sc1)) 
        #sys.stderr.write("sc2: %s\n" % str(sc2)) 
        #sys.stderr.write("sc1.shape: %s\n" % str(sc1.shape)) 
        #sys.stderr.write("sc2.shape: %s\n" % str(sc2.shape)) 
        seg_votes = np.zeros_like(hits).astype('uint16')
        obj_votes = np.zeros(self._nsource, dtype='uint16')
        #sys.stderr.write("matches: %s\n" % str(matches))
        for i1,i2 in zip(sc1, sc2):
            c1seg, c2seg = sidx1[i1], sidx2[i2]
            seg_votes[i1, i2] += 1
            obj_votes[c1seg[0], c2seg[0]] += 1
            obj_votes[c1seg[1], c2seg[1]] += 1
            #sys.stderr.write("c1seg,c2seg: %s,%s\n" % (str(c1seg), str(c2seg)))
        #return hits, matches[0], matches[1]
        return (seg_votes, obj_votes)

    # Tally segment and object votes from existing all_diffs array:
    def _make_seg_obj_vtest(self, params, tols):
        """Voting results reported using RAW catalog indexes."""
        this_func = sys._getframe().f_code.co_name
        sidx1, sidx2 = [sss['idx'] for sss in self._seg_data]
        sc1, sc2 = self._segments_from_diffs(params, tols)
        #seg_votes = np.zeros(self._nsegment).astype('uint16')
        #obj_votes = np.zeros(self._nsource, dtype='uint16')
        seg_votes, obj_votes = self._empty_seg_obj_vote_arrays()
        for i1,i2 in zip(sc1, sc2):
            c1seg, c2seg = sidx1[i1], sidx2[i2]
            seg_votes[i1, i2] += 1
            obj_votes[c1seg[0], c2seg[0]] += 1
            obj_votes[c1seg[1], c2seg[1]] += 1
        return (seg_votes, obj_votes)

    # --------------------------------------------------------
    # Get best (most popular) parameters from segment differences for the
    # specified histogram binning scheme:
    @staticmethod
    def _params_from_hbins(diffs, ranges, n_bins, iname=None):
        #hh, edges = np.histogramdd(diffs.reshape(-1, 3),
        #                                bins=n_bins, range=ranges)
        hh, edges = np.histogram(diffs.flatten(), 
                bins=n_bins[0], range=ranges[0])
        #sys.stderr.write("hh.shape: %s\n" % str(hh.shape))
        peak = np.unravel_index(hh.argmax(), hh.shape)
        idx10 = np.argsort(hh.flatten())[::-1][:10]     # top 10 voted indexes
        top10 = hh.flatten()[idx10]                     # top 10 vote counts
        nvotes = top10[0]                               # votes for hist peak
        margin = (top10[0] - top10[1]) if (top10.size > 1) else nvotes
        #sys.stderr.write("nvotes:  %4d\n" % nvotes)
        #sys.stderr.write("top10:  %s\n" % str(top10))
        #sys.stderr.write("margin:  %4d\n" % margin)

        # optionally save historgram as FITS image:
        if iname:
            qsave(iname, hh, clobber=True)

        # values at bin centers:
        #vparams = [0.5*(ee[ii]+ee[ii+1]) for (ee,ii) in zip(edges, peak)]
        vparams = [0.5 * (edges[peak[0]] + edges[peak[0]+1])]
        return (np.array(vparams), nvotes, margin)
        #return (vparams, nvotes, margin)


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
            sys.stderr.write("dh_pars: %s\n" % str(dh_pars))
            sys.stderr.write("htols: %s\n" % str(htols))

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

    # --------------------------------------------------------- #
    #       Efficient Segment Differencing (all vs all):        #
    # --------------------------------------------------------- #

    # --------------------------------------------------------
    def _calc_all_diffs(self):
        segments1 = self._seg_data[0]['seg']
        segments2 = self._seg_data[1]['seg']
        ndata = self._get_ndata()
        if (ndata > self._MAXSEG):
            sys.stderr.write("Probably a terrible idea ...\n")
            raise MemoryError

        tseg1 = segments1[:, 1]
        tseg2 = segments2[:, 1]

        # Catalog2 difference from catalog1:
        diffs = np.atleast_3d(tseg1[:, None] - tseg2[None, :])
        #diffs = tseg1[:, None] - tseg2[None, :]

        ## Ensure angle in range [-180, 180]:
        #d_ang = diffs[:, :, 1]
        #d_ang[(d_ang >  180.0)] -= 360.0
        #d_ang[(d_ang < -180.0)] += 360.0
        return diffs

    # Check for existing all-vs-all differences:
    def _have_all_diffs(self):
        if not isinstance(self.all_diffs, np.ndarray):
            return False
        else:
            return True

    # --------------------------------------------------------- #
    #             Segment/Object Matching Routines:             #
    # --------------------------------------------------------- #

    def _pairs_from_votes(self, votes, vmin=5, vmfrac=0.6, tally=False):
        # note vote-space dimensions:
        n1, n2 = votes.shape
        maxmatch = min(votes.shape)
        which, swap = (1, False) if (n1 <= n2) else (0, True)

        # Determine vote threshold for pair acceptance. Arg vmin is an absolute
        # floor on the number of acceptable votes. Also set a relative floor
        # based on the peak vote count seen. Use whichever is stricter:
        vmin_rel = np.ceil(vmfrac * np.max(votes))
        minkeep = min(vmin, vmin_rel)

        # Note votes/index of most popular pairing for each row/col:
        max_val = np.max(votes, axis=which)
        max_idx = np.argmax(votes, axis=which)
        matches = []
        pv_list = []
        #sys.stderr.write("swap: %s\n" % str(swap))
        for i,(peakval, peakidx) in enumerate(zip(max_val, max_idx)):
            #sys.stderr.write("i: %2d, %4d@%2d\n" % (i, peakval, peakidx))
            if (peakval >= minkeep):
                matches.append((i, peakidx))
                pv_list.append(peakval)
        if swap:
            matches = [(y, x) for (x, y) in matches]
        return (matches, pv_list)

    def _pairs_from_votes_both(self, votes, vmin=5, vmfrac=0.6):
        # note vote-space dimensions:
        n1, n2 = votes.shape
        maxmatch = min(votes.shape)
        #which, swap = (1, False) if (n1 <= n2) else (0, True)
        axiswap = [(0, True), (1, False)]

        # Determine vote threshold for pair acceptance. Arg vmin is an absolute
        # floor on the number of acceptable votes. Also set a relative floor
        # based on the peak vote count seen. Use whichever is stricter:
        vmin_rel = np.ceil(vmfrac * np.max(votes))
        #minkeep = min(vmin, vmin_rel)
        minkeep = max(vmin, vmin_rel)

        # Note votes/index of most popular pairing for each row/col:
        match_sets = []
        tally_data = []
        for which,swap in axiswap:
            max_val = np.max(votes, axis=which)
            max_idx = np.argmax(votes, axis=which)
            matches = []
            pv_list = []
            #sys.stderr.write("swap: %s\n" % str(swap))
            for i,(peakval, peakidx) in enumerate(zip(max_val, max_idx)):
                #sys.stderr.write("i: %2d, %4d@%2d\n" % (i, peakval, peakidx))
                if (peakval >= minkeep):
                    matches.append((i, peakidx))
                    pv_list.append(peakval)
            if swap:
                matches = [(y, x) for (x, y) in matches]
            match_sets.append(set(matches))
            #tally_data.append(zip(matches, pv_list))
            tally_data.append(dict(zip(matches, pv_list)))
            pass
        #npairs = len(common_pairs)
        #sys.stderr.write("common_pairs: %s\n" % str(common_pairs))
        #sys.stderr.write("common_tally: %s\n" % str(common_tally))
        pairs = list(match_sets[0] & match_sets[1])
        tally = [np.average([x[op] for x in tally_data]) for op in pairs]
        npairs = len(pairs)
        if (npairs > 1):
            pairs, tally = zip(*sorted(zip(pairs, tally)))
        return (pairs, tally)


    # --------------------------------------------------------- #
    #               Extract Matches from Votes                  #
    # --------------------------------------------------------- #

    # --------------------------------------------------------
    def matched_source_indexes(self, vmin=5, vmfrac=0.6):
        if not isinstance(self.obj_votes, np.ndarray):
            sys.stderr.write("No segment votes to analyze!\n")
            return None
        self.obj_pairs, self.obj_tally = \
                self._pairs_from_votes_both(self.obj_votes,
                        vmin=vmin, vmfrac=vmfrac)
        return self.obj_pairs

    def matched_segment_indexes(self, vmin=5, vmfrac=0.6):
        if not isinstance(self.seg_votes, np.ndarray):
            sys.stderr.write("No segment votes to analyze!\n")
            return None
        self.seg_pairs, self.seg_tally = \
                self._pairs_from_votes_both(self.seg_votes,
                        vmin=vmin, vmfrac=vmfrac)
        return self.seg_pairs

    # Get raw data for matched sources:
    def get_matched_coords(self, vmin=5, vmfrac=0.6):
        """Get coordinates of stars that matched between catalogs.

        Returns: (coords1, coords2)
            
        coords1 and coords2 are numpy arrays of shape (Nmatch, 2).
        """
        pairs = self.matched_source_indexes(vmin, vmfrac)
        if (len(pairs) == 0):
            sys.stderr.write("WARNING: no matches!\n")
            return (None, None)
        return [coo[x,] for x,coo in zip(zip(*pairs), self._raw_data)]
        #return [coo[x, :2] for x,coo in zip(zip(*pairs), self._raw_data)]


##--------------------------------------------------------------------------##




######################################################################
# CHANGELOG (wl_solve_test.py):
#---------------------------------------------------------------------
#
#  2019-02-12:
#     -- Increased __version__ to 0.1.0.
#     -- First created wl_solve_test.py.
#
