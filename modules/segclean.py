#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Identify and exclude likely false positives from matched line lists.
#
# Rob Siverd
# Created:       2019-03-06
# Last modified: 2019-03-06
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
#import datetime as dt
#from dateutil import parser as dtp
#import scipy.linalg as sla
#import scipy.signal as ssig
#import scipy.ndimage as ndi
#import scipy.optimize as opti
#import scipy.interpolate as stp
#import scipy.spatial.distance as ssd

##--------------------------------------------------------------------------##

## Home-brew robust statistics:
try:
    import robust_stats
    reload(robust_stats)
    rs = robust_stats
except ImportError:
    sys.stderr.write("\nError!  robust_stats module not found!\n"
           "Please install and try again ...\n\n")
    sys.exit(1)

##--------------------------------------------------------------------------##
## Simple polynomial fitting in numpy:
def polyfit(x, y, deg):
   if (deg < 1):
      return np.average(y)
   nmat = np.ones_like(y)
   for pow in range(1, deg+1, 1):
      nmat = np.column_stack((nmat, x**pow))
   return np.linalg.lstsq(nmat, y)[0]

## Evaluation of best fit:
def polyval(x, mod):
   z = np.zeros_like(x)
   for i in range(mod.size):
      z += mod[i] * x**i
   return z

##--------------------------------------------------------------------------##
##------------------     Segment-Match Results Clean-Up     ----------------##
##--------------------------------------------------------------------------##

class SegCleanUp(object):

    def __init__(self, min_keep=5, vlevel=0,
            dirty_thresh=0.25, improve_thresh=0.5):
        self.minpts = min_keep
        self.vlevel = vlevel
        self.dirty_thresh = dirty_thresh
        self.improve_thresh = improve_thresh
        return

    # ------------------------------------
    # Set up coordinates for analysis:
    def setup(self, coords1, coords2, params):
        self.coord1 = np.copy(coords1)
        self.coord2 = np.copy(coords2)
        self.ndata  = self.coord1.shape[0]
        self.index  = np.arange(self.ndata, dtype='uint16')
        self.keeper = np.ones(self.ndata, dtype='bool')
        #self.xy_shift = np.zeros(2, dtype='float32')       # not yet known
        self.start_params = np.copy(params)
        self.tuned_params = np.copy(params)
        return

    # ------------------------------------
    # Set up coordinates for analysis:
    def drop_crap(self, order, thresh=2.0):
        npts1 = np.sum(self.keeper)
        while (np.sum(self.keeper) > self.minpts):
            tcoo1, tcoo2, ikept = self._remaining_data()
            rms_0 = self._calc_rms(tcoo1, tcoo2, order)
            too_rms = self._calc_too_rms(tcoo1, tcoo2, order)
            rmsleft = rms_0 - too_rms
            rms_med, rms_iqrn = rs.calc_ls_med_IQR(rmsleft)
            #sys.stderr.write("rms_med: %10.5f\n" % rms_med)
            #sys.stderr.write("rms_iqrn: %10.5f\n" % rms_iqrn)
            #sys.stderr.write("rmsleft: %s\n" % str(rmsleft))
            npoints = rmsleft.size
            worst, w_idx = rmsleft.max(), rmsleft.argmax()
            #sys.stderr.write("worst:   %10.5f\n" % worst)
            #sys.stderr.write("badfrac: %10.5f\n" % badfrac)
            avg_bad = rms_0 / float(npoints)
            #sys.stderr.write("avg_bad: %10.5f\n" % avg_bad)
            if (worst >= thresh * avg_bad):
                self.keeper[ikept[w_idx]] = False
                sys.stderr.write("Dropped match %d\n" % ikept[w_idx])
            else:
                sys.stderr.write("No clear improvement, process converged!\n")
                break
            pass
        npts2 = np.sum(self.keeper)
        if (self.vlevel >= 0):
            sys.stderr.write("Outlier purge complete!\n")
            sys.stderr.write("Kept %d of %d data points.\n" % (npts2, npts1))
        return self._remaining_data()[2]
        #return self._remaining_data()[:2]

    # ------------------------------------ #
    #           Helper Routines            #
    # ------------------------------------ #

    # Fetch whatever is left:
    def _remaining_data(self):
        keep = self.keeper
        return (self.coord1[keep], self.coord2[keep], self.index[keep])

    # Views into coordinate arrays without the specified point:
    @staticmethod
    def _exclude_point(idx, coords1, coords2):
        keep = (np.arange(coords1.shape[0]) != idx)
        return (coords1[keep], coords2[keep])

    # Get residual of fit from points:
    @staticmethod
    def _calc_resid(coords1, coords2, order):
        model = polyfit(coords1, coords2, order)
        return coords2 - polyval(coords1, model)

    @staticmethod
    def _calc_rms(coords1, coords2, order):
        model = polyfit(coords1, coords2, order)
        resid = coords2 - polyval(coords1, model)
        return np.sqrt(np.sum(resid * resid))

    # ----------------------------------------
    # Calculate throw-one-out RMS improvement:
    def _calc_too_rms(self, coo1, coo2, order):
        subset_rms = []
        for i in range(coo1.shape[0]):
            tc1, tc2 = self._exclude_point(i, coo1, coo2)
            subset_rms.append(self._calc_rms(tc1, tc2, order))
        return np.array(subset_rms)



##--------------------------------------------------------------------------##




######################################################################
# CHANGELOG (segclean.py):
#---------------------------------------------------------------------
#
#  2019-03-06:
#     -- Increased __version__ to 0.1.0.
#     -- First created segclean.py.
#
