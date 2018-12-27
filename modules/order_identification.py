#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Useful routines to assist in order identification. Primary tasks
# include 
# 1) identifying which traces correspond to the same orders
# 2) approximate wavelength identification from Argon lines
#
# Rob Siverd
# Created:       2018-12-27
# Last modified: 2018-12-27
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.0.5"

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
import copy
import time
import numpy as np
#from numpy.lib.recfunctions import append_fields
#import scipy.linalg as sla
#import scipy.signal as ssig
#import scipy.ndimage as ndi
#import scipy.optimize as opti
#import scipy.interpolate as stp
#import scipy.spatial.distance as ssd
#from functools import partial
#from collections import OrderedDict
#import multiprocessing as mp
#import pandas as pd
#import statsmodels.api as sm
#import statsmodels.formula.api as smf
#from statsmodels.regression.quantile_regression import QuantReg
#import theil_sen as ts
#import window_filter as wf
#import itertools as itt

##--------------------------------------------------------------------------##

## Home-brew robust statistics:
#try:
#    import robust_stats
#    reload(robust_stats)
#    rs = robust_stats
#except ImportError:
#    sys.stderr.write("\nError!  robust_stats module not found!\n"
#           "Please install and try again ...\n\n")
#    sys.exit(1)

## Home-brew KDE:
#try:
#    import my_kde
#    reload(my_kde)
#    mk = my_kde
#except ImportError:
#    sys.stderr.write("\nError!  my_kde module not found!\n"
#           "Please install and try again ...\n\n")
#    sys.exit(1)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Utility class to pair orders by ThAr content correlation:
class AdjacentDoubleSimilarity(object):

    def __init__(self):
        self._nearly_zero = 0.1
        return

    # Normalize spectrum by median prior to checking correlation:
    def normalize_spectrum(self, spec_data):
        norm_data = copy.deepcopy(spec_data)
        for odata in norm_data:
            odata['spec'] /= np.median(odata['spec'])
            too_low = (odata['spec'] < self._nearly_zero)
            odata['spec'][too_low] = self._nearly_zero
            pass
        return norm_data

    # ---------------------------------------------------

    #def correl_sim(flx1, flx2, pctcut=75):
    #    return np.sum(flx1 * flx2)

    # Identify overlapping X coordinates:
    @staticmethod
    def _get_overlap(xpix1, xpix2):
        xleft  = max(xpix1[ 0], xpix2[ 0])
        xright = min(xpix1[-1], xpix2[-1])
        which1 = (xleft <= xpix1) & (xpix1 <= xright)
        which2 = (xleft <= xpix2) & (xpix2 <= xright)
        return (xleft, xright, which1, which2)

    @staticmethod
    def _best_correl_offset_and_score(data1, data2, shiftmin, shiftmax,
            pctcut):
        cor_scores = []
        for nn in range(shiftmin, shiftmax+1):
            cor_scores.append(np.sum(np.roll(data1, nn) * data2))
        best_offset = np.argmax(cor_scores)
        return (best_offset, cor_scores[best_offset])

    # Correlation-based similarity of two spectra:
    def correl_calc_similarity(self, spec1, spec2,
                        pctcut=75, nudgemin=0, nudgemax=20):
        xl, xr, msk1, msk2 = self._get_overlap(spec1['xpix'], spec2['xpix'])
        overlap_flx1 = spec1['spec'][msk1]
        overlap_flx2 = spec2['spec'][msk2]
    
        return self._best_correl_offset_and_score(overlap_flx1, overlap_flx2,
                                    nudgemin, nudgemax, pctcut)

    # Compute similarity of all adjacent traced spectra:
    def adj_similarity_scores(self, thar_norm):
        ncomparisons = len(thar_norm) - 1
        #match_scores = np.zeros(ncomparisons, dtype='float')
        #match_shifts = np.zeros(ncomparisons, dtype='int')
        #trace_indexes = []
        match_summary = []
        for i in range(len(thar_norm) - 1):
            sys.stderr.write("order: %3d\n" % i)
            #trace_indexes.append((i, i+1))
            mshift, mscore = \
                    self.correl_calc_similarity(thar_norm[i], thar_norm[i+1])
            tridx = (i, i+1)
            result = {'tridx':(i, i+1), 'mshift':mshift, 'mscore':mscore}
            #match_summary.append((mshift, mscore))
            #match_summary.append((tridx, result))
            match_summary.append(result)
           #match_shifts[i], match_scores[i] = \
           #        self.correl_calc_similarity(thar_norm[i], thar_norm[i+1])
        return match_summary

    # resolve pairs in score order:
    #def resolve_trace_pairs(self, match_shifts, match_scores):
    #def resolve_trace_pairs(self, trace_indexes, match_summary):
    def resolve_trace_pairs(self, match_summary):
        n_traces       = len(match_summary) + 1
        unpaired_trace = list(range(n_traces))
        possible_pairs = {x['tridx']:x for x in match_summary}
        #possible_pairs = {(i,i+1):x for i,x in enumerate(match_summary)}
        #possible_pairs = {x:y for x,y in zip(trace_indexes, match_summary)}
        #possible_pairs = {x:y for x,y in match_summary}
        detected_pairs = {}
        #mshifts, match_scores = zip(*match_summary)
        match_scores   = np.array([x['mscore'] for x in match_summary])
        for botidx in np.argsort(match_scores)[::-1]:
            sys.stderr.write("botidx: %d\n" % botidx)
            prev_pair = (botidx - 1, botidx)
            this_pair = (botidx, botidx + 1)
            next_pair = (botidx + 1, botidx + 2)
            all_pairs = prev_pair, this_pair, next_pair
            all_items = (botidx - 1, botidx, botidx + 1)
        
            sys.stderr.write("this_pair: %s\n" % str(this_pair))
            if this_pair in possible_pairs.keys():
                #sys.stderr.write("pair still live!\n")
                # Grab pair information, rule out overlapping matches:
                detected_pairs[this_pair] = possible_pairs[this_pair]
                for toss in all_pairs:
                    possible_pairs.pop(toss, None)
                    pass
        
                # Update list of unpaired orders:
                for tnum in botidx, botidx+1:
                    if tnum in unpaired_trace:
                        unpaired_trace.remove(tnum)
                        pass
                    pass
        
                pass
            else:
                sys.stderr.write("PAIR NOT AVAILABLE!!\n")
                pass
            pass
        return (detected_pairs, unpaired_trace)



##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##



##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##








######################################################################
# CHANGELOG (order_identification.py):
#---------------------------------------------------------------------
#
#  2018-12-27:
#     -- Increased __version__ to 0.0.5.
#     -- First created order_identification.py.
#
