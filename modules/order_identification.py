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
