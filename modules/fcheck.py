#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Fit comparison/checking
#
# Rob Siverd
# Created:       2018-11-30
# Last modified: 2018-11-30
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.0.1"

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


##--------------------------------------------------------------------------##
## Fit for Y/X using numpy linalg.lstsq:
def fit_yxratio_numpy_lstsq(xvals, yvals, rcond=-1, full=False):
    result = np.linalg.lstsq(xvals[:, None], yvals, rcond=rcond)
    if full:
        return result       # everything we know
    else:
        return result[0]    # parameters only

##--------------------------------------------------------------------------##
## Fit for Y/X using direct summation:
def fit_yxratio_direct_sums(xvals, yvals):
    """Fit Y = a*X. Returns `a'."""
    return np.sum(xvals * yvals) / np.sum(xvals * xvals)



##--------------------------------------------------------------------------##




######################################################################
# CHANGELOG (modules/fcheck.py):
#---------------------------------------------------------------------
#
#  2018-11-30:
#     -- Increased __version__ to 0.0.1.
#     -- First created modules/fcheck.py.
#
