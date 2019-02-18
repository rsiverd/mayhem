#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Some useful routines to aid in vector raytracing. These are helpful when
# figuring out 3-D positioning of light rays through a spectrograph.
#
# Background information on the mathematics can be found at:
# https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
# https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
#
# Rob Siverd
# Created:       2019-02-18
# Last modified: 2019-02-18
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
#from functools import partial
#from collections import OrderedDict
#import multiprocessing as mp
#np.set_printoptions(suppress=True, linewidth=160)
#import pandas as pd
#import statsmodels.api as sm
#import statsmodels.formula.api as smf
#from statsmodels.regression.quantile_regression import QuantReg
#import PIL.Image as pli
#import seaborn as sns
#import cmocean
#import theil_sen as ts
#import window_filter as wf
#import itertools as itt

##--------------------------------------------------------------------------##
## Reflection and refraction in 3 dimensions (vector forms):
def calc_surface_vectors(v_incident, surf_norm, n1_n2_ratio):
    cti = -1.0 * np.dot(v_incident, surf_norm)     # cos(theta_i)
    nnr = n1_n2_ratio
    v_reflect = v_incident + 2. * cti * surf_norm
    smult = nnr*cti - np.sqrt(1.0 - nnr**2 * (1.0 - cti**2))
    v_refract = nnr * v_incident + smult * surf_norm
    return v_reflect, v_refract

##--------------------------------------------------------------------------##







######################################################################
# CHANGELOG (raytrace_tools.py):
#---------------------------------------------------------------------
#
#  2019-02-18:
#     -- Increased __version__ to 0.1.0.
#     -- First created raytrace_tools.py.
#
