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
# Last modified: 2019-02-24
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
## Intersection of a line and plane. In this model, both lines and planes are
## described by (point, vector pairs). In the case of a line, the vector points
## "along" the line. For the plane, the vector is surface normal. 
#def line_plane_intersection(line, plane):
    #lpoint, lvector =  line['point'],  line['vector']
    #ppoint, pvector = plane['point'], plane['vector']

def line_plane_intersection(lpoint, lvector, ppoint, pnormal):
    """
    All inputs should be 3-element numpy ndarray type. The case of parallel
    line and plane should be handled ~correctly.

    lpoint  -  any point (x, y, z) on the line
    lvector -  any vector (dx, dy, dz) along the line
    ppoint  -  any point (x, y, z) on the plane
    pnormal -  any vector (dx, dy, dz) normal to plane surface

    Returns:
    intersection  -  (x, y, z) point of intersection
    """

    pl_sep = np.dot(ppoint - lpoint, pnormal)
    angsep = np.dot(lvector, pnormal)
    # Parallel line/plane is handled separately:
    if (angsep == 0.0):
        sys.stderr.write("WARNING: line and plane are PARALLEL!\n")
        if (pl_sep == 0.0):
            sys.stderr.write("Line lies within plane!\n")
            return lpoint
            #return (0.0, lpoint)
        else:
            sys.stderr.write("Line and plane do not intersect!\n")
            return None
            #return None, None
    # If not parallel, get distance and intersection: 
    distance = pl_sep / angsep
    return lpoint + distance * lvector
    #isect = lpoint + distance * lvector
    #truedist = np.sqrt(np.sum((isect - lpoint)**2))
    #sys.stderr.write("truedist: %10.5f\n" % truedist)
    #return distance, isect

##--------------------------------------------------------------------------##
## Point-in-polygon test routine. This will be useful to check whether a
## calculated line-plane intersection point resides within the boundaries of
## a face as defined by its vertices.
def point_in_polygon(point, vtx_list):
    return


######################################################################
# CHANGELOG (raytrace_tools.py):
#---------------------------------------------------------------------
#
#  2019-02-18:
#     -- Increased __version__ to 0.1.0.
#     -- First created raytrace_tools.py.
#
