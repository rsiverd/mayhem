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
# Last modified: 2019-03-29
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.3.0"

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

def refracted_ray(v_incident, surf_norm, n1_n2_ratio):
    cti = -1.0 * np.dot(v_incident, surf_norm)     # cos(theta_i)
    nnr = n1_n2_ratio
    smult = nnr*cti - np.sqrt(1.0 - nnr**2 * (1.0 - cti**2))
    v_refract = nnr * v_incident + smult * surf_norm
    return v_refract

##--------------------------------------------------------------------------##
## Intersection of a line and plane. In this model, both lines and planes are
## described by (point, vector pairs). In the case of a line, the vector points
## "along" the line. For the plane, the vector is surface normal. 
#def line_plane_intersection(line, plane):
    #lpoint, lvector =  line['point'],  line['vector']
    #ppoint, pvector = plane['point'], plane['vector']
#def line_plane_intersection(lpoint, lvector, ppoint, pnormal):
#    """
#    All inputs should be 3-element numpy ndarray type. The case of parallel
#    line and plane should be handled ~correctly.
#
#    lpoint  -  any point (x, y, z) on the line
#    lvector -  any vector (dx, dy, dz) along the line
#    ppoint  -  any point (x, y, z) on the plane
#    pnormal -  any vector (dx, dy, dz) normal to plane surface
#
#    Returns:
#    intersection  -  (x, y, z) point of intersection
#    """
#
#    pl_sep = np.dot(ppoint - lpoint, pnormal)
#    angsep = np.dot(lvector, pnormal)
#    # Parallel line/plane is handled separately:
#    if (angsep == 0.0):
#        sys.stderr.write("WARNING: line and plane are PARALLEL!\n")
#        if (pl_sep == 0.0):
#            sys.stderr.write("Line lies within plane!\n")
#            return lpoint
#            #return (0.0, lpoint)
#        else:
#            sys.stderr.write("Line and plane do not intersect!\n")
#            return None
#            #return None, None
#    # If not parallel, get distance and intersection: 
#    distance = pl_sep / angsep
#    return lpoint + distance * lvector
#    #isect = lpoint + distance * lvector
#    #truedist = np.sqrt(np.sum((isect - lpoint)**2))
#    #sys.stderr.write("truedist: %10.5f\n" % truedist)
#    #return distance, isect

##--------------------------------------------------------------------------##
## Point-in-polygon test routine. This will be useful to check whether a
## calculated line-plane intersection point resides within the boundaries of
## a face as defined by its vertices.
#def point_in_polygon(point, vtx_list):
#    return

##--------------------------------------------------------------------------##
##------------------         Grating Diffraction            ----------------##
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
##------------------     End-to-end Spectrograph Raytrace   ----------------##
##--------------------------------------------------------------------------##

class E2ERT(object):

    def __init__(self, prism, grating, ccd):
        #self._faces = [x for x in face_list]
        #self._faces = copy.deepcopy(face_list)
        self._probj = prism
        self._grobj = grating
        self._cmobj = ccd
        self._prf1  = self._probj.get_face('face1')
        self._prf2  = self._probj.get_face('face2')
        return

    def follow(self, xyz0, traj0, wl_um, spec_order):
        #sys.stderr.write("--------------------------------------------\n")
        prf1   = self._probj.get_face('face1')
        prf2   = self._probj.get_face('face2')
        p_n1n2 = self._probj._n1n2_ratio(wl_um) # n_glass / n_air
        grbot  = self._grobj.get_face('bot')
        sensor = self._cmobj.get_face('front')
        light_path = [(xyz0, traj0)]
        #sys.stderr.write("v_initial (NEW): %s\n" % str(traj0))
        #sys.stderr.write("prf1_norm (NEW): %s\n" % str(prf1['normal']))
        #sys.stderr.write("n1n2ratio (NEW): %s\n" % str(p_n1n2))

        # -------------------------------------------------------
        # Enter prism through face1 (point + trajectory):
        valid, f1_isect = prf1.get_intersection(xyz0, traj0)
        if not valid:
            sys.stderr.write("Input beam missed prism!\n")
            return light_path
        new_traj = calc_surface_vectors(traj0, prf1['normal'], p_n1n2)[1]
        #sys.stderr.write("new_traj after face1: %s\n" % str(new_traj))
        light_path.append((f1_isect, new_traj))

        # Exit through prism face2 (point + trajectory):
        valid, f2_isect = prf2.get_intersection(f1_isect, new_traj)
        if not valid:
            sys.stderr.write("Beam failed to exit prism through face2!\n")
            return light_path
        new_traj = calc_surface_vectors(new_traj, 
                        -prf2['normal'], 1. / p_n1n2)[1]

        light_path.append((f2_isect, new_traj))
        #sys.stderr.write("--------------------------------------------\n")

        # -------------------------------------------------------
        # Intersect grating and change direction (diffract):
        hits_grating, gr_isect = grbot.get_intersection(f2_isect, new_traj)
        if not valid:
            sys.stderr.write("Light ray misses grating!\n")
            return light_path
        valid, diffr_vec = \
                self._grobj.diffracted_ray(new_traj, wl_um, spec_order)
        if not valid:
            light_path.append((gr_isect, None))
            sys.stderr.write("No valid diffracted beam from grating!\n")
            return light_path
        light_path.append((gr_isect, diffr_vec))
 
        # -------------------------------------------------------
        # Re-enter prism through face2 (point + trajectory):
        hits_prism, f2_isect = prf2.get_intersection(gr_isect, diffr_vec)
        if not hits_prism:
            sys.stderr.write("Diffracted ray does not return to prism!\n")
            return light_path
        new_traj = calc_surface_vectors(diffr_vec, prf2['normal'], p_n1n2)[1]
        light_path.append((f2_isect, new_traj))

        # Exit through prism face1 (point + trajectory):
        valid, f1_isect = prf1.get_intersection(f2_isect, new_traj)
        if not valid:
            sys.stderr.write("Beam failed to exit prism through face1!\n")
            return light_path
        new_traj = calc_surface_vectors(new_traj, 
                        -prf1['normal'], 1. / p_n1n2)[1]
        light_path.append((f1_isect, new_traj))

        # -------------------------------------------------------
        # Check for intersection with CCD plane:
        valid, cam_isect = sensor.get_intersection(f1_isect, new_traj)
        if not valid:
            sys.stderr.write("Beam misses CCD sensor!\n")
            return light_path
        light_path.append((cam_isect, None))

        return light_path


######################################################################
# CHANGELOG (raytrace_tools.py):
#---------------------------------------------------------------------
#
#  2019-02-18:
#     -- Increased __version__ to 0.1.0.
#     -- First created raytrace_tools.py.
#
