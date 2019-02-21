#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Describe NRES optical components, primarily the prism and grating, as a set
# of polygons. This method should make it easier to (a) produce an internally
# consistent spectrograph layout model and (b) perform raytracing operations on
# those components once placed.
#
# A basic OpticsPolygon object provides methods such as:
# * shift_xyz()
# * xrotate()
# * yrotate()
# * zrotate()
#
# These methods will apply the selected transformation to all of the internal
# pieces, including vertex points, center-of-face points, and surface normals.
#
#
# Rob Siverd
# Created:       2019-02-20
# Last modified: 2019-02-20
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

## Rotations in 3D:
import fov_rotation
reload(fov_rotation)
r3d = fov_rotation.Rotate3D()

##--------------------------------------------------------------------------##

_rot_map = {
        'x' :   r3d.xrotate_xyz,
        'y' :   r3d.yrotate_xyz,
        'z' :   r3d.zrotate_xyz,
        }


##--------------------------------------------------------------------------##
##------------------         Basic Optical Polygon          ----------------##
##--------------------------------------------------------------------------##

class OpticalPolygon(object):

    def __init__(self):
        self._vtx   = {}
        self._faces = {}
        return

    # ------------------------------------
    # Property getters:
    def get_center(self):
        return np.average(self._vtx['all'], axis=0)

    def get_vertices(self, which='all'):
        return self._vtx[which]

    # ------------------------------------
    # Polygon face initializer:
    def _calc_normal(self, face_vtx_list, reverse=False):
        v1, v2, v3 = face_vtx_list[:3, :]
        d21 = v2 - v1
        d31 = v3 - v1
        normvec = np.cross(d21, d31)
        return normvec / self._vec_length(normvec)

    def _face_center(self, face_vtx_list):
        midpoint = np.average(face_vtx_list, axis=0)
        return midpoint

    def _make_face(self, face_vtx_list):
        face = {}
        face['vertices'] = np.copy(face_vtx_list)
        face['center'] = self._face_center(face_vtx_list)
        face['normal'] = self._calc_normal(face_vtx_list)
        return face

    @staticmethod
    def _vec_length(vector):
        return np.sqrt(np.sum(vector**2))

    # ------------------------------------
    # Polygon movements:
    def recenter_origin(self):
        centroid = self.get_center()
        return self.shift_vec(-centroid)

    def shift_xyz(self, dx, dy, dz):
        return self.shift_vec(np.array([dx, dy, dz]))

    def shift_vec(self, displacement):
        # All vertices are displaced:
        for kk,vv in self._vtx.items():
            self._vtx[kk] = vv + displacement
 
        # For faces, shift center and vertices (normal unchanged):
        for ff in self._faces.keys():
            for kk in ('center', 'vertices',):
                self._faces[ff][kk] += displacement
        return

    def _rotate(self, ang_rad, axname):
        # NOTE: 2-D point sets need transposition but vectors do not
        _rfunc = _rot_map[axname]

        # Rotate vertices (point sets):
        for kk,vv in self._vtx.items():
            temp = _rfunc(ang_rad, vv.T)
            self._vtx[kk] = np.array(temp.T)

        # Rotate faces:
        for ff in self._faces.keys():
            for kk in ('center', 'normal'):
                self._faces[ff][kk] = _rfunc(ang_rad, self._faces[ff][kk]).A1
            for kk in ('vertices',):
                temp = _rfunc(ang_rad, self._faces[ff][kk].T)
                self._faces[ff][kk] = np.array(temp.T)
        return

    def xrotate(self, ang_rad):
        return self._rotate(ang_rad, 'x')

    def yrotate(self, ang_rad):
        return self._rotate(ang_rad, 'y')

    def zrotate(self, ang_rad):
        return self._rotate(ang_rad, 'z')

##--------------------------------------------------------------------------##
##------------------       Specific Optical Components      ----------------##
##--------------------------------------------------------------------------##

class PolygonPrism(OpticalPolygon):

    def __init__(self, apex_angle_deg, apex_edge_mm, height_mm, unit='mm'):
        super(PolygonPrism, self).__init__()
        self._unit       = unit
        #self._apex_deg   = apex_angle_deg
        self._apex_rad   = np.radians(apex_angle_deg)
        self._height_mm  = height_mm
        self._ap_edge_mm = apex_edge_mm
        #self._o_edges_mm = 0.5 * self._ap_edge_mm / np.sin(self._apex_rad / 2.)
        self._symaxis_mm = 0.5 * self._ap_edge_mm / np.tan(self._apex_rad / 2.)
        self._vtx = {}
        self._vtx['bot'] = self._bottom_vertices()
        self._vtx['top'] = self._vtx['bot'] \
                            + np.array([0.0, 0.0, self._height_mm])
        self._vtx['all'] = np.vstack((self._vtx['bot'], self._vtx['top']))
        #self.barycenter  = self.get_center()
        self.recenter_origin()
        self._faces = {}
        self._faces['top'] = self._make_face(self._vtx['top'])
        self._faces['bot'] = self._make_face(self._vtx['bot'][::-1, :])
        #pr_norm_top = self._calc_normal(self.get_vertices('top'))
        #pr_norm_bot = self._calc_normal(self.get_vertices('bot')[::-1, :])
        #sys.stderr.write("pr_norm_top: %s\n" % str(pr_norm_top))
        #sys.stderr.write("pr_norm_bot: %s\n" % str(pr_norm_bot))
        return

    def _bottom_vertices(self):
        vtx1 = np.array([ 0.5 * self._ap_edge_mm,               0.0, 0.0])
        vtx2 = np.array([                    0.0,  self._symaxis_mm, 0.0])
        vtx3 = np.array([-0.5 * self._ap_edge_mm,               0.0, 0.0])
        return np.array([vtx1, vtx2, vtx3])


class PolygonGrating(OpticalPolygon):

    def __init__(self, width_mm, length_mm, height_mm, unit='mm'):
        super(PolygonGrating, self).__init__()
        self._unit      = unit
        self._width_mm  = width_mm
        self._length_mm = length_mm
        self._height_mm = height_mm
        #g_normal = np.array([    0.0,   0.0, -1.0])  # grating normal
        self._vtx = {}
        self._vtx['bot'] = self._bottom_vertices()
        self._vtx['top'] = self._vtx['bot'] \
                            + np.array([0.0, 0.0, self._height_mm])
        self._vtx['all'] = np.vstack((self._vtx['bot'], self._vtx['top']))
        self.recenter_origin()
        #sys.stderr.write("Initial 'bot' vtx:\n%s\n" % str(self._vtx['bot']))
        self._faces = {}
        self._faces['top'] = self._make_face(self._vtx['top'])
        self._faces['bot'] = self._make_face(self._vtx['bot'][::-1, :])
        #gr_norm_top = self._calc_normal(self.get_vertices('top'))
        #gr_norm_bot = self._calc_normal(self.get_vertices('bot')[::-1, :])
        #sys.stderr.write("gr_norm_top: %s\n" % str(gr_norm_top))
        #sys.stderr.write("gr_norm_bot: %s\n" % str(gr_norm_bot))
        return

    def _bottom_vertices(self):
        vtx1 = np.array([           0.0,             0.0,  0.0])
        vtx2 = np.array([self._width_mm,             0.0,  0.0])
        vtx3 = np.array([self._width_mm, self._length_mm,  0.0])
        vtx4 = np.array([           0.0, self._length_mm,  0.0])
        return np.array([vtx1, vtx2, vtx3, vtx4])


######################################################################
# CHANGELOG (polygon_optics.py):
#---------------------------------------------------------------------
#
#  2019-02-20:
#     -- Increased __version__ to 0.1.0.
#     -- First created polygon_optics.py.
#
