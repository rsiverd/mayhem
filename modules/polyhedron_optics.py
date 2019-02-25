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
# Last modified: 2019-02-25
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.2.5"

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
##------------------       Polygon Face for Polyhedra       ----------------##
##--------------------------------------------------------------------------##

class PolygonFace(object):

    """The PolygonFace class represents a regular polygon face of a larger
    polyhedron. The face is initialized using a list of vertices and internally
    tracks numerous properties to simplify raytracing calculations such as the
    intersection of lines and planes. If the polygon center is known, this
    class internally ensures that the normal vector points outward."""

    def __init__(self, face_vtx_list):
        self._verts     = np.atleast_2d(face_vtx_list).copy()
        self._basis     = self._calc_basis_vectors(self._verts)
        self._perim     = self._calc_perimeter(self._verts)
        self._center    = self._calc_center(self._verts)
        self._normal    = np.cross(*self._basis)
        self._uv_origin = self._verts[0]     # origin of UV coordinate system

        # Lookup-table for dictionary-like access:
        self._mapping  = {'vertices':self._verts, 'basis':self._basis,
                'perim':self._perim, 'center':self._center,
                'normal':self._normal, 'uv_origin':self._uv_origin}
        return

    # -----------------------------------
    # Dictionary type emulation:
    # -----------------------------------

    def __contains__(self, key):
        return (key in self._mapping.keys())

    def __getitem__(self, key):
        if not key in self._mapping.keys():
            raise KeyError
        return self._mapping[key]

    def __setitem__(self, key, value):
        if not key in self._mapping.keys():
            raise KeyError
        self._mapping[key] = value

    def keys(self):
        return self._mapping.keys()

    def values(self):
        return self._mapping.values()

    def items(self):
        return self._mapping.items()

    # -----------------------------------
    # Face parameter calculations:
    # -----------------------------------

    # Vector length calculation:
    @staticmethod
    def _vec_length(vector):
        return np.sqrt(np.sum(vector**2))

    # Calculate perimeter from array of vertices:
    @staticmethod
    def _calc_perimeter(face_vtx_list):
        diffs = np.roll(face_vtx_list, -1, axis=0) - face_vtx_list
        return np.sum(np.sqrt(np.sum(diffs**2, axis=1)))

    # Calculate midpoint of face:
    @staticmethod
    def _calc_center(face_vtx_list):
        return np.average(face_vtx_list, axis=0)

    # Calculate and store orthogonal unit vectors that span the face:
    def _calc_basis_vectors(self, face_vtx_list):
        v1, v2, v3 = face_vtx_list[:3, :]
        d21 = v2 - v1
        d31 = v3 - v1
        basis1 = d21 / self._vec_length(d21)
        basis2 = d31 - basis1 * np.dot(basis1, d21) # non-unit vector
        basis2 /= self._vec_length(basis2)          # now unit vector!
        return np.vstack((basis1, basis2))

    # -----------------------------------
    # Natural face coordinates:
    # -----------------------------------

    # Convert a single XYZ point to 'natural' UV coordinates:
    def _xyz2uv_s(self, point):
        return [np.dot(bb, point) for bb in self._basis]

    # Convert an array of XYZ points to 'natural' UV coordinates:
    def _xyz2uv_m(self, xyz_list):
        bu, bv = self._basis
        uu = np.array([np.dot(bu, pnt) for pnt in xyz_list])
        vv = np.array([np.dot(bv, pnt) for pnt in xyz_list])
        return np.array((uu, vv))

##--------------------------------------------------------------------------##
##------------------      Polyhedral Optic Base Class       ----------------##
##--------------------------------------------------------------------------##

class PolyhedralOptic(object):

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

    def get_face(self, fname):
        return self._faces.get(fname, None)

    # ------------------------------------
    # Polygon face initializer:
    def _make_face(self, face_vtx_list):
        face = PolygonFace(face_vtx_list)
        return face

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
            #for kk in ('basis',):
            #    temp = _rfunc(ang_rad, np.array(self._faces[ff][kk]).T
            for kk in ('vertices', 'basis'):
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
##------------------       Isosceles Triangular Prism       ----------------##
##--------------------------------------------------------------------------##


## Isosceles triangular prism:
class IsosPrismPolyhedron(PolyhedralOptic):

    def __init__(self, apex_angle_deg, apex_edge_mm, height_mm, unit='mm'):
        super(IsosPrismPolyhedron, self).__init__()
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

        self._faces['top'] = self._make_face(self._vtx['top'])
        self._faces['bot'] = self._make_face(self._vtx['bot'][::-1, :])
        self._faces['face1'] = self._prism_face((0, 1))
        self._faces['face2'] = self._prism_face((1, 2))
        self._faces['face3'] = self._prism_face((2, 0))
        return

    def _bottom_vertices(self):
        vtx1 = np.array([ 0.5 * self._ap_edge_mm,               0.0, 0.0])
        vtx2 = np.array([                    0.0,  self._symaxis_mm, 0.0])
        vtx3 = np.array([-0.5 * self._ap_edge_mm,               0.0, 0.0])
        return np.array([vtx1, vtx2, vtx3])

    def _prism_face(self, bvlist):
        tmpvtx = [self._vtx['bot'][x] for x in bvlist]      # bottom vertices
        tmpvtx += [self._vtx['top'][x] for x in reversed(bvlist)]   # add tops
        return self._make_face(np.array(tmpvtx))
        #return PolygonFace(tmpvtx)

##--------------------------------------------------------------------------##
##------------------      Diffraction Grating Polyhedron    ----------------##
##--------------------------------------------------------------------------##


## Grating represented by rectangular prism:
class GratingPolyhedron(PolyhedralOptic):

    def __init__(self, width_mm, length_mm, height_mm, unit='mm'):
        super(GratingPolyhedron, self).__init__()
        self._unit      = unit
        self._width_mm  = width_mm
        self._length_mm = length_mm
        self._height_mm = height_mm
        self._vtx['bot'] = self._bottom_vertices()
        self._vtx['top'] = self._vtx['bot'] \
                            + np.array([0.0, 0.0, self._height_mm])
        self._vtx['all'] = np.vstack((self._vtx['bot'], self._vtx['top']))
        self.recenter_origin()
        self._faces['top'] = self._make_face(self._vtx['top'])
        self._faces['bot'] = self._make_face(self._vtx['bot'][::-1, :])
        return

    def _bottom_vertices(self):
        vtx1 = np.array([           0.0,             0.0,  0.0])
        vtx2 = np.array([self._width_mm,             0.0,  0.0])
        vtx3 = np.array([self._width_mm, self._length_mm,  0.0])
        vtx4 = np.array([           0.0, self._length_mm,  0.0])
        return np.array([vtx1, vtx2, vtx3, vtx4])


######################################################################
# CHANGELOG (polyhedron_optics.py):
#---------------------------------------------------------------------
#
#  2019-02-24:
#     -- Increased __version__ to 0.2.0.
#     -- Renamed module to polyhedron_optics.py (also renamed classes).
#
#  2019-02-20:
#     -- Increased __version__ to 0.1.0.
#     -- First created polygon_optics.py.
#
