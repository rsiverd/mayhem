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
# Last modified: 2019-06-09
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.3.2"

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
#from functools import partial
#from collections import OrderedDict
#import multiprocessing as mp
#np.set_printoptions(suppress=True, linewidth=160)
#import itertools as itt

## Indexes of refraction for glass types:
from glass import Glass

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

    def __init__(self, face_vtx_list, name='unnamed', debug=False):
        self._name      = name
        self._debug     = debug
        self._verts     = np.atleast_2d(face_vtx_list).copy()
        self._basis     = self._calc_basis_vectors(self._verts)
        self._perim     = self._calc_perimeter(self._verts)
        self._center    = self._calc_center(self._verts)
        self._normal    = np.cross(*self._basis)
        self._uv_origin = np.copy(self._verts[0])       # UV coordinate origin
        self._uv_verts  = self._xyz2uv_m(self._verts)   # UV vertex coo (fixed)
        self._vtxcount  = self._verts.shape[0]          # number of vertices

        # Lookup-table for dictionary-like access:
        self._mapping  = {'vertices':'_verts', 'basis':'_basis',
                'perim':'_perim', 'center':'_center',
                'normal':'_normal', 'uv_origin':'_uv_origin'}

        # Displacement configuration:
        self._shift_these = ['vertices', 'center', 'uv_origin']

        # Rotation configuration:
        self._rotate_vecs = ['center', 'normal', 'uv_origin']
        self._rotate_arrs = ['vertices', 'basis']

        if self._debug:
            self._sanity_check()

        # Point-in-polygon configuration:
        if self._vtxcount == 4:
            self._contains_point = self._rect_contains
        else:
            self._contains_point = self._wind_contains
        return

    # -----------------------------------
    # Internal consistency checker:
    # -----------------------------------

    def _sanity_check(self):
        caller_name = sys._getframe(1).f_code.co_name
        this_func = sys._getframe().f_code.co_name
        b1, b2 = self._basis
        b1size = self._vec_length(b1)
        b2size = self._vec_length(b2)
        fn     = self._normal
        fnsize = self._vec_length(fn)
        sys.stderr.write("Face %s: running %s called by %s ... "
                % (self._name, this_func, caller_name))
        for nn,ll,vv in zip(('b1', 'b2', 'fn'),
                            (b1size, b2size, fnsize), 
                            (b1, b2, fn)):
            if (np.abs(ll - 1.0) > 1e-5):
                sys.stderr.write("FAILURE!!\n\n")
                sys.stderr.write("Bogus %s length: %10.5f\n" % (nn, ll))
                sys.stderr.write("%s vector: %s\n\n" % (nn, str(vv)))
                raise ValueError
        sys.stderr.write("passed!\n")
        return

    # -----------------------------------
    # Dictionary type emulation:
    # -----------------------------------

    def __contains__(self, key):
        return (key in self._mapping.keys())

    def __getitem__(self, key):
        if not key in self._mapping.keys():
            raise KeyError
        return getattr(self, self._mapping[key])

    def __setitem__(self, key, value):
        if not key in self._mapping.keys():
            raise KeyError
        setattr(self, self._mapping[key], value)
        return

    def keys(self):
        return self._mapping.keys()

    def values(self):
        return [getattr(self, x) for x in self._mapping.values()]

    def items(self):
        return [(kk, getattr(self, vv)) for kk,vv in self._mapping.items()]

    # -----------------------------------
    # Post-initialization configuration:
    # -----------------------------------

    def set_name(self, face_name):
        self._name = face_name
        return

    # Ensure outward-pointing normal using polyhedron center:
    def _ensure_outward(self, ph_center):
        #bc_dist = self.get_intersect_distance(ph_center, self._normal)
        bc_dist = self._calc_isect_distance(ph_center, self._normal)
        #sys.stderr.write("bc_dist: %10.5f, normal: %s\n" 
        #        % (bc_dist, str(self._normal)))
        if (bc_dist < 0.0):
            sys.stderr.write("Inward-facing normal detected! FAIL!\n")
            sys.exit(1)
            #sys.stderr.write("Inversion detected and corrected!\n")
            #self._basis[1, :] *= -1.0   # reverse 2nd basis vector
            #self._normal      *= -1.0   # reverse normal vector
        return

    # -----------------------------------
    # Polygon shifts and rotations:
    # -----------------------------------

    def _apply_shift(self, displacement):
        for kk in self._shift_these:
            thing = getattr(self, self._mapping[kk])
            thing += displacement
        if self._debug:
            self._sanity_check()
        return

    def _apply_rotation(self, ang_rad, axname):
        # NOTE: 2-D point sets need transposition but vectors do not
        _rfunc = _rot_map[axname]

        # Rotate vertices en masse:
        for kk in self._rotate_arrs:
            which = self._mapping[kk]
            temp = _rfunc(ang_rad, getattr(self, which).T)
            setattr(self, which, np.array(temp.T))

        # Rotate individual vectors separately:
        for kk in self._rotate_vecs:
            which = self._mapping[kk]
            setattr(self, which, _rfunc(ang_rad, getattr(self, which)).A1)

        if self._debug:
            self._sanity_check()
        return

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
        v1, v2, v3 = np.copy(face_vtx_list[:3, :])
        d21 = v2 - v1
        d31 = v3 - v1
        basis1 = d21 / self._vec_length(d21)
        basis2 = d31 - basis1 * np.dot(basis1, d31) # non-unit vector
        basis2 /= self._vec_length(basis2)          # now unit vector!
        return np.vstack((basis1, basis2))

    # -----------------------------------
    # Natural face coordinates:
    # -----------------------------------

    # Convert a single XYZ point to 'natural' UV coordinates:
    def _xyz2uv_s(self, point):
        return np.array([np.dot(bb, point - self._uv_origin) \
                                            for bb in self._basis])

    # Convert an array of XYZ points to 'natural' UV coordinates:
    def _xyz2uv_m(self, xyz_list):
        bu, bv = self._basis
        uu = np.array([np.dot(bu, pnt - self._uv_origin) for pnt in xyz_list])
        vv = np.array([np.dot(bv, pnt - self._uv_origin) for pnt in xyz_list])
        return np.array((uu, vv))

    # -----------------------------------
    # Point-in-polygon face bounds tests:
    # -----------------------------------

    # Efficient test for rectangular faces:
    def _rect_contains(self, point, rtol=1e-6, emergency=False):
        """
        This method computes the distance in the UV plane between a test 
        point (assumed to be on that plane!) and the four corners of the face 
        in UV coordinates (i.e., using basis vectors parallel to the sides of
        the rectangle). The absolute sum of the U- and V-coordinate differences
        (i.e., abs_sum = Σ (|ΔU_i| + |ΔV_i|) for i = 1...4) will equal the 
        rectangle perimeter when the test point lies inside the polygon. If the
        point is outside, this sum exceeds the perimeter. The test point is
        considered inside when ((abs_sum - perim) / perim) < rtol.

        As noted above, the test point is assumed to be in the UV plane (which
        contains the rectangle). More precisely, the xyz2uv conversion gives 
        the projected UV position (any distance from the UV plane is dropped).
        As a result, this routine is not by itself appropriate for finding 
        whether a point lies *ON* a face or not. Proper use entails first 
        computing the intersection of a ray with the UV plane. The intersection
        point is then tested for 'insideness' using this method.

        NOTE TO FUTURE SELF: When writing this description I wondered for a bit
        why I don't compute full uvw coordinates from xyz instead. A routine
        to find the W distance, get_intersect_distance(), exists already. The
        answer is that this routine is intended for in-plane use only, hence
        the name _rect_contains (and NOT _face_contains).

        """
        point_uv = self._xyz2uv_s(point)
        diffs_uv = self._uv_verts - point_uv[:, None]
        totalsep = np.sum(np.abs(diffs_uv))
        fracdiff = np.abs(totalsep - self._perim) / self._perim
        if emergency:
            sys.stderr.write("point_uv:  %s\n" % str(point_uv))
            sys.stderr.write("diffs_uv:  %s\n" % str(diffs_uv))
            sys.stderr.write("totalsep:  %s\n" % str(totalsep))
            sys.stderr.write("perimeter: %s\n" % str(self._perim))
            sys.stderr.write("fracdiff:  %s\n" % str(fracdiff))
        return (fracdiff < rtol)

    # Orientation difference calculator for vectors. This finds
    # the angle subtended by the line segment connecting the "ends"
    # of two input vectors. Use for winding number calculation.
    @staticmethod
    def _calc_rel_angle(delta_1, delta_2):
        twopi = 2.0 * np.pi
        dx1, dy1 = delta_1
        dx2, dy2 = delta_2
        theta1 = np.arctan2(dy1, dx1)
        theta2 = np.arctan2(dy2, dx2)
        dtheta = (theta2 - theta1) % twopi
        dtheta[(dtheta > np.pi)] -= twopi
        return dtheta

    # Calculate winding number of vertices relative to test point:
    def _count_windings(self, point):
        twopi   = 2.0 * np.pi
        uvpoint = self._xyz2uv_s(point)
        diffs_1 = self._uv_verts - uvpoint[:, None]
        diffs_2 = np.roll(diffs_1, 1, axis=1)
        n_winds = np.sum(self._calc_rel_angle(diffs_1, diffs_2)) / twopi
        return np.round(np.abs(n_winds), 5)

    # General winding-number 'insideness' check:
    def _wind_contains(self, point):
        return (self._count_windings(point) == 1.0)

    # -----------------------------------
    # Intersections of lines and planes:
    # -----------------------------------

    def _line_intersection(self, lpoint, lvector):
        ppoint, pnormal = self._center, self._normal
        pl_sep = np.dot(ppoint - lpoint, pnormal)   # plane-line separation
        angsep = np.dot(lvector, pnormal)           # angular separation
        if (angsep == 0.0):
            sys.stderr.write("WARNING: line and plane are PARALLEL!\n")
            if (pl_sep == 0.0):
                sys.stderr.write("Line lies within plane!\n")
                return lpoint
            else:
                sys.stderr.write("Line and plane do not intersect!\n")
                return np.array([np.nan, np.nan, np.nan])

        # If not parallel, get distance and intersection:
        distance = pl_sep / angsep
        return lpoint + distance * lvector

    # Distance to intersection point (no bounds checking):
    def _calc_isect_distance(self, lpoint, lvector):
        isect = self._line_intersection(lpoint, lvector)
        return np.dot(isect - lpoint, lvector) / self._vec_length(lvector)

    # Line intersection WITH face bounds check:
    def get_intersection(self, lpoint, lvector):
        """Check whether specified line intersects this face.
        Returns:
            valid  -- True if intersection occurs, otherwise False
            isect  -- point of intersection if possible, otherwise None
        """
        isect = self._line_intersection(lpoint, lvector)
        if self._contains_point(isect):
            return True, isect
        else:
            return False, None

    # User-friendly distance-to-intersection calculator (handles bounds):
    def get_intersect_distance(self, lpoint, lvector):
        valid, isect = self.get_intersection(lpoint, lvector)
        if not valid:
            return np.nan
        else:
            return np.dot(isect - lpoint, lvector) / self._vec_length(lvector)


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

    def list_faces(self):
        """A list of faces for use with get_face()."""
        return list(self._faces.keys())

    def _prism_face(self, bvlist):
        tmpvtx = [self._vtx['bot'][x] for x in bvlist]      # bottom vertices
        tmpvtx += [self._vtx['top'][x] for x in reversed(bvlist)]   # add tops
        return self._make_face(np.array(tmpvtx))

    # ------------------------------------
    # Polygon face initializer:
    def _make_face(self, face_vtx_list):
        face = PolygonFace(face_vtx_list)
        return face

    # Post-initialization face naming:
    def _update_face_names(self):
        for kk,face in self._faces.items():
            face.set_name(kk)
        return

    # Post-initialization outward normals:
    def _update_face_normals(self):
        barycenter = self.get_center()
        for face in self._faces.values():
            face._ensure_outward(barycenter)
        return

    # Set names and update normals:
    def _update_face_names_norms(self):
        barycenter = self.get_center()
        for kk,face in self._faces.items():
            #sys.stderr.write("Updating face %s ... \n" % kk)
            face.set_name(kk)
            face._ensure_outward(barycenter)
            face._debug = self._debug
        return

    # ------------------------------------
    # Polygon movements:
    def recenter_origin(self):
        centroid = self.get_center()
        return self.shift_vec(-centroid)

    def shift_xyz(self, dx, dy, dz):
        return self.shift_vec(np.array([dx, dy, dz]))

    def shift_vec(self, displacement):
        # All vertices are displaced:
        have_top = 'top' in self._vtx.keys()
        for kk,vv in self._vtx.items():
            #self._vtx[kk] = vv + displacement
            self._vtx[kk] += displacement
 
        # Faces handle displacement internally:
        for face in self._faces.values():
            face._apply_shift(displacement)

        # Sanity check:
        if self._debug:
            self._brute_consistency_check()
        return

    def _rotate(self, ang_rad, axname):
        # NOTE: 2-D point sets need transposition but vectors do not
        _rfunc = _rot_map[axname]

        # Rotate vertices (point sets):
        for kk,vv in self._vtx.items():
            temp = _rfunc(ang_rad, vv.T)
            self._vtx[kk] = np.array(temp.T)

        # Rotate faces:
        for face in self._faces.values():
            face._apply_rotation(ang_rad, axname)

        # Sanity check:
        if self._debug:
            self._brute_consistency_check()
        return

    def xrotate(self, ang_rad):
        return self._rotate(ang_rad, 'x')

    def yrotate(self, ang_rad):
        return self._rotate(ang_rad, 'y')

    def zrotate(self, ang_rad):
        return self._rotate(ang_rad, 'z')

    # A brute-force run-time sanity check on polyhedron-face consistency:
    def _brute_consistency_check(self):
        caller_name = sys._getframe(1).f_code.co_name
        this_func = sys._getframe().f_code.co_name
        sys.stderr.write("Executing consistency check (from %s) ... "
                % caller_name)
        common_keys = [x for x in self._vtx.keys() if x in self._faces.keys()]

        # ignore bottom for now (order switched):
        #if 'bot' in common_keys:
        #    common_keys.remove('bot') # FIXME: need consistent vertex ordering
        #sys.stderr.write("common_keys: %s\n" % str(common_keys))
        #sys.stderr.write("vtx keys: %s\n" % str(self._vtx.keys()))
        #sys.stderr.write("facekeys: %s\n" % str(self._faces.keys()))
        for kk in common_keys:
            diffs = self._faces[kk]['vertices'] - self._vtx[kk]
            success = np.all(diffs == 0.0)
            if not success:
                sys.stderr.write("FAILURE!!!\n\n")
                sys.stderr.write("Inconsistency in %s vertices!!!!\n" % kk)
                sys.stderr.write("%s diffs: %s\n" % (kk, str(diffs)))
                sys.exit(1)
        sys.stderr.write("passed!\n")
        return

##--------------------------------------------------------------------------##
##------------------       Isosceles Triangular Prism       ----------------##
##--------------------------------------------------------------------------##


## Isosceles triangular prism:
class IsosPrismPolyhedron(PolyhedralOptic):

    def __init__(self, apex_angle_deg, apex_edge_mm, height_mm,
            glass_type, n_air=1.0, unit='mm', debug=False):
        super(IsosPrismPolyhedron, self).__init__()
        self._unit       = unit
        self._debug      = debug
        self._apex_rad   = np.radians(apex_angle_deg)
        self._height_mm  = height_mm
        self._ap_edge_mm = apex_edge_mm
        self._symaxis_mm = 0.5 * self._ap_edge_mm / np.tan(self._apex_rad / 2.)
        self._vtx = {}
        self._vtx['bot'] = self._bottom_vertices()
        self._vtx['top'] = self._vtx['bot'] \
                            + np.array([0.0, 0.0, self._height_mm])
        self._vtx['all'] = np.vstack((self._vtx['bot'], self._vtx['top']))
        #self.barycenter  = self.get_center()
        self.recenter_origin()

        self._faces['top'] = self._make_face(self._vtx['top'])
        # FIXME:
        self._faces['bot'] = self._make_face(self._vtx['bot'][::-1, :])
        #self._faces['bot'] = self._make_face(self._vtx['bot'])
        self._faces['face1'] = self._prism_face((0, 1))
        self._faces['face2'] = self._prism_face((1, 2))
        self._faces['face3'] = self._prism_face((2, 0))
        #j = calc_surface_vectors(traj0, prf1['normal'], p_n1n2)[1]
        self._update_face_names()
        #self._update_face_normals()
        self._update_face_names_norms()
        
        # Indexes of refraction:
        self._n_air      = n_air
        self._glass      = Glass(glass_type)
        return

    def _bottom_vertices(self):
        vtx1 = np.array([ 0.5 * self._ap_edge_mm,               0.0, 0.0])
        vtx2 = np.array([                    0.0,  self._symaxis_mm, 0.0])
        vtx3 = np.array([-0.5 * self._ap_edge_mm,               0.0, 0.0])
        return np.array([vtx1, vtx2, vtx3])

    #def _prism_face(self, bvlist):
    #    tmpvtx = [self._vtx['bot'][x] for x in bvlist]      # bottom vertices
    #    tmpvtx += [self._vtx['top'][x] for x in reversed(bvlist)]   # add tops
    #    return self._make_face(np.array(tmpvtx))

    def _n1n2_ratio(self, wl_um):
        return self._n_air / self._glass.refraction_index(wl_um)
        #return self._glass.refraction_index(wl_um) / self._n_air

##--------------------------------------------------------------------------##
##------------------          CCD Camera Polyhedron         ----------------##
##--------------------------------------------------------------------------##

class CameraPolyhedron(PolyhedralOptic):
    
    def __init__(self, width_mm, length_mm, height_mm,
            unit='mm', debug=False):
        super(CameraPolyhedron, self).__init__()
        self._ccd_roll  = 0.0
        self._unit      = unit
        self._debug     = debug
        self._width_mm  = width_mm
        self._length_mm = length_mm
        self._height_mm = height_mm
        self._vtx['bot'] = self._bottom_vertices()
        self._vtx['top'] = self._vtx['bot'] \
                            + np.array([0.0, 0.0, self._height_mm])
        self._vtx['all'] = np.vstack((self._vtx['bot'], self._vtx['top']))
        self.recenter_origin()
        self._faces['front'] = self._prism_face((2, 3))
        #self._faces['face2'] = self._prism_face((1, 2))
        #self._faces['face3'] = self._prism_face((2, 0))
        #self.set_sensor_roll_deg(roll_deg)
        return

    def _bottom_vertices(self):
        vtx1 = np.array([           0.0,             0.0,  0.0])
        vtx2 = np.array([self._width_mm,             0.0,  0.0])
        vtx3 = np.array([self._width_mm, self._length_mm,  0.0])
        vtx4 = np.array([           0.0, self._length_mm,  0.0])
        return np.array([vtx1, vtx2, vtx3, vtx4])

    #def set_sensor_roll_deg(self, roll_degrees):
    #    """This rotation routine is used to set and track a 'roll' angle of
    #    the sensor about the optical axis (the 'front' face normal vector).
    #    Importantly, the roll position is tracked internally to simplify
    #    transformation from spectrograph to CCD coordinates.
    #    
    #    Effectively, this is a sensor adjustment INSIDE the polygon face
    #    This roll serves two purposes:
    #    1) Within the polygon, it orients the sensor in the spectrograph 
    #    coordinate system so that
    #    intersections can be identified in raytracing.
    #    2) that 'rolls' the CCD sensor about the (potentially arbitrary)
    #    optical axis (the 'front' face normal).  and tracks the
    #    rotation. rather than
    #    """

##--------------------------------------------------------------------------##
##------------------      Diffraction Grating Polyhedron    ----------------##
##--------------------------------------------------------------------------##


## Grating represented by rectangular prism:
class GratingPolyhedron(PolyhedralOptic):

    def __init__(self, width_mm, length_mm, height_mm, ruling_lmm, 
            unit='mm', debug=False):
        super(GratingPolyhedron, self).__init__()
        self._unit      = unit
        self._debug     = debug
        self._width_mm  = width_mm
        self._length_mm = length_mm
        self._height_mm = height_mm
        self._ruling_lmm = ruling_lmm
        self._spacing_um = 1e3 / self._ruling_lmm
        self._vtx['bot'] = self._bottom_vertices()
        self._vtx['top'] = self._vtx['bot'] \
                            + np.array([0.0, 0.0, self._height_mm])
        self._vtx['all'] = np.vstack((self._vtx['bot'], self._vtx['top']))
        self.recenter_origin()
        self._faces['top'] = self._make_face(self._vtx['top'])
        self._faces['bot'] = self._make_face(self._vtx['bot'][::-1, :])
        #self._faces['bot'] = self._make_face(self._vtx['bot'])
        #self._update_face_names()
        #self._update_face_normals()
        self._update_face_names_norms()
        return

    def _bottom_vertices(self):
        vtx1 = np.array([           0.0,             0.0,  0.0])
        vtx2 = np.array([self._width_mm,             0.0,  0.0])
        vtx3 = np.array([self._width_mm, self._length_mm,  0.0])
        vtx4 = np.array([           0.0, self._length_mm,  0.0])
        return np.array([vtx1, vtx2, vtx3, vtx4])

    # 3-D diffraction from grating bottom. Basis vectors are:
    # * b1 ~ parallel to grooves
    # * b2 ~ perpendicular to grooves
    def diffracted_ray(self, u_incident, wlen_um, spec_ord):
        """
        Diffracted ray direction is produced by superposition of the components
        parallel and perpendicular to grating grooves. The groove-parallel 
        component undergoes specular reflection while the groove-perpendicular
        component direction change is determined by path difference. Once those
        in-grating-plane components are found, the grating-normal component 
        follows from the normalization condition (total length == 1).
    
        b_para -- basis vector parallel to grooves
        b_perp -- basis vector perpendicular to grooves
    
        Returns:
        valid     --> True/False to indicate if real-valued solution exists
        diffr_hat --> unit vector in direction of diffracted ray
        """
        bottom = self.get_face('bot')
        b_para, b_perp = bottom['basis']
        g_norm = bottom['normal']
        v_para = np.dot(u_incident, b_para)
        v_perp = np.dot(u_incident, b_perp) \
                    + (spec_ord * wlen_um / self._spacing_um)
        normsq = 1.0 - v_para**2 - v_perp**2
        #sys.stderr.write("NEW v_para, v_perp, normsq: %10.5f, %10.5f, %10.5f\n"
        #        % (v_para, v_perp, normsq))
        if (normsq < 0.0):
            sys.stderr.write("Imaginary diffracted ray!\n")
            sys.stderr.write("normsq: %s\n" % str(normsq))
            return False, np.array([np.nan, np.nan, np.nan])
        v_norm = np.sqrt(normsq)
        #sys.stderr.write("NEW v_para, v_perp, v_norm: %10.5f, %10.5f, %10.5f\n"
        #        % (v_para, v_perp, v_norm))
        diffr_vec = v_para * b_para + v_perp * b_perp + v_norm * g_norm
        diffr_vec /= bottom._vec_length(diffr_vec)
        return True, diffr_vec


######################################################################
# CHANGELOG (polyhedron_optics.py):
#---------------------------------------------------------------------
#
#  2019-03-31:
#     -- Increased __version__ to 0.3.0.
#     -- Added diffraction to grating object and refraction to prism object.
#
#  2019-02-25:
#     -- Increased __version__ to 0.2.5.
#     -- Added new PolygonFace() class to abstract away lots of useful stuff.
#
#  2019-02-24:
#     -- Increased __version__ to 0.2.0.
#     -- Renamed module to polyhedron_optics.py (also renamed classes).
#
#  2019-02-20:
#     -- Increased __version__ to 0.1.0.
#     -- First created polygon_optics.py.
#
