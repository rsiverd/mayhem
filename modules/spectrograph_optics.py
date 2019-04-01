#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Routines related to spectrograph optical properties and calculations.
#
# Rob Siverd
# Created:       2018-12-29
# Last modified: 2019-03-30
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.2.0"

## Python version-agnostic module reloading:
#try:
#    reload                              # Python 2.7
#except NameError:
#    try:
#        from importlib import reload    # Python 3.4+
#    except ImportError:
#        from imp import reload          # Python 3.0 - 3.3

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
import scipy.optimize as opti
#import scipy.interpolate as stp
#import scipy.spatial.distance as ssd
from functools import partial
#from collections import OrderedDict
#import multiprocessing as mp
#np.set_printoptions(suppress=True, linewidth=160)
#import pandas as pd
#import statsmodels.api as sm
#import statsmodels.formula.api as smf
#from statsmodels.regression.quantile_regression import QuantReg
#import theil_sen as ts
#import window_filter as wf
#import itertools as itt

## Glass indexes of refraction:
from glass import Glass

###--------------------------------------------------------------------------##
### Calculating index of refraction for several glass types:
#class Glass(object):
#
#    def __init__(self, glasstype):
#        self._bcoeffs = {
#            'SiO2':np.array([0.67071081e0, 0.433322857e0, 0.877379057e0]),
#            'LLF1':np.array([1.21640125e0, 1.33664540e-1, 8.83399468e-1]),
#            'PBM2':np.array([1.39446503e0, 1.59230985e-1, 2.45470216e-1]),
#             'LF5':np.array([1.28035628e0, 1.6350597e-1,  8.93930112e-1]), }
#        self._ccoeffs = {
#            'SiO2':np.array([0.00449192312e0, 0.0132812976e0, 95.8899878e0]),
#            'LLF1':np.array([8.57807248e-3,   4.20143003e-2,   1.07593060e+2]),
#            'PBM2':np.array([1.10571872e-2,   5.07194882e-2,   3.14440142e1]),
#             'LF5':np.array([9.29854416e-3,   4.49135769e-2,   1.10493685e2]), }
#        if self._unknown_glass(glasstype):
#            raise
#        self._gtype = glasstype
#        self._coeffs = zip(self._bcoeffs[self._gtype],
#                           self._ccoeffs[self._gtype])
#        return
#
#    def _unknown_glass(self, glasstype):
#        if not glasstype in self._bcoeffs.keys():
#            sys.stderr.write("Unknown glass type: %s\n" % glasstype)
#            return True
#        else:
#            return False
#
#    # Squared index of refraction for specified wavelengths:
#    def refraction_index_squared(self, wlen_um):
#        lam_um_sq = wlen_um**2
#        n_squared = np.ones_like(wlen_um, dtype='float')
#        for bb,cc in self._coeffs:
#            n_squared += (lam_um_sq * bb) / (lam_um_sq - cc)
#        return n_squared
#
#    def refraction_index(self, wlen_um):
#        return np.sqrt(self.refraction_index_squared(wlen_um))
#    
#    def glass_dn_dlambda_easy(self, wlen_um, glasstype, epsfrac=1e-5):
#        #if self._unknown_glass(glasstype):
#        #    raise
#
#        wlen_lower = wlen_um * (1.0 - epsfrac)
#        wlen_upper = wlen_um * (1.0 + epsfrac)
#
#        nn_lower = self.refraction_index(wlen_lower) #, glasstype)
#        nn_upper = self.refraction_index(wlen_upper) #, glasstype)
#
#        return (nn_upper - nn_lower) / (wlen_upper - wlen_lower)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Prism deflection:
#def prism_deflection_n(incid_r, apex_r, n):
#    ptemp = np.sqrt(n**2 - np.sin(incid_r)**2) * np.sin(apex_r) \
#                - np.cos(apex_r) * np.sin(incid_r)
#    return incid_r - apex_r + np.arcsin(ptemp)

#def prism_deflection_glass_wl(incid_r, apex_r, gtype, wlen):

#def wiki_prism_deflection_n(i, A, n):
#    return i - A + np.arcsin(n * np.sin(A - np.arcsin(np.sin(i) / n)))

##--------------------------------------------------------------------------##
## Prism object to calculate deflections for the specified material and shape:
class Prism(object):

    def __init__(self, glasstype, apex_deg, n_air=1.0):
        #self._apex_deg = apex_deg
        #self._apex_rad = np.radians(apex_deg)
        self._material = Glass(glasstype)
        self.set_apex_deg(apex_deg)
        self.set_n_air(n_air)    # index of refraction of air in spectrograph
        return

    # -------------------------------------
    # Getters and setters:
    def get_apex_rad(self):
        return self._apex_rad

    def get_apex_deg(self):
        return self._apex_deg

    def set_apex_deg(self, apex_deg):
        self._apex_deg = apex_deg
        self._apex_rad = np.radians(apex_deg)
        return

    def set_n_air(self, index_refr):
        self._n_air = index_refr
        return

    # -------------------------------------

    @staticmethod
    def _wiki_deflection_rad_nr(i, A, nr):
        """Deflection angle using formula from Wikipedia (which seems
        to be identical to others used here. Inputs are:
        i  -- incidence angle (RADIANS)
        A  -- prism apex angle (RADIANS)
        nr -- index of refraction RATIO at desired wavelength(s): 
                    n_glass / n_air
        """
        return i - A + np.arcsin(nr * np.sin(A - np.arcsin(np.sin(i) / nr)))

    def deflection_rad_nr(self, incidence_r, nr):
        """Deflection angle as a function of incidence angle and 
        index of refraction RATIO (n_glass / n_air)."""
        return self._wiki_deflection_rad_nr(incidence_r, self._apex_rad, nr)

    def deflection_rad_nr2(self, incidence_r, nr2):
        """Calculate deflection angle from incidence and SQUARED index of
        refraction RATIO (glass / air). Units of RADIANS used throughout."""
        ptemp = np.sqrt(nr2 - np.sin(incidence_r)**2) * np.sin(self._apex_rad) \
                - np.cos(self._apex_rad) * np.sin(incidence_r)
        return incidence_r - self._apex_rad + np.arcsin(ptemp)

    def deflection_rad_wl(self, incidence_r, wavelength_um):
        """Calculate deflection angle given incidence angle and wavelength
        in microns. Units of RADIANS used throughout."""
        n2_glass = self._material.refraction_index_squared(wavelength_um)
        n2_air   = self._n_air**2
        return self.deflection_rad_nr2(incidence_r, n2_glass/n2_air)
        #ptemp = np.sqrt(n2 - np.sin(incidence_r)**2) * np.sin(self._apex_rad) \
        #        - np.cos(self._apex_rad) * np.sin(incidence_r)
        #return incidence_r - self._apex_rad + np.arcsin(ptemp)

    def deflection_deg_wl(self, incidence_d, wavelength_um):
        """Calculate deflection angle given incidence angle and wavelength
        in microns. Units of RADIANS used throughout."""
        incidence_r = np.radians(incidence_deg)
        return np.degrees(self.deflection_rad_wl(incidence_r, wavelength_um))


##--------------------------------------------------------------------------##
## Notes on notation, relations, identities, etc.:

# theta_B  --> blaze angle
# lambda_B --> blaze wavelength
# theta    --> facet illumination angle
# alpha    --> angle of incidence
# beta     --> angle of diffraction
# gamma    --> angle of incidence relative to facet normal SPECIFICALLY
#                   in the plane that is parallel to grooves

# Related by:
# alpha = theta_B + theta
# beta  = theta_B - theta

# Note that:
# * theta is SMALL in Littrow configuration, alpha ~= beta

## Grating methods implemented:
class GratingTools(object):

    #groove_spacing_units = {'mm':1.0, 'um':1e3, 'nm':1e6}

    def __init__(self, gratio=None, lines_per_mm=None):
        # Placeholders:
        self._blaze_angle_rad = 0.0
        self._g_spacing_units = {'mm':1.0, 'um':1e3, 'nm':1e6}
        self._groove_spacing = {x:0.0 for x in self._g_spacing_units.keys()}

        # Handle useful inputs:
        if gratio:
            self._blaze_angle_rad = self._calc_blaze_angle(gratio)
        if lines_per_mm:
            for kk,mult in self._g_spacing_units.items():
                self._groove_spacing[kk] = mult / lines_per_mm

        return

    # --------------------------------------
    # Blaze angle and wavelength:

    @staticmethod
    def _calc_blaze_angle(ratio):
        """For R4 grating, use ratio = 4.  Blaze angle returns in RADIANS."""
        return np.arctan(ratio)

    @staticmethod
    def _calc_blaze_wavelength(spec_orders, groove_spacing,
            theta_B_rad, facet_rad, gamma_rad):
        """
        Compute blaze wavelength for an array of order numbers.
        NOTE: wavelength computed in same units as groove spacing.
        """
        ord_vec = np.atleast_1d(spec_orders)
        if not np.all(ord_vec == np.int_(ord_vec)):
            sys.stderr.write("Error: non-integer spec_orders provided!\n")
            return np.nan * spec_orders

        ang_tmp = np.sin(theta_B_rad) * np.cos(facet_rad) * np.cos(gamma_rad)
        return 2.0 * groove_spacing * ang_tmp / np.float_(spec_orders)

    def get_blaze_wavelengths(self, spec_orders, facet_rad=0.0, gamma_rad=0.0,
            units='nm'):
        return self._calc_blaze_wavelength(spec_orders,
                self._groove_spacing[units], self._blaze_angle_rad, 
                facet_rad=facet_rad, gamma_rad=gamma_rad)

    @staticmethod
    def _calc_line_tilt_ctr(blaze_rad, gamma_rad):
        """Calculate line tilt at order center (at blaze wavelength).
        Formulae from Barnes (2004) thesis, page 11:
        tan(φ) = (sin(α) + sin(β)) / cos(β) * (sin(γ) / cos(γ))
               = λ * dβ/dλ * tan(γ)
        tan(φ) = 2 * tan(θ_B) * tan(γ)      # at blaze wavelength λ_B
        """
        return np.arctan(2.0 * np.tan(blaze_rad) * np.tan(gamma_rad))

    @staticmethod
    def _calc_line_tilt_any(alpha_rad, beta_rad, gamma_rad):
        """Calculate line tilt at order center (at blaze wavelength).
        Formulae from Barnes (2004) thesis, page 11:
        tan(φ) = (sin(α) + sin(β)) / cos(β) * (sin(γ) / cos(γ))
               = λ * dβ/dλ * tan(γ)
        tan(φ) = 2 * tan(θ_B) * tan(γ)      # at blaze wavelength λ_B
        """
        numer = (np.sin(alpha_rad) + np.sin(beta_rad)) * np.sin(gamma_rad)
        denom = np.cos(beta_rad) * np.cos(gamma_rad)
        return np.arctan(numer / denom)

    @staticmethod
    def _calc_order_FSR(wlen_cen, spec_orders):
        return wlen_cen / spec_orders
    
    @staticmethod
    def _calc_order_width(wlen_cen, groove_spacing, 
            theta_B_rad, facet_rad, gamma_rad):
        """
        Compute angular span of the order (free spectral range) with central
        wavelength wlen_cen. 
        
        NOTES: 
        * groove_spacing and wlen_cen need same units for correct result!
        * diff_angle_rad should correspond to central wavelength (wlen_cen)
        * wlen_cen supports numpy arrays
        """
        diff_angle_rad = theta_B_rad - facet_rad
        cos_diff_gamma = np.cos(diff_angle_rad) * np.cos(gamma_rad)
        return wlen_cen / (groove_spacing * cos_diff_gamma)

    def get_order_params(self, spec_orders, facet_rad=0.0, gamma_rad=0.0,
            units='nm'):
        """
        Computes central (blaze) wavelength and angular extent of specified
        orders using chosen geometry/angles.

        Outputs:
        central_wlen -- corresponding central wavelengths in requested units
        order_FSR    -- order free spectral ranges (units of central_wlen)
        order_angwid -- order angular width (RADIANS)
        """
        use_spacing = self._groove_spacing[units]
        central_wlen = self._calc_blaze_wavelength(spec_orders, use_spacing,
                self._blaze_angle_rad, facet_rad, gamma_rad)
        order_FSR = self._calc_order_FSR(central_wlen, spec_orders)
        order_angwid = self._calc_order_width(central_wlen, use_spacing,
                self._blaze_angle_rad, facet_rad, gamma_rad)
        return (central_wlen, order_FSR, order_angwid)


##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Grating class (to be overloaded for double pass:

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## More sophisticated grating+prism class (closer to real design):
class DoublePassPrismGratingPrism(object):

    def __init__(self, config={}):
        # PARAMETERS LIST:
        # 0. coordinate system:
        #       * x-axis points "right" towards input fiber (grating at -X).
        #       * y-axis points "ahead" from camera to grating (optical axis)
        #       * z-axis points "up" towards sky from bench surface
        # 1. prism:
        #       * glass type
        #       * apex angle (degrees)
        #
        # 2. grating:
        #       * blaze ratio (i.e., R4)
        #       * groove spacing in lines/mm
        #
        # 3. layout:
        #       * angle of initial beam w.r.t. optical axis (assume 0?)
        #       * rotation of prism apex w.r.t. optical axis. Let 0 degrees
        #           represent apex perpendicular to optical axis.
        #       * direction to grating w.r.t. optical axis. This is really
        #           the prism deflection angle that "points" to the grating
        #       * alpha or sine(alpha), tilt of grating w.r.t. optical bench.
        #           --> could also be expressed as facet angle
        #       * rotation of grating w.r.t. optical bench normal. This is
        #           effectively the nominal gamma angle of the ray pointing
        #           directly at the grating. EQUIVALENTLY can use orientation
        #           of grating on optical bench (same units as prism base)

        # For prism apex 'perpendicular' to optical axis, incidence ~ apex/2
        self.prism_glass    = config.get('prism_glass', "PBM2")
        self.apex_angle_deg = config.get('apex_angle_deg', 55.0)

        # How far prism is "turned" w.r.t. main optical axis (def from CAD):
        self.prism_turn_deg = config.get('prism_turn_deg', 23.507)

        # Turn of incident beam (from fiber/collimator) 
        # w.r.t. main optical axis (positive angles towards grating):
        # MAY VARY WITH WAVELENGTH???
        self.input_turn_deg = config.get('input_turn_deg', 2.0)

        # Allow non-unity index of refraction for air:
        self.air_index_refr = config.get('air_index_refr', 1.0)

        # Front and rear prism face orientations w.r.t. z-axis:
        self.prism_face1_deg = self.prism_turn_deg + 0.5 * self.apex_angle_deg
        self.prism_face2_deg = self.prism_turn_deg - 0.5 * self.apex_angle_deg

        # Angle of incidence of input beam w.r.t. prism face1 normal:
        self.prism_front_incid_deg = self.prism_face1_deg - self.input_turn_deg

        # Create prism object with given settings:
        self.pr_obj = Prism(self.prism_glass, self.apex_angle_deg)
        self.pr_obj.set_n_air(self.air_index_refr)

        # ---------------------------------------------------

        # Grating ruling in lines / millimeter:
        self.grating_ruling_lmm = config.get('grating_ruling_lmm', 41.59)
        self.grating_spacing_um = 1e3 / self.grating_ruling_lmm

        # Grating face angle w.r.t. optical axis (turn about bench normal):
        self.grating_turn_deg = config.get('grating_turn_deg', 44.827)

        # Grating face tilt w.r.t. optical bench:
        self.grating_tilt_deg = config.get('grating_tilt_deg', 13.786)

        # Derived quantities:
        self.input_turn_rad   = np.radians(self.input_turn_deg)
        self.grating_turn_rad = np.radians(self.grating_turn_deg)
        self.grating_tilt_rad = np.radians(self.grating_tilt_deg)

        self.alpha_angle_rad  = np.radians(90.0 - self.grating_tilt_deg)
        self.blaze_angle_rad  = config.get('blaze_angle_rad', np.arctan(4.))
        self.facet_angle_rad  = self.alpha_angle_rad - self.blaze_angle_rad

        self.fixed_geometry   = 2.0 * self.grating_spacing_um \
                * np.cos(self.facet_angle_rad) * np.sin(self.blaze_angle_rad)

        ## deal with angle changes between prism/grating for second prism pass:
        #pass1_out_heading = self.input_turn_deg + np.degrees(deflect_r)
        #gamma_eff_grating = pass1_out_heading - self.grating_turn_rad
        #pass1_out_incid_d = pass1_out_heading - self.prism_face2_deg    # emergent
        #pass2_new_incid_d = pass1_out_incid_d - 2.0 * gamma_eff_grating

        # Lastly include some distances:
        self.coll_focallen_mm = config.get('coll_focallen_mm', 380.0)
        #self.coll_focallen_mm = 375.0
        #self.prism_grating_mm = 284.0       # approximately
        #self.prism_grating_mm = 275.0       # approximately
        #self.prism_grating_mm = 250.0       # approximately
        self.prism_grating_mm = config.get('prism_grating_mm', 100.0)
        self.lens_compression = config.get('lens_compression',   2.0)
        return

    # ---------------------------------------------------------
    # Joint prism+grating central wavelength and deflection:
    def _prism_pass1_deflection(self, wavelength_um):
        """
        Calculate the deflection angle and new direction of travel for the
        specified wavelength after first pass through prism (using configured
        prism and layout). Assumes all rays incident on prism are parallel.
        Output angles are in RADIANS.
        """
        incid_1_r = np.radians(self.prism_front_incid_deg)
        deflect_r = self.pr_obj.deflection_rad_wl(incid_1_r, wavelength_um)
        return deflect_r

    def _deflection_gamma(self, wavelength_um):
        """
        Calculate the effective gamma angle at the grating for the specified
        wavelength using configured prism and layout. Assumes all rays
        incident on prism are parallel. Output angle in RADIANS.
        """
        deflect_r = self._prism_pass1_deflection(wavelength_um)
        #incid_1_r = np.radians(self.prism_front_incid_deg)
        #deflect_r = self.pr_obj.deflection_rad_wl(incid_1_r, wavelength_um)
        #heading_r = self.input_turn_deg + deflect_r     # outbound direction
        return self.input_turn_rad + deflect_r - self.grating_turn_rad

    def _lamcen_residual(self, wavelength_um, order=0):
        ls = wavelength_um * order
        rs = self.fixed_geometry * np.cos(self._deflection_gamma(wavelength_um))
        return ls - rs

    def _iter_calc_lamcen(self, order):
        kw = {'order':order}
        runme = partial(self._lamcen_residual, **kw)
        return opti.bisect(runme, 0.0, 10.0)

    def calc_central_wlen_um(self, spec_order_list):
        return np.array([self._iter_calc_lamcen(x) for x in spec_order_list])

    def gamma_from_wlen_um(self, wlen_list_um):
        return np.array([self._deflection_gamma(x) for x in wlen_list_um])

    def fancy_deflections(self, spec_order_list):
        order_ctr_wlen_um = self.calc_central_wlen_um(spec_order_list)
        order_ctr_gamma_r = self.gamma_from_wlen_um(order_ctr_wlen_um)
        return order_ctr_wlen_um, order_ctr_gamma_r

    def two_pass_deflection(self, wavelength_um):
        """
        incid_1_r   --> angle of incidence of light onto prism (first pass)
        deflect_1_r --> deflection angle produced by first pass
        heading_1_r --> direction of deflected ray in benchtop units
        gamma_eff_r --> effective gamma angle at grating for this wavelength
        excid_1_r   --> angle of deflected ray to exit face normal

        incid_2_r   --> angle of incidence on prism face 2 (second pass)
        deflect_2_r --> deflection angle produced by second pass
        """
        incid_1_r   = np.radians(self.prism_front_incid_deg)
        deflect_1_r = self.pr_obj.deflection_rad_wl(incid_1_r, wavelength_um)
        heading_1_r = self.input_turn_rad + deflect_1_r
        gamma_eff_r = heading_1_r - self.grating_turn_rad
        excid_1_r   = heading_1_r - np.radians(self.prism_face2_deg)

        # Non-zero gamma causes change in incidence for second pass:
        incid_2_r = excid_1_r - 2.0 * gamma_eff_r
        refl_beam_r = heading_1_r + np.pi - 2 * gamma_eff_r
        deflect_2_r = self.pr_obj.deflection_rad_wl(incid_2_r, wavelength_um)
        heading_2_r = refl_beam_r - deflect_2_r

        # Total Y-shift produced by a combination of prism<-->grating offset
        # and deflection over the collimator focal length:
        pg_yshift =  2.0 * self.prism_grating_mm * np.tan(gamma_eff_r)
        pc_yshift = -1.0 * self.coll_focallen_mm * np.tan(heading_2_r - np.pi)
        wlen_nm = 1e3 * wavelength_um
        gamma_deg = np.degrees(gamma_eff_r)
        sys.stderr.write("λ=%6.1f nm, γ= %6.3f: PG,PC shifts: %6.2f, %6.2f | i2=%5.2f\n"
                % (wlen_nm, gamma_deg, pg_yshift, pc_yshift, np.degrees(incid_2_r)))
        return heading_2_r, pg_yshift, pc_yshift
        #return deflect_1_r, refl_beam_r, deflect_2_r, heading_2_r

##--------------------------------------------------------------------------##




######################################################################
# CHANGELOG (spectrograph_optics.py):
#---------------------------------------------------------------------
#
#  2018-12-29:
#     -- Increased __version__ to 0.1.0.
#     -- First created spectrograph_optics.py.
#
