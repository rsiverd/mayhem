#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Module to simplify access to line lists and wavelength-related tools.
#
# Rob Siverd
# Created:       2019-02-08
# Last modified: 2019-03-06
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
import time
import numpy as np
from numpy.lib.recfunctions import append_fields
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
import pandas as pd
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

## Fast FITS I/O:
#try:
#    import fitsio
#except ImportError:
#    sys.stderr.write("\nError: fitsio module not found!\n")
#    sys.exit(1)

## FITS I/O:
try:
    import astropy.io.fits as pf
except ImportError:
    try:
       import pyfits as pf
    except ImportError:
        sys.stderr.write("\nError!  No FITS I/O module found!\n"
               "Install either astropy.io.fits or pyfits and try again!\n\n")
        sys.exit(1)

##--------------------------------------------------------------------------##
## Save FITS image with clobber (astropy / pyfits):
#def qsave(iname, idata, header=None, padkeys=1000, **kwargs):
#    this_func = sys._getframe().f_code.co_name
#    sys.stderr.write("Writing to '%s' ... " % iname)
#    if header:
#        while (len(header) < padkeys):
#            header.append() # pad header
#    if os.path.isfile(iname):
#        os.remove(iname)
#    pf.writeto(iname, idata, header=header, **kwargs)
#    sys.stderr.write("done.\n")

##--------------------------------------------------------------------------##
## Save FITS image with clobber (fitsio):
#def qsave(iname, idata, header=None, **kwargs):
#    this_func = sys._getframe().f_code.co_name
#    sys.stderr.write("Writing to '%s' ... " % iname)
#    #if os.path.isfile(iname):
#    #    os.remove(iname)
#    fitsio.write(iname, idata, clobber=True, header=header, **kwargs)
#    sys.stderr.write("done.\n")


##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Refractive index of air calculation (from Bacigalupo thesis):
def n_air(wlen_um, t=18., p=101325.):
    nminus1 = 0.0472326 / (173.3 - (1.0 / wlen_um)**2)
    return 1.0 + nminus1

## IAU vacuum -> air wavelength conversion (Angstroms):
## https://classic.sdss.org/dr7/products/spectra/vacwavelength.html
##   AIR = VAC / (1.0 + 2.735182E-4 + 131.4182 / VAC^2 + 2.76249E8 / VAC^4)
def iau_vac2air(lam_vac_aa):
    denom = (1.0 + 2.735182e-4 + 131.4182 / lam_vac_aa**2 \
                                + 2.76249e8 / lam_vac_aa**4)
    return lam_vac_aa / denom

## Tim's code to convert air --> vacuum wavelength:
#function airlam,lam,z0
#; This routine returns lam(vacuum), with input lam(air) in microns,
#; with pressure such that n(NaD = 0.589 micron)=1.+z0.
#; For dry air at STP, set z0 = 0.000295
#; result is returned in a double-precision vector
#
#lamm2=1./lam^2
#nam2=1./0.589d0^2
#zs=0.05792105d0/(238.0185d0-lamm2) + 0.00167917d0/(57.362d0-lamm2) 
#zna=0.05792105d0/(238.0185d0-nam2) + 0.00167917d0/(57.362d0-nam2) 
#
#zair=zs*z0/zna
#
#return,lam*(1.d0+zair)


##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Track down line lists:
mayhem_root = os.path.join(os.getenv('HOME'), 'NRES', 'mayhem')
if not os.path.isdir(mayhem_root):
    sys.stderr.write("Folder not found: %s\n" % mayhem_root)
    mayhem_root = None
    raise
wavelen_dir = os.path.join(mayhem_root, 'wavelength')
#if not os.path.isdir(wavelen_dir):
#    sys.stderr.write("Folder not found: %s\n" % wavelen_dir)
#    raise

## NIST Argon list (cleaned up a bit):
nist_dir = os.path.join(wavelen_dir, 'NIST')
lope_dir = os.path.join(wavelen_dir, 'lovis_pepe_2007')
nist_thar_path  = os.path.join(nist_dir, 'NIST_ThAr_lines.fits')
nist_argon_path = os.path.join(nist_dir, 'clean_argon.csv')
lovis_pepe_path = os.path.join(lope_dir, 'lovis_pepe.csv')
nist_fft_thar_path = os.path.join(nist_dir, 'NIST_spectrum.all.fits')


##--------------------------------------------------------------------------##
## Load NIST Argon line list. Data file columns include:
##  lam_obs_vac_nm, err_lam_obs_vac_nm
##  lam_ritz_vac_nm, err_lam_ritz_vac_nm
##  rel_intensity, einstein_aki, ion_name
def load_nist_argon_pd(data_file=nist_argon_path):
    ndata = pd.read_csv(data_file)
    ndata = ndata.assign(lam_obs_nm=lambda x: x.lam_obs_vac_nm)
    ndata = ndata.assign(lam_obs_aa=lambda x: x.lam_obs_vac_nm * 10.)
    ndata = ndata.assign(lam_obs_um=lambda x: x.lam_obs_vac_nm / 1e3)
    return ndata

def load_nist_argon(data_file=nist_argon_path):
    nist_data = np.genfromtxt(data_file, dtype=None, names=True)
    lam_vac_um = nist_data['lam_obs'] / 1e3
    nist_data = append_fields(nist_data, 'lam_obs_um', lam_obs_um)
    nist_data = append_fields(nist_data, 'lam_obs_nm', nist_data['lam_obs'])
    return nist_data

def load_nist_thar(data_file=nist_thar_path):
    ndata = pd.DataFrame(data_file)
    lam_vac_um = nist_data['lam_obs'] / 1e3
    nist_data = append_fields(nist_data, 'lam_obs_um', lam_obs_um)
    nist_data = append_fields(nist_data, 'lam_obs_nm', nist_data['lam_obs'])
    return nist_data

##--------------------------------------------------------------------------##
## Load Lovis & Pepe line list:
def load_lovis_pepe_thar(data_file=lovis_pepe_path):
    lope_data = np.genfromtxt(lovis_pepe_path, dtype=None, names=True,
            delimiter=',')
    lam_vac_nm = lope_data['lam_vac'] / 10.0
    lope_data = append_fields(lope_data, 'lam_vac_nm', lam_vac_nm)
    lope_data = append_fields(lope_data, 'lam_vac_um', lam_vac_nm / 1e3)
    return lope_data


##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## On-demand delivery of line lists for wavelength calibration:
class WLFetcher(object):

    def __init__(self, vlevel=0,
            nist_argon_path=nist_argon_path,
            nist_thar_path=nist_thar_path,
            lope_thar_path=lovis_pepe_path):
        self._vlevel    = vlevel
        self._nist_argon_path = nist_argon_path
        self._nist_thar_path  = nist_thar_path
        self._lope_thar_path  = lope_thar_path
        self._nist_argon = None
        self._nist_thar  = None
        self._lope_data  = None
        self.load_lines()
        return

    # --------------------------------------------
    # Verbosity-specific messages:
    def vlwrite(self, msg_vlevel, msg_text):
        if (self._vlevel >= msg_vlevel):
            sys.stderr.write(msg_text)
        return

    # --------------------------------------------
    # Load the NIST Argon line list:
    @staticmethod
    def _load_nist_argon(data_file):
        ndata = pd.read_csv(data_file)
        ndata = ndata.assign(lam_obs_nm=lambda x: x.lam_obs_vac_nm)
        ndata = ndata.assign(lam_obs_aa=lambda x: x.lam_obs_vac_nm * 10.)
        ndata = ndata.assign(lam_obs_um=lambda x: x.lam_obs_vac_nm / 1e3)
        return ndata
 
    # Load lines I extracted from the NIST ThAr reference spectrum:
    @staticmethod
    def _load_nist_thar(data_file):
        ndata = pd.DataFrame(pf.getdata(data_file))
        ndata = ndata.assign(lam_obs_vac_nm=lambda x: x.lam_obs_vac_aa / 10.)
        ndata = ndata.assign(lam_obs_vac_um=lambda x: x.lam_obs_vac_nm / 1e3)
        ndata = ndata.assign(lam_obs_aa=lambda x: x.lam_obs_vac_aa)
        ndata = ndata.assign(lam_obs_nm=lambda x: x.lam_obs_vac_nm)
        ndata = ndata.assign(lam_obs_um=lambda x: x.lam_obs_vac_um)
        return ndata

    # Load Lovis & Pepe line list:
    @staticmethod
    def _load_lovis_pepe_thar(data_file):
        ldata = pd.read_csv(data_file)
        ldata = ldata.assign(lam_vac_nm=lambda x: x.lam_vac / 10.)
        ldata = ldata.assign(lam_vac_um=lambda x: x.lam_vac_nm / 1e3)
        return ldata

    ## Load Lovis & Pepe line list:
    #@staticmethod
    #def _load_lovis_pepe_thar_old(data_file):
    #    lope_data = np.genfromtxt(data_file, 
    #            dtype=None, names=True, delimiter=',')
    #    lam_vac_nm = lope_data['lam_vac'] / 10.0
    #    lope_data = append_fields(lope_data, 'lam_vac_nm', lam_vac_nm)
    #    lope_data = append_fields(lope_data, 'lam_vac_um', lam_vac_nm / 1e3)
    #    return lope_data

    # Load all the line lists:
    def load_lines(self):
        self._nist_argon = self._load_nist_argon(self._nist_argon_path)
        self._nist_thar  = self._load_nist_thar(self._nist_thar_path)
        #self._nist_fcol = 'rel_intensity'
        self._lope_data  = self._load_lovis_pepe_thar(self._lope_thar_path)
        #self._lope_fcol = 'flux'
        return

    # --------------------------------------------
 
    # Line selection with some smarts:
    def get_nist_argon_lines(self, wl1, wl2, reltol=0.001, minflx=100.,
            lamcol='lam_obs_nm', flxcol='rel_intensity'):
        twlen = self._nist_argon[lamcol]
        tflux = self._nist_argon[flxcol]
        which = (wl1 < twlen) & (twlen < wl2) & (tflux > minflx)
        neato = self._nist_argon[which]
        nkept = which.sum()
        self.vlwrite(1, "NIST Argon: %d lines found.\n" % nkept)
        if self._vlevel >= 2:
            for ww,ff in neato[[lamcol, flxcol]].values:
                sys.stderr.write("%10.5f --- %10.3f\n" % (ww, ff))
        if not nkept:
            return np.array([]), np.array([])

        # stick to brightest of lines found:
        thresh = neato[flxcol].max() * reltol
        #sys.stderr.write("thresh: %10.5f\n" % thresh)
        smart = (neato[flxcol] >= thresh)   # relative to high peak
        bright = neato[smart]
        self.vlwrite(1, "After peak-rel-cut, have %d lines.\n" % smart.sum())
        return (bright[lamcol].values, bright[flxcol].values)

    # Line selection with some smarts:
    def get_nist_thar_lines(self, wl1, wl2, reltol=0.001, minflx=0.0,
            lamcol='lam_obs_nm', flxcol='rel_intensity'):
        twlen = self._nist_thar[lamcol]
        tflux = self._nist_thar[flxcol]
        which = (wl1 < twlen) & (twlen < wl2) & (tflux > minflx)
        neato = self._nist_thar[which]
        nkept = which.sum()
        self.vlwrite(1, "NIST ThAr: %d lines found.\n" % nkept)
        if self._vlevel >= 2:
            for ww,ff in neato[[lamcol, flxcol]].values:
                sys.stderr.write("%10.5f --- %10.3f\n" % (ww, ff))
        if not nkept:
            return np.array([])
    
        # stick to brightest of lines found:
        thresh = neato[flxcol].max() * reltol
        #sys.stderr.write("thresh: %10.5f\n" % thresh)
        smart = (neato[flxcol] >= thresh)   # relative to high peak
        bright = neato[smart]
        self.vlwrite(1, "After peak-rel-cut, have %d lines.\n" % smart.sum())
        return bright[lamcol].values, bright[flxcol].values

    # Retrieve lines by wavelength range from Lovis+Pepe (2007):
    def get_lope_thar_lines(self, wl1, wl2, nmax=30, reltol=0.1, minflx=100.,
            lamcol='lam_vac_nm', flxcol='flux'):
        twlen = self._lope_data[lamcol]
        tflux = self._lope_data[flxcol]
        flcut = np.percentile(tflux, 75)
        #which = (wl1 < twlen) & (twlen < wl2) & (tflux > flcut)
        which = (wl1 < twlen) & (twlen < wl2) & (tflux > minflx)
        neato = self._lope_data[which]
        nkept = which.sum()
        self.vlwrite(1, "Lovis+Pepe: %d lines found.\n" % nkept)
        if self._vlevel >= 2:
            for ww,ff in zip(neato[lamcol], neato[flxcol]):
                sys.stderr.write("%10.5f --- %10.3f\n" % (ww, ff))
        if not nkept:
            return np.array([]), np.array([])
    
        # stick to brightest of lines found:
        thresh = neato[flxcol].max() * reltol
        smart = (neato[flxcol] >= thresh)   # relative to high peak
        bright = neato[smart]
        self.vlwrite(1, "After peak-rel-cut, have %d lines.\n" % smart.sum())
    
        #top_few_idx = np.argsort(neato[flxcol])[-nmax:]
        top_few_idx = np.argsort(bright[flxcol])[-nmax:]
        self.vlwrite(1, "Selecting top %d with highest flux ...\n" % nmax)
        return (bright[lamcol].values[top_few_idx],
                bright[flxcol].values[top_few_idx])

        #return neato[lamcol][top_few_idx].values

    def get_combined_lines(self, wl1, wl2, mtol=0.001):
        lam_nist, flx_nist = self.get_nist_argon_lines(wl1, wl2)
        lam_lope, flx_lope = self.get_lope_thar_lines(wl1, wl2)
        lam_comb = np.array(lam_nist.tolist() + lam_lope.tolist())
        flx_comb = np.array(flx_nist.tolist() + flx_lope.tolist())
        wl_order = np.argsort(lam_comb)
        return (lam_comb[wl_order], flx_comb[wl_order])
        #ll_comb = np.array(sorted(ll_lope.tolist() + ll_nist.tolist()))
        #return ll_comb


    # --------------------------------------------





######################################################################
# CHANGELOG (wavelength_reference.py):
#---------------------------------------------------------------------
#
#  2019-02-08:
#     -- Increased __version__ to 0.0.1.
#     -- First created wavelength_reference.py.
#
