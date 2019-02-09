#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Module to simplify access to line lists and wavelength-related tools.
#
# Rob Siverd
# Created:       2019-02-08
# Last modified: 2019-02-09
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
nist_argon_path = os.path.join(wavelen_dir, 'NIST', 'clean_argon.csv')
lovis_pepe_path = os.path.join(wavelen_dir, 'lovis_pepe_2007', 'lovis_pepe.csv')


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
## Quick ASCII I/O:
#data_file = 'data.txt'
#all_data = np.genfromtxt(data_file, dtype=None)
#all_data = np.genfromtxt(data_file, dtype=None, names=True, autostrip=True)
#all_data = np.genfromtxt(data_file, dtype=None, names=True, autostrip=True,
#                 delimiter='|', comments='%0%0%0%0')
#                 loose=True, invalid_raise=False)
#all_data = aia.read(data_file)
#all_data = pd.read_csv(data_file)
#all_data = pd.read_table(data_file, delim_whitespace=True)
#all_data = pd.read_table(data_file, sep='|')
#fields = all_data.dtype.names
#if not fields:
#    x = all_data[:, 0]
#    y = all_data[:, 1]
#else:
#    x = all_data[fields[0]]
#    y = all_data[fields[1]]

##--------------------------------------------------------------------------##
## Quick FITS I/O:
#data_file = 'image.fits'
#img_vals = pf.getdata(data_file)
#hdr_keys = pf.getheader(data_file)
#img_vals, hdr_keys = pf.getdata(data_file, header=True)
#img_vals, hdr_keys = pf.getdata(data_file, header=True, uint=True) # USHORT
#img_vals, hdr_keys = fitsio.read(data_file, header=True)

#date_obs = hdr_keys['DATE-OBS']
#site_lat = hdr_keys['LATITUDE']
#site_lon = hdr_keys['LONGITUD']

## Initialize time:
#img_time = astt.Time(hdr_keys['DATE-OBS'], scale='utc', format='isot')
#img_time += astt.TimeDelta(0.5 * hdr_keys['EXPTIME'], format='sec')
#jd_image = img_time.jd

## Initialize location:
#observer = ephem.Observer()
#observer.lat = np.radians(site_lat)
#observer.lon = np.radians(site_lon)
#observer.date = img_time.datetime

#pf.writeto('new.fits', img_vals)
#qsave('new.fits', img_vals)
#qsave('new.fits', img_vals, header=hdr_keys)

## Star extraction:
#pse.set_image(img_vals, gain=3.6)
#objlist = pse.analyze(sigthresh=5.0)

##--------------------------------------------------------------------------##
## Solve prep:
#ny, nx = img_vals.shape
#x_list = (0.5 + np.arange(nx)) / nx - 0.5            # relative (centered)
#y_list = (0.5 + np.arange(ny)) / ny - 0.5            # relative (centered)
#xx, yy = np.meshgrid(x_list, y_list)                 # relative (centered)
#xx, yy = np.meshgrid(nx*x_list, ny*y_list)           # absolute (centered)
#xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))   # absolute
#yy, xx = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij') # absolute
#yy, xx = np.nonzero(np.ones_like(img_vals))          # absolute
#yy, xx = np.mgrid[0:ny,   0:nx].astype('uint16')     # absolute (array)
#yy, xx = np.mgrid[1:ny+1, 1:nx+1].astype('uint16')   # absolute (pixel)

## 1-D vectors:
#x_pix, y_pix, ivals = xx.flatten(), yy.flatten(), img_vals.flatten()
#w_vec = np.ones_like(ivals)            # start with uniform weights
#design_matrix = np.column_stack((np.ones(x_pix.size), x_pix, y_pix))

## Image fitting (statsmodels etc.):
#data = sm.datasets.stackloss.load()
#ols_res = sm.OLS(ivals, design_matrix).fit()
#rlm_res = sm.RLM(ivals, design_matrix).fit()
#rlm_model = sm.RLM(ivals, design_matrix, M=sm.robust.norms.HuberT())
#rlm_res = rlm_model.fit()
#data = pd.DataFrame({'xpix':x_pix, 'ypix':y_pix})
#rlm_model = sm.RLM.from_formula("ivals ~ xpix + ypix", data)




######################################################################
# CHANGELOG (wavelength_reference.py):
#---------------------------------------------------------------------
#
#  2019-02-08:
#     -- Increased __version__ to 0.0.1.
#     -- First created wavelength_reference.py.
#
