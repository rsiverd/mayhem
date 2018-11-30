#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Uncomplicated spectrum extraction.
#
# Rob Siverd
# Created:       2018-08-03
# Last modified: 2018-11-30
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
import argparse
import signal
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
#from collections import OrderedDict
import matplotlib.pyplot as plt
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

import nres_extraction
reload(nres_extraction)
trio = nres_extraction.TraceIO()
nrex = nres_extraction

import my_kde
reload(my_kde)
mk = my_kde

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

## Time conversion:
#try:
#    import astropy.time as astt
#except ImportError:
#    sys.stderr.write("\nError: astropy module not found!\n"
#           "Please install and try again.\n\n")
#    sys.exit(1)

##--------------------------------------------------------------------------##
## Colors for fancy terminal output:
NRED    = '\033[0;31m'   ;  BRED    = '\033[1;31m'
NGREEN  = '\033[0;32m'   ;  BGREEN  = '\033[1;32m'
NYELLOW = '\033[0;33m'   ;  BYELLOW = '\033[1;33m'
NBLUE   = '\033[0;34m'   ;  BBLUE   = '\033[1;34m'
NMAG    = '\033[0;35m'   ;  BMAG    = '\033[1;35m'
NCYAN   = '\033[0;36m'   ;  BCYAN   = '\033[1;36m'
NWHITE  = '\033[0;37m'   ;  BWHITE  = '\033[1;37m'
ENDC    = '\033[0m'

## Suppress colors in cron jobs:
if (os.getenv('FUNCDEF') == '--nocolors'):
    NRED    = ''   ;  BRED    = ''
    NGREEN  = ''   ;  BGREEN  = ''
    NYELLOW = ''   ;  BYELLOW = ''
    NBLUE   = ''   ;  BBLUE   = ''
    NMAG    = ''   ;  BMAG    = ''
    NCYAN   = ''   ;  BCYAN   = ''
    NWHITE  = ''   ;  BWHITE  = ''
    ENDC    = ''

## Fancy text:
degree_sign = u'\N{DEGREE SIGN}'

## Dividers:
halfdiv = "----------------------------------------"
fulldiv = halfdiv + halfdiv

##--------------------------------------------------------------------------##
## Catch interruption cleanly:
def signal_handler(signum, frame):
    sys.stderr.write("\nInterrupted!\n\n")
    sys.exit(1)

signal.signal(signal.SIGINT, signal_handler)

##--------------------------------------------------------------------------##
## Save FITS image with clobber (astropy / pyfits):
def qsave(iname, idata, header=None, padkeys=1000, **kwargs):
    this_func = sys._getframe().f_code.co_name
    sys.stderr.write("Writing to '%s' ... " % iname)
    if header:
        while (len(header) < padkeys):
            header.append() # pad header
    if os.path.isfile(iname):
        os.remove(iname)
    pf.writeto(iname, idata, header=header, **kwargs)
    sys.stderr.write("done.\n")

##--------------------------------------------------------------------------##
def ldmap(things):
    return dict(zip(things, range(len(things))))

def argnear(vec, val):
    return (np.abs(vec - val)).argmin()

## Robust location/scale estimate using median/MAD:
def calc_ls_med_MAD(a, axis=None):
    """Return median and median absolute deviation of *a* (scaled to normal)."""
    med_val = np.median(a, axis=axis)
    sig_hat = (1.482602218 * np.median(np.abs(a - med_val), axis=axis))
    return (med_val, sig_hat)

## Robust location/scale estimate using median/IQR:
def calc_ls_med_IQR(a, axis=None):
    """Return median and inter-quartile range of *a* (scaled to normal)."""
    pctiles = np.percentile(a, [25, 50, 75], axis=axis)
    med_val = pctiles[1]
    sig_hat = (0.741301109 * (pctiles[2] - pctiles[0]))
    return (med_val, sig_hat)

## Select inliners given specified sigma threshold:
def pick_inliers(data, sig_thresh):
    med, sig = calc_ls_med_IQR(data)
    return ((np.abs(data - med) / sig) <= sig_thresh)

##--------------------------------------------------------------------------##
## Find distribution modal value using KDE:
def calc_kde_mode(data, save_kde=False, **kwargs):
    """
    Find the peak of the 'data' distribution using a KDE. 

    Returns:
        peak                            if save_kde=False (default)
        peak, (kde_grid, kde_vals)      if save_kde=True
    """
    kde_grid, kde_vals = mk.go(data, **kwargs)
    peak_val = kde_grid[kde_vals.argmax()]
    if save_kde:
        return peak_val, (kde_grid, kde_vals)
    else:
        return peak_val

##--------------------------------------------------------------------------##
## Parse arguments and run script:
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)
if __name__ == '__main__':

    # ------------------------------------------------------------------
    descr_txt = """
    Uncomplicated spectrum extraction.
    
    Version: %s
    """ % __version__
    parser = argparse.ArgumentParser(
            prog=os.path.basename(__file__),
            description=descr_txt)
    parser = MyParser(prog=os.path.basename(__file__), description=descr_txt)
    # ------------------------------------------------------------------
    parser.set_defaults(extr_method="dumb")
    # ------------------------------------------------------------------
    #parser.add_argument('-d', '--dayshift', required=False, default=0,
    #        help='Switch between days (1=tom, 0=today, -1=yest', type=int)
    parser.add_argument('-q', '--quiet', action='count', default=0,
            help='less progress/status reporting')
    parser.add_argument('-v', '--verbose', action='count', default=0,
            help='more progress/status reporting')
    parser.add_argument('--debug', dest='debug', default=False,
            help='Enable extra debugging messages', action='store_true')
    #parser.add_argument('remainder', help='other stuff', nargs='*')
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    ifgroup = parser.add_argument_group('Required I/O')
    ifgroup.add_argument('-S', '--spectrum', default=None, required=True,
            help='image with spectrum to extract')
    ifgroup.add_argument('-L', '--lampflat', default=None, required=True,
            help='lampflat spectrum image (for normalization)')
    ifgroup.add_argument('-T', '--traces', default=None, required=True,
            help='FITS file with trace position parameters')
    ifgroup.add_argument('-o', '--output_file', required=True,
            default=None, help='output FITS file for extracted data')
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    extgroup = parser.add_argument_group('Extraction Method')
    extgroup = extgroup.add_mutually_exclusive_group()
    extgroup.add_argument('--dumb', dest='extr_method', action='store_const',
            const='dumb', help='dumb/obvious method [def]')
    extgroup.add_argument('--fox', dest='extr_method', action='store_const',
            const='fox', help='flat-relative extraction')
    # ------------------------------------------------------------------

    context = parser.parse_args()
    #context.vlevel = context.verbose - context.quiet
    context.vlevel = 99 if context.debug else (context.verbose-context.quiet)

##--------------------------------------------------------------------------##
## Load input spectrum:
if context.spectrum:
    sys.stderr.write("Loading spectrum ... ")
    spec_data, spec_hdrs = pf.getdata(context.spectrum, header=True)
    sys.stderr.write("done.\n")

## Load input lampflat:
if context.lampflat:
    sys.stderr.write("Loading lampflat ... ")
    lamp_data, lamp_hdrs = pf.getdata(context.lampflat, header=True)
    sys.stderr.write("done.\n")

## Load input trace list:
if context.traces:
    sys.stderr.write("Loading trace list ... ")
    tdata = trio.load_traces(context.traces)
    sys.stderr.write("done.\n")

##--------------------------------------------------------------------------##
##------------------         Dimensionality Checking        ----------------##
##--------------------------------------------------------------------------##

if (spec_data.shape != lamp_data.shape):
    sys.stderr.write("Spectrum and lampflat dimensions differ:\n")
    sys.stderr.write("  --> spectrum.shape: %s\n" % str(spec_data.shape))
    sys.stderr.write("  --> lampflat.shape: %s\n" % str(lamp_data.shape))
    sys.exit(1)

##--------------------------------------------------------------------------##
##------------------         Spectrum Data Extraction       ----------------##
##--------------------------------------------------------------------------##

## Generate pixel masks from trace parameters:
sys.stderr.write("Calculating pixels from trace parameters ... ")
trace_pixel_pos = nrex.mask_from_traces(spec_data.shape, tdata)
sys.stderr.write("done.\n")

## Extract data at specified positions:
n_traces = len(trace_pixel_pos)
spec_chunks = []
fox_results = []
for ii,trace_pos in enumerate(trace_pixel_pos, 1):
    sys.stderr.write("\rExtracting blob %d of %d ... " % (ii, n_traces))
    ycoo, xcoo = trace_pos

    # Get CCD Y- and X-positions for each lambda in the blob:
    spec_rows = np.average(ycoo, axis=0) + 1.0
    spec_cols = np.average(xcoo, axis=0) + 1.0
    #spec_cols = xcoo[0] + 1.0                   # (all rows equal, take first)

    # My dumb method:
    if (context.extr_method == 'dumb'):
        spec_blob = spec_data[ycoo, xcoo]
        lamp_blob = lamp_data[ycoo, xcoo]
        lamp_vsum = np.sum(lamp_blob, axis=0)
        lamp_rflx = lamp_blob / lamp_vsum

        spec_vsum = np.sum(spec_blob, axis=0)
        spec_wsum = np.sum(spec_blob * lamp_rflx, axis=0)
        #normalize = calc_kde_mode(spec_vsum / spec_wsum, **kde_opts)
        #normalize = np.median(spec_vsum / spec_wsum)
        #spec_wsum *= normalize  # to preserve total counts
        #norm_valu = np.sum(
        #spec_wsum = np.sum(spec_blob * lamp_fwei, axis=0)
        #lamp_vals = np.sum(lamp_blob, axis=0)
        #spec_chunks.append((spec_cols, spec_rows, spec_vsum, lamp_vsum))
        spec_chunks.append((spec_cols, spec_rows, 
                            spec_vsum, spec_wsum, lamp_vsum))

    # FOX method:
    if (context.extr_method == 'fox'):
        schunk = spec_data[trace_pos]
        fchunk = lamp_data[trace_pos]
        fox_spec = nrex.flat_rel_solver(schunk, fchunk)
        fox_results.append((spec_cols, spec_rows, fox_spec))
sys.stderr.write("done.\n")


#kde_opts = {'xmin':0.0, 'xmax':40.0, 'npts':500} #'bw':0.5}
#calc_kde_mode(spec_vsum / spec_wsum, **kde_opts)
##--------------------------------------------------------------------------##
kde_opts = {'xmin':0.0, 'xmax':40.0, 'npts':500, 'bw':0.5}
fig_dims = (12, 10)
fig = plt.figure(1, figsize=fig_dims)

def qplot(blob):
    col_lims = (-5, 4100)
    xpix, ypix, sflux, wflux, lflux = blob
    kde_mode = calc_kde_mode(sflux / wflux, **kde_opts)
    wei_flux = kde_mode * wflux
    #use_flux = sflux / lflux
    use_flux = sflux
    high_flux = np.percentile(use_flux, 98.0)
    lamp_relcounts = lflux / lflux.mean()
    clean_specflux = wei_flux / lamp_relcounts
    fig.clf()

    # Raw and weighted counts from spectrum:
    ax1 = fig.add_subplot(221)
    ax1.grid(True)
    #ax1.plot(xpix, sflux, lw=1.0) 
    ax1.plot(xpix, use_flux, lw=1.50, label='spec_vsum') 
    ax1.plot(xpix, wei_flux, lw=0.25, c='r', label='spec_wsum') 
    ax1.set_xlabel('X pixel')
    ax1.set_ylabel('Counts (e-)')
    ax1.set_xlim(col_lims)
    ax1.set_ylim(-0.07*high_flux, 1.25*high_flux)
    #ax1.set_yscale('log')
    ax1.legend(loc='best')
    sys.stderr.write("Trace X range: %d to %d\n" % (xpix.min(), xpix.max()))
    sys.stderr.write("Trace Y middle: %10.5f\n" % np.average(ypix)) 

    # Corresponding lampflat counts, summed across order:
    ax2 = fig.add_subplot(222, sharex=ax1)
    ax2.grid(True)
    ax2.plot(xpix, lflux, lw=1.0, label='lamp_vsum') 
    ax2.legend(loc='best')
    ax2.set_xlabel('X pixel')
    ax2.set_ylabel('Counts (e-)')

    # Attempt to normalize using lampflat ...
    ax3 = fig.add_subplot(223, sharex=ax1)
    ax3.grid(True)
    ax3.plot(xpix, clean_specflux, lw=0.5, c='g', label='clean_specflux')
    ax3.set_xlabel('X pixel')
    ax3.set_ylabel('Counts (e-)')
    ax3.legend(loc='best')

    # Illustrate CCD position of blob:
    ax4 = fig.add_subplot(224, aspect='equal')
    ax4.grid(True)
    ax4.scatter(xpix, ypix, lw=0, s=5, label='CCD position')
    ax4.set_xlim(col_lims)
    ax4.set_ylim(col_lims)
    ax4.set_xlabel("X pixel")
    ax4.set_ylabel("Y pixel")
    ax4.legend(loc='best')

    fig.tight_layout() # adjust boundaries sensibly, matplotlib v1.1+
    plt.draw()

def _ann_image_name(tag):
    return "lookie_%s.fits" % tag

def analyze(tnum, show_spec=True, show_flat=False):
    home_dir = os.getenv('HOME')
    func_cmd = "source %s/config/functions/DS9.sh" % home_dir
    view_cmd = "fztf --lin --slock" 
    show_cmd = func_cmd + ";" + view_cmd

    # Make plot:
    qplot(spec_chunks[tnum])
    if show_spec:
        annotated = _ann_image_name("spec")
        speccheck = np.copy(spec_data).astype('float32')
        speccheck[trace_pixel_pos[tnum]] = np.nan
        qsave(annotated, speccheck)
        os.system("%s %s %s &" % (show_cmd, annotated, context.spectrum))
    if show_flat:
        annotated = _ann_image_name("lamp")
        lampcheck = np.copy(lamp_data).astype('float32')
        lampcheck[trace_pixel_pos[tnum]] = np.nan
        qsave(annotated, lampcheck)
        os.system("%s %s %s &" % (show_cmd, annotated, context.lampflat))


######################################################################
# CHANGELOG (mayhem-dumb-extract.py):
#---------------------------------------------------------------------
#
#  2018-08-03:
#     -- Increased __version__ to 0.1.0.
#     -- First created mayhem-dumb-extract.py.
#
