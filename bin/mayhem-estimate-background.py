#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Trace-assisted background fitting for NRES spectra.
#
# Rob Siverd
# Created:       2018-06-18
# Last modified: 2018-12-24
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.1.6"

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
#import getopt
#import resource
#import signal
#import glob
#import gc
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

## Fancy background estimation (wrappers around C programs):
import fastbg
reload(fastbg)

## NRES extraction helper routines:
try:
    import nres_extraction
    reload(nres_extraction)
except ImportError:
    sys.stderr.write("\nError!  nres_extraction module not found!\n")
    sys.exit(1)
trio = nres_extraction.TraceIO()
nrex = nres_extraction

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

## ASCII I/O:
#try:
#    import astropy.io.ascii as aia
#except ImportError:
#    sys.stderr.write("\nError: astropy module not found!\n")
#    sys.exit(1)

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
# Save FITS image with clobber (astropy / pyfits):
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
## Save FITS image with clobber (fitsio):
#def qsave(iname, idata, header=None, **kwargs):
#    this_func = sys._getframe().f_code.co_name
#    sys.stderr.write("Writing to '%s' ... " % iname)
#    #if os.path.isfile(iname):
#    #    os.remove(iname)
#    fitsio.write(iname, idata, clobber=True, header=header, **kwargs)
#    sys.stderr.write("done.\n")

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

## Settings:
debug = False
timer = False
vlevel = 0
#prog_name = 'mayhem-fit-background.py'
prog_name = os.path.basename(__file__)
full_prog = sys.argv[0]
base_prog = os.path.basename(full_prog)
num_todo = 0

## hsmooth options:
#bkg_quantile = 0.45
bkg_quantile = 0.40
bkg_halfkern = 37

## Options:
save_file = None

##--------------------------------------------------------------------------##
## Argument type-checking:
def is_integer(asdf):
    try:
        int(asdf)
        return True
    except ValueError:
        return False

def is_float(asdf):
    try:
        float(asdf)
        return True
    except ValueError:
        return False

def parse_limit_string(lstring, ordered=True):
    """Returns lower,upper,valid tuple."""
    pieces = lstring.split(',')

    # Limit string must contain two elements:
    if len(pieces) != 2:
        return (None, None, False)      # exactly two elements required

    # Limits must be numbers:
    if not all([is_float(x) for x in pieces]):
        return (None, None, False)      # limits must be numbers

    # Extract:
    lower, upper = [float(x) for x in pieces]
    if ordered:
        lower, upper = (lower, upper) if (lower < upper) else (upper, lower)
    return lower, upper, True

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
    Trace-assisted background fitting for NRES spectra.
    
    Version: %s
    """ % __version__
    parser = MyParser(prog=os.path.basename(__file__), description=descr_txt)
    parser.add_argument('data_file', help='input NRES spectrum')
    parser.add_argument('-A', '--apron', required=False, default=None,
            help='override trace apron value')
    parser.add_argument('-B', '--baffle', required=False, default=None,
            help='baffle mask from FITS image')
    parser.add_argument('-Y', '--yprofile', required=False, default=None,
            help='load and use flat-field image for assistance')
    parser.add_argument('-T', '--trace', required=True,
            default=None, help='FITS table of trace parameters')
    #parser.add_argument('-d', '--dayshift', required=False, default=0,
    #        help='Switch between days (1=tom, 0=today, -1=yest', type=int)
    parser.add_argument('-q', '--quiet', action='count', default=0)
    parser.add_argument('-v', '--verbose', action='count', default=0)
    parser.add_argument('--debug', dest='debug', default=False,
             help='Enable extra debugging messages', action='store_true')
    #parser.add_argument('remainder', help='other stuff', nargs='*')
    # ------------------------------------------------------------------
    ofgroup = parser.add_argument_group('Output files (choose at least one)')
    ofgroup.add_argument('-o', '--save_bkg', required=False, default=None,
            help='Save background to filename (FITS image)')
    ofgroup.add_argument('-s', '--save_sub', required=False, default=None,
            help='Save subtracted image to filename (FITS image)')
    #ofgroup = parser.add_argument_group('Output format')
    #fmtparse = ofgroup.add_mutually_exclusive_group()
    #fmtparse.add_argument('--python', required=False, dest='output_mode',
    #        help='Return Python dictionary with results [default]',
    #        default='pydict', action='store_const', const='pydict')
    #bash_var = 'ARRAY_NAME'
    #bash_msg = 'output Bash code snippet (use with eval) to declare '
    #bash_msg += 'an associative array %s containing results' % bash_var
    #fmtparse.add_argument('--bash', required=False, default=None,
    #        help=bash_msg, dest='bash_array', metavar=bash_var)
    # ------------------------------------------------------------------

    context = parser.parse_args()
    #context.vlevel = context.verbose - context.quiet
    context.vlevel = 99 if context.debug else (context.verbose-context.quiet)

    # Some output image is required:
    if (not context.save_bkg) and (not context.save_sub):
        sys.stderr.write("Error: at least one output is required!\n")
        sys.exit(1)


### Output file is required:
#if not save_file:
#    sys.stderr.write(BRED + "\nOutput file required!" + ENDC + "\n")
#    usage(sys.stderr)
#    sys.exit(1)
#
### Input file is required:
#if (len(remainder) < 1):
#    sys.stderr.write(BRED + "\nInput file required!" + ENDC + "\n")
#    usage(sys.stderr)
#    sys.exit(1)
#
### Check for input file:
#data_file = remainder[0]
#if not os.path.isfile(data_file):
#    msg = "%s error:  file not found: %s" % (prog_name, data_file)
#    sys.stderr.write("\n" + BRED + msg + ENDC + "\n")
#    sys.exit(1)

## The clock starts now!
tstart = time.time()

##--------------------------------------------------------------------------##
## Load spectrum image:
img_vals, hdr_keys = pf.getdata(context.data_file, header=True)

##--------------------------------------------------------------------------##
## Load auxiliary files:
baff_mask = None
#prof_data = None
prof_data = np.ones_like(img_vals, dtype='float32')
#yprofile_path = "./spec_yprofile_med.fits.fz"
#baffmask_path = "./synthetic_baffle_mask_border.fits.fz"
#baffmask_path = "./vtx-baffle-lsc.b01.fits"
#yprofile_path = "/home/rsiverd/NRES/handy/spec_yprofile_med.fits.fz"
#baffmask_path = "/home/rsiverd/NRES/handy/synthetic_baffle_mask_border.fits.fz"
if os.path.isfile(context.yprofile):
    prof_data = pf.getdata(context.yprofile)
    sys.stderr.write("Successfully loaded Y-profile '%s'\n" % context.yprofile)

if os.path.isfile(context.baffle):
    baff_mask = np.bool_(pf.getdata(context.baffle))
    sys.stderr.write("Successfully loaded baffle mask '%s'\n" % context.baffle)

## Load trace data:
trace_params_list = None
if os.path.isfile(context.trace):
    trdata = trio.load_traces(context.trace)
    trace_params_list = trdata.get_trace_list()

##--------------------------------------------------------------------------##
## If requested, flatten background with Y-profile:
look_pix = (img_vals / prof_data).astype('float32')
#if isinstance(prof_data, np.ndarray):
#    look_pix = (img_vals / prof_data).astype('float32')
#else:
#    look_pix = img_vals.astype('float32')

##--------------------------------------------------------------------------##
## Reconstitute an order from its trace:
def trace2ridge(trace):
    xpixels = np.arange(trace['xmin'], trace['xmax'] + 1)
    ypixels = nrex.ridge_eval(trace['params'], xpixels)
    return (xpixels, ypixels)

#def trace2mask(trace, custom_apron=None):
#    return None


# test trace rebuild:
#trace_pixel_positions = [trace2ridge(x) for x in trace_params_list]
img_orders_masked = np.copy(look_pix)
for trace in trace_params_list:
    xpos, ypos = trace2ridge(trace)
    #apron = trace['apron']
    ymid = np.int_(np.round(ypos))
    apron = trace['apron']
    apron = 6
    for yshift in range(-apron, apron + 1):
        img_orders_masked[ymid + yshift, xpos] = np.nan

#qsave('masked_orders_pass1.fits', img_orders_masked, overwrite=True)

## Brute-force background estimate:
smoothed_bkg = fastbg.hsmooth(img_orders_masked, bkg_halfkern,
                                    bkg_quantile, autolims=True)

##--------------------------------------------------------------------------##
##------------------         Second Pass Improvement        ----------------##
##--------------------------------------------------------------------------##

## Measure excess light using initial background estimate:
excess_light = img_orders_masked - smoothed_bkg
excess_light[baff_mask] = np.nan
#qsave('derp.fits', excess_light, overwrite=True)

## Note additional bright stuff in the image:
excess_thresh = 100.0
#remaining_idx = (~np.isnan(excess_light)).nonzero()[0]
#remaining_val = excess_light[remaining_idx]
#exclude_which = (remaining_val > excess_thresh).nonzero()[0]
#exclude_these = tuple([x[exclude_which] for x in remaining_idx])
exclude_these = ((~np.isnan(excess_light)) & \
                    (excess_light > excess_thresh)).nonzero()
img_orders_masked[exclude_these] = np.nan
qsave('masked_orders_pass2.fits', img_orders_masked, overwrite=True)

## Re-smooth with better masking:
smoothed_bkg = fastbg.hsmooth(img_orders_masked, bkg_halfkern, 
                                    bkg_quantile, autolims=True)

##--------------------------------------------------------------------------##
##------------------       Save Final Background Image      ----------------##
##--------------------------------------------------------------------------##

## De-flatten the result after smoothing (if necessary):
bkg_estimate = (smoothed_bkg * prof_data).astype('float32')

## Save background estimate, if requested:
if context.save_bkg:
    qsave(context.save_bkg, bkg_estimate, overwrite=True)

## Save background-subtracted image, if requested:
if context.save_sub:
    sub_image = img_vals - bkg_estimate
    qsave(context.save_sub, sub_image, overwrite=True, header=hdr_keys)

tstop = time.time()
sys.stderr.write("Elapsed time: %.3f seconds\n" % (tstop - tstart))




######################################################################
# CHANGELOG (mayhem-fit-background.py):
#---------------------------------------------------------------------
#
#  2018-08-06:
#     -- Increased __version__ to 0.1.5.
#     -- Added -s, --save_sub option to save background-subtracted image.
#
#  2018-06-23:
#     -- Increased __version__ to 0.1.0.
#     -- Basic functionality now works as intended.
#     -- Now initialize prof_data as all 1s so that an image is always used.
#     -- Continued development. 
#
#  2018-06-18:
#     -- Increased __version__ to 0.0.1.
#     -- First created mayhem-estimate-background.py.
#
