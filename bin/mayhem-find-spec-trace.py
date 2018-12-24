#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
#    Identify ridges of light from 2D spectrum image. Record for future use.
#
# Rob Siverd
# Created:       2017-10-03
# Last modified: 2018-08-06
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.4.5"

## Modules:
import argparse
from math import floor
import getopt
#import shutil
import signal
#import glob
import copy
import os
import sys
import time
#import ephem
import numpy as np
#import datetime as dt
import scipy.signal as ssig
#import scipy.ndimage as ndi
#import scipy.optimize as opti
#import scipy.interpolate as stp
#import scipy.spatial.distance as ssd
#from functools import partial
#from collections import OrderedDict
#import multiprocessing as mp
#np.set_printoptions(suppress=True)

#import pandas as pd
#import statsmodels.api as sm
#import statsmodels.formula.api as smf
#from statsmodels.regression.quantile_regression import QuantReg

#import PIL.Image as pli
#import seaborn as sns
#import theil_sen as ts
import window_filter as wf

### fastbg:
#import fastbg
#reload(fastbg)

## NRES extraction helper routines:
try:
    import nres_extraction
    reload(nres_extraction)
except ImportError:
    sys.stderr.write("\nError!  nres_extraction module not found!\n")
    sys.exit(1)
trio = nres_extraction.TraceIO()
nrex = nres_extraction

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
    sys.stderr.write("Writing to '%s' ... " % iname)
    if header:
        while (len(header) < padkeys):
            header.append() # pad header
    if os.path.isfile(iname):
        os.remove(iname)
    pf.writeto(iname, idata, header=header, **kwargs)
    sys.stderr.write("done.\n")

##--------------------------------------------------------------------------##
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

## Settings:
debug = False
timer = False
vlevel = 0
prog_name = 'rextract.py'
full_prog = sys.argv[0]
base_prog = os.path.basename(full_prog)
num_todo = 0

## Output files and options:
clobber      = True         # by default, allow output clobbering
save_traces  = None         # output FITS file for detected traces
clean_imsave = None         # if specified, write pre-trace (clean) image here
save_div_img = None         # optional output image with order divisions shown

## Input files:
baffmask_path = None        # path to baffle mask image (really good idea)
yprofile_path = None        

## Trace-finding options:
extract_apron = 7           # 8 pixels on either side of ridge
guess_norders = 66          # approx. number of spectrograph orders
trbox_xcenter = 2350        # where to place the search box for trace-finding
#trbox_xpixels = 20
#trbox_xpixels = 40
trbox_halfpix = 30          # width of search area in pixels
trbox_smethod = 'median'    # smoothing method for trace search box
trace_polyord = 2           # order of polynomial fit to trace

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Default tracing process metadata:
proc_metadata = {
        'EXTRVERS'      : [__version__, 'extraction script version'],
        'BAFFMASK'      : ['none', 'baffle mask used in extraction'],
        'YPROFILE'      : ['none', 'Y-profile nomralization image used'],
        }

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

##--------------------------------------------------------------------------##
##*********************       Parse command line:      *********************##
##--------------------------------------------------------------------------##

## Parse arguments and run script:
#if __name__ == '__main__':
#
#    descr  = """
#    Find (horizontal) spectrum traces in input image. Procedure:
#    1) pre-treat image as specified to highlight details
#    2) smooth search area along rows, find local maxima (traces)
#    3) follow light ridge until stopping condition met
#
#    Version: %s
#    """ % __version__
#
#    parser = argparse.ArgumentParser(
#            prog=os.path.basename(__file__),
#            description=descr)
#    parser.add_argument('firstpos', help='first positional argument')
#    # ------------------------------------------------------------------
#    #parser.add_argument('-s', '--site',
#    #        help='Site to retrieve data for', required=True)
#    parser.add_argument('-n', '--number_of_days', default=1,
#            help='Number of days of data to retrieve.')
#    parser.add_argument('-o', '--output_file', 
#            default='observations.csv', help='Output filename.')
#    parser.add_argument("--start", type=str, default=None, 
#            help="Start time for date range query.")
#    parser.add_argument("--end", type=str, default=None,
#            help="End time for date range query.")
#    parser.add_argument('-d', '--dayshift', required=False, default=0,
#            help='Switch between days (1=tom, 0=today, -1=yest', type=int)
#    parser.add_argument('-s', '--site', nargs=1, required=False,
#            help='Site to make URL for', choices=all_sites, default=all_sites)
#    parser.add_argument('--debug', dest='debug', default=False,
#             help='Enable extra debugging messages', action='store_true')
#    parser.add_argument('-q', '--quiet', action='count', default=0)
#    parser.add_argument('-v', '--verbose', action='count', default=0)
#    parser.add_argument('remainder', help='other stuff', nargs='*')
#    # ------------------------------------------------------------------
#    ofgroup = parser.add_argument_group('Output format')
#    fmtparse = ofgroup.add_mutually_exclusive_group()
#    fmtparse.add_argument('--python', required=False, dest='output_mode',
#            help='Return Python dictionary with results [default]',
#            default='pydict', action='store_const', const='pydict')
#    bash_var = 'ARRAY_NAME'
#    bash_msg = 'output Bash code snippet (use with eval) to declare '
#    bash_msg += 'an associative array %s containing results' % bash_var
#    fmtparse.add_argument('--bash', required=False, default=None,
#            help=bash_msg, dest='bash_array', metavar=bash_var)
#
#
#    context = parser.parse_args()
#    #context.vlevel = context.verbose - context.quiet
#    context.vlevel = 99 if context.debug else (context.verbose-context.quiet)

##--------------------------------------------------------------------------##
##*********************     Help and options menu:     *********************##
##--------------------------------------------------------------------------##

## Syntax / how to run:
def usage(stream):
    stream.write("\n"
        + "Usage: %s [options] nres_flat.fits\n" % base_prog
        + "Find (horizontal) spectrum traces in input image. Procedure:\n"
        + "1) pre-treat image as specified to highlight details\n"
        + "2) horizontally smooth search box and find local maxima (traces)\n"
        + "3) follow light ridge until stopping condition(s) met\n"
        + "\n"
        + "Version: %s\n" % __version__
        + "\n"
        + "Input files:\n"
        + "   -B, --baffle=IMAGE   load baffle mask from FITS image\n"
        + "   -Y, --yprof=IMAGE    load Y-profile (flat-field) from IMAGE\n"
        + "\n"
        + "Output files:\n"
        + "   -c, --clean=FILE     save cleaned image (pre-trace) to FILE\n"
        + "   -o, --output=FILE    save results to FILE\n"
        + "\n"
        + "Trace search options:\n"
        + "   -A, --apron=NPIX     ridge extraction box half-height "
        +                                      "[def: %d]\n" % extract_apron
        + "   -W, --hxwidth=NPIX   half-width of trace search region "
        +                                      "[def: %d]\n" % trbox_halfpix
        + "   -X, --xcenter=PIXEL  place search region center at X=PIX "
        +                                      "[def: %d]\n" % trbox_xcenter
        + "       --trmedian       trace box median smoothing [default]\n"
        + "       --traverage      trace box average smoothing (for ThAr)\n"
        #+ "       --tsmooth=FUNC   trace box smoothing method (avg or med) "
        #+                                      "[def: %s]\n" % trbox_smethod
        + "\n"
        + "Other options:\n"
        + "       --debug          extra debugging info\n"
        + "   -h, --help           print this page\n"
        + "   -q, --quiet          suppress unnecessary output\n"
        + "   -t, --timer          report program run-time\n"
        + "   -v, --verbose        more status updates\n"
        + "\n")
        #+ "   -n, --numtodo=N     stop after N iterations\n"
        #+ "   -s, --sigcut=N      clip data beyond N*sigma\n"

##--------------------------------------------------------------------------##
##*********************       Parse command line:      *********************##
##--------------------------------------------------------------------------##

## Options:
short_opts = 'B:Y:c:o:A:W:X:hqtv' # n:s:
long_opts = ['baffle=', 'yprof=', 'clean=', 'output=', 
                'apron=', 'hxwidth=', 'xcenter=',
                'debug', 'help', 'quiet', 'timer', 'verbose']
# 'numtodo=', 'sigcut='

## GNU-style parsing (with exception handling):
try:
    options, remainder = getopt.gnu_getopt(sys.argv[1:], short_opts, long_opts)
except getopt.GetoptError, err:
    sys.stderr.write("%s\n" % str(err))
    usage(sys.stderr)
    sys.exit(2)

## Handle selected options:
for opt, arg in options:
    # ------------------------------------------------
    if (opt == '--debug'):
        debug = True
        sys.stderr.write(BRED + "\nDebugging output enabled!" + ENDC + "\n")
    # ------------------------------------------------
    #elif ((opt == '-n') or (opt == '--numtodo')):
    #    if not is_integer(arg):
    #        sys.stderr.write("Error!  Non-integer argument: %s\n\n" % arg)
    #        sys.exit(1)
    #    num_todo = int(arg)
    #    if (vlevel >= 0):
    #        msg = "Stopping after %d items." % num_todo
    #        sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    # ------------------------------------------------
    elif ((opt == '-B') or (opt == '--baffle')):
        if not os.path.isfile(arg):
            msg = """
            %s error: baffle image not found:\n
            --> %s
            """ % (prog_name, arg)
            sys.stderr.write("%s\n\n" % msg)
            sys.exit(1)
        baffmask_path = arg
        if (vlevel >= 0):
            msg = "Loading baffle mask from: " + arg
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    elif ((opt == '-Y') or (opt == '--yprof')):
        if not os.path.isfile(arg):
            msg = """
            %s error: y-profile image not found:\n
            --> %s
            """ % (prog_name, arg)
            sys.stderr.write("%s\n\n" % msg)
            sys.exit(1)
        yprofile_path = arg
        if (vlevel >= 0):
            msg = "Loading Y-profile image from: " + arg
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    # ------------------------------------------------
    elif (opt == '--clobber'):
        clobber = True
        if (vlevel >= 0):
            sys.stderr.write("Output clobber enabled!\n")
    elif (opt == '--noclobber'):
        clobber = False
        if (vlevel >= 0):
            sys.stderr.write("Output clobber disabled!\n")
    # ------------------------------------------------
    elif ((opt == '-c') or (opt == '--clean')):
        if (not clobber) and os.path.isfile(arg):
            msg = """%s error: clobber=True and file exists:
            --> '%s'""" % (prog_name, arg)
            sys.stderr.write("%s\n\n" % msg)
            sys.exit(1)
        clean_imsave = arg
        if (vlevel >= 0):
            msg = "Saving pre-trace (clean) image to: " + arg
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    elif ((opt == '-o') or (opt == '--output')):
        if (not clobber) and os.path.isfile(arg):
            sys.stderr.write("Error: clobber=True and file exists:\n")
            sys.stderr.write("--> '%s'\n" % arg)
            sys.exit(1)
        save_traces = arg
        if (vlevel >= 0):
            msg = "Saving traces to: " + arg
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    # ------------------------------------------------
    elif ((opt == '-A') or (opt == '--apron')):
        if not is_integer(arg):
            sys.stderr.write("Error!  Non-integer apron: %s\n\n" % arg)
            sys.exit(1)
        extract_apron = int(arg)
        if (vlevel >= 0):
            msg = "Using extraction apron: %d pixels" % extract_apron
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    elif ((opt == '-W') or (opt == '--hxwidth')):
        if not is_integer(arg):
            sys.stderr.write("Error!  Non-integer hxwidth: %s\n\n" % arg)
            sys.exit(1)
        trbox_halfpix = int(arg)
        if (vlevel >= 0):
            msg = "Trace search box width: %d pixels" % trbox_halfpix
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    elif ((opt == '-X') or (opt == '--xcenter')):
        if not is_integer(arg):
            sys.stderr.write("Error!  Non-integer X position: %s\n\n" % arg)
            sys.exit(1)
        trbox_xcenter = int(arg)
        if (vlevel >= 0):
            msg = "Centering trace search at X=%d." % trbox_xcenter
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    elif ((opt == '--traverage')):
        trbox_smethod = 'average'
        if (vlevel >= 0):
            msg = "Using 'average' smooth in trace box!"
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    elif ((opt == '--trmedian')):
        trbox_smethod = 'median'
        if (vlevel >= 0):
            msg = "Using 'average' smooth in trace box!"
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    # ------------------------------------------------
    #elif ((opt == '-s') or (opt == '--sigcut')):
    #    if not is_float(arg):
    #        sys.stderr.write("Error!  Non-numeric argument: %s\n\n" % arg)
    #        sys.exit(1)
    #    sigcut = float(arg)
    #    if (vlevel >= 0):
    #        msg = "Using %.2f-sigma outlier threshold." % sigcut
    #        sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    # ------------------------------------------------
    
    # ------------------------------------------------
    # ------------------------------------------------
    elif ((opt == '-h') or (opt == '--help')):
        usage(sys.stdout)
        sys.exit(0)
    elif ((opt == '-q') or (opt == '--quiet')):
        vlevel -= 1
    elif ((opt == '-t') or (opt == '--timer')):
        timer = True
    elif ((opt == '-v') | (opt == '--verbose')):
        vlevel += 1
        sys.stderr.write(NYELLOW + "Increasing verbosity." + ENDC + "\n")
    # ------------------------------------------------
    else:
        msg = "Unhandled option: %s" % opt
        sys.stderr.write(BRED + "\n" + msg + ENDC + "\n\n")
        sys.exit(1)
    pass

## Verbosity:
if (vlevel >= 1):
    sys.stderr.write("%sVerbosity level: %d%s\n" % (NYELLOW, vlevel, ENDC))

## Full command line if highly verbose:
if (vlevel >= 2):
    sys.stderr.write("%s\nFull command line:%s\n" % (NCYAN, ENDC))
    sys.stderr.write("   %s\n" % sys.argv)

##--------------------------------------------------------------------------##
## Requirements check:
requirements = [(baffmask_path, 'baffle mask image'),
                #(yprofile_path, 'yprofile (BG flattener) image'),
                (save_traces, 'output traces file'),]
for var,info in requirements:
    if var == None:
        sys.stderr.write(BRED + "No %s provided!" % info + ENDC + "\n")
        usage(sys.stderr)
        sys.exit(1)

##--------------------------------------------------------------------------##
## Output file is required:
#if not save_file:
#    sys.stderr.write(BRED + "\nOutput file required!" + ENDC + "\n")
#    usage(sys.stderr)
#    sys.exit(1)

## Input file is required:
if (len(remainder) < 1):
    sys.stderr.write(BRED + "\nInput file required!" + ENDC + "\n")
    usage(sys.stderr)
    sys.exit(1)

## Check for input file:
data_file = remainder[0]
if not os.path.isfile(data_file):
    msg = "%s error:  file not found: %s" % (prog_name, data_file)
    sys.stderr.write("\n" + BRED + msg + ENDC + "\n")
    sys.exit(1)

##--------------------------------------------------------------------------##
## Annotate an image by filling in a list of sections:
def save_annotated_image(filename, idata, sections, 
        drawval=np.nan, clobber=True):
    tmp_image = np.copy(idata)
    for area in sections:
        tmp_image[area] = drawval
    qsave(filename, tmp_image, overwrite=clobber)
    return

##--------------------------------------------------------------------------##
## Load auxiliary files:
baff_mask = None
prof_data = None
#yprofile_path = "./spec_yprofile_med.fits.fz"
#baffmask_path = "./synthetic_baffle_mask_border.fits.fz"
#baffmask_path = "./vtx-baffle-lsc.b01.fits"
#yprofile_path = "/home/rsiverd/NRES/handy/spec_yprofile_med.fits.fz"
#baffmask_path = "/home/rsiverd/NRES/handy/synthetic_baffle_mask_border.fits.fz"
if yprofile_path and os.path.isfile(yprofile_path):
    prof_data = pf.getdata(yprofile_path)
    sys.stderr.write("Successfully loaded Y-profile '%s'\n" % yprofile_path)
    proc_metadata["YPROFILE"][0] = yprofile_path

if baffmask_path and os.path.isfile(baffmask_path):
    baff_mask = pf.getdata(baffmask_path)
    sys.stderr.write("Successfully loaded baffle mask '%s'\n" % baffmask_path)
    proc_metadata["BAFFMASK"][0] = baffmask_path

##--------------------------------------------------------------------------##
## Simple smoother:
def boxcar_smooth(values, box_size):
    kernel = np.ones(int(floor(box_size))) / float(box_size)
    return np.convolve(values, kernel, mode='same')

##--------------------------------------------------------------------------##
## Load spectrum image (usually a flat):
#data_file = 'image.fits'
#img_vals = pf.getdata(data_file)
#hdr_keys = pf.getheader(data_file)
img_vals, hdr_keys = pf.getdata(data_file, header=True)
proc_metadata['TR_IMAGE'] = [data_file, 'source image of traces']
proc_metadata['SRC_XPIX'] = [img_vals.shape[1], 'source image X-pixels']
proc_metadata['SRC_YPIX'] = [img_vals.shape[0], 'source image Y-pixels']

## If requested, flatten background with Y-profile:
if isinstance(prof_data, np.ndarray):
    look_pix = (img_vals / prof_data).astype('float32')
else:
    look_pix = img_vals.astype('float32')

#interior = ~np.bool_(baff_mask)
#inner_bg = np.percentile(look_pix[interior], 25.0)
#look_pix -= inner_bg

## Save pre-trace image, if requested:
if clean_imsave != None:
    qsave(clean_imsave, look_pix, overwrite=True)

#look_pix = np.copy(img_vals)

##--------------------------------------------------------------------------##
## Make/announce trace search box:

#xwidth = 20
#xwidth = 40
#xwidth = 30
#xwidth = 100

#xcenter = (float(look_pix.shape[1]) + 1.0) / 2.0
#xcenter = 2350
sys.stdout.write("trace box halfsize: %.1f\n" % trbox_halfpix)
sys.stdout.write("trace box X-center: %.1f\n" % trbox_xcenter)

trbox_xlower = int(trbox_xcenter - trbox_halfpix - 1.0)
trbox_xupper = int(trbox_xcenter + trbox_halfpix + 2.0)
trbox_columns = slice(trbox_xlower, trbox_xupper)

proc_metadata['HALFBPIX'] = [trbox_halfpix, 'half-size of tracing box']
proc_metadata['TRMETHOD'] = [trbox_smethod, 'trace box smoothing method']

##--------------------------------------------------------------------------##

## Horizontal (dispersion direction) smoothing:
if trbox_smethod == 'average':
    row_avgs = np.average(look_pix[:, trbox_columns], axis=1)
elif trbox_smethod == 'median':
    row_avgs = np.median(look_pix[:, trbox_columns], axis=1)
else:
    sys.stderr.write("FIXME: unsupported trbox_smethod: %s\n" % trbox_smethod)
    sys.exit(1)
use_avgs = boxcar_smooth(row_avgs, 7)
row_indx = np.arange(row_avgs.size)

## Normalize to background:
#n_orders = 66
#n_fibers = 2
approx_order_size = int(float(use_avgs.size) / float(guess_norders))
decent_kernel_size = 5 * approx_order_size
#running_mean = boxcar_smooth(use_avgs, decent_kernel_size)

running_25th = wf.window_filter(row_indx, row_indx, use_avgs,
        xwidth=decent_kernel_size, pctile=25)

running_50th = wf.window_filter(row_indx, row_indx, use_avgs,
        xwidth=decent_kernel_size, pctile=50)

#rel_avgs_q25 = use_avgs / running_25th

nsig = 5.0
nsig = 10.0
counts_above_bg = use_avgs - running_25th
ridge_thresh = nsig * np.std(running_25th)
clean_ridges = np.clip(counts_above_bg, ridge_thresh, None)

#maxima = ssig.argrelmax(use_avgs, axis=0, order=3)[0]
#sys.stderr.write("Found %d potential ridges.\n" % maxima.size)

#smarter = ssig.argrelmax(clean_ridges, axis=0, order=3)[0]
maxima = ssig.argrelmax(clean_ridges, axis=0, order=3)[0]
#maxima = ssig.argrelmax(counts_above_bg, axis=0, order=3)[0]
sys.stderr.write("Found %d potential ridges.\n" % maxima.size)

##--------------------------------------------------------------------------##
## Draw ridges on input image for inspection:
draw_ridges = False
ridge_apron = 3
if draw_ridges:
    tmp_image = np.copy(look_pix)
    for yshift in range(-ridge_apron, ridge_apron+1):
        tmp_image[maxima+yshift, trbox_columns] = np.nan
    qsave('ridge_inspection.fits', tmp_image, overwrite=True)    


do_extraction = False
do_extraction = True

asdf = nrex.Ridge(look_pix, bmask=baff_mask)
rlist = []
#apron = 9
apron = extract_apron
#params = {'vlevel':1, 'mincounts':25}
#params = {'vlevel':0, 'mincounts':25}
mincounts =  50
vlevel = 0
ntofit = 0
mstart = 0 #136
mdepth = 0 #25

if do_extraction:
    params = {'vlevel':vlevel, 'maxdepth':mdepth, 'mincounts':mincounts}
    for i,ypos in enumerate(maxima[mstart:]):
        rlist.append(asdf.extract(ypos, trbox_columns, apron, **params))
        if (ntofit > 0) and (i >= ntofit):
            break

##--------------------------------------------------------------------------##
## Fit ridges with polynomial:
sys.stderr.write("\n")
trace_data = []
n_ridges = len(rlist)
for i,(xlist,ylist) in enumerate(rlist, 1):
    sys.stderr.write("\rFitting order %d of %d ... " % (i, n_ridges))
    trace_fit = {'xmin':xlist.min(), 'xmax':xlist.max(), 'apron':apron}
    trace_fit['params'] = \
            nrex.fit_polynomial(xlist, ylist, trace_polyord)['params']
    #trace_fit = nrex.fit_polynomial(xlist, ylist, trace_polyord)['params']
    #keep_data = [xlist.min(), xlist.max(), 
    #xmin, xmax = xlist.min(), xlist.max()
    trace_data.append(trace_fit)
sys.stderr.write("done.\n")

## Write traces to file:
if save_traces:
    trio.store_traces(save_traces, trace_data, hdata=proc_metadata)


sys.exit(0)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##


##--------------------------------------------------------------------------##
## Group fibers according to separation:
spacings = np.diff(maxima)
ymiddles = maxima[:-1] + 0.5 * spacings
ordernum = np.arange(spacings.size)
if ordernum.size < 3:
    sys.stderr.write("Too few orders detected for grouping!\n")
    sys.exit(0)
params_1 = nrex.fit_polynomial(ordernum, spacings, poly=1)['params']
params_2 = nrex.fit_polynomial(ordernum, spacings, poly=2)['params']

## Evaluate for plotting:
line_fit = nrex.ridge_eval(params_1, ordernum)
quad_fit = nrex.ridge_eval(params_2, ordernum)

## Average spacing of orders as a function of y-pixel:
ysep_pix_model = nrex.fit_polynomial(ymiddles, spacings, poly=1)['params']
avg_ysep = nrex.ridge_eval(ysep_pix_model, ymiddles)

## Select 'above' subset:
above = (spacings > line_fit)
u_ordernum = ordernum[above]
u_spacings = spacings[above]
u_ymiddles = ymiddles[above]
uonum_params_1 = nrex.fit_polynomial(u_ordernum, u_spacings, poly=1)['params']
uonum_params_2 = nrex.fit_polynomial(u_ordernum, u_spacings, poly=2)['params']
uypix_params_2 = nrex.fit_polynomial(u_ymiddles, u_spacings, poly=2)['params']
uonum_line = nrex.ridge_eval(uonum_params_1, ordernum)
uonum_quad = nrex.ridge_eval(uonum_params_2, ordernum)
uypix_quad = nrex.ridge_eval(uypix_params_2, ymiddles)


### Quantile regression:
#useful = (ordernum > 20)
#df = pd.DataFrame({'ordernum':ordernum[useful], 'spacings':spacings[useful]})
##mod = smf.quantreg('spacings ~ ordernum', df)
#mod = smf.quantreg('spacings ~ ordernum^2', df)
#res10 = mod.fit(q=0.10)
#res50 = mod.fit(q=0.50)
#res75 = mod.fit(q=0.75)
#res90 = mod.fit(q=0.90)
#pct10_fit = nrex.ridge_eval(res10.params, ordernum)
#pct50_fit = nrex.ridge_eval(res50.params, ordernum)
#pct75_fit = nrex.ridge_eval(res75.params, ordernum)
#pct90_fit = nrex.ridge_eval(res90.params, ordernum)

##--------------------------------------------------------------------------##
## Group a list of positions according to a list of boundaries:
def organize_orders(data, boundaries):
    # blue-to-red sorted order dividers (with endcaps added):
    lo_cap, hi_cap = data.min() - 1, data.max() + 1
    use_bounds = [hi_cap] + np.sort(boundaries)[::-1].tolist() + [lo_cap]

    grouped_idxs = []
    grouped_vals = []
    for top,bot in zip(use_bounds, use_bounds[1:]):
        which = (top > data) & (data >= bot)
        grouped_vals.append(data[which])
        grouped_idxs.append(which.nonzero()[0])
    return (grouped_idxs, grouped_vals)

##--------------------------------------------------------------------------##
## Measure frequency of values in an array:
def find_repeat_offenders(values):
    uvals = np.unique(values)
    nhits = np.int_([np.sum(x == values) for x in uvals])
    hfrac = np.float_(nhits) / float(np.sum(nhits))
    order = np.int_(np.argsort(nhits)[::-1])
    return (uvals[order], hfrac[order])

##--------------------------------------------------------------------------##
## Draw separations between orders on original image:
ydiv_pix = ymiddles[(spacings > avg_ysep)]
ann_list = [(slice(y, y+2), trbox_columns) for y in np.int_(ydiv_pix)]
if save_div_img:
    save_annotated_image(save_div_img, look_pix, ann_list)

###--------------------------------------------------------------------------##
### Cross-correlate found maxima with midpoints of short separations:
#yorders = np.int_(np.cumsum(avg_ysep))
#annotations = [(slice(y, y+2), trbox_columns) for y in yorders]
#save_annotated_image('easydiv.fits', look_pix, annotations)



##--------------------------------------------------------------------------##
## Group the identified orders using separators:
midx, mgroups = organize_orders(maxima, ydiv_pix)
ridges_per_group = np.array([x.size for x in midx])
vals, frac = find_repeat_offenders(ridges_per_group)

## The number of fibers per order must be reached by large majority:
fminimum = 0.8
if (frac[0] < fminimum):
    sys.stderr.write("No consensus in order groups!\n")
    sys.stderr.write("%s\n" % str(np.column_stack((vals, frac))))
    sys.exit(0)

gsize, gfrac = vals[0], frac[0]
sys.stderr.write("Detected group size: %d (%5.1f%%)\n" % (gsize, 1e2*gfrac))

##--------------------------------------------------------------------------##
## Check whether bluemost order has all fibers:
keep_midx = copy.deepcopy(midx)
#keep_mgroups = copy.deepcopy(mgroups)
if (keep_midx[0].size != gsize):
    sys.stderr.write("Bluemost order is incomplete and will be dropped!\n")
    keep_midx.pop(0)
keep_midx = np.hstack(keep_midx)


## Exclude bogus maxima, re-sort remaining traces (blue-to-red):
good_maxima = np.sort(maxima[keep_midx])[::-1]

## Assuming completeness, eliminate stragglers on red end:
nstragglers = good_maxima.size % gsize
if (nstragglers > 0):
    sys.stderr.write("Dropping %d stragglers on red end ...\n")
    good_maxima = good_maxima[:-nstragglers]

grouped_maxima = good_maxima.reshape(-1, gsize)

group_ycenters = np.average(grouped_maxima, axis=1)

## Illustrate order midpoints:
annotations = [(slice(y, y+2), trbox_columns) for y in np.int_(group_ycenters)]
save_annotated_image('order_midpoints.fits', look_pix, annotations)


#sys.exit(0)


##--------------------------------------------------------------------------##
## Draw this ridge on a test image:
#tmp_image = np.copy(look_pix)
#for rx,ry in segments:
#    ry_low = np.int_(np.floor(ry))
#    tmp_image[ry_low+0, rx] = np.nan
#    tmp_image[ry_low+1, rx] = np.nan
#
#qsave('asdf.fits', tmp_image, overwrite=True)


##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

###--------------------------------------------------------------------------##
### Trimmer for 'drawing' ridges:
#def trim_to_image_dims(xcoo, ycoo, image):
#    ny, nx = image.shape
#    useful = (0 <= xcoo) & (xcoo < nx) & (0 <= ycoo) & (ycoo < ny)
#    return (ycoo[useful], xcoo[useful])
#
##tmp_image = np.copy(look_pix)
### Draw all ridges:
#tmp_image = np.copy(look_pix)
##n_ridges = len(rlist)
##fit_degree = 2
##for i,(xlist,ylist) in enumerate(rlist):
#for i,trace_fit in enumerate(trace_data, 1):
#    sys.stderr.write("\rPainting order %d of %d ... " % (i+1, len(trace_data)))
#    #ordfit_params = nrex.fit_polynomial(xlist, ylist, fit_degree)['params']
#    xlist = np.arange(trace_fit['xmin'], trace_fit['xmax']).astype('uint16')
#    ordfit_ycoord = nrex.ridge_eval(trace_fit['params'], xlist)
#    ylower = np.int_(np.floor(ordfit_ycoord))
#    ylower_safe = trim_to_image_dims(xlist, ylower + 0, tmp_image)
#    yupper_safe = trim_to_image_dims(xlist, ylower + 1, tmp_image)
#    tmp_image[trim_to_image_dims(xlist, ylower + 0, tmp_image)] = np.nan
#    tmp_image[trim_to_image_dims(xlist, ylower + 1, tmp_image)] = np.nan
#sys.stderr.write("done.\n")
#qsave('allthree.fits', tmp_image, overwrite=True)

##--------------------------------------------------------------------------##
## Misc:

## Convenient, percentile-based plot limits:
def nice_limits(vec, pctiles=[1,99], pad=1.2):
    ends = np.percentile(vec, pctiles)
    middle = np.average(ends)
    return (middle + pad * (ends - middle))

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Imports for plotting:
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#import matplotlib.ticker as mt
#import matplotlib._pylab_helpers as hlp
#from matplotlib.colors import LogNorm
#from matplotlib import colors
#import matplotlib.colors as mplcolors
#import matplotlib.gridspec as gridspec

## Diagnostic plots:
fig_dims = (12, 10)
fig = plt.figure(1, figsize=fig_dims)
plt.gcf().clf()
fg2 = plt.figure(2, figsize=fig_dims)
#fig, axs = plt.subplots(2, 2, sharex=True, figsize=fig_dims, num=1)
# sharex='col' | sharex='row'
#fig.frameon = False # disable figure frame drawing
#fig.subplots_adjust(left=0.07, right=0.95)
ax1 = fig.add_subplot(211)
ax1.grid(True)

ax1.plot(use_avgs, c='b', label='use_avgs')
#ax1.plot(running_mean, c='g', label='running_mean')
ax1.plot(running_25th, c='g', label='running_25th')
ax1.plot(running_50th, c='m', label='running_50th')

for xpos in maxima:
    ax1.axvline(xpos, c='r')

ax1.legend(loc='upper right')

ax2 = fig.add_subplot(212)
ax2.grid(True)
ax2.plot(counts_above_bg, c='b', label='shazam')
ax2.plot(clean_ridges, c='g', label='cleaned')
ax2.legend(loc='best')

## Disable axis offsets:
#ax1.xaxis.get_major_formatter().set_useOffset(False)
#ax1.yaxis.get_major_formatter().set_useOffset(False)

#ax1.plot(kde_pnts, kde_vals)

#blurb = "some text"
#ax1.text(0.5, 0.5, blurb, transform=ax1.transAxes)
#ax1.text(0.5, 0.5, blurb, transform=ax1.transAxes,
#      va='top', ha='left', bbox=dict(facecolor='white', pad=10.0))

#colors = cm.rainbow(np.linspace(0, 1, len(plot_list)))
#for camid, c in zip(plot_list, colors):
#    cam_data = subsets[camid]
#    xvalue = cam_data['CCDATEMP']
#    yvalue = cam_data['PIX_MED']
#    yvalue = cam_data['IMEAN']
#    ax1.scatter(xvalue, yvalue, color=c, lw=0, label=camid)

#mtickpos = [2,5,7]
#ndecades = 1.0   # for symlog, set width of linear portion in units of dex
#nonposx='mask' | nonposx='clip' | nonposy='mask' | nonposy='clip'
#ax1.set_xscale('log', basex=10, nonposx='mask', subsx=mtickpos)
#ax1.set_xscale('log', nonposx='clip', subsx=[3])
#ax1.set_yscale('symlog', basey=10, linthreshy=0.1, linscaley=ndecades)
#ax1.xaxis.set_major_formatter(formatter) # re-format x ticks
#ax1.set_ylim(ax1.get_ylim()[::-1])
#ax1.set_xlabel('whatever', labelpad=30)  # push X label down 

#ax1.set_xticks([1.0, 3.0, 10.0, 30.0, 100.0])
#ax1.set_xticks([1, 2, 3], ['Jan', 'Feb', 'Mar'])

#ax1.set_xlim(nice_limits(xvec, pctiles=[1,99], pad=1.2))
#ax1.set_ylim(nice_limits(yvec, pctiles=[1,99], pad=1.2))

#spts = ax1.scatter(x, y, lw=0, s=5)
#cbar = fig.colorbar(spts, orientation='vertical')
#cbar.formatter.set_useOffset(False)
#cbar.update_ticks()

fig.tight_layout() # adjust boundaries sensibly, matplotlib v1.1+
plt.draw()

##--------------------------------------------------------------------------##
## Investigate order spacing:

fg2.clf()
ax3 = fg2.add_subplot(211)
ax3.grid(True)

#ax3.scatter(ordernum, spacings, c='b', s=3, lw=0)
ax3.plot(ordernum, spacings, c='b', label='spacings', lw=0.5)
ax3.plot(ordernum, line_fit, c='g', label='linear fit')
ax3.plot(ordernum, quad_fit, c='r', label='quadratic fit')

ax3.plot(ordernum, uonum_line, c='m', label='upper line')
ax3.plot(ordernum, uonum_quad, c='k', label='upper quad')

#ax3.plot(ordernum, pct10_fit, c='k', label='10th percentile')
#ax3.plot(ordernum, pct50_fit, c='c', label='50th percentile')
#ax3.plot(ordernum, pct75_fit, c='y', label='75th percentile')
#ax3.plot(ordernum, pct90_fit, c='y', label='90th percentile')
ax3.set_xlim(nice_limits(ordernum, pctiles=[1,99], pad=1.2))
ax3.set_ylim(nice_limits(spacings, pctiles=[1,99], pad=1.2))

ax3.legend(loc='best')

ax4 = fg2.add_subplot(212)
ax4.grid(True)

ax4.plot(ymiddles, spacings, c='b', label='spacings', lw=0.5)
#ax4.plot(ymiddles, line_fit, c='g', label='linear fit')
#ax4.plot(ymiddles, quad_fit, c='r', label='quadratic fit')
ax4.plot(ymiddles, avg_ysep, c='r', label='average ysep')

#ax4.plot(ymiddles, uypix_quad, c='k', label='upper quad')
#ax4.plot(ymiddles, uypix_quad, c='k', label='upper quad')
ax4.legend(loc='best')

fg2.tight_layout()
plt.draw()



######################################################################
# CHANGELOG (mayhem-find-spec-trace.py):
#---------------------------------------------------------------------
#
#  2018-08-06:
#     -- Increased __version__ to 0.4.5.
#     -- Decreased default extraction apron size from 9 to 7 pixels. The
#           slimmer box prevents jumping between channels in the redder
#           channels of NRES spectra (where ridges get close together).
#     -- Added -A, --apron=NPIX option to manually specify apron size
#           for ridge extraction.
#
#  2018-08-03:
#     -- Increased __version__ to 0.4.3.
#     -- Fixed error raised when no Y-profile provided.
#
#  2018-06-19:
#     -- Increased __version__ to 0.4.2.
#     -- Now only save 'divided.fits' equivalent if output is specified. Still
#           need to implement a command line option for this.
#
#  2018-06-18:
#     -- Increased __version__ to 0.4.1.
#     -- Removed 'fastbg' import (not used, moved to separate script).
#
#  2018-05-07:
#     -- Increased __version__ to 0.4.0.
#     -- Replaced 'clobber' with 'overwrite' in astropy writeto() calls.
#     -- Replaced numerous instances of 'xcolumns' with trbox_columns.
#     -- Implemented -Y, --yprof=IMAGE command line option to specify
#           image to use for background normalization.
#
#  2018-04-15:
#     -- Increased __version__ to 0.3.5.
#     -- Moved matplotlib import to end of script (can be avoided if no plots).
#     -- Command line parsing updates:
#           * added global 'clobber' setting
#           * added CL options --clobber and --noclobber
#           * added -c, --clean=FILE to save cleaned image (respects clobber)
#           * added -W, --hxwidth=NPIX option to specify trbox half-size
#           * added -X, --xcenter=PIXEL option to specify trbox X-center
#
#  2018-04-14:
#     -- Increased __version__ to 0.3.1.
#     -- Renamed script from rextract to mayhem-find-spec-trace, moved to
#           new location.
#
#  2017-11-01:
#     -- Increased __version__ to 0.3.0.
#     -- Tons of updates ... this morphed into a feature-tester.
#
#  2017-09-10:
#     -- Increased __version__ to 0.1.0.
#     -- First created spread-from-center.py.
#
