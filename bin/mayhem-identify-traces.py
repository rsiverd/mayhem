#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Using ThAr DOUBLE spectra, idenitfy which 'raw' traces are paired
# with one another. Save updated trace file to disk.
#
# Rob Siverd
# Created:       2018-12-26
# Last modified: 2018-12-26
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.0.5"

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
import os
import sys
import copy
import time
import signal
import numpy as np
import scipy.signal as ssig
#from numpy.lib.recfunctions import append_fields
#from functools import partial
#from collections import OrderedDict
#import multiprocessing as mp
#import pandas as pd
#import statsmodels.api as sm
#import statsmodels.formula.api as smf
#from statsmodels.regression.quantile_regression import QuantReg

## Mayhem extraction tools:
import nres_extraction
reload(nres_extraction)
nrex = nres_extraction
trio = nres_extraction.TraceIO()
frox = nres_extraction.FlatRelativeOptimalExtraction()

## Order identification tools:
import order_identification
reload(order_identification)
ads = order_identification.AdjacentDoubleSimilarity()

##--------------------------------------------------------------------------##

## Home-brew robust statistics:
#try:
#    import robust_stats
#    reload(robust_stats)
#    rs = robust_stats
#except ImportError:
#    sys.stderr.write("\nError!  robust_stats module not found!\n"
#           "Please install and try again ...\n\n")
#    sys.exit(1)

## Home-brew KDE:
#try:
#    import my_kde
#    reload(my_kde)
#    mk = my_kde
#except ImportError:
#    sys.stderr.write("\nError!  my_kde module not found!\n"
#           "Please install and try again ...\n\n")
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
halfdiv = '-' * 40
fulldiv = '-' * 80

##--------------------------------------------------------------------------##
## Catch interruption cleanly:
def signal_handler(signum, frame):
    sys.stderr.write("\nInterrupted!\n\n")
    sys.exit(1)

signal.signal(signal.SIGINT, signal_handler)

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
##------------------         Parse Command Line             ----------------##
##--------------------------------------------------------------------------##

## Parse arguments and run script:
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

## Enable raw text AND display of defaults:
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                        argparse.RawDescriptionHelpFormatter):
    pass

## Parse the command line:
if __name__ == '__main__':

    # ------------------------------------------------------------------
    descr_txt = """
    Identify corresponding traces using various methods, including:
    * cross-correlation / similarity of ThAr spectrum in adjacent traces
    * comparison of inter-trace spacings (same-order traces close together)
    
    Version: %s
    """ % __version__
    parser = MyParser(prog=os.path.basename(__file__), description=descr_txt,
                          formatter_class=argparse.RawTextHelpFormatter)
    # ------------------------------------------------------------------
    #parser.set_defaults(thing1='value1', thing2='value2')
    # ------------------------------------------------------------------
    #parser.add_argument('firstpos', help='first positional argument')
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
    ifgroup.add_argument('-D', '--double', default=None, required=True,
            help='matching DOUBLE image (ThAr Spectrum)', dest='tharfile')
    ifgroup.add_argument('-L', '--lampflat', default=None, required=True,
            help='lampflat spectrum image (for normalization)')
    ifgroup.add_argument('-T', '--traces', default=None, required=True,
            help='FITS file with trace position parameters')
    ifgroup.add_argument('-o', '--output_file', required=True,
            default=None, help='output FITS file for extracted data')
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------

    context = parser.parse_args()
    context.vlevel = 99 if context.debug else (context.verbose-context.quiet)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
## Load ThAr spectrum:
if context.tharfile:
    sys.stderr.write("Loading ThAr spectrum ... ")
    thar_data, thar_hdrs = pf.getdata(context.tharfile, header=True)
    thar_fobj = thar_hdrs['OBJECTS']
    sys.stderr.write("done.\n")

## Load input lampflat:
if context.lampflat:
    sys.stderr.write("Loading lampflat ... ")
    lamp_data, lamp_hdrs = pf.getdata(context.lampflat, header=True)
    lamp_fobj = lamp_hdrs['OBJECTS']
    sys.stderr.write("done.\n")

## Load input trace list:
if context.traces:
    sys.stderr.write("Loading trace list ... ")
    trdata = trio.load_traces(context.traces)
    sys.stderr.write("done.\n")

## Ensure corresponding channels on DOUBLE and 

##--------------------------------------------------------------------------##
##------------------         Dimensionality Checking        ----------------##
##--------------------------------------------------------------------------##

if (thar_data.shape != lamp_data.shape):
    sys.stderr.write("Spectrum and lampflat dimensions differ:\n")
    sys.stderr.write("  --> spectrum.shape: %s\n" % str(thar_data.shape))
    sys.stderr.write("  --> lampflat.shape: %s\n" % str(lamp_data.shape))
    sys.exit(1)

##--------------------------------------------------------------------------##
##------------------          ThAr DOUBLE Extraction        ----------------##
##--------------------------------------------------------------------------##

thar_data = frox.extract(thar_data, lamp_data, trdata)
thar_norm = ads.normalize_spectrum(thar_data)

##--------------------------------------------------------------------------##
##------------------      Spectrum Order Identification     ----------------##
##--------------------------------------------------------------------------##

sys.stderr.write("Comparing adjacent fibers/traces ... ")
match_summary = ads.adj_similarity_scores(thar_norm)

sys.stderr.write("resolving pairs ... ")
detected_pairs, unpaired_traces = ads.resolve_trace_pairs(match_summary)
pairs_list = [(b,a) for a,b in detected_pairs.keys()]
sys.stderr.write("done.\n")

## TODO: need a way to consistently identify a specific order as
## starting point for downstream work. This should safely exclude
## orders that are likely to be inconsistently identified by these
## procedures.

### Look for lines (potential Argon maxima):
#peak_pctg = 99.5
##thar_smoo = copy.deepcopy(thar_data)
#line_centers = []
#for odata in thar_data:
#    peak_thresh = np.percentile(odata['spec'], peak_pctg)
#    #odata['spec'] = nrex.boxcar_smooth(odata['spec'], 3)
#    maxima_idx = ssig.argrelmax(odata['spec'], order=3)[0]
#    #maxima_val = odata['spec'][maxima_idx] 
#    which_high = (odata['spec'][maxima_idx] >= peak_thresh)
#    maxima_use = maxima_idx[which_high]
#    line_centers.append(maxima_use)
#    sys.stderr.write("peak_thresh: %10.5f\n" % peak_thresh) 
#    #sys.stderr.write("maxima_idx: %s\n" % str(maxima_idx))
#    #sys.stderr.write("maxima_val: %s\n" % str(maxima_val))
#    sys.stderr.write("maxima_use: %s\n" % str(maxima_use))
#    sys.stderr.write("\n")


#fig = plt.figure()
#def illustrate(tnum):

## Deduce fibers in use:
n_update = nrex.spec_assign_fibnum(thar_fobj, thar_data, pairs_list)

n_update = nrex.traces_update_fibnum(thar_fobj, 
        trdata.get_trace_list(), pairs_list)

## Updated trace info includes fiber number/position:
trio.store_TraceData('derp.fits', trdata)

## Build a single, concatenated spectrum for fiddling:



######################################################################
# CHANGELOG (mayhem-identify-traces.py):
#---------------------------------------------------------------------
#
#  2018-12-26:
#     -- Increased __version__ to 0.0.5.
#     -- First created mayhem-identify-traces.py.
#
