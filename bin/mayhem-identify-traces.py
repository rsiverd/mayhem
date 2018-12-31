#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Using ThAr DOUBLE spectra, idenitfy which 'raw' traces are paired
# with one another. Save updated trace file to disk.
#
# Rob Siverd
# Created:       2018-12-26
# Last modified: 2018-12-28
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.2.0"

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
import matplotlib.pyplot as plt
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
## NRES specs:
nres_gratio = 4.0           # NRES uses R4 grating
nres_ruling_lmm = 41.59     # lines per mm ruled
#nres_spacing_lmm = 24.0442  # grating lines per mm

nres_prism_glass = "PBM2"   # glass type used in cross-dispersing prism
nres_prism_apex_deg = 55.0  # apex angle of cross-dispersing prism

nres_focallen_mm = 375.15   # approximate camera focal length
nres_center_wl_um = 0.479   # [I THINK] light wavelength nearest CCD center
nres_pix_size_mm = 0.015

useful_orders = 52.0 + np.arange(67.0)

## Spectrograph/optics brilliance:
import spectrograph_optics
reload(spectrograph_optics)
ogt = spectrograph_optics.GratingTools(nres_gratio, 
        lines_per_mm=nres_ruling_lmm)
#spec_order_list = np.arange(10, 151)
spec_order_list = np.copy(useful_orders)
spec_order_wlmid = ogt.get_blaze_wavelengths(spec_order_list, units='um')
spec_order_table = {kk:vv for kk,vv in zip(spec_order_list, spec_order_wlmid)}
for ii,ww in enumerate(spec_order_wlmid):
    sys.stderr.write("oid %3d --> %10.5f nm\n" % (ii, 1e3 * ww))
spec_order_FSR = spec_order_wlmid / spec_order_list

## Prism index of refraction for each order:
sog = spectrograph_optics.Glass()
spec_order_nn = sog.refraction_index(spec_order_wlmid, nres_prism_glass)

## Minimum deviation angle of prism at central wavelength:
## D = 2 * arcsin[n * sin(a / 2)] - a
center_wl_nn = sog.refraction_index(nres_center_wl_um, nres_prism_glass)
apex_rad = np.radians(nres_prism_apex_deg)      # NRES prism apex angle in RADIANS
incident_ang_rad = np.arcsin(center_wl_nn * np.sin(0.5 * apex_rad))
min_dev_rad = 2.0 * incident_ang_rad - apex_rad
#halfdev_rad = 0.5 * min_dev_rad

minimum_dev_r = 2.0 * np.arcsin(center_wl_nn * np.sin(0.5 * apex_rad)) - apex_rad
incident_ang_1_r = 0.5 * (minimum_dev_r + apex_rad)

## Calculate prism angular dispersion (vs wavelength):
bB_ratio = np.tan(incident_ang_rad) / center_wl_nn / np.tan(0.5 * apex_rad)
dn_dlambda = sog.glass_dn_dlambda_easy(spec_order_wlmid, nres_prism_glass)
prism_dispersion = bB_ratio * dn_dlambda
spec_order_sep_mm = nres_focallen_mm * prism_dispersion * spec_order_FSR

## Calculate order Y-position differences:
order_nn_diff = spec_order_nn - center_wl_nn
order_ypos_mm = nres_focallen_mm * 2.0 * min_dev_rad * order_nn_diff
order_ypos_pix = order_ypos_mm / nres_pix_size_mm
yspan_pixels = order_ypos_pix.max() - order_ypos_pix.min()
sys.stderr.write("Y-span (pixels): %8.2f\n" % yspan_pixels)

shifts = (np.roll(order_ypos_pix, -1) - order_ypos_pix)
zip(spec_order_list, shifts)  

## TESTING brute-force prism deflection:
incidence_1_r = incident_ang_rad * np.ones_like(spec_order_nn)
deflections_1_r = spectrograph_optics.prism_deflection_n(incidence_1_r,
                        apex_rad, spec_order_nn)

inc_change_r = deflections_1_r - min_dev_rad
incidence_2_r = incidence_1_r + inc_change_r
deflections_2_r = spectrograph_optics.prism_deflection_n(incidence_2_r,
                        apex_rad, spec_order_nn)

wdefl1 = spectrograph_optics.wiki_prism_deflection_n(incidence_1_r,
                        apex_rad, spec_order_nn)

ychange_mm = (2.0 * inc_change_r) * nres_focallen_mm
ychange_pix = ychange_mm / nres_pix_size_mm

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
            default=None, help='output FITS file for updated traces')
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
##------------------    Pair Traces and Identify Channels   ----------------##
##--------------------------------------------------------------------------##

sys.stderr.write("Comparing adjacent fibers/traces ... ")
match_summary = ads.adj_similarity_scores(thar_norm)

sys.stderr.write("resolving pairs ... ")
detected_pairs, unpaired_traces = ads.resolve_trace_pairs(match_summary)
pairs_list = [(b,a) for a,b in detected_pairs.keys()]
sys.stderr.write("done.\n")

##--------------------------------------------------------------------------##
##------------------    Store Updated d Identify Channels   ----------------##
##--------------------------------------------------------------------------##

## Add fiber/channel number to trace metadata:
n_update = nrex.traces_update_fibnum(thar_fobj, 
        trdata.get_trace_list(), pairs_list)

## TODO: identify spectroscopic orders using inter-order spacing:

## Updated trace info includes fiber number/position:
if context.output_file:
    trio.store_TraceData(context.output_file, trdata)

## -----------------------------------------------------------------------
## Brute-force Y-position of traces in central columns:

## 
all_traces = trdata.get_trace_list()
#rx, ry = trdata._ridge_from_trace(all_traces[5])

def select_xpix_range(xpix, x1, x2):
    return np.where((x1 <= xpix) & (xpix <= x2))[0]

def get_ridge_midpoints(ridges, x1, x2):
    midpoints = []
    for rx,ry in ridges:
        which = np.where((x1 <= rx) & (rx <= x2))[0]
        midpoints.append((np.average(rx[which]), np.average(ry[which])))
    return midpoints

## List of available fibers/channels:
have_fibers = list(set([x['fnum'] for x in all_traces]))
fib0_traces = [x for x in all_traces if x['fnum']==0]
fib1_traces = [x for x in all_traces if x['fnum']==1]
fib0_ridges = [trdata._ridge_from_trace(x) for x in fib0_traces]
fib1_ridges = [trdata._ridge_from_trace(x) for x in fib1_traces]

#fib0_y_central = []
#xlower, xupper = 2300, 2350
#for rx,ry in fib0_ridges:
#    which = select_xpix_range(rx, xlower, xupper)
#    ycentral = ry[which]
#    fib0_y_central.append(np.average(ycentral))

xlower, xupper = 2300, 2350
fib0_midpts = get_ridge_midpoints(fib0_ridges, xlower, xupper)
fib1_midpts = get_ridge_midpoints(fib1_ridges, xlower, xupper)
f0_xmid, f0_ymid = zip(*fib0_midpts)
f1_xmid, f1_ymid = zip(*fib1_midpts)

ydeltas = np.diff(f1_ymid)[1:]
norm_ydelta = ydeltas / ydeltas.max()

inv_blaze_wlen = 1.0 / spec_order_wlmid
norm_inv_blaze_wlen = inv_blaze_wlen / inv_blaze_wlen.max()


## -----------------------------------------------------------------------
## Compare estimate with measurement:

norm_ychange_pix = ychange_pix - np.average(ychange_pix)
norm_ychange_pix /= ychange_pix.max() - ychange_pix.min()
comp_f0_ymid_pix = np.array(f0_ymid) - np.average(f0_ymid)
comp_f0_ymid_pix /= comp_f0_ymid_pix.max() - comp_f0_ymid_pix.min()


### Load NIST data for fiddling:
#nist_data, nist_hdrs = pf.getdata('NIST_spectrum.top.fits', header=True)
#nist_data, nist_hdrs = pf.getdata('NIST_spectrum.all.fits', header=True)
#
#
### TODO: need a way to consistently identify a specific order as
### starting point for downstream work. This should safely exclude
### orders that are likely to be inconsistently identified by these
### procedures.
#
#### Look for lines (potential Argon maxima):
#peak_pctg = 99.5
#peak_pctg = 95.0
###thar_smoo = copy.deepcopy(thar_data)
#line_centers = []
#peak_fluxes  = []
#for odata in thar_data:
#    peak_thresh = np.percentile(odata['spec'], peak_pctg)
#    #odata['spec'] = nrex.boxcar_smooth(odata['spec'], 3)
#    maxima_idx = ssig.argrelmax(odata['spec'], order=3)[0]
#    #maxima_val = odata['spec'][maxima_idx] 
#    which_high = (odata['spec'][maxima_idx] >= peak_thresh)
#    maxima_use = maxima_idx[which_high]
#    line_centers.append(maxima_use)
#    peak_fluxes.append(odata['spec'][maxima_use])
#    sys.stderr.write("peak_thresh: %10.5f\n" % peak_thresh) 
#    #sys.stderr.write("maxima_idx: %s\n" % str(maxima_idx))
#    #sys.stderr.write("maxima_val: %s\n" % str(maxima_val))
#    sys.stderr.write("maxima_use: %s\n" % str(maxima_use))
#    sys.stderr.write("\n")
#intensities = np.concatenate(peak_fluxes)

## Build a single, concatenated spectrum for fiddling:

## Inspection routine:
#def qcomb(wlen, flux):
#    plt.clf()
#    for 

#SINALP  =       0.971747764900 / Sine of echelle incidence angle
#FL      =        375.146862776 / [mm] Camera focal length
#Y0      =       -22.1621828561 / [mm] y-position on CCD where gamma=0
#Z0      =    0.000267784405245 / Air (n-1) refractive index in spectrograph




######################################################################
# CHANGELOG (mayhem-identify-traces.py):
#---------------------------------------------------------------------
#
#  2018-12-28:
#     -- Increased __version__ to 0.2.0.
#     -- Basic functionality achieved!
#
#  2018-12-26:
#     -- Increased __version__ to 0.0.5.
#     -- First created mayhem-identify-traces.py.
#
