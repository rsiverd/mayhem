#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Using ThAr DOUBLE spectra, idenitfy which 'raw' traces are paired
# with one another. Save updated trace file to disk.
#
# Rob Siverd
# Created:       2018-12-26
# Last modified: 2019-02-14
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
import argparse
import os
import sys
import copy
import time
import signal
import numpy as np
import scipy.signal as ssig
import scipy.optimize as opti
import matplotlib.pyplot as plt
from numpy.lib.recfunctions import append_fields
from functools import partial
#from collections import OrderedDict
#import multiprocessing as mp
#import pandas as pd
#import statsmodels.api as sm
#import statsmodels.formula.api as smf
#from statsmodels.regression.quantile_regression import QuantReg
import itertools as itt

## Rotation matries and more:
import fov_rotation
reload(fov_rotation)
r3d = fov_rotation.Rotate3D()

## Theil-Sen fitting:
import theil_sen as ts

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
nres_grating_tilt_deg = 13.786  # grating tilt w.r.t. optical bench surface
nres_alpha_angle_rad = np.radians(90.0 - nres_grating_tilt_deg)

nres_prism_glass = "PBM2"   # glass type used in cross-dispersing prism
nres_prism_apex_deg = 55.0  # apex angle of cross-dispersing prism

bluemost_order = 119        # spectroscopic order of 'upper' order

#nres_focallen_mm = 375.15   # approximate camera focal length
nres_focallen_mm = 400.00   # TESTING
nres_focallen_mm = 390.00   # TESTING
nres_focallen_mm = 385.00   # TESTING
nres_focallen_mm = 380.00   # TESTING
nres_center_wl_um = 0.479   # [I THINK] light wavelength nearest CCD center
nres_pix_size_mm = 0.015

useful_orders = 52.0 + np.arange(67.0)
useful_orders = 51.0 + np.arange(69.0)
#useful_orders = 54.0 + np.arange(67.0)

## Spectrograph/optics brilliance:
import spectrograph_optics
reload(spectrograph_optics)
nrp = spectrograph_optics.Prism(nres_prism_glass, nres_prism_apex_deg)
ogt = spectrograph_optics.GratingTools(nres_gratio, 
        lines_per_mm=nres_ruling_lmm)
dppgp = spectrograph_optics.DoublePassPrismGratingPrism()

#spec_order_list = np.arange(10, 151)
spec_order_list = np.copy(useful_orders)
#spec_order_wlmid = ogt.get_blaze_wavelengths(spec_order_list, units='um')
spec_order_wlmid, spec_order_FSR, spec_order_angsize = \
        ogt.get_order_params(spec_order_list, units='um')
spec_order_table = {kk:vv for kk,vv in zip(spec_order_list, spec_order_wlmid)}
for ii,ww in enumerate(spec_order_wlmid):
    sys.stderr.write("oid %3d --> %10.5f nm\n" % (ii, 1e3 * ww))
#spec_order_FSR = spec_order_wlmid / spec_order_list

## Prism index of refraction for each order:
#sog = spectrograph_optics.Glass()
sog = spectrograph_optics.Glass(nres_prism_glass)
spec_order_nn = sog.refraction_index(spec_order_wlmid) #, nres_prism_glass)

## Minimum deviation angle of prism at central wavelength:
## D = 2 * arcsin[n * sin(a / 2)] - a
center_wl_nn = sog.refraction_index(nres_center_wl_um) #, nres_prism_glass)
apex_rad = np.radians(nres_prism_apex_deg)      # NRES prism apex angle in RADIANS
incident_ang_rad = np.arcsin(center_wl_nn * np.sin(0.5 * apex_rad))
min_dev_rad = 2.0 * incident_ang_rad - apex_rad
#halfdev_rad = 0.5 * min_dev_rad

minimum_dev_r = 2.0 * np.arcsin(center_wl_nn * np.sin(0.5 * apex_rad)) - apex_rad
incident_ang_1_r = 0.5 * (minimum_dev_r + apex_rad)

## TESTING brute-force prism deflection:
incidence_1_r = incident_ang_rad * np.ones_like(spec_order_nn)
deflections_1_r = nrp.deflection_rad_n2(incidence_1_r, spec_order_nn**2)

inc_change_r = deflections_1_r - min_dev_rad
incidence_2_r = incidence_1_r + inc_change_r
deflections_2_r = nrp.deflection_rad_n2(incidence_2_r, spec_order_nn**2)

ychange_mm = (2.0 * inc_change_r) * nres_focallen_mm
ychange_pix = ychange_mm / nres_pix_size_mm


##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Solve for ACTUAL undeflected wavelength ...
## --> from CAD, prism is inclined by 51.007 deg w.r.t. camera face
## --> from CAD, deflection direct to grating is 44.827 degrees

all_wavelengths = np.linspace(0.375, 0.925, 1000)
#all_indx_refrac = sog.refraction_index(all_wavelengths)
cad_incidence_1_r = np.radians(51.007) * np.ones_like(all_wavelengths)
cad_deflections_r = nrp.deflection_rad_wl(cad_incidence_1_r, all_wavelengths)
def deflect_resid(guess_lam_um):
    _cad_incid_r = np.radians(51.007)
    _cad_gturn_r = np.radians(44.827)
    this_deflect_r = nrp.deflection_rad_wl(_cad_incid_r, guess_lam_um)
    return this_deflect_r - _cad_gturn_r

## Solve for central wavelength via bisection:
answer = opti.bisect(deflect_resid, 0.3, 0.5)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Reproduce Tim's wavelength model:
y0_mm = -22.1621828561      # CCD Y-coordinate where gamma angle is 0.0

## ADDITIONAL PARAMETERS:
# * define alpha angle (incidence onto grating). This should be the set by the
# orientation of the grating w.r.t. optical bench and (I think) is the same for
# all wavelengths. Set once. Could also define this as facet angle.
# * Need to know/set the gamma angle that corresponds to the minimally
# deflected wavelength through the prism. This is effectively the "rotation" of
# the grating on the optical bench surface away from The gamma angles of individual
# wavelengths reaching grating are then (after one pass thru prism):
#       -- gamma_mindef + inc_change_r

##--------------------------------------------------------------------------##

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

## Separate the ThAr channels too:
f0_thar_data = [y for x,y in zip(all_traces, thar_norm) if x['fnum']==0]
f1_thar_data = [y for x,y in zip(all_traces, thar_norm) if x['fnum']==1]
thar_specord = bluemost_order - np.arange(len(f0_thar_data))[::-1]

#def nearest_xy(x1, y1, x2, y2, xmax=110):
#    tdiff = []
#    for tx,ty in zip(x1, y1):
#        which = (np.abs(tx - x2) <= xmax)
#        if (np.sum(which) > 0):
#            tdiff.append(np.min(np.hypot(tx - x2[which], ty - y2[which])))
#    return np.min(tdiff)
#
### Closest separations:
#def get_minseps(ridges):
#    ndiffs = len(ridges) - 1
#    separations = []
#    for ii in range(ndiffs):
#        x1, y1 = ridges[ii]
#        x2, y2 = ridges[ii+1]
#        tsep = [np.min(np.hypot(x2 - tx, y2 - ty)) for tx,ty in zip(x1,y1)]
#        separations.append(np.min(tsep))
#    return np.array(separations)
#
### Closest separations:
#def get_minseps2(ridges):
#    ndiffs = len(ridges) - 1
#    separations = []
#    for ii in range(ndiffs):
#        x1, y1 = ridges[ii]
#        x2, y2 = ridges[ii+1]
#        sys.stderr.write("matching orders %d,%d\n" % (ii, ii+1))
#        separations.append(nearest_xy(x1, y1, x2, y2))
#    return np.array(separations)
#
##sys.stderr.write("Timing approach 1 ... ")
##tik = time.time()
##yay_seps_1 = get_minseps(fib0_ridges)
##tok = time.time()
##sys.stderr.write("took %.3f seconds.\n" % (tok-tik))
## TAKES 20 SECONDS
#
#sys.stderr.write("Timing approach 2 ... ")
#tik = time.time()
#yay_seps_2 = get_minseps2(fib0_ridges)
#tok = time.time()
#sys.stderr.write("took %.3f seconds.\n" % (tok-tik))
## TAKES 14 SECONDS

## -----------------------------------------------------------------------
## -----------------------------------------------------------------------
## -----------------------------------------------------------------------
## Expected order angular size vs size in pixels:

fib0_xwidth = np.float_([(tt['xmax'] - tt['xmin']) for tt in fib0_traces])
fib1_xwidth = np.float_([(tt['xmax'] - tt['xmin']) for tt in fib1_traces])
tmp_f0_xpix = fib0_xwidth[2:]
tmp_f1_xpix = fib1_xwidth[2:]

full_fib0_xwidth = fib0_xwidth[(fib0_xwidth <= 4090)]
full_fib1_xwidth = fib1_xwidth[(fib1_xwidth <= 4090)]
full_fib0_consec = full_fib0_xwidth / np.roll(full_fib0_xwidth, 1)
full_fib1_consec = full_fib1_xwidth / np.roll(full_fib1_xwidth, 1)


theory_consec = spec_order_angsize / np.roll(spec_order_angsize, 1)

## -----------------------------------------------------------------------
## -----------------------------------------------------------------------
## -----------------------------------------------------------------------

#fib0_y_central = []
#xlower, xupper = 2300, 2350
#for rx,ry in fib0_ridges:
#    which = select_xpix_range(rx, xlower, xupper)
#    ycentral = ry[which]
#    fib0_y_central.append(np.average(ycentral))

xlower, xupper = 2300, 2350
fib0_midpts = get_ridge_midpoints(fib0_ridges, xlower, xupper)
fib1_midpts = get_ridge_midpoints(fib1_ridges, xlower, xupper)
f0_xmid, f0_ymid = np.array(fib0_midpts).T
f1_xmid, f1_ymid = np.array(fib1_midpts).T

ydeltas = np.diff(f1_ymid)[1:]
norm_ydelta = ydeltas / ydeltas.max()

#inv_blaze_wlen = 1.0 / spec_order_wlmid
#norm_inv_blaze_wlen = inv_blaze_wlen / inv_blaze_wlen.max()

## Scale to match data:
#shift, scale = ts.linefit(ychange_pix, np.array(f0_ymid))


## -----------------------------------------------------------------------
## Compare estimate with measurement:

norm_ychange_pix = ychange_pix - np.min(ychange_pix)
norm_ychange_pix /= ychange_pix.max() - ychange_pix.min()
comp_f0_ymid_pix = np.array(f0_ymid) - np.min(f0_ymid)
comp_f0_ymid_pix /= comp_f0_ymid_pix.max() - comp_f0_ymid_pix.min()

#norm_ychange_rng = norm_ychange_pix.max() - norm_ychange_pix.min()
#comp_f0_ymid_rng = comp_f0_ymid_pix.max() - comp_f0_ymid_pix.min()

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

## -----------------------------------------------------------------------
## -----------------------------------------------------------------------
## -----------------------------------------------------------------------
## -----------------------------------------------------------------------

nres_nominal_gamma = 0.0    # [radians] gamma angle for nres_center_wl_um
nres_prism_apex_deg = 55.0  # apex angle of cross-dispersing prism
nres_prism_apex_rad = np.radians(nres_prism_apex_deg)
nres_spacing_um = 1e3 / nres_ruling_lmm     # grating spacing in microns


### -----------------------------------------------------------------------
### Effective gamma angle is the nominal gamma angle plus the changes in prism
### deflection angle due to difference from nominal wavelength.
#def gamma_eff(gamma_nom, wlen_um):
#    # delta_gamma = deflection_rad - min_dev_rad
#    #deflect_r = spectrograph_optics.prism_deflection_n(incident_ang_rad,
#    #        nres_prism_apex_rad, sog.refraction_index(wlen_um))
#    deflect_r = nrp.deflection_rad_wl(incident_ang_rad, wlen_um)
#    return gamma_nom + deflect_r - min_dev_rad

#use_gamma_eff = partial(gamma_eff, nres_nominal_gamma)

## From existing solutions ...
nres_sine_alpha = 0.971747764900
blaze_r = np.arctan(4.0)
#facet_r = np.radians(0.5)
facet_r = np.arcsin(nres_sine_alpha) - blaze_r
#fixed_geometry = 2.0 * nres_spacing_um * np.cos(facet_r) * np.sin(blaze_r)

#def lamcen_residual(wlen, order=0):
#    ls = wlen * order
#    rs = fixed_geometry * np.cos(use_gamma_eff(wlen))
#    return ls - rs
#
#def iter_calc_lamcen(order):
#    kw = {'order':order}
#    runme = partial(lamcen_residual, **kw)
#    return opti.bisect(runme, 0.0, 10.0)

## FIBER CHOICE:
fib_which = 0

##-----------------------------------------------------------------------
## Quick test of central wavelength code:
lam_cen_dppgp = dppgp.calc_central_wlen_um(spec_order_list)

ctr_wlen, ctr_gamma = dppgp.fancy_deflections(spec_order_list)

ctr_headings, pg_yshifts, pc_yshifts = \
        np.array([dppgp.two_pass_deflection(x) for x in ctr_wlen]).T
dp_yshifts_mm = pg_yshifts + pc_yshifts
dp_yshifts_pix = dp_yshifts_mm / nres_pix_size_mm
#dp_yshifts_range = dp_yshifts_pix.max() - dp_yshifts_pix.min()
#normed_dp_yshifts = (dp_yshifts_pix - dp_yshifts_pix.min()) / dp_yshifts_range

## Compute corresponding line tilts (TESTING):
ctr_line_tilts = ogt._calc_line_tilt_ctr(blaze_r, ctr_gamma)
ctr_tilts_deg  = np.degrees(ctr_line_tilts)

## -----------------------------------------------------------------------
## Transform CCD -> spectrograph coordinates:
def ccd2spec_xy(ccdx, ccdy, rot_deg, xnudge=0, ynudge=0,
        xcenter=0, ycenter=0):
    if ccdx.shape != ccdy.shape:
        raise ValueError("CCD coordinate arrays have mismatched shape\n:"
                + "%s != %s\n" % (str(ccdx.shape), str(ccdy.shape)))
    if len(ccdx.shape) != 1:
        raise ValueError("Expected 1-D input, have shape %s" % str(ccdx.shape))
    # MEMORY INTENSIVE!!
    old_dim = ccdx.shape
    #ccd_xyz = np.vstack((ccdx.flatten(), 
    #                     ccdy.flatten(),
    #                     np.zeros(ccdx.size)))
    ccd_xyz = np.vstack((ccdx - xcenter, ccdy - ycenter, np.zeros(ccdx.size)))
    sx, sy, _ = r3d.zrot(np.radians(rot_deg), ccd_xyz)
    return sx.A1 + xcenter + xnudge, sy.A1 + ycenter + ynudge
    #return np.squeeze(np.asarray(sx)), np.squeeze(np.asarray(sy))
    #return np.array(sx), np.array(sy)
    #return sx.reshape(old_dim), sy.reshape(old_dim)
    #return np.array(sx).reshape(old_dim), np.array(sy).reshape(old_dim)

## -----------------------------------------------------------------------
## Coordinate rotation time!
spec_rotation = 13.091

## -----------------------------------------------------------------------
## Initial crack at wavelength solution:
nres_focallen_mm = 391.0
rlist = fib0_ridges if fib_which==0 else fib1_ridges
xpix_beta_c = 2048.5        # X-pixel where beta=beta_c
xpix_beta_c =    0.0        # X-pixel where beta=beta_c
xpix_beta_c = 4080.0        # X-pixel where beta=beta_c
xpix_beta_c = 2100.0        # X-pixel where beta=beta_c
xpix_beta_c = 2000.0        # X-pixel where beta=beta_c
xpix_beta_c = 2100.0        # X-pixel where beta=beta_c
xpix_beta_c = 2050.0        # X-pixel where beta=beta_c
xpix_beta_c = 2010.0        # X-pixel where beta=beta_c
xpix_beta_c = 2150.0        # X-pixel where beta=beta_c
xpix_beta_c = 2200.0        # X-pixel where beta=beta_c
xpix_beta_c = 2300.0        # X-pixel where beta=beta_c
xpix_beta_c = 2615.0        # X-pixel where beta=beta_c
wavelengths = {}
wavelength2 = {}
some_xpix = []
some_mmsx = []
#xpix_bet2_c = 1850.0
xpix_bet2_c = 1800.0
#xpix_bet2_c = 2000.0
for ii,spord in enumerate(spec_order_list):
    sys.stderr.write("\rOrder %3d ... " % spord)
    sys.stderr.write("\n")
    rx, ry = rlist[ii]
    center_wl = ctr_wlen[ii]
    cos_gamma = np.cos(ctr_gamma[ii])
    #mmrx = (rx - xpix_beta_c) * nres_pix_size_mm
    mmrx = (xpix_beta_c - rx) * nres_pix_size_mm
    #beta = np.arcsin(nres_sine_alpha) - np.arcsin(mmrx / nres_focallen_mm)
    beta = np.arcsin(nres_sine_alpha) - np.arctan(mmrx / nres_focallen_mm)
    sys.stderr.write("--> beta min,max: %8.5f, %8.5f\n" % 
            (np.degrees(np.min(beta)), np.degrees(np.max(beta))))
    #sine_beta = np.sin(np.arctan(mmrx / nres_focallen_mm))
    #center_wl = spec_order_wlmid[ii]
    #center_wl = iter_calc_lamcen(spord)
    #cos_gamma = np.cos(use_gamma_eff(center_wl))
    tlam = nres_spacing_um / float(spord) * cos_gamma \
            * (nres_sine_alpha + np.sin(beta))
    wavelengths[int(spord)] = tlam
    #wavelengths.append(tlam)

    sxx, syy = ccd2spec_xy(rx, ry, spec_rotation, xnudge=-100)
    #sxx, syy = ccd2spec_xy(rx, ry, spec_rotation, 
    #        xcenter=2048.5, ycenter=2048.5, xnudge=-100)
    mmsx = (xpix_bet2_c - sxx) * nres_pix_size_mm
    bet2 = np.arcsin(nres_sine_alpha) - np.arctan(mmsx / nres_focallen_mm)
    sys.stderr.write("--> bet2 min,max: %8.5f, %8.5f\n" % 
            (np.degrees(np.min(bet2)), np.degrees(np.max(bet2))))
    slam = nres_spacing_um / float(spord) * cos_gamma \
            * (nres_sine_alpha + np.sin(bet2))
    wavelength2[int(spord)] = slam
    #wavelength2.append(slam)
    some_xpix.append(sxx)
    some_mmsx.append(mmsx)
sys.stderr.write("done.\n")

mean_xpix = [x.mean() for x in some_xpix]
lower, upper = np.min(mean_xpix), np.max(mean_xpix)
spread = upper - lower
sys.stderr.write("lower:  %9.3f\n" % lower)
sys.stderr.write("upper:  %9.3f\n" % upper)
sys.stderr.write("spread: %9.3f\n" % spread)

##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
##-----------------------------------------------------------------------
## Some plotting ...

def shift_normalizer(ypos):
    yrange = ypos.max() - ypos.min()
    ynudge = ypos - ypos.min()
    return ynudge / yrange

fig = plt.figure(8, figsize=(12,7))
fig.clf()
ax1 = fig.add_subplot(121)
ax1.grid(True)
ax1.plot(shift_normalizer(ychange_mm), label='ychange_mm')
ax1.plot(shift_normalizer(dp_yshifts_mm), label='dp_yshifts_mm')
ax1.plot(shift_normalizer(ydeltas), label='YDELTAS')
ax1.legend(loc='upper left')
ax2 = fig.add_subplot(122)
ax2.grid(True)
ax2.plot(shift_normalizer(pg_yshifts), label='prism-GRATING shifts')
ax2.plot(shift_normalizer(pc_yshifts), label='prism-CAMERA  shifts')
ax2.plot(shift_normalizer(ydeltas), label='DATA')
ax2.legend(loc='best')
fig.tight_layout()
plt.draw()

## ----------------------------------------------------------------------- ##
## ----------------------------------------------------------------------- ##
## Wavelength references:
import wavelength_reference
reload(wavelength_reference)
#wlr = wavelength_reference
wlf = wavelength_reference.WLFetcher()

## Wavelength fitting helpers:
import wl_solve_test
reload(wl_solve_test)
tlf = wl_solve_test.LineFinder()
afsr = wl_solve_test.ApproxFSR()
wlsm = wl_solve_test.WLSegMatch()

## Compute line positions for every order:
sys.stderr.write("Computing line positions for fib_which=%d ...\n" % fib_which)
corresponding_thar = f0_thar_data if fib_which==0 else f1_thar_data
measured_line_xpix = []
slow_args = {'pctile':True, 'shallow':0.01}
#fast_args = {'pctile':False, 'shallow':0.1}
fast_args = {'pctile':False, 'shallow':0.01}
for i,tdata in enumerate(corresponding_thar, 1):
    sys.stderr.write("\rScanning order %d of %d ... "
            % (i, len(corresponding_thar)))
    linepos = tlf.extract_lines_xpix(tdata['xpix'], tdata['spec'], **fast_args)
    measured_line_xpix.append(linepos)
sys.stderr.write("done.\n")

## Approximate (greedy) wavelength limits for the specified order:
spec_order_wl_lims = []
spec_order_line_sets = []
for sord,ctrwl_nm in zip(spec_order_list, 1e3 * ctr_wlen):
    wl_lims_nm = afsr.greedy_wl_limits(ctrwl_nm, sord, nFSR=1.5)
    comb_lines = wlf.get_combined_lines(*wl_lims_nm)
    spec_order_wl_lims.append(wl_lims_nm)
    spec_order_line_sets.append(comb_lines)

## ----------------------------------------------------------------------- ##
## ----------------------------------------------------------------------- ##
## Segment matching for lines?

#def list_segments(xcoords):
#    pairs = itt.combinations(range(len(xcoords)), 2)
#    ii,jj = zip(*pairs)
#    x_sep = np.log10(xcoords[jj,] - xcoords[ii,])
#    params = np.column_stack((xcoords[ii], x_sep))
#    indices = np.column_stack((ii, jj)).astype('uint16')
#    return {'idx':indices, 'seg':params}
#    #return np.log10(diffs)
#

## ----------------------------------------------------------------------- ##
## ----------------------------------------------------------------------- ##
## Brute force comparison ...
soi = lambda x: int(spec_order_list[x])
wl_refpoints_nm = {}

## 11th order (spec order 62) has two booming lines in it:
#tord = int(spec_order_list[11])
wl_refpoints_nm[soi( 2)] = np.array([871.162590, 875.043326, 876.064871, 
                877.798276, 884.361065, 887.126858])
wl_refpoints_nm[soi( 3)] = np.array([857.547551, 866.786549, 
                                            867.032496, 871.162590])
wl_refpoints_nm[soi( 4)] = np.array([841.052601, 841.904021, 842.031154,
            842.353957, 842.697974, 844.780782,
            844.883236, 848.068736, 848.086191, 852.378507])



wl_refpoints_nm[soi(11)] = np.array([750.59341792, 751.67241877])
wl_refpoints_nm[soi(12)] = np.array([738.60150497])
wl_refpoints_nm[soi(13)] = np.array([727.49377621])
wl_refpoints_nm[soi(14)] = np.array([717.0870892])
wl_refpoints_nm[soi(15)] = np.array([706.9167041])
wl_refpoints_nm[soi(16)] = np.array([696.735506])
wl_refpoints_nm[soi(20)] = np.array([651.416368, 653.314665, 655.597097,
                657.903089, 658.572414, 659.035969, 659.330542, 659.575998])
wl_refpoints_nm[soi(21)] = np.array([641.367089, 641.538735, 645.906730,
                646.439855, 649.253065, 651.416368])



## ----------------------------------------------------------------------- ##
## ----------------------------------------------------------------------- ##

## Quick Gaussian evaluator:
def eval_gauss(xgrid, mu, sigma, height=1.0):
    zz = (xgrid - mu) / sigma
    ggrid = height * np.exp(-0.5 * zz * zz)
    return ggrid

def grid_tally_gauss(xgrid, centers, gsigma):
    return np.sum([eval_gauss(xgrid, x, gsigma) for x in centers], axis=0)

## Cross-reference X-coordinates (determine degree of overlap) by
## replacing lines with profiles of some width and multiplying.
def crude_crossref(xvals1, xvals2, gsigma, pad=0.05, sfactor=1.0):
    xmin = min(xvals1.min(), xvals2.min())
    xmax = max(xvals1.max(), xvals2.max())
    x_lo = np.floor(xmin - pad * (xmax - xmin))
    x_hi = np.ceil(xmax + pad * (xmax - xmin))
    npix = x_hi - x_lo + 1
    xgrid = np.arange(npix) + x_lo
    sys.stderr.write("xmin,xmax: %8.2f,%8.2f\n" % (xmin, xmax))
    #sys.stderr.write("x_lo: %8.2f\n" % x_lo)
    #sys.stderr.write("x_hi: %8.2f\n" % x_hi)
    gvals1 = grid_tally_gauss(xgrid, xvals1, gsigma)
    gvals2 = grid_tally_gauss(xgrid, xvals2, gsigma)
    gs1, gs2 = np.sum(gvals1), np.sum(gvals2)
    gcc = np.sum(gvals1 * gvals2)
    g12 = np.sqrt(gs1 * gs2)
    return gcc / g12
    #return gvals1, gvals2

## How to check the fit for a specific order:
fancy = True
using_wlmod = wavelength2 if fancy else wavelengths
def wlcheck(oidx):
    sord = int(spec_order_list[oidx])
    tdata = corresponding_thar[oidx]
    #wlref_wl_nm = spec_order_line_sets[oidx]
    line_ref_nm = spec_order_line_sets[oidx]
    #line_wl_um = line_wl_nm / 1e3   # convert to um
    model_wl_nm = using_wlmod[sord] * 1e3   # convert to nm
    lcenter_pix = measured_line_xpix[oidx]
    #lcenter_lam = np.interp(lcenter_pix, tdata['xpix'], model_wl_nm)
        #ar_xpixels = np.interp(ar_lines, wlen, thar['xpix'])
    #sys.stderr.write("lcenter_lam: %s\n" % str(lcenter_lam))
    line_ref_xx = np.interp(line_ref_nm, model_wl_nm, tdata['xpix'])
    #sys.stderr.write("\nDetected line centers (X):\n")
    #[sys.stderr.write("--> %8.2f\n" % x) for x in lcenter_pix]
    #sys.stderr.write("\nExpected line positions:")
    #[sys.stderr.write("--> %8.2f\n" % x) for x in line_ref_xx]
    result = crude_crossref(lcenter_pix, line_ref_xx, 5.0)
    sys.stderr.write("result: %10.5f\n" % result)
    return result
    #return tdata, model_wl_nm, line_wl_nm


tidx = 45
tord = int(spec_order_list[tidx])
tdata = corresponding_thar[tidx]
line_ref_nm = spec_order_line_sets[tidx]
model_wl_nm = using_wlmod[tord] * 1e3
line_xpix = measured_line_xpix[tidx]
line_refx = np.interp(line_ref_nm, model_wl_nm, tdata['xpix'])
#segs_meas_data = wlsm._list_segments(line_xpix)
#segs_lref_data = wlsm._list_segments(line_refx)
#segs_meas = segs_meas_data['seg']
#segs_lref = segs_lref_data['seg']
#diffs = segs_meas[:, None] - segs_lref[None, :]
#nseg_dims = (len(segs_meas), len(segs_lref))
#nobj_dims = (len(line_xpix), len(line_refx))

wlsm.set_measured_peaks(line_xpix)
wlsm.set_reference_lines(line_refx)

#len_range = (-0.2, 0.2)
len_tol   = np.log10(1.1)
len_bins  = 11
len_range = wlsm.bintol_range(len_bins, len_tol)
tdivs = (3,)

use_ranges = (len_range,)
use_nbins  = (len_bins,)
best_pars = wlsm.dither_hist_best_fit(use_ranges, use_nbins,
        tdivs, mode='weighted')

line_pairs = wlsm.matched_source_indexes()
midx, ridx = zip(*line_pairs) 

print(line_xpix[midx,])
print(line_refx[ridx,])

#scale_nbins = 50
#scale_range = (-0.2, 0.2)
#nb_attempts = [30, 35, 40, 45, 50, 55, 60]

#mh_results = multi_search(diffs, scale_range, nb_attempts)

#hh, edges = np.histogram(diffs.flatten(), 
#        bins=scale_nbins, range=scale_range)
#peak = np.unravel_index(hh.argmax(), hh.shape)
##vparams = [0.5*(ee[ii]+ee[ii+1]) for (ee,ii) in zip(edges, peak)]
#vpeak = 0.5 * (edges[peak[0]] + edges[peak[0]+1])
#bsize = (max(scale_range) - min(scale_range)) / float(scale_nbins)

## Extract pairings:
#tol = 0.75 * bsize
#hits = np.abs(diffs - vpeak) < 0.5*bsize



#sys.exit(0)

## ----------------------------------------------------------------------- ##
## ----------------------------------------------------------------------- ##

###-----------------------------------------------------------------------

max_lines_per_order = 30

## Visual inspection of ThAr data vs wavelength solution:
corresponding_thar = f0_thar_data if fib_which==0 else f1_thar_data
def oinspect(oidx, ww2=False, wlmode=False,
        sdata=corresponding_thar, pad=0.1):
    thar = sdata[oidx]
    sord = thar_specord[oidx]
    wlen = 1e3 * wavelength2[sord] if ww2 else 1e3 * wavelengths[sord]
    wl2pix = lambda x: np.interp(x, wlen, thar['xpix'])
    pix2wl = lambda x: np.interp(x, thar['xpix'], wlen)
    #wlen *= 1e3     # switch to nm
    #wlen = wavelengths[oidx] * 1e3  # switch to nm
    #sys.stderr.write("wlen.size: %d\n" % wlen.size)
    #sys.stderr.write("xpix.size: %d\n" % thar['xpix'].size)
    wlrange = wlen.max() - wlen.min()
    wl1 = wlen.min() - pad * wlrange
    wl2 = wlen.max() + pad * wlrange
    sys.stderr.write("oidx %d covers wavelength range: %.3f to %.3f\n" 
            % (oidx, wlen.min(), wlen.max()))

    fig = plt.figure(1, figsize=(10,5))
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.grid(True)

    #ax1.plot(thar['xpix'], thar['spec'])
    ax1.plot(wlen, thar['spec'])
    #ax1.set_yscale('log')
    ax1.set_yscale('linear')
    ax1.set_xlabel('Wavelength (nm)')
    ax1.set_xlim(wl1, wl2)

    # Overplot NIST Argon lines:
    ar_lines = wlf.get_nist_lines(wl1, wl2, reltol=0.05) #reltol=0.001)
    if (ar_lines.size > 0):
        #ar_xpixels = np.interp(ar_lines, wlen, thar['xpix'])
        #ar_show = wl2pix(ar_lines)     # pixels
        ar_show = ar_lines              # wavelength
        for line in ar_show[:-1]:
            ax1.axvline(line, ls=':', c='r')
        ax1.axvline(ar_show[-1], ls=':', c='r', label='NIST Argon')

    # Overplot Lovis & Pepe (2007) lines:
    thar_lines = wlf.get_lope_lines(wl1, wl2, reltol=0.05)
    if (thar_lines.size > 0):
        #thar_xpix = wl2pix(thar_lines)
        #thar_show = wl2pix(thar_lines)  # pixels
        thar_show = thar_lines          # wavelength
        for line in thar_show[:-1]:
            ax1.axvline(line, ls=':', c='g')
        ax1.axvline(thar_show[-1], ls=':', c='g', label='Lovis_Pepe_2007')

    #return
    ax3 = ax1.twiny()
    ax3.set_xlim(ax1.get_xlim())
    ax3.set_xlabel("X Pixel")
    xpix_ticks_xpix = 1e3 * np.arange(5)
    xpix_ticks_wlen = pix2wl(xpix_ticks_xpix)
    #ax3.set_xticks(xpix_ticks_wlen, xpix_ticks_xpix)
    ax3.set_xticks(xpix_ticks_wlen)
    ax3.set_xticklabels(xpix_ticks_xpix)


    ax1.legend(loc='upper right')
    fig.tight_layout()
    plt.draw()
    return



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
