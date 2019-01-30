#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Compute order spacing with and without accounting for half the cross-
# dispersion occurring prior to grating. Effectively, the "accounting" case
# allows the gamma angle to vary with wavelength whereas the naive case uses
# a fixed gamma angle for all wavelengths.
#
# Rob Siverd
# Created:       2019-01-20
# Last modified: 2019-01-30
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.1.0"

## Optional matplotlib control:
#from matplotlib import use, rc, rcParams
#from matplotlib import use
#from matplotlib import rc
#from matplotlib import rcParams
#use('GTKAgg')  # use GTK with Anti-Grain Geometry engine
#use('agg')     # use Anti-Grain Geometry engine (file only)
#use('ps')      # use PostScript engine for graphics (file only)
#use('cairo')   # use Cairo (pretty, file only)
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('font',**{'sans-serif':'Arial','family':'sans-serif'})
#rc('text', usetex=True) # enables text rendering with LaTeX (slow!)
#rcParams['axes.formatter.useoffset'] = False   # v. 1.4 and later

## Python version-agnostic module reloading:
try:
    reload                              # Python 2.7
except NameError:
    try:
        from importlib import reload    # Python 3.4+
    except ImportError:
        from imp import reload          # Python 3.0 - 3.3

## Modules:
#import argparse
#import mimetypes
#import linecache
#import getopt
#import shutil
#import resource
#import signal
#import glob
#import gc
import os
import sys
import time
#import vaex
#import calendar
#import ephem
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
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#import matplotlib.ticker as mt
#import matplotlib._pylab_helpers as hlp
#from matplotlib.colors import LogNorm
#from matplotlib import colors
#import matplotlib.colors as mplcolors
#import matplotlib.gridspec as gridspec
from functools import partial
#from collections import OrderedDict
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

## Fast FITS I/O:
#try:
#    import fitsio
#except ImportError:
#    sys.stderr.write("\nError: fitsio module not found!\n")
#    sys.exit(1)

## FITS I/O:
#try:
#    import astropy.io.fits as pf
#except ImportError:
#    try:
#       import pyfits as pf
#    except ImportError:
#        sys.stderr.write("\nError!  No FITS I/O module found!\n"
#               "Install either astropy.io.fits or pyfits and try again!\n\n")
#        sys.exit(1)

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

## Mayhem extraction tools:
import nres_extraction
reload(nres_extraction)
nrex = nres_extraction
trio = nres_extraction.TraceIO()
frox = nres_extraction.FlatRelativeOptimalExtraction()
#test_traces = "./scratch/trace1_med_nres01_tung01_20181125_07d_06b_00a.fits"
test_traces = "derp.fits"
double_path = "images/med_nres01_thar01_20181125_07d_06b_00a.fits.fz"
for item in (test_traces, double_path):
    if not os.path.isfile(item):
        sys.stderr.write("Can't find file: '%s'\n" % item)
        sys.exit(1)

## Load traces for fiddling:
trdata = trio.load_traces(test_traces)

## Spectrograph/optics brilliance:
import spectrograph_optics
reload(spectrograph_optics)

##--------------------------------------------------------------------------##
## Spectrograph hardware config items (constants):
nres_pixel_size_mm = 0.015  # 15-micron pixels
nres_gratio = 4.0           # NRES uses R4 grating
nres_ruling_lmm = 41.59     # lines per mm ruled
nres_blaze_angle_rad = np.arctan(nres_gratio)
nres_spacing_um = 1e3 / nres_ruling_lmm     # grating spacing in microns

nres_prism_glass = "PBM2"   # glass type used in cross-dispersing prism
nres_prism_apex_deg = 55.0  # apex angle of cross-dispersing prism
nres_prism_apex_rad = np.radians(nres_prism_apex_deg)

nres_center_wl_um = 0.479   # [I THINK] light wavelength nearest CCD center

nres_nominal_gamma = 0.0    # [radians] gamma angle for nres_center_wl_um
#nres_nominal_gamma = np.radians(1.0)

nres_focallen_mm = 375.15   # approximate camera focal length



##--------------------------------------------------------------------------##
sog = spectrograph_optics.Glass(nres_prism_glass)
ogt = spectrograph_optics.GratingTools(nres_gratio,
        lines_per_mm=nres_ruling_lmm)

## Refractive index of central wavelength:
center_wl_nn = sog.refraction_index(nres_center_wl_um)

## Nominal incidence angle onto prism (~same for all wavelengths):
incident_ang_rad = np.arcsin(center_wl_nn * np.sin(0.5 * nres_prism_apex_rad))

## Minimum deflection angle at nominal wavelength:
min_dev_rad = 2.0 * incident_ang_rad - nres_prism_apex_rad



## Effective gamma angle is the nominal gamma angle plus the changes in prism
## deflection angle due to difference from nominal wavelength.
def gamma_eff(gamma_nom, wlen_um):
    # delta_gamma = deflection_rad - min_dev_rad
    deflect_r = spectrograph_optics.prism_deflection_n(incident_ang_rad,
            nres_prism_apex_rad, sog.refraction_index(wlen_um))
    return gamma_nom + deflect_r - min_dev_rad

use_gamma_eff = partial(gamma_eff, nres_nominal_gamma)

## From existing solutions ...
nres_sine_alpha = 0.971747764900

## To get central wavelengths, want to solve:
## m * λ = 2 * d * cos(γ(λ)) * cos(θ) * sin(θ_B)
## where
##  * θ   = facet illumination angle
##  * θ_B = blaze angle
#blaze_angle_rad = np.arctan(4.0)    # R4 grating
#gamma_angle_rad = 0.0               # grating gamma angle
#facet_angle_rad = np.radians(2.0)   # how close to blaze angle
#alpha_angle_rad = blaze_angle_rad + facet_angle_rad
blaze_r = np.arctan(4.0)
#facet_r = np.radians(0.5)
facet_r = np.arcsin(nres_sine_alpha) - blaze_r
fixed_geometry = 2.0 * nres_spacing_um * np.cos(facet_r) * np.sin(blaze_r)

def lamcen_residual(wlen, order=0):
    ls = wlen * order
    rs = fixed_geometry * np.cos(use_gamma_eff(wlen))
    return ls - rs

def iter_calc_lamcen(order):
    kw = {'order':order}
    runme = partial(lamcen_residual, **kw)
    return opti.bisect(runme, 0.0, 10.0)

## NOTES:
## alpha = blaze_angle + facet_angle
##  beta = blaze_angle - facet_angle    (center)

##--------------------------------------------------------------------------##
## Orders to evaluate:
useful_orders = 51.0 + np.arange(69.0)

## Old way:
spec_order_list = np.copy(useful_orders)
spec_order_wlmid, spec_order_FSR, spec_order_angsize = \
        ogt.get_order_params(spec_order_list, units='um')
spec_order_table = {kk:vv for kk,vv in zip(spec_order_list, spec_order_wlmid)}

## New way:
iter_order_wlmid = np.array([iter_calc_lamcen(x) for x in spec_order_list])
iter_order_table = {kk:vv for kk,vv in zip(spec_order_list, iter_order_wlmid)}

## Comparison of methods:
for ii,ww in enumerate(spec_order_wlmid):
    w2 = iter_order_wlmid[ii]
    sys.stderr.write("oid %3d --> %10.5f nm || %10.5f\n" % (ii, 1e3 * ww, 1e3 * w2))

## Store results for external use/comparison:
sys.stderr.write("Storing central wavelengths ... ")
with open('lam_ctr_data.csv', 'w') as f:
    f.write("oidx,spec_ord,lam_dumb,lam_iter\n")
    for ii,w1 in enumerate(spec_order_wlmid):
        so = spec_order_list[ii]
        w2 = iter_order_wlmid[ii]
        f.write("%2d,%3d,%8.6f,%8.6f\n" % (ii, so, w1, w2))
        pass
    pass
sys.stderr.write("done.\n")



##--------------------------------------------------------------------------##
## Given a set of central wavelengths, compute order separation:
def calc_order_sep(central_wlen_um):
    refr_idx = sog.refraction_index(central_wlen_um)
    deflections = spectrograph_optics.prism_deflection_n(incident_ang_rad,
                nres_prism_apex_rad, refr_idx)
    angle_change = deflections - min_dev_rad
    return 2*angle_change

nres_flen_pix = nres_focallen_mm / nres_pixel_size_mm
old_sep = calc_order_sep(spec_order_wlmid) * nres_flen_pix
new_sep = calc_order_sep(iter_order_wlmid) * nres_flen_pix




##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Dispersion is a function of angles and wavelength:
## dβ/dλ = n / (d * cosγ * cosβ)
## dβ/dλ = (sinα + sinβ) / (λ * cosβ)
#dbeta_dlambda = (nres_sine_alpha + 

##--------------------------------------------------------------------------##
## Arrays of spectrograph coordinates:
ny, nx = 4096, 4096
x_list = (0.5 + np.arange(nx)) / nx - 0.5            # relative (centered)
y_list = (0.5 + np.arange(ny)) / ny - 0.5            # relative (centered)
#xx, yy = np.meshgrid(x_list, y_list)                 # relative (centered)
xx, yy = np.meshgrid(nx*x_list, ny*y_list)           # absolute (centered)
#xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))   # absolute
#yy, xx = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij') # absolute
#yy, xx = np.nonzero(np.ones_like(img_vals))          # absolute
#yy, xx = np.mgrid[0:ny,   0:nx].astype('uint16')     # absolute (array)
#yy, xx = np.mgrid[1:ny+1, 1:nx+1].astype('uint16')   # absolute (pixel)
#xx *= nres_pixel_size_mm
#yy *= nres_pixel_size_mm

gamma2d = np.arcsin(yy * nres_pixel_size_mm / nres_focallen_mm)
##gamma2d = np.arctan(yy * nres_pixel_size_mm / nres_focallen_mm)
beta_2d = np.arcsin(xx * nres_pixel_size_mm / nres_focallen_mm) # beta - beta_c
del xx, yy, x_list, y_list

#ordernum = 100
#wlen2d_test = nres_spacing_um * 

##--------------------------------------------------------------------------##
## Plot config:

# gridspec examples:
# https://matplotlib.org/users/gridspec.html

#gs1 = gridspec.GridSpec(4, 4)
#gs1.update(wspace=0.025, hspace=0.05)  # set axis spacing

#ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=3) # top-left + center + right
#ax2 = plt.subplot2grid((3, 3), (1, 0), colspan=2) # mid-left + mid-center
#ax3 = plt.subplot2grid((3, 3), (1, 2), rowspan=2) # mid-right + bot-right
#ax4 = plt.subplot2grid((3, 3), (2, 0))            # bot-left
#ax5 = plt.subplot2grid((3, 3), (2, 1))            # bot-center


##--------------------------------------------------------------------------##
fig_dims = (12, 10)
fig = plt.figure(1, figsize=fig_dims)
plt.gcf().clf()
#fig, axs = plt.subplots(2, 2, sharex=True, figsize=fig_dims, num=1)
# sharex='col' | sharex='row'
#fig.frameon = False # disable figure frame drawing
#fig.subplots_adjust(left=0.07, right=0.95)
#ax1 = plt.subplot(gs[0, 0])
#ax1 = fig.add_subplot(111)
#ax1 = fig.add_axes([0, 0, 1, 1])
#ax1.patch.set_facecolor((0.8, 0.8, 0.8))
#ax1.grid(True)
#ax1.axis('off')

## Disable axis offsets:
#ax1.xaxis.get_major_formatter().set_useOffset(False)
#ax1.yaxis.get_major_formatter().set_useOffset(False)

#ax1.plot(kde_pnts, kde_vals)

#blurb = "some text"
#ax1.text(0.5, 0.5, blurb, transform=ax1.transAxes)
#ax1.text(0.5, 0.5, blurb, transform=ax1.transAxes,
#      va='top', ha='left', bbox=dict(facecolor='white', pad=10.0))
#      fontdict={'family':'monospace'}) # fixed-width

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
#for label in ax1.get_xticklabels():
#    label.set_rotation(30)
#    label.set_fontsize(14) 

#ax1.set_xlim(nice_limits(xvec, pctiles=[1,99], pad=1.2))
#ax1.set_ylim(nice_limits(yvec, pctiles=[1,99], pad=1.2))

#spts = ax1.scatter(x, y, lw=0, s=5)
#cbar = fig.colorbar(spts, orientation='vertical')
#cbar.formatter.set_useOffset(False)
#cbar.update_ticks()

#fig.tight_layout() # adjust boundaries sensibly, matplotlib v1.1+
#plt.draw()
#fig.savefig(plot_name, bbox_inches='tight')

# cyclical colormap ... cmocean.cm.phase
# cmocean: https://matplotlib.org/cmocean/




######################################################################
# CHANGELOG (compare_order_Y_positions.py):
#---------------------------------------------------------------------
#
#  2019-01-20:
#     -- Increased __version__ to 0.0.1.
#     -- First created compare_order_Y_positions.py.
#
