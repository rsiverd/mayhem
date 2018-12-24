#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# This script is written to investigate what (if any) mathematical
# relationship exists among the orders that could be used to simplify
# fitting and reduce jitter.
#
# Rob Siverd
# Created:       2018-05-08
# Last modified: 2018-12-24
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.0.2"

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
import resource
import signal
#import glob
import gc
import os
import sys
import time
import copy
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
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#import matplotlib.ticker as mt
#import matplotlib._pylab_helpers as hlp
#from matplotlib.colors import LogNorm
#from matplotlib import colors
#import matplotlib.colors as mplcolors
#import matplotlib.gridspec as gridspec
#from functools import partial
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
import theil_sen as ts
#import window_filter as wf
#import itertools as itt

import nres_extraction
reload(nres_extraction)
trio = nres_extraction.TraceIO()

##--------------------------------------------------------------------------##

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

## ASCII I/O:
#try:
#    import astropy.io.ascii as aia
#except ImportError:
#    sys.stderr.write("\nError: astropy module not found!\n")
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
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Get trace file from command line:
if (len(sys.argv) != 2):
    sys.stderr.write("Syntax: %s trace_file.fits\n"
            % os.path.basename(__FILE__))
    sys.exit(1)

trace_file = sys.argv[1]
if not os.path.isfile(trace_file):
    sys.stderr.write("File not found: %s\n" % trace_file)
    sys.exit(1)

##--------------------------------------------------------------------------##
## Load traces for analysis:

trdata = trio.load_traces(trace_file)
traces = trdata.get_traces()
pars = np.array([x['params'].tolist() for x in traces])

c0, c1, c2 = pars.T
fiber1 = slice(0, None, 2)
fiber2 = slice(1, None, 2)

##--------------------------------------------------------------------------##
## Theil-Sen line-fitting:
#model = ts.linefit(xvals, yvals)
#icept, slope = ts.linefit(xvals, yvals)

icept2, slope2 = ts.linefit(c0, np.sqrt(c0) * c1)
line2 = icept2 + c0 * slope2

icept3, slope3 = ts.linefit(c0, c0 * c2)
line3 = icept3 + c0 * slope3

f1_icept3, f1_slope3 = ts.linefit(c0[fiber1], (c0 * c2)[fiber1])
f2_icept3, f2_slope3 = ts.linefit(c0[fiber2], (c0 * c2)[fiber2])
f1_line3 = f1_icept3 + c0[fiber1] * f1_slope3
f2_line3 = f2_icept3 + c0[fiber2] * f2_slope3


f1_resid3 = f1_line3 - (c0 * c2)[fiber1]
f2_resid3 = f2_line3 - (c0 * c2)[fiber2]

##--------------------------------------------------------------------------##
## Force c2 onto best-fit line:
linearized_pars = np.copy(pars)
linearized_pars[fiber1, 2] = f1_line3 / c0[fiber1]
linearized_pars[fiber2, 2] = f2_line3 / c0[fiber2]

linearized_data = copy.deepcopy(data)
for tfit,ldata in zip(linearized_data, linearized_pars):
    tfit['params'] = ldata

mod_save = "modified_trace.fits"
sys.stderr.write("Saving adjusted trace file '%s' ... " % mod_save)
trio.store_traces(mod_save, linearized_data)
sys.stderr.write("done.\n")

##--------------------------------------------------------------------------##
## Misc:
#def log_10_product(x, pos):
#   """The two args are the value and tick position.
#   Label ticks with the product of the exponentiation."""
#   return '%.2f' % (x)  # floating-point
#
#formatter = plt.FuncFormatter(log_10_product) # wrap function for use

## Convenient, percentile-based plot limits:
#def nice_limits(vec, pctiles=[1,99], pad=1.2):
#    ends = np.percentile(vec, pctiles)
#    middle = np.average(ends)
#    return (middle + pad * (ends - middle))

## Convenient plot limits for datetime/astropy.Time content:
#def nice_time_limits(tvec, buffer=0.05):
#    lower = tvec.min()
#    upper = tvec.max()
#    ndays = upper - lower
#    return ((lower - 0.05*ndays).datetime, (upper + 0.05*ndays).datetime)

## Convenient limits for datetime objects:
#def dt_limits(vec, pad=0.1):
#    tstart, tstop = vec.min(), vec.max()
#    trange = (tstop - tstart).total_seconds()
#    tpad = dt.timedelta(seconds=pad*trange)
#    return (tstart - tpad, tstop + tpad)

##--------------------------------------------------------------------------##
## Plot config:

# gridspec examples:
# https://matplotlib.org/users/gridspec.html

#gs1 = gridspec.GridSpec(4, 4)
#gs1.update(wspace=0.025, hspace=0.05)  # set axis spacing

##--------------------------------------------------------------------------##
fig_dims = (12, 10)
fig = plt.figure(1, figsize=fig_dims)
plt.gcf().clf()
#fig, axs = plt.subplots(3, 2, sharex=True, figsize=fig_dims, num=1)
fig, axs = plt.subplots(3, 3, sharex='col', figsize=fig_dims, num=1)
# sharex='col' | sharex='row'
#fig.frameon = False # disable figure frame drawing
#fig.subplots_adjust(left=0.07, right=0.95)
#ax1 = plt.subplot(gs[0, 0])
#ax1 = fig.add_subplot(111)
#ax1 = fig.add_axes([0, 0, 1, 1])
#ax1.patch.set_facecolor((0.8, 0.8, 0.8))
#ax1.grid(True)
#ax1.axis('off')

pkwargs = {'ls':'', 'ms':2, 'marker':'.'}
ypos = pars[:, 0]
yrng = np.arange(ypos.size)
for i in range(3):
    #ax = axs[i, 0]
    ppar0 = np.copy(pars)
    ppar1 = pars * np.sqrt(ypos[:, None])
    ppar2 = pars * ypos[:, None]
    #rpar = pars * np.sqrt(ypos[:, None])
    #axs[i, 0].plot(yrng[fiber1], ppar0[fiber1, i], **pkwargs)
    #axs[i, 0].plot(yrng[fiber2], ppar0[fiber2, i], **pkwargs)
    axs[i, 0].plot(ypos[fiber1], ppar0[fiber1, i], label='p0', **pkwargs)
    axs[i, 0].plot(ypos[fiber2], ppar0[fiber2, i], **pkwargs)
    axs[i, 1].plot(ypos[fiber1], ppar1[fiber1, i], **pkwargs)
    axs[i, 1].plot(ypos[fiber2], ppar1[fiber2, i], **pkwargs)
    axs[i, 2].plot(ypos[fiber1], ppar2[fiber1, i], **pkwargs)
    axs[i, 2].plot(ypos[fiber2], ppar2[fiber2, i], **pkwargs)

axs[1, 1].plot(ypos, line2, c='r', lw=0.5)
axs[2, 2].plot(ypos, line3, c='r', lw=0.5)
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

#ax1.set_xlim(nice_limits(xvec, pctiles=[1,99], pad=1.2))
#ax1.set_ylim(nice_limits(yvec, pctiles=[1,99], pad=1.2))

#spts = ax1.scatter(x, y, lw=0, s=5)
#cbar = fig.colorbar(spts, orientation='vertical')
#cbar.formatter.set_useOffset(False)
#cbar.update_ticks()

fig.tight_layout() # adjust boundaries sensibly, matplotlib v1.1+
plt.draw()

# cyclical colormap ... cmocean.cm.phase
# cmocean: https://matplotlib.org/cmocean/




######################################################################
# CHANGELOG (plot_traces.py):
#---------------------------------------------------------------------
#
#  2018-05-08:
#     -- Increased __version__ to 0.0.1.
#     -- First created plot_traces.py.
#
