#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Measure noise vs. NUMSTACK in stacked NRES biases.
#
# Rob Siverd
# Created:       2018-11-19
# Last modified: 2018-11-19
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.0.1"

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
import pandas as pd
#import statsmodels.api as sm
#import statsmodels.formula.api as smf
#from statsmodels.regression.quantile_regression import QuantReg
#import PIL.Image as pli
#import seaborn as sns
#import cmocean
import theil_sen as ts
#import window_filter as wf
#import itertools as itt

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Quick ASCII I/O:
data_file = "img_noise.txt"
#data_file = 'data.txt'
#all_data = np.genfromtxt(data_file, dtype=None)
#all_data = np.genfromtxt(data_file, dtype=None, names=True, autostrip=True)
#all_data = np.genfromtxt(data_file, dtype=None, names=True, autostrip=True,
#                 delimiter='|', comments='%0%0%0%0')
#                 loose=True, invalid_raise=False)
#all_data = aia.read(data_file)
#all_data = pd.read_csv(data_file)
all_data = pd.read_table(data_file, delim_whitespace=True)
#all_data = pd.read_table(data_file, sep='|')
#fields = all_data.dtype.names
#if not fields:
#    x = all_data[:, 0]
#    y = all_data[:, 1]
#else:
#    x = all_data[fields[0]]
#    y = all_data[fields[1]]

##--------------------------------------------------------------------------##
## Theil-Sen line-fitting:
#model = ts.linefit(xvals, yvals)
#icept, slope = ts.linefit(xvals, yvals)
xvals = np.log10(all_data['NUMSTACK'].values)
yvals = np.log10(all_data['ISIGMA'].values)
icept, slope = ts.linefit(xvals, yvals)
fit_exponent = slope
fit_multiplier = 10**icept

def fitted_x2y(xvec, icept, slope):
    return icept + slope*xvec

def fitted_y2x(yvec, icept, slope):
    return (yvec - icept) / slope

#def loglog_eval(xvals, icept, slope):
#    return 
def noisecalc(nstack, icept, slope):
    return 10.0**(fitted_x2y(np.log10(nstack), icept, slope))

def stackcalc(bias_noise, icept, slope):
    return 10.0**(fitted_y2x(np.log10(bias_noise), icept, slope))

#def noisecalc(nstack, fit_mult, fit_expo):
#    return fitted_x2y(np.log10(nstack), 

#def stackcalc(bias_noise, fit_mult, fit_expo):
#    return (bias_noise / fit_mult)**(1.0 / fit_expo)

fstack = np.arange(100) + 1.0
fnoise = noisecalc(fstack, icept, slope)

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
fig_dims = (8, 7)
fig = plt.figure(1, figsize=fig_dims)
plt.gcf().clf()
#fig, axs = plt.subplots(2, 2, sharex=True, figsize=fig_dims, num=1)
# sharex='col' | sharex='row'
#fig.frameon = False # disable figure frame drawing
#fig.subplots_adjust(left=0.07, right=0.95)
#ax1 = plt.subplot(gs[0, 0])
ax1 = fig.add_subplot(111)
#ax1 = fig.add_axes([0, 0, 1, 1])
#ax1.patch.set_facecolor((0.8, 0.8, 0.8))
ax1.grid(True)
#ax1.axis('off')

## Disable axis offsets:
#ax1.xaxis.get_major_formatter().set_useOffset(False)
#ax1.yaxis.get_major_formatter().set_useOffset(False)

#ax1.plot(kde_pnts, kde_vals)
ax1.scatter(all_data['NUMSTACK'], all_data['ISIGMA'])
ax1.plot(fstack, fnoise, c='r')
ax1.set_xlabel("N stacked")
ax1.set_ylabel("StdDev (electrons)")
ax1.set_xlim(-1.0, 100.0)

plot_name = 'bias_noise_curve.png'
fig.tight_layout() # adjust boundaries sensibly, matplotlib v1.1+
plt.draw()
fig.savefig(plot_name, bbox_inches='tight')

# cyclical colormap ... cmocean.cm.phase
# cmocean: https://matplotlib.org/cmocean/




######################################################################
# CHANGELOG (noise_check.py):
#---------------------------------------------------------------------
#
#  2018-11-19:
#     -- Increased __version__ to 0.0.1.
#     -- First created noise_check.py.
#
