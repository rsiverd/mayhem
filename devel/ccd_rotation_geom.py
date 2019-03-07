#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Illustrate CCD rotation geometry.
#
# Rob Siverd
# Created:       2019-03-07
# Last modified: 2019-03-07
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
import resource
#from numpy.lib.recfunctions import append_fields
#import datetime as dt
#from dateutil import parser as dtp
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

## Rotation matries and more:
import fov_rotation
reload(fov_rotation)
r3d = fov_rotation.Rotate3D()


##--------------------------------------------------------------------------##

unlimited = (resource.RLIM_INFINITY, resource.RLIM_INFINITY)
if (resource.getrlimit(resource.RLIMIT_DATA) == unlimited):
    resource.setrlimit(resource.RLIMIT_DATA,  (3e9, 6e9))
if (resource.getrlimit(resource.RLIMIT_AS) == unlimited):
    resource.setrlimit(resource.RLIMIT_AS, (3e9, 6e9))

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

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




##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## CCD config:
spec_rotation = 13.091
#ny, nx = 1024, 1024
ny, nx =  512,  512

## Make unrotated CCD coordinates:
x_list = (0.5 + np.arange(nx)) / nx - 0.5            # relative (centered)
y_list = (0.5 + np.arange(ny)) / ny - 0.5            # relative (centered)
#xx, yy = np.meshgrid(x_list, y_list)                 # relative (centered)
xx, yy = np.meshgrid(nx*x_list, ny*y_list)           # absolute (centered)

## Make rotated CCD coordinates:
xxs, yys = [], []
for rowx,rowy in zip(xx, yy):
    tx, ty = ccd2spec_xy(rowx, rowy, spec_rotation)
    xxs.append(tx)
    yys.append(ty)
xxs = np.array(xxs)
yys = np.array(yys)

## Note spectrum Y-axis:
midline = (np.abs(xxs) < 1.0)



##--------------------------------------------------------------------------##
## Plot config:

##--------------------------------------------------------------------------##
fig_dims = (12, 10)
fig1 = plt.figure(1, figsize=fig_dims)
fig2 = plt.figure(2, figsize=fig_dims)
plt.gcf().clf()
#fig, axs = plt.subplots(2, 2, sharex=True, figsize=fig_dims, num=1)
# sharex='col' | sharex='row'
#fig.frameon = False # disable figure frame drawing
#fig.subplots_adjust(left=0.07, right=0.95)
#ax1 = plt.subplot(gs[0, 0])
ax1 = fig1.add_subplot(111, aspect='equal')
ax2 = fig2.add_subplot(111, aspect='equal')
#ax1 = fig.add_axes([0, 0, 1, 1])
#ax1.patch.set_facecolor((0.8, 0.8, 0.8))
#ax1.grid(True)
ax1.grid(False)
#ax1.axis('off')

## Show rotated CCD coordinates with X=0 highlighted:
ax1.scatter(xxs.flatten(), yys.flatten(), lw=0, s=3, c='b')
ax1.scatter(xxs[midline], yys[midline], lw=0, s=4, c='r')
ax1.set_xlabel("Rotated CCD-Y / native spectrograph")
ax1.set_ylabel("Rotated CCD-X / native spectrograph")

## -----------------------------------------------------------------------
## Unrotated CCD coordinates:
ax2.scatter(xx.flatten(), yy.flatten(), lw=0, s=3, c='b')
ax2.scatter(xx[midline], yy[midline], lw=0, s=4, c='r')
ax2.set_xlabel("Unrotated CCD-X / Rotated Spectrograph")
ax2.set_ylabel("Unrotated CCD-Y / Rotated Spectrograph")

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

fig1.tight_layout() # adjust boundaries sensibly, matplotlib v1.1+
fig2.tight_layout() # adjust boundaries sensibly, matplotlib v1.1+
plt.draw()
#fig.savefig(plot_name, bbox_inches='tight')

# cyclical colormap ... cmocean.cm.phase
# cmocean: https://matplotlib.org/cmocean/




######################################################################
# CHANGELOG (ccd_rotation_geom.py):
#---------------------------------------------------------------------
#
#  2019-03-07:
#     -- Increased __version__ to 0.0.1.
#     -- First created ccd_rotation_geom.py.
#
