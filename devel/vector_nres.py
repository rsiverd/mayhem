#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Make a vector representation of the spectrograph (approximate).
#
# Rob Siverd
# Created:       2019-02-19
# Last modified: 2019-02-19
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
np.set_printoptions(suppress=True, linewidth=160)
#import pandas as pd
#import statsmodels.api as sm
#import statsmodels.formula.api as smf
#from statsmodels.regression.quantile_regression import QuantReg
#import PIL.Image as pli
#import seaborn as sns
#import cmocean
#import theil_sen as ts
#import window_filter as wf
import itertools as itt

## Rotation kit:
import fov_rotation
reload(fov_rotation)
r3d = fov_rotation.Rotate3D()

## Polygon kit:
import polygon_optics
reload(polygon_optics)

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
def ldmap(things):
    return dict(zip(things, range(len(things))))

def argnear(vec, val):
    return (np.abs(vec - val)).argmin()




##--------------------------------------------------------------------------##
##------------------         Parse Command Line             ----------------##
##--------------------------------------------------------------------------##

## Parse arguments and run script:
#class MyParser(argparse.ArgumentParser):
#    def error(self, message):
#        sys.stderr.write('error: %s\n' % message)
#        self.print_help()
#        sys.exit(2)
#
### Enable raw text AND display of defaults:
#class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
#                        argparse.RawDescriptionHelpFormatter):
#    pass
#
### Parse the command line:
#if __name__ == '__main__':
#
#    # ------------------------------------------------------------------
#    descr_txt = """
#    PUT DESCRIPTION HERE.
#    
#    Version: %s
#    """ % __version__
#    parser = argparse.ArgumentParser(
#            prog='PROGRAM_NAME_HERE',
#            prog=os.path.basename(__file__),
#            #formatter_class=argparse.RawTextHelpFormatter)
#            description='PUT DESCRIPTION HERE.')
#            #description=descr_txt)
#    parser = MyParser(prog=os.path.basename(__file__), description=descr_txt)
#                          #formatter_class=argparse.RawTextHelpFormatter)
#    # ------------------------------------------------------------------
#    parser.set_defaults(thing1='value1', thing2='value2')
#    # ------------------------------------------------------------------
#    parser.add_argument('firstpos', help='first positional argument')
#    parser.add_argument('-s', '--site',
#            help='Site to retrieve data for', required=True)
#    parser.add_argument('-n', '--number_of_days', default=1,
#            help='Number of days of data to retrieve.')
#    parser.add_argument('-o', '--output_file', 
#            default='observations.csv', help='Output filename.')
#    parser.add_argument('--start', type=str, default=None, 
#            help="Start time for date range query.")
#    parser.add_argument('--end', type=str, default=None,
#            help="End time for date range query.")
#    parser.add_argument('-d', '--dayshift', required=False, default=0,
#            help='Switch between days (1=tom, 0=today, -1=yest', type=int)
#    parser.add_argument('-e', '--encl', nargs=1, required=False,
#            help='Encl to make URL for', choices=all_encls, default=all_encls)
#    parser.add_argument('-s', '--site', nargs=1, required=False,
#            help='Site to make URL for', choices=all_sites, default=all_sites)
#    parser.add_argument('-q', '--quiet', action='count', default=0,
#            help='less progress/status reporting')
#    parser.add_argument('-v', '--verbose', action='count', default=0,
#            help='more progress/status reporting')
#    parser.add_argument('--debug', dest='debug', default=False,
#            help='Enable extra debugging messages', action='store_true')
#    parser.add_argument('remainder', help='other stuff', nargs='*')
#    # ------------------------------------------------------------------
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
#    fmtparse.set_defaults(output_mode='pydict')
#    # ------------------------------------------------------------------
#    # Miscellany:
#    miscgroup = parser.add_argument_group('Miscellany')
#    miscgroup.add_argument('-q', '--quiet', action='count', default=0,
#            help='less progress/status reporting')
#    miscgroup.add_argument('-v', '--verbose', action='count', default=0,
#            help='more progress/status reporting')
#    miscgroup.add_argument('--debug', dest='debug', default=False,
#            help='Enable extra debugging messages', action='store_true')
#    # ------------------------------------------------------------------
#
#    context = parser.parse_args()
#    context.vlevel = 99 if context.debug else (context.verbose-context.quiet)
#
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
## COORDINATE SYSTEM SEEN FROM ABOVE (from +Z axis):

##  +Y
##   ^
##   |
##   |
##   |
##   ------> +X

## Z-axis "up" towards viewer

##--------------------------------------------------------------------------##
## Prism:
apex_angle_deg =  55.0
apex_angle_rad = np.radians(apex_angle_deg)
#apex_length_mm = 174.0
#long_length_mm = 205.7
#half_short_mm  = long_length_mm * np.sin(0.5 * apex_angle_rad)
prism_turn_deg = 23.507 # prism base 'turn' w.r.t. optical axis
short_edge_mm = 190.0
#large_edge_mm = 0.5 * short_edge_mm / np.sin(0.5 * apex_angle_rad)
#symmetryax_mm = 0.5 * short_edge_mm / np.tan(0.5 * apex_angle_rad)
height_mm     = 130.0

prpoly = polygon_optics.PolygonPrism(apex_angle_deg, short_edge_mm, height_mm)
prpoly.zrotate(np.radians(-90.0))
prpoly.zrotate(np.radians(prism_turn_deg))

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Grating geometry class:
class GratingPolygon(object):
    def __init__(self):
        self.data = {}
        return

## Grating:
gr_length_mm = 320.0
gr_width_mm  = 165.0
gr_height_mm =  50.5
grating_turn_deg = 44.827
grating_tilt_deg = 13.786

#vtx1_bot = np.array([        0.0,          0.0,  0.0])
#vtx2_bot = np.array([gr_width_mm,          0.0,  0.0])
#vtx3_bot = np.array([gr_width_mm, gr_length_mm,  0.0])
#vtx4_bot = np.array([        0.0, gr_length_mm,  0.0])
g_normal = np.array([        0.0,          0.0, -1.0])  # grating normal
#grating_vtx_bot = np.array([vtx1_bot, vtx2_bot, vtx3_bot, vtx4_bot])
#grating_vtx = np.vstack((grating_vtx_bot, 
#                         grating_vtx_bot + np.array([0.0, 0.0, gr_height_mm])))
#grating_vtx -= np.average(grating_vtx, axis=0)

#grating_vtx = grating_vtx.T     # transpose for rotations
#grating_vtx = r3d.xrotate_xyz(np.radians(-grating_tilt_deg), grating_vtx)
#grating_vtx = r3d.zrotate_xyz(np.radians(grating_turn_deg), grating_vtx)
#grating_vtx = np.array(grating_vtx.T)
g_normal = r3d.xrotate_xyz(np.radians(-grating_tilt_deg), g_normal).A1
g_normal = r3d.zrotate_xyz(np.radians(grating_turn_deg), g_normal).A1
#grating_vtx += np.array([-165.0, 235.0, 0.0])

grpoly = polygon_optics.PolygonGrating(gr_width_mm, gr_length_mm, gr_height_mm)
grpoly.xrotate(np.radians(-grating_tilt_deg))
grpoly.zrotate(np.radians(grating_turn_deg))
grpoly.shift_xyz(-165.0, 235.0, 0.0)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Edges from vertices:
def edges_from_vertices(vtx_array):
    nhalf = vtx_array.shape[0] / 2
    vtx_bot = vtx_array[:nhalf]
    vtx_top = vtx_array[-nhalf:]
    idx_list = np.arange(nhalf)
    idxpairs = zip(idx_list, np.roll(idx_list, -1))
    #bot_pairs = zip(idx_list, np.roll(idx_list, -1))
    #top_pairs = zip(idx_list + nhalf, np.roll(idx_list + nhalf, -1))
    #sys.stderr.write("nhalf: %d\n" % nhalf)
    #sys.stderr.write("top:\n%s\n" % str(vtx_top))
    #sys.stderr.write("bot:\n%s\n" % str(vtx_bot))
    #sys.stderr.write("idx_list: %s\n" % str(idx_list))
    #sys.stderr.write("idxpairs: %s\n" % str(idxpairs))
    edges = []
    edges.extend([vtx_bot[:, :2][x,] for x in idxpairs])
    edges.extend([vtx_top[:, :2][x,] for x in idxpairs])
    #sys.stderr.write("edges: %s\n" % str(edges))
    return edges


##--------------------------------------------------------------------------##
fig_dims = (12, 10)
fig = plt.figure(1, figsize=fig_dims)
plt.gcf().clf()
#fig, axs = plt.subplots(2, 2, sharex=True, figsize=fig_dims, num=1)
# sharex='col' | sharex='row'
#fig.frameon = False # disable figure frame drawing
#fig.subplots_adjust(left=0.07, right=0.95)
#ax1 = plt.subplot(gs[0, 0])
ax1 = fig.add_subplot(111, aspect='equal')
#ax1 = fig.add_axes([0, 0, 1, 1])
#ax1.patch.set_facecolor((0.8, 0.8, 0.8))
ax1.grid(True)
#ax1.axis('off')

## Disable axis offsets:
#ax1.xaxis.get_major_formatter().set_useOffset(False)
#ax1.yaxis.get_major_formatter().set_useOffset(False)

# -----------------------------------------------------------------------
# Prism vertices:
#for x,y,z in prism_vtx:
for x,y,z in prpoly.get_vertices():
    ax1.scatter(x, y, lw=0, s=30, c='b')

# Prism edges:
vx, vy, vz = zip(*prpoly.get_vertices())
for p1,p2 in itt.combinations(range(3), 2):
    ax1.plot([vx[p1], vx[p2]], [vy[p1], vy[p2]], c='b')
    pass

# -----------------------------------------------------------------------
# Grating vertices:
for x,y,z in grpoly.get_vertices():
    ax1.scatter(x, y, lw=0, s=30, c='g')

# Grating edges:
for pair in edges_from_vertices(grpoly.get_vertices()):
    vx, vy = pair.T
    ax1.plot(vx, vy, c='g')

ax1.set_xlim(-400, 200)
ax1.set_ylim(-500, 500)
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

fig.tight_layout() # adjust boundaries sensibly, matplotlib v1.1+
plt.draw()
#fig.savefig(plot_name, bbox_inches='tight')

# cyclical colormap ... cmocean.cm.phase
# cmocean: https://matplotlib.org/cmocean/




######################################################################
# CHANGELOG (vector_nres.py):
#---------------------------------------------------------------------
#
#  2019-02-19:
#     -- Increased __version__ to 0.0.1.
#     -- First created vector_nres.py.
#
