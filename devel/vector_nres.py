#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Make a vector representation of the spectrograph (approximate).
#
# Rob Siverd
# Created:       2019-02-19
# Last modified: 2019-02-24
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
import polyhedron_optics
reload(polyhedron_optics)
po = polyhedron_optics

## Raytracing stuff:
import raytrace_tools
reload(raytrace_tools)
rt = raytrace_tools

## Some NRES hardware specs:
nres_prism_glass = "PBM2"       # glass type used in cross-dispersing prism
nres_prism_apex_deg = 55.0

## More optics goodies:
import spectrograph_optics
reload(spectrograph_optics)
sog = spectrograph_optics.Glass(nres_prism_glass)
#nrp = spectrograph_optics.Prism(nres_prism_glass, nres_prism_apex_deg)
dppgp = spectrograph_optics.DoublePassPrismGratingPrism()


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
#apex_angle_deg =  55.0
apex_angle_deg = nres_prism_apex_deg
apex_angle_rad = np.radians(apex_angle_deg)
#apex_length_mm = 174.0
#long_length_mm = 205.7
#half_short_mm  = long_length_mm * np.sin(0.5 * apex_angle_rad)
prism_turn_deg = 23.507 # prism base 'turn' w.r.t. optical axis
short_edge_mm = 190.0
#large_edge_mm = 0.5 * short_edge_mm / np.sin(0.5 * apex_angle_rad)
#symmetryax_mm = 0.5 * short_edge_mm / np.tan(0.5 * apex_angle_rad)
height_mm     = 130.0

prpoly = po.IsosPrismPolyhedron(apex_angle_deg, short_edge_mm, height_mm)
prpoly.zrotate(np.radians(-90.0))
prpoly.zrotate(np.radians(prism_turn_deg))

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Grating:
gr_length_mm = 320.0
gr_width_mm  = 165.0
gr_height_mm =  50.5
grating_turn_deg = 44.827
grating_tilt_deg = 13.786

grpoly = po.GratingPolyhedron(gr_width_mm, gr_length_mm, gr_height_mm)
grpoly.xrotate(np.radians(-grating_tilt_deg))
grpoly.zrotate(np.radians(grating_turn_deg))
grpoly.shift_xyz(-165.0, 235.0, 0.0)
#grpoly.shift_xyz(-165.0, 235.0, 0.0)
#grpoly.shift_vec(np.array([-165.0, 235.0, 0.0]))

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Center prism face1 on origin (move grating and prism in tandem):
nudge = -1.0 * np.copy(prpoly.get_face('face1')['center'])
prpoly.shift_vec(nudge)
grpoly.shift_vec(nudge)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## An initial input beam:
wl_initial = 0.8                        # ray wavelength (microns)
v_initial = np.array([0.0, 1.0, 0.0])   # headed in Y-direction
n_spec_air = 1.0                        # spectrograph air index of refraction
n_pris_glass = sog.refraction_index(wl_initial) # prism index of refraction
n1_n2_ratio = n_spec_air / n_pris_glass
v_reflect, v_refract = rt.calc_surface_vectors(v_initial,
                        prpoly.get_face('face1')['normal'], n1_n2_ratio)

## Record directions for plotting:
path1 = v_initial.copy()
path2 = v_refract.copy()

## Exit point from prism face2:
f1_isect = np.zeros(3)
prf2 = prpoly.get_face('face2')
f2_isect = rt.line_plane_intersection(np.zeros(3), path2, prf2['center'], prf2['normal'])

## Exit direction from prism face2:
v_reflect, v_refract = rt.calc_surface_vectors(path2,
                        -prf2['normal'], n_pris_glass / n_spec_air)
path3 = v_refract.copy()

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Normalized prf2 vertices:
b1, b2 = prf2['basis']
use_origin = prf2['vertices'][0]
#use_origin = prf2['center'][0]
rvtx_list = np.array([x-use_origin for x in prf2['vertices']])
for item in rvtx_list:
    sys.stderr.write("item: %s\n" % str(item))
    ucoord = np.dot(b1, item)
    vcoord = np.dot(b2, item - ucoord * b1)
    #sys.stderr.write("ucoord: %s\n" % str(ucoord))
    #sys.stderr.write("vcoord: %s\n" % str(vcoord))
    sys.stderr.write("u,v = %.3f, %.3f\n" % (ucoord, vcoord))
    rebuild = ucoord * b1 + vcoord * b2
    sys.stderr.write("rebuild: %s\n" % str(rebuild))
    #sys.stderr.write("u, v: %.3f\n" % (ucoord))

def xyz2uv_m(xyz_list, basis_vecs):
    b1, b2 = basis_vecs
    uu = np.array([np.dot(b1, pnt) for pnt in xyz_list])
    vv = np.array([np.dot(b2, pnt) for pnt in xyz_list])
    #sys.stderr.write("uu: %s\n" % str(uu))
    #sys.stderr.write("vv: %s\n" % str(vv))
    return np.array((uu, vv))

def xyz2uv_s(xyz_point, basis_vecs):
    return [np.dot(bb, xyz_point) for bb in basis_vecs]

def total_insep(point, face):
    origin  = face['vertices'][0, :]
    use_vtx = np.array([x-origin for x in face['vertices']])
    use_pnt = point - origin
    #sys.stderr.write("use_vtx: %s\n" % str(use_vtx))
    #sys.stderr.write("use_pnt: %s\n" % str(use_pnt))
    pntu, pntv = xyz2uv_s(use_pnt, face['basis'])
    vtxu, vtxv = xyz2uv_m(use_vtx, face['basis'])
    udiff = vtxu - pntu
    vdiff = vtxv - pntv
    edgetot = np.sum(np.abs(udiff)) + np.sum(np.abs(vdiff))
    return edgetot

#herpderp = po.PolygonFace(prf2['vertices'])

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
##--------------------------------------------------------------------------##
## Plot layout and annotations as sanity check:
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
    pass

# Prism edges:
vx, vy, vz = zip(*prpoly.get_vertices())
for p1,p2 in itt.combinations(range(3), 2):
    ax1.plot([vx[p1], vx[p2]], [vy[p1], vy[p2]], c='b')
    pass

# Draw prism face normals:
arr_len = 50.0
for which in ['face1', 'face2']: #, 'face3']:
    fdata = prpoly.get_face(which)
    arrx1, arry1, _ = fdata['center']
    ndx, ndy, ndz = arr_len * fdata['normal']
    ax1.arrow(arrx1, arry1, ndx, ndy, head_width=5, head_length=5, color='b')
    pass

# -----------------------------------------------------------------------
# Grating vertices:
gcolor = 'k'
for x,y,z in grpoly.get_vertices():
    ax1.scatter(x, y, lw=0, s=30, c=gcolor)

# Grating edges:
for pair in edges_from_vertices(grpoly.get_vertices()):
    vx, vy = pair.T
    ax1.plot(vx, vy, c=gcolor)

ax1.set_xlim(-400, 200)
ax1.set_ylim(-500, 500)

# -----------------------------------------------------------------------
# Routine to connect two xyz points:
def qconnect(xyz1, xyz2, **kwargs):
    x1, y1, z1 = xyz1
    x2, y2, z2 = xyz2
    ax1.plot([x1, x2], [y1, y2], **kwargs)
    return

# Light paths (green):
grkw = {'ls':'--', 'c':'g'}
beam_start = np.array([0.0, -400.0, 0.0])
qconnect(beam_start, f1_isect, **grkw)
qconnect(f1_isect, f2_isect, **grkw)

path3_len = 300.0
path3_end = f2_isect + path3_len * path3
qconnect(f2_isect, path3_end, **grkw)


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
