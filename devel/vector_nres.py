#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Make a vector representation of the spectrograph (approximate).
#
# Rob Siverd
# Created:       2019-02-19
# Last modified: 2019-04-01
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.3.0"

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

## Colorize-by-wavelength:
import wavelength_colors
reload(wavelength_colors)

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
dppgp = spectrograph_optics.DoublePassPrismGratingPrism()

##--------------------------------------------------------------------------##
## Spectroscopic orders and a lookup routine:
useful_orders = 51 + np.arange(69)
#useful_orders = 54 + np.arange(66)
center_wlen_um = dppgp.calc_central_wlen_um(useful_orders.astype('float'))
spec_order_FSR = center_wlen_um / useful_orders.astype('float')
spec_order_wlmin = center_wlen_um - spec_order_FSR
spec_order_wlmax = center_wlen_um + spec_order_FSR

def argnear(vec, val):
    return (np.abs(vec - val)).argmin()

def wlen2orders(wlen_um):
    which = (spec_order_wlmin <= wlen_um) & (wlen_um <= spec_order_wlmax)
    return useful_orders[which]

def wlen2order(wlen_um):
    return useful_orders[argnear(center_wlen_um, wlen_um)]

## Select test wavelengths from each order to test raytracing:
rtwl_samples = {}
#samp_per_order = 11
samp_per_order = 21
#FSRs_per_order = 1.3
FSRs_per_order = 1.0
sample_offsets = FSRs_per_order * np.linspace(-0.5, 0.5, samp_per_order)
for nord,wlcen in zip(useful_orders, center_wlen_um):
    this_FSR = wlcen / float(nord)
    rtwl_samples[nord] = wlcen + this_FSR * sample_offsets


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
def veclen(vector):
    return np.sqrt(np.sum(vector**2))

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

prpoly = po.IsosPrismPolyhedron(apex_angle_deg, short_edge_mm, height_mm,
                                            nres_prism_glass)
prpoly.zrotate(np.radians(-90.0))
prpoly.zrotate(np.radians(prism_turn_deg))

#ft = prpoly.get_face('top')

# for experimentation:
#f1 = prpoly.get_face('face1')


##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Grating:
gr_length_mm = 320.0
gr_width_mm  = 165.0
gr_height_mm =  50.5
gr_ruling_lmm = 41.59                   # lines per mm
gr_spacing_um = 1e3 / gr_ruling_lmm     # groove spacing in microns
grating_turn_deg = 44.827
grating_tilt_deg = 13.786

grpoly = po.GratingPolyhedron(gr_width_mm, gr_length_mm, gr_height_mm,
                                gr_ruling_lmm)
grpoly.xrotate(np.radians(-grating_tilt_deg))
grpoly.zrotate(np.radians(grating_turn_deg))
grpoly.shift_xyz(-165.0, 235.0, 0.0)

## Adjust grating so center of bottom face lies on z=0:
gr_bot_ctr = np.copy(grpoly.get_face('bot')['center'])
#sys.stderr.write("grating bottom-center: %s\n" % str(gr_bot_ctr))
grpoly.shift_vec(gr_bot_ctr * np.array([0.0, 0.0, -1.0]))

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Center prism face1 on origin (move grating and prism in tandem):
nudge = -1.0 * np.copy(prpoly.get_face('face1')['center'])
prpoly.shift_vec(nudge)
grpoly.shift_vec(nudge)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Keep a list of traced vertexes for plotting:
#light_path = []

## Input beam direction:
input_turn_deg = 5.0
#input_turn_deg = 7.0
input_uptilt_deg = 0.0
input_phi = np.radians(90.0 + input_turn_deg)
input_theta = np.radians(90.0 - input_uptilt_deg)
v_initial = np.array([np.sin(input_theta) * np.cos(input_phi),
                      np.sin(input_theta) * np.sin(input_phi),
                      np.cos(input_theta)])

## An initial input beam (center fiber):
wl_initial = 0.80                        # ray wavelength (microns)
prvtx_xmax = prpoly.get_vertices('top')[:, 0].max()     # prism +X extremum
fiber_xpos = 0.85 * prvtx_xmax
ferrule_fsep = 0.250                                    # mm between fibers

## Fibers 0+2 are offset from fiber1 in X-Y plane:
frac_offset = np.array([-np.sin(input_phi), np.cos(input_phi), 0.0])
#fib2_offset = np.array([ , , 0.0])

fiber_exit = np.array([fiber_xpos, -400.0, 0.0])

## Starting positions for all three fibers:
fiber_ends = [fiber_exit.copy() + ferrule_fsep * frac_offset,
              fiber_exit.copy(),
              fiber_exit.copy() - ferrule_fsep * frac_offset]
#light_path.append(fiber_exit)
#keep_going = True

## Make CCD polyhedron:
ccdpoly = po.CameraPolyhedron(200.0, 0.1, 200.0)
ccdpoly.shift_xyz(0.0, -450.0, 0.0)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Initialize end-to-end raytracer with polygons:
spectr = rt.E2ERT(prpoly, grpoly, ccdpoly, vlevel=1)
owhich = wlen2order(wl_initial)
test_path = spectr.follow(fiber_exit, v_initial, wl_initial, owhich)
test_path = spectr.follow(fiber_exit, v_initial, wl_initial, owhich+1)
#test_path = spectr.follow(fiber_exit, v_initial, 0.799, 59)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Raytrace numerous orders:
sys.stderr.write("Tracing spectrum ...\n")
path_details = {}
spec_mod_pts = {}
trace_fibers = [0, 1, 2]
#for ff in range(3):
#    path_details[ff] = {}
#    t_spec_mod_pts = {}
#all_path_details = {}
#all_spec_mod_pts = {}
#for ff in trace_fibers:
#    path_details = {}
#    spec_mod_pts = {}
for i,(sord, wsamp) in enumerate(rtwl_samples.items(), 1):
    sys.stderr.write("\rsord: %3d (order %d of %d) ... " 
            % (sord, i, useful_orders.size))
    path_details[sord] = {}
    spec_mod_pts[sord] = {}
    for ff in trace_fibers:
        launch = fiber_ends[ff]
        paths = [spectr.follow(launch, v_initial, ww, sord) for ww in wsamp]
        path_details[sord][ff] = paths
        endpoints = [x[-1] for x in paths]
        kept = [(ww, cx, cz) for ww,((cx,cy,cz),tt) \
                            in zip(wsamp, endpoints) if (tt is None)]
        #spec_mod_pts[ff][sord] = np.array(kept)
        spec_mod_pts[sord][ff] = np.array(kept)
    
    #all_path_details[ff] = path_details
    #all_spec_mod_pts[ff] = spec_mod_pts
sys.stderr.write("done.\n")

#pickone = 0
#path_details = all_path_details[pickone]
#spec_mod_pts = all_spec_mod_pts[pickone]

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Edges from vertices:
def edges_from_vertices(vtx_array):
    nhalf = int(vtx_array.shape[0] / 2)
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
fig_dims = (9, 10)
fig = plt.figure(1, figsize=fig_dims)
#plt.gcf().clf()
fig.clf()
#fig, axs = plt.subplots(2, 2, sharex=True, figsize=fig_dims, num=1)
# sharex='col' | sharex='row'
#fig.frameon = False # disable figure frame drawing
#fig.subplots_adjust(left=0.07, right=0.95)
#ax1 = plt.subplot(gs[0, 0])
#ax1 = fig.add_subplot(111, aspect='equal')
ax1 = fig.add_axes([0.05, 0.05, 0.93, 0.93], aspect='equal')
#ax1.patch.set_facecolor((0.8, 0.8, 0.8))
ax1.grid(True)
#ax1.axis('off')

fig2 = plt.figure(2, figsize=(6, 6))
fig2.clf()
ax2 = fig2.add_subplot(111, aspect='equal')
ax2.grid(True)
sskw = {'lw':0, 's':15}
flip_plot = True
everything = []
for onum,odata in spec_mod_pts.items():
    for ff in trace_fibers:
        wlen, xcoo, ycoo = odata[ff].T
        everything.extend([x for x in odata.values()])
        #ax2.scatter(xcoo, ycoo, **sskw)
        wlnm = 1e3 * wlen.mean()
        color = wavelength_colors.wave2rgb(wlnm)
        if flip_plot:
            ax2.plot(ycoo, -xcoo, color=color)
        else:
            ax2.plot(xcoo, ycoo, color=color)

## Combined data set:
#mod_wlen, mod_xcoo, mod_ycoo = \
#        np.concatenate([x for x in spec_mod_pts.values()]).T
mod_wlen, mod_xcoo, mod_ycoo = np.concatenate(everything).T
if flip_plot:
    mod_xcoo, mod_ycoo = mod_ycoo, -mod_xcoo
xmin, xmax = mod_xcoo.min(), mod_xcoo.max()
ymin, ymax = mod_ycoo.min(), mod_ycoo.max()
span = 1.05 * max(xmax - xmin, ymax - ymin)
span = 1.05 * (ymax - ymin)
xmid, ymid = 0.50 * (xmin + xmax), 0.5 * (ymin + ymax)
#xavg = np.average(mod_xcoo)
#yavg = np.average(mod_ycoo)
#ax2.set_xlim(xavg - 65., xavg + 65.)
#ax2.set_ylim(yavg - 53., yavg + 65.)
ax2.set_xlim(xmid - 0.5*span, xmid + 0.5*span)
ax2.set_ylim(ymid - 0.5*span, ymid + 0.5*span)
fig2.tight_layout()

## Disable axis offsets:
#ax1.xaxis.get_major_formatter().set_useOffset(False)
#ax1.yaxis.get_major_formatter().set_useOffset(False)

# -----------------------------------------------------------------------
# Mark CCD position:
#ax1.scatter(ccd_xseg, ccd_yseg, lw=0, s=30, c='m')
#ax1.plot(ccd_xseg, ccd_yseg, c='m')

# CCD vertices:
for x,y,z in ccdpoly.get_vertices():
    ax1.scatter(x, y, lw=0, s=30, c='m')
    pass

## CCD edges:
vx, vy, vz = zip(*ccdpoly.get_vertices())
for p1,p2 in itt.combinations(range(3), 2):
    ax1.plot([vx[p1], vx[p2]], [vy[p1], vy[p2]], c='m')
    pass


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

## Draw prism face normals:
draw_normals = False
if draw_normals:
    arr_len = 50.0
    for which in ['face1', 'face2']: #, 'face3']:
        fdata = prpoly.get_face(which)
        arrx1, arry1, _ = fdata['center']
        ndx, ndy, ndz = arr_len * fdata['normal']
        ax1.arrow(arrx1, arry1, ndx, ndy,
                head_width=5, head_length=5, color='b')
        pass

## -----------------------------------------------------------------------
## Grating vertices:
gcolor = 'k'
for x,y,z in grpoly.get_vertices():
    ax1.scatter(x, y, lw=0, s=30, c=gcolor)

## Grating edges:
for pair in edges_from_vertices(grpoly.get_vertices()):
    vx, vy = pair.T
    ax1.plot(vx, vy, c=gcolor)

ax1.set_xlim(-500, 300)
ax1.set_ylim(-500, 500)
ax1.set_xlabel("X axis (mm)")
ax1.set_ylabel("Y axis (mm)")

## -----------------------------------------------------------------------
## Routine to connect two xyz points:
def qconnect(xyz1, xyz2, **kwargs):
    x1, y1, z1 = xyz1
    x2, y2, z2 = xyz2
    ax1.plot([x1, x2], [y1, y2], **kwargs)
    return

### Light paths (green):
#grkw = {'ls':'--', 'c':'g'}
#for pstart,pstop in zip(light_path[:-1], light_path[1:]):
#    qconnect(pstart, pstop, **grkw)

## Make plottable series of points:
test_verts = [x[0] for x in test_path]

## Add another point to show final trajectory if destination not reached:
if isinstance(test_path[-1][1], np.ndarray):
    draw_dist = 500.
    last_posn, last_traj = test_path[-1]
    after_posn = last_posn + draw_dist * last_traj
    test_verts.append(after_posn)

## Light paths (green):
grkw = {'ls':'--', 'c':'g'}
for pstart,pstop in zip(test_verts[:-1], test_verts[1:]):
    qconnect(pstart, pstop, **grkw)

#beam_start = np.array([0.0, -400.0, 0.0])
#qconnect(beam_start, f1_isect, **grkw)
#qconnect(f1_isect, f2_isect, **grkw)
#
##path3_len = 300.0
##path3_end = f2_isect + path3_len * path3
#qconnect(f2_isect, gr_isect, **grkw)
#
#qconnect(gr_isect, next_verts[58.0], **grkw)
##path4_len = 300.0


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
