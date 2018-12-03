#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Look for bloomed spectral lines to identify with Argon wavelengths.
#
# Rob Siverd
# Created:       2018-12-03
# Last modified: 2018-12-03
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
import argparse
import os
import sys
import time
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
#import theil_sen as ts
#import window_filter as wf
#import itertools as itt

import kelt_satur_tools
reload(kelt_satur_tools)
kstbc = kelt_satur_tools.BloomCentroid()

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
try:
    import astropy.io.fits as pf
except ImportError:
    try:
       import pyfits as pf
    except ImportError:
        sys.stderr.write("\nError!  No FITS I/O module found!\n"
               "Install either astropy.io.fits or pyfits and try again!\n\n")
        sys.exit(1)

## WCS handling:
#try:
#    from astropy.wcs import WCS
#    import astropy.wcs as awcs
#except ImportError:
#    sys.stderr.write("\nError: astropy module not found!\n"
#           "Please install and try again.\n\n")
#    sys.exit(1)

## Star extraction:
#try:
#    import easy_sep
#    reload(easy_sep)
#except ImportError:
#    sys.stderr.write("Error: easy_sep module not found!\n\n")
#    sys.exit(1)
#pse = easy_sep.EasySEP()

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
#halfdiv = "----------------------------------------"
#fulldiv = halfdiv + halfdiv
halfdiv = '-' * 40
fulldiv = '-' * 80

##--------------------------------------------------------------------------##
## Save FITS image with clobber (astropy / pyfits):
def qsave(iname, idata, header=None, padkeys=1000, **kwargs):
    this_func = sys._getframe().f_code.co_name
    sys.stderr.write("Writing to '%s' ... " % iname)
    if header:
        while (len(header) < padkeys):
            header.append() # pad header
    if os.path.isfile(iname):
        os.remove(iname)
    pf.writeto(iname, idata, header=header, **kwargs)
    sys.stderr.write("done.\n")

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
    Look for bloomed Argon lines in an NRES spectrum.
    
    Version: %s
    """ % __version__
    parser = MyParser(prog=os.path.basename(__file__), description=descr_txt,
                          formatter_class=argparse.RawTextHelpFormatter)
    # ------------------------------------------------------------------
    parser.set_defaults(rthresh=0.9)
    # ------------------------------------------------------------------
    parser.add_argument('spec_file', help='NRES spectrum with ThAr lines')
    parser.add_argument('-r', '--rthresh', required=False, type=float,
            help='use fraction of image max as saturation threshold')
    parser.add_argument('-t', '--sthresh', required=False, type=float,
            help='use absolute saturation threshold')
    parser.add_argument('-o', '--output', default=None, required=False,
            help='Output FITS filename')
    #parser.add_argument('-d', '--dayshift', required=False, default=0,
    #        help='Switch between days (1=tom, 0=today, -1=yest', type=int)
    #parser.add_argument('-e', '--encl', nargs=1, required=False,
    #        help='Encl to make URL for', choices=all_encls, default=all_encls)
    #parser.add_argument('-s', '--site', nargs=1, required=False,
    #        help='Site to make URL for', choices=all_sites, default=all_sites)
    parser.add_argument('-q', '--quiet', action='count', default=0,
            help='less progress/status reporting')
    parser.add_argument('-v', '--verbose', action='count', default=0,
            help='more progress/status reporting')
    parser.add_argument('--debug', dest='debug', default=False,
            help='Enable extra debugging messages', action='store_true')
    #parser.add_argument('remainder', help='other stuff', nargs='*')
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    #ofgroup = parser.add_argument_group('Output format')
    #fmtparse = ofgroup.add_mutually_exclusive_group()
    #fmtparse.add_argument('--python', required=False, dest='output_mode',
    #        help='Return Python dictionary with results [default]',
    #        default='pydict', action='store_const', const='pydict')
    #bash_var = 'ARRAY_NAME'
    #bash_msg = 'output Bash code snippet (use with eval) to declare '
    #bash_msg += 'an associative array %s containing results' % bash_var
    #fmtparse.add_argument('--bash', required=False, default=None,
    #        help=bash_msg, dest='bash_array', metavar=bash_var)
    #fmtparse.set_defaults(output_mode='pydict')
    # ------------------------------------------------------------------
    # Miscellany:
    #miscgroup = parser.add_argument_group('Miscellany')
    #miscgroup.add_argument('-q', '--quiet', action='count', default=0,
    #        help='less progress/status reporting')
    #miscgroup.add_argument('-v', '--verbose', action='count', default=0,
    #        help='more progress/status reporting')
    #miscgroup.add_argument('--debug', dest='debug', default=False,
    #        help='Enable extra debugging messages', action='store_true')
    # ------------------------------------------------------------------

    context = parser.parse_args()
    context.vlevel = 99 if context.debug else (context.verbose-context.quiet)

##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
## Quick ASCII I/O:
#data_file = 'data.txt'
#all_data = np.genfromtxt(data_file, dtype=None)
#all_data = np.genfromtxt(data_file, dtype=None, names=True, autostrip=True)
#all_data = np.genfromtxt(data_file, dtype=None, names=True, autostrip=True,
#                 delimiter='|', comments='%0%0%0%0')
#                 loose=True, invalid_raise=False)
#all_data = aia.read(data_file)
#all_data = pd.read_csv(data_file)
#all_data = pd.read_table(data_file, delim_whitespace=True)
#all_data = pd.read_table(data_file, sep='|')
#fields = all_data.dtype.names
#if not fields:
#    x = all_data[:, 0]
#    y = all_data[:, 1]
#else:
#    x = all_data[fields[0]]
#    y = all_data[fields[1]]

##--------------------------------------------------------------------------##
## Timestamp modification:
#def time_warp(jdutc, jd_offset, scale):
#    return (jdutc - jd_offset) * scale

## Self-consistent time-modification for plotting:
#tfudge = partial(time_warp, jd_offset=tstart.jd, scale=24.0)    # relative hrs
#tfudge = partial(time_warp, jd_offset=tstart.jd, scale=1440.0)  # relative min

##--------------------------------------------------------------------------##
## Quick FITS I/O:
#data_file = 'image.fits'
#img_vals = pf.getdata(data_file)
#hdr_keys = pf.getheader(data_file)
img_vals, hdr_keys = pf.getdata(context.spec_file, header=True)
#img_vals, hdr_keys = pf.getdata(data_file, header=True, uint=True) # USHORT
#img_vals, hdr_keys = fitsio.read(data_file, header=True)

## Working copy of image with some tweaks:
wrk_vals = img_vals.astype('float32')
wrk_vals[wrk_vals < 0.0] = 0.0
wrk_vals = np.sqrt(wrk_vals)

#maxval = img_vals.max()
sat_thresh = context.rthresh * wrk_vals.max()
satur_mask = (wrk_vals >= sat_thresh).nonzero()

## Save masked image copy:
if context.output:
    tdata = wrk_vals.astype('float32')
    tdata[satur_mask] = np.nan
    qsave(context.output, tdata, overwrite=True)
    del tdata

sys.exit(0)
## Try existing tools:
pix_origin = 1
kstbc.create_sat_mask(wrk_vals, sat_thresh)
ccd_xx, ccd_yy, ccd_npix = kstbc.find_stars(minpix=5, origin=pix_origin)


## Mark up an image:
tdata = wrk_vals.astype('float32')
yrel, xrel = np.mgrid[-1:2, -1:2]
yrel = yrel.flatten()
xrel = xrel.flatten()
ixpos = np.int_(np.floor(ccd_xx))
iypos = np.int_(np.floor(ccd_yy))
#for xmid,ymid in zip(ixpos, iypos):
for xmid,ymid in zip(ccd_xx, ccd_yy):
    ycoo = int(np.floor(ymid)) + yrel
    xcoo = int(np.floor(xmid)) + xrel
    tdata[ycoo, xcoo] = np.nan
if context.output:
    qsave(context.output, tdata, overwrite=True)

## Star extraction:
#pse.set_image(img_vals, gain=3.6)
#objlist = pse.analyze(sigthresh=5.0)

##--------------------------------------------------------------------------##
## Solve prep:
#ny, nx = img_vals.shape
#x_list = (0.5 + np.arange(nx)) / nx - 0.5            # relative (centered)
#y_list = (0.5 + np.arange(ny)) / ny - 0.5            # relative (centered)
#xx, yy = np.meshgrid(x_list, y_list)                 # relative (centered)
#xx, yy = np.meshgrid(nx*x_list, ny*y_list)           # absolute (centered)
#xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))   # absolute
#yy, xx = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij') # absolute
#yy, xx = np.nonzero(np.ones_like(img_vals))          # absolute
#yy, xx = np.mgrid[0:ny,   0:nx].astype('uint16')     # absolute (array)
#yy, xx = np.mgrid[1:ny+1, 1:nx+1].astype('uint16')   # absolute (pixel)

## 1-D vectors:
#x_pix, y_pix, ivals = xx.flatten(), yy.flatten(), img_vals.flatten()
#w_vec = np.ones_like(ivals)            # start with uniform weights
#design_matrix = np.column_stack((np.ones(x_pix.size), x_pix, y_pix))

## Image fitting (statsmodels etc.):
#data = sm.datasets.stackloss.load()
#ols_res = sm.OLS(ivals, design_matrix).fit()
#rlm_res = sm.RLM(ivals, design_matrix).fit()
#rlm_model = sm.RLM(ivals, design_matrix, M=sm.robust.norms.HuberT())
#rlm_res = rlm_model.fit()
#data = pd.DataFrame({'xpix':x_pix, 'ypix':y_pix})
#rlm_model = sm.RLM.from_formula("ivals ~ xpix + ypix", data)

##--------------------------------------------------------------------------##
## KDE:
#kde_pnts, kde_vals = mk.go(data_vec)

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
# CHANGELOG (find_bloomed_argon.py):
#---------------------------------------------------------------------
#
#  2018-12-03:
#     -- Increased __version__ to 0.0.1.
#     -- First created find_bloomed_argon.py.
#
