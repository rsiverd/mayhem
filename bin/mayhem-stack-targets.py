#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Combine multiple TARGET frames (images) into a single file to reduce noise
# and reject cosmic rays. Output file metadata (dates etc.) are sensibly
# (if not quite correctly) estimated.
#
# Rob Siverd
# Created:       2018-08-07
# Last modified: 2018-08-07
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
import argparse
import os
import sys
import time
import numpy as np

##--------------------------------------------------------------------------##

## Fast FITS I/O:
try:
    import fitsio
except ImportError:
    sys.stderr.write("\nError: fitsio module not found!\n")
    sys.exit(1)

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

## Time conversion:
try:
    import astropy.time as astt
except ImportError:
    sys.stderr.write("\nError: astropy module not found!\n"
           "Please install and try again.\n\n")
    sys.exit(1)

## WCS handling:
#try:
#    from astropy.wcs import WCS
#    import astropy.wcs as awcs
#except ImportError:
#    sys.stderr.write("\nError: astropy module not found!\n"
#           "Please install and try again.\n\n")
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

## Header dividers:
divider_string = \
    "-----------------------------------------------------------------------"
divider_card = pf.Card('COMMENT', divider_string)

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

## Configuration:
min_stack_size        = 3
#require_objskey_match = True

## Keywords that should match among all input images:
hkeys_requiring_match = ['OBJECTS', 'BLKUID', 'TRACKNUM', 'MOLUID', 'EXPTIME']

## Keywords that must differ among all input images (ensures uniqueness):
hkeys_requiring_diffs = ['DATE-OBS', 'ORIGNAME']

##--------------------------------------------------------------------------##
## Parse arguments and run script:
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

if __name__ == '__main__':

    # ------------------------------------------------------------------
    descr_txt = """
    Quick stacking of target frames. Does something halfway-sensible with
    dates, timestamps, and other metadata.

    Version: %s
    """ % __version__
    parser = MyParser(prog=os.path.basename(__file__), description=descr_txt)
    #parser.add_argument('firstpos', help='first positional argument')
    parser.add_argument('-n', '--chunks', default=5, required=False, type=int,
            help='number of chunks to split images into')
    parser.add_argument('-o', '--output_file', default=None, required=True,
            help='Output filename for combined spectrum')
    parser.add_argument('-q', '--quiet', action='count', default=0,
            help='less progress/status reporting')
    parser.add_argument('-v', '--verbose', action='count', default=0,
            help='more progress/status reporting')
    parser.add_argument('--debug', dest='debug', default=False,
            help='Enable extra debugging messages', action='store_true')
    parser.add_argument('input_images', nargs='*',
            help='spectrum images to stack')
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
    # ------------------------------------------------------------------

    context = parser.parse_args()
    #context.vlevel = context.verbose - context.quiet
    context.vlevel = 99 if context.debug else (context.verbose-context.quiet)

    context.input_images.sort()     # enforce chronological order

##--------------------------------------------------------------------------##
## Input image sanity checks:
#def file_contains_extname(filename, extname):
#    hdu_list = fitsio.FITS(filename)

##--------------------------------------------------------------------------##
## Assert minimum stack size:
stack_size = len(context.input_images)
if (stack_size < min_stack_size):
    sys.stderr.write("\nError: minimum stack size is %d images.\n\n"
            % min_stack_size)
    sys.exit(1)

## Ensure all images are present:
for item in context.input_images:
    if not os.path.isfile(item):
        sys.stderr.write("\nError: file not found: %s\n\n" % item)
        sys.exit(1)

## Ensure all images contain a SPECTRUM extension and collect headers:
spec_HDUs = []
for item in context.input_images:
    # Look for SPECTRUM extension in each image:
    try:
        keep_HDUs = [x for x in fitsio.FITS(item) \
                            if (x.get_extname() == 'SPECTRUM')]
    except:
        sys.stderr.write("\nError: failed to load FITS file: %s\n" % item)
        sys.exit(1)

    # Each file must have ONE extension named SPECTRUM:
    if (len(keep_HDUs) != 1):
        sys.stderr.write("Detected bogus content in %s ...\n" % item)
        sys.exit(1)

    # Collect header and dimensions:
    spec_HDUs.append(keep_HDUs[0])
    pass

## Verify that dimensions match:
spec_dims = [x.get_dims() for x in spec_HDUs]
if not all([(x == spec_dims[0]) for x in spec_dims]):
    sys.stderr.write("Dimensions are not identical!\n")
    sys.exit(1)
spec_dims = spec_dims[0]

## Only 2-D images are supported at the moment:
if (len(spec_dims) != 2):
    sys.stderr.write("Error: currently only 2-D images are supported!\n\n")
    sys.exit(1)

##--------------------------------------------------------------------------##
## Extract headers from each SPECTRUM HDU:
spec_hdrs = [x.read_header() for x in spec_HDUs]

## Require exact matches for several keywords:
for hkey in hkeys_requiring_match:
    meta_vals = [x[hkey] for x in spec_hdrs]
    if (len(set(meta_vals)) != 1):
        sys.stderr.write("\nError: %s keywords do not match:\n" % hkey)
        for things in zip(context.input_images, meta_vals):
            sys.stderr.write("%s --> %s\n" % things)
        sys.stderr.write("\n") 
        sys.exit(1)

## Also require differences for certain keywords:
for hkey in hkeys_requiring_diffs:
    meta_vals = [x[hkey] for x in spec_hdrs]
    if (len(set(meta_vals)) != len(meta_vals)):
        sys.stderr.write("\nError: %s keywords are not unique:\n" % hkey)
        for things in zip(context.input_images, meta_vals):
            sys.stderr.write("%s --> %s\n" % things)
        sys.stderr.write("\n") 
        sys.exit(1)

## Could verify that DATE-OBS ordering matches filename ordering ...

##--------------------------------------------------------------------------##



##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
## Image midpoints:
obs_dates = astt.Time([x['DATE-OBS'] for x in spec_hdrs], 
        scale='utc', format='isot')
exp_times = astt.TimeDelta([x['EXPTIME'] for x in spec_hdrs], format='sec')
exp_stops = obs_dates + exp_times

## Effective stacked exposure runs from start of first thru end of last:
spec_tstart   = obs_dates.min()
spec_tstop    = exp_stops.max()
total_exp_sec = np.sum(exp_times).sec
total_dur_sec = (spec_tstop - spec_tstart).sec
down_time_sec = total_dur_sec - total_exp_sec

## Start output header using data from first spectrum:
first_spec_file = context.input_images[obs_dates.argmin()]
final_spec_file = context.input_images[obs_dates.argmax()]
_, first_header = pf.getdata(first_spec_file, extname='SPECTRUM', header=True)
_, final_header = pf.getdata(final_spec_file, extname='SPECTRUM', header=True)

output_header = first_header.copy(strip=True)
output_header['UTSTOP']  = final_header['UTSTOP']
output_header['EXPTIME'] = total_dur_sec    # for UTSTART/UTSTOP consistency

##--------------------------------------------------------------------------##
## Real / special data goes in a separate header section:
output_header.append(divider_card, bottom=True)
output_header.append(pf.Card("EXPTRUE", total_exp_sec, 
                "[s] true sum of input image exposure times"), bottom=True)

## Record names of all input images:
orig_names = [x['ORIGNAME'] for x in spec_hdrs]
for i,oname in enumerate(orig_names, 1):
    new_card = pf.Card("IMSRC%03d" % i, oname, "input image %d" % i)
    output_header.append(new_card, bottom=True)

output_header.append(divider_card, bottom=True)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

def slice_from_arange(data, origin):
    return slice(data[0] - origin, data[-1] + 1 - origin, None)

## Median-combine input images to create output image:
tik = time.time()
r_origin = 1
row_list = np.arange(spec_dims[0]) + r_origin
n_chunks = context.chunks

med_chunks = []
for i,some_rows in enumerate(np.array_split(row_list, n_chunks), 1):
    sys.stderr.write("\rStacking chunk %d of %d ... " % (i, n_chunks))
    section = slice_from_arange(some_rows, r_origin)
    med_chunks.append(np.median([x[section, :] for x in spec_HDUs], axis=0))
    pass
med_stack = np.concatenate(med_chunks).astype('float32')
tok = time.time()
sys.stderr.write("done.\nTotal stacking time: %7.3f sec (%d chunks)\n"
        % (tok - tik, n_chunks))

## Save result:
qsave(context.output_file, med_stack, header=output_header)




######################################################################
# CHANGELOG (mayhem-stack-targets.py):
#---------------------------------------------------------------------
#
#  2018-08-07:
#     -- Increased __version__ to 0.1.0.
#     -- First created mayhem-stack-targets.py.
#
