#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
#    Check the location of hot pixels in NRES images. This simple test will
# help identify readout problems with nres01.
#
# Rob Siverd
# Created:       2017-08-11
# Last modified: 2017-08-13
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.1.5"

## Modules:
import getopt
import signal
import os
import sys
import time
import numpy as np
np.set_printoptions(suppress=True)

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
#try:
#    import astropy.time as astt
#except ImportError:
#    sys.stderr.write("\nError: astropy module not installed!\n"
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

##--------------------------------------------------------------------------##
## Catch interruption cleanly:
def signal_handler(signum, frame):
    sys.stderr.write("\nInterrupted!\n\n")
    sys.exit(1)

signal.signal(signal.SIGINT, signal_handler)

##--------------------------------------------------------------------------##
## Save FITS image with clobber (astropy / pyfits):
#def qsave(iname, idata, header=None, padkeys=1000, **kwargs):
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
#    sys.stderr.write("Writing to '%s' ... " % iname)
#    #if os.path.isfile(iname):
#    #    os.remove(iname)
#    fitsio.write(iname, idata, clobber=True, header=header, **kwargs)
#    sys.stderr.write("done.\n")

##--------------------------------------------------------------------------##
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

## Settings:
debug = False
timer = False
vlevel = 0
prog_name = 'verify-hot-pixels.py'
full_prog = sys.argv[0]
base_prog = os.path.basename(full_prog)
num_todo = 0

## Options:
save_file = None
camera_id = None

##--------------------------------------------------------------------------##
## Argument type-checking:
def is_integer(asdf):
    try:
        int(asdf)
        return True
    except ValueError:
        return False

def is_float(asdf):
    try:
        float(asdf)
        return True
    except ValueError:
        return False

##--------------------------------------------------------------------------##
##*********************     Help and options menu:     *********************##
##--------------------------------------------------------------------------##

## Syntax / how to run:
def usage(stream):
    stream.write("\n"
        + "Usage: %s [options] nres_spectrum.fits\n" % base_prog
        + "Verify the location of hot pixels in NRES spectra.\n"
        + "Version: %s\n" % __version__
        + "\n"
        + "Camera choice:\n"
        + "       --nres01         verify nres01 hot pixels\n"
        + "\n"
        + "Output file (REQUIRED):\n"
        + "   -o, --output=FILE    save results to FILE\n"
        + "\n"
        + "Available options:\n"
        + "       --debug          extra debugging info\n"
        + "   -h, --help           print this page\n"
        + "   -q, --quiet          suppress unnecessary output\n"
        + "   -t, --timer          report program run-time\n"
        + "   -v, --verbose        more status updates\n"
        + "\n")
        #+ "   -n, --numtodo=N     stop after N iterations\n"
        #+ "   -s, --sigcut=N      clip data beyond N*sigma\n"

##--------------------------------------------------------------------------##
##*********************       Parse command line:      *********************##
##--------------------------------------------------------------------------##

## Options:
short_opts = 'o:hqtv' # n:s:
long_opts = ['nres01', 'output=',
                'debug', 'help', 'quiet', 'timer', 'verbose']
# 'numtodo=', 'sigcut='

## GNU-style parsing (with exception handling):
try:
    options, remainder = getopt.gnu_getopt(sys.argv[1:], short_opts, long_opts)
except getopt.GetoptError, err:
    sys.stderr.write("%s\n" % str(err))
    usage(sys.stderr)
    sys.exit(2)

## Handle selected options:
for opt, arg in options:
    # ------------------------------------------------
    if (opt == '--debug'):
        debug = True
        sys.stderr.write(BRED + "\nDebugging output enabled!" + ENDC + "\n")
    # ------------------------------------------------
    #elif ((opt == '-n') or (opt == '--numtodo')):
    #    if not is_integer(arg):
    #        sys.stderr.write("Error!  Non-integer argument: %s\n\n" % arg)
    #        sys.exit(1)
    #    num_todo = int(arg)
    #    if (vlevel >= 0):
    #        msg = "Stopping after %d items." % num_todo
    #        sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    elif (opt == '--nres01'):
        camera_id = 'nres01'
        if (vlevel >= 0):
            msg = "Selected camera: " + camera_id
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    elif ((opt == '-o') or (opt == '--output')):
        save_file = arg
        if (vlevel >= 0):
            msg = "Saving results to: " + arg
            sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    #elif ((opt == '-s') or (opt == '--sigcut')):
    #    if not is_float(arg):
    #        sys.stderr.write("Error!  Non-numeric argument: %s\n\n" % arg)
    #        sys.exit(1)
    #    sigcut = float(arg)
    #    if (vlevel >= 0):
    #        msg = "Using %.2f-sigma outlier threshold." % sigcut
    #        sys.stderr.write(NYELLOW + msg + ENDC + "\n")
    # ------------------------------------------------
    elif ((opt == '-h') or (opt == '--help')):
        usage(sys.stdout)
        sys.exit(0)
    elif ((opt == '-q') or (opt == '--quiet')):
        vlevel -= 1
    elif ((opt == '-t') or (opt == '--timer')):
        timer = True
    elif ((opt == '-v') | (opt == '--verbose')):
        vlevel += 1
        sys.stderr.write(NYELLOW + "Increasing verbosity." + ENDC + "\n")
    # ------------------------------------------------
    else:
        msg = "Unhandled option: %s" % opt
        sys.stderr.write(BRED + "\n" + msg + ENDC + "\n\n")
        sys.exit(1)
    pass

## Verbosity:
if (vlevel >= 1):
    sys.stderr.write("%sVerbosity level: %d%s\n" % (NYELLOW, vlevel, ENDC))

## Full command line if highly verbose:
if (vlevel >= 2):
    sys.stderr.write("%s\nFull command line:%s\n" % (NCYAN, ENDC))
    sys.stderr.write("   %s\n" % sys.argv)

##--------------------------------------------------------------------------##
## Output file is required:
if not save_file:
    sys.stderr.write(BRED + "\nOutput file required!" + ENDC + "\n")
    usage(sys.stderr)
    sys.exit(1)

## Input file is required:
if (len(remainder) < 1):
    sys.stderr.write(BRED + "\nInput file required!" + ENDC + "\n")
    usage(sys.stderr)
    sys.exit(1)

## Check for input file:
data_file = remainder[0]
if not os.path.isfile(data_file):
    msg = "%s error:  file not found: %s" % (prog_name, data_file)
    sys.stderr.write("\n" + BRED + msg + ENDC + "\n")
    sys.exit(1)

##--------------------------------------------------------------------------##
## Known hot pixels (PIXEL COORDINATES [1-indexed]).
hotpix_coords = {}
hotpix_coords['nres01'] = [(3304, 584)]

## Hot pixel detection criterion:
min_sigs = 20.0             # hot pixel must be >= N*sigma above median
sig_buffer = 10.0           # hot pixel must be >= M*sigma above other pixels

## Make sure camera is specified and recognized:
if not camera_id:
    sys.stderr.write("No camera selected!\n")
    usage(sys.stderr)
    sys.exit(1)
if not camera_id in hotpix_coords.keys():
    sys.stderr.write("Unrecognized camera: %s\n" % camera_id)
    usage(sys.stderr)
    sys.exit(1)

##--------------------------------------------------------------------------##
## Shift config:
deltas = np.arange(-1, 2)
xx, yy = np.meshgrid(deltas, deltas)
dx_list = xx.flatten()
dy_list = yy.flatten()

##--------------------------------------------------------------------------##
## Quick FITS I/O:
#data_file = 'image.fits'
#img_vals = pf.getdata(data_file)
#hdr_keys = pf.getheader(data_file)
#img_vals, hdr_keys = pf.getdata(data_file, header=True)
img_vals, hdr_keys = pf.getdata(data_file, header=True, uint=True) # USHORT
#img_vals, hdr_keys = fitsio.read(data_file, header=True)

## Bulk stats:
img_med, img_iqrn = calc_ls_med_IQR(img_vals)

## Begin results:
saving = {}
saving['results'] = [ data_file,   img_med,   img_iqrn]
saving['formats'] = [      '%s',   '%9.2f',    '%9.2f']
saving['columns'] = ['FILENAME', 'PIX_MED', 'PIX_IQRN']

## Check area around each known hot pixel:
shift_list = []
for (xpix, ypix) in hotpix_coords[camera_id]:
    if (vlevel >= 1):
        sys.stderr.write("xpix: %s\n" % str(xpix))
        sys.stderr.write("ypix: %s\n" % str(ypix))

    # select pixels and sort by value:
    values = img_vals[dy_list + (ypix - 1), dx_list + (xpix - 1)]
    order  = np.argsort(values)[::-1]
    values = values[order]
    maxdex = order[0]

    # check for identical values:
    if np.unique(values).size == 1:
        sys.stderr.write("WARNING: identical pixel values!\n")
        sys.stderr.write("    -->  %s\n" % str(values))
        continue

    # stats to ensure confidence:
    pix_med, pix_iqrn = calc_ls_med_IQR(values)
    sigmas = (values - pix_med) / pix_iqrn
    if (vlevel >= 1):
        sys.stderr.write("pix_med:  %9.3f\n" % pix_med)
        sys.stderr.write("pix_iqrn: %9.3f\n" % pix_iqrn)
        sys.stderr.write("sigmas: %s\n" % str(sigmas))
    if (sigmas[0] > min_sigs) and (sigmas[0] - sigmas[1] > sig_buffer):
        delta = (dx_list[maxdex], dy_list[maxdex])
        shift_list.append(delta)
        if (vlevel >= 1):
            sys.stderr.write("delta (x,y):  (%d, %d)\n" % delta)

## Typical shift:
n_shifts = len(shift_list)
if n_shifts > 0:
    avg_shift = np.average(shift_list, axis=0)
else:
    sys.stderr.write("WARNING: shift measurement inconclusive!\n") 
    avg_shift = np.array([0.0, 0.0])

## Update results:
saving['results'].extend([avg_shift[0], avg_shift[1],  n_shifts])
saving['formats'].extend([     '%6.2f',      '%6.2f',     '%2d'])
saving['columns'].extend([    'XSHIFT',     'YSHIFT', 'NSHIFTS'])

##--------------------------------------------------------------------------##
## Create output file with headers (if not present):
if not os.path.isfile(save_file):
    with open(save_file, 'w') as f:
        f.write(' '.join(saving['columns']) + '\n')

## Append results to file:
linefmt = ' '.join(saving['formats']) + '\n'
savetxt = linefmt % tuple(saving['results'])
with open(save_file, 'a') as f:
    f.write(savetxt)
if vlevel >= -1:
    sys.stderr.write(savetxt)
    #f.write("%s %9.2f %9.2f %6.2f %6.2f %2d\n" % 
    #    (data_file, img_med, img_iqrn, avg_shift[0], avg_shift[1], n_shifts))




######################################################################
# CHANGELOG (verify-hot-pixels.py):
#---------------------------------------------------------------------
#
#  2017-08-13:
#     -- Increased __version__ to 0.1.5.
#     -- Output file is now required.
#     -- Removed unused matplotlib import.
#     -- Now also calculate and report median/IQRN pixel value.
#     -- Increased __version__ to 0.1.1.
#     -- Added check for identically-valued pixels (bad image).
#
#  2017-08-11:
#     -- Increased __version__ to 0.1.0.
#     -- First created verify-hot-pixels.py.
#
