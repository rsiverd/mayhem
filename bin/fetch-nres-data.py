#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Retrieve NRES spectra from the LCO archive.
#
# Rob Siverd
# Created:       2018-09-05
# Last modified: 2018-11-27
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.4.0"

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
import signal
import shutil
import errno
import os
import sys
import copy
import time
import numpy as np
import requests
import traceback

##--------------------------------------------------------------------------##
## Catch interruption cleanly:
#def signal_handler(signum, frame):
#    sys.stderr.write("\nInterrupted!\n\n")
#    sys.exit(1)
#
#signal.signal(signal.SIGINT, signal_handler)

##--------------------------------------------------------------------------##
## LCO site information:
try:
    import lco_site_info
    reload(lco_site_info)
except ImportError:
    sys.stderr.write("\nError: lco_site_info not found! Install and retry.\n")
    sys.exit(1)
lsi = lco_site_info.info
nres_sites = [x for x in lsi if lsi[x]['nres_spec']]
nres_cam_map = dict([(x, lsi[x]['nres_spec'][0]['cam']) for x in nres_sites])

## LCO API interaction:
obs_api_token = 'aa63829a45739210d124c6ec1c7edb0d62e5e860'
#archive_token = '2de3ffb5590fe7411e426d1d28d04376e77d05d1' # OLD, lolwut?
archive_token = 'f05b4a4f098a009b3008460b2949339a398662ed'
import lco_api_tools
reload(lco_api_tools)
lcofrm = lco_api_tools.LCO_Frames(archive_token)
#lcoapi = lco_api_tools.LCO_API(obs_api_token)

## Fancy/pretty file downloading:
#_fancy_downloading = False
try:
    import downloading
    reload(downloading)
    #_fancy_downloading = True
    fdl = downloading.Downloader()
except ImportError:
    sys.stderr.write("\nRequired 'downloading' module not found!\n"
                    "Progress messages and rate limiting will be disabled.\n")
    sys.exit(1)

## Time conversion:
try:
    import astropy.time as astt
except ImportError:
    sys.stderr.write("\nError: astropy module not found!\n"
           "Please install and try again.\n\n")
    sys.exit(1)

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
## Recursive directory creation (emulates `mkdir -p`):
#def mkdir_p(path):
#    try:
#        os.makedirs(path)
#    except OSError as exc:  # Python >2.5
#        if exc.errno == errno.EEXIST and os.path.isdir(path):
#            pass
#        else:
#            raise

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## ISO date format checker:
def date_is_iso(datestr):
    try:
        when = astt.Time(datestr, scale='utc', format='iso')
        return True
    except:
        return False

## Promote to full ISO date (including hh:mm:ss):
def promote_full_iso(datestr):
    return astt.Time(datestr).iso

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
    #allowed_sites = ['lsc', 'elp', 'cpt', 'tlv']

    # ------------------------------------------------------------------
    descr_txt = """
    Retrieve NRES data from LCO archive.
    
    Version: %s
    """ % __version__
    parser = MyParser(prog=os.path.basename(__file__), description=descr_txt,
                          formatter_class=argparse.RawTextHelpFormatter)
    #parser.add_argument('firstpos', help='first positional argument')
    parser.add_argument('-q', '--quiet', action='count', default=0,
            help='less progress/status reporting')
    parser.add_argument('-v', '--verbose', action='count', default=0,
            help='more progress/status reporting')
    parser.add_argument('--cronjob', default=False, action='store_true',
            help='use cron-friendly output messages')
    parser.add_argument('--debug', dest='debug', default=False,
            help='Enable extra debugging messages', action='store_true')
    #parser.add_argument('remainder', help='other stuff', nargs='*')
    #parser.set_defaults(thing='value')
    parser.set_defaults(data_rlevel=None, do_download=True)
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    typegroup = parser.add_argument_group('NRES Data Types')
    calib_sci = typegroup.add_mutually_exclusive_group()
    calib_sci.add_argument('--cal_only', required=False, default=False,
            dest='cal_only', action='store_true',
            help='only fetch NRES calibration data')
    calib_sci.add_argument('--sci_only', required=False, default=False,
            dest='sci_only', action='store_true',
            help='only fetch NRES science data')
    rdx_level = typegroup.add_mutually_exclusive_group()
    rdx_level.add_argument('--raw', required=False, dest='data_rlevel',
            action='store_const', const=0,
            help='only retrieve raw data (RLEVEL 00)')
    rdx_level.add_argument('--proc', required=False, dest='data_rlevel',
            action='store_const', const=91,
            help='only retrieve processed data (RLEVEL 91)')

    timegroup = parser.add_argument_group('Data Time Range')
    timegroup.add_argument('-N', '--ndays', required=False, default=0.0,
            type=float, help='restrict search to data from past N days')
    timegroup.add_argument('--start', required=False, default=None, type=str,
            help='start of time window (YYYY-MM-DD [hh:mm:ss])')
    timegroup.add_argument('--end', required=False, default=None, type=str,
            help='end of time window (YYYY-MM-DD [hh:mm:ss])')


    ordergroup = parser.add_argument_group('Download Order')
    ordergroup = ordergroup.add_mutually_exclusive_group()
    ordergroup.add_argument('--new-first', required=False, default=False,
            action='store_false', dest='oldest_first',
            help='download in reverse chrono order (newest files first)')
    ordergroup.add_argument('--old-first', required=False, default=False,
            action='store_true', dest='oldest_first',
            help='download in chrono order (oldest files first)')

    datagroup = parser.add_argument_group('Data Type and Site Choice')
    #datagroup.add_argument('--token', required=True, default=None,
    #        help='LCO API/Archive authentication token')

   #authgroup = parser.add_argument_group('Accounts and Authentication')
   #authgroup.add_argument('--obs_api_token', required=False,
   #        default=obs_api_token,
   #        help='LCO observing API (Valhalla) authentication token')
   #authgroup.add_argument('--archive_token', required=False,
   #        default=archive_token,
   #        help='LCO Archive authentication token')

    filegroup = parser.add_argument_group('Local File I/O')
    filegroup.add_argument('-o', '--save_root', required=True, default=None,
            help='data storage root (i.e., /archive/engineering)')

    sitegroup = parser.add_argument_group('Site Selection')
    sitegroup = sitegroup.add_mutually_exclusive_group()
    #sitegroup.add_argument('--skip-site', required=False, dest='skip_sites',
    #        choices=nres_sites, action='append', default=[],
    #        help='do not download files from SITE')
    sitegroup.add_argument('-a', '--all-sites', required=False, default=False,
            dest='all_sites', action='store_true',
            help='retrieve data for all NRES sites')
    sitegroup.add_argument('-s', '--site', required=False, dest='one_site',
            choices=nres_sites, default=None,
            help='retrieve data for specific site')
            #help='add site to retrieval list (works multiple times)')

    miscgroup = parser.add_argument_group('Miscellany')
    miscgroup.add_argument('--max_depth', required=False, default=0,
            type=int, help='max search iterations')
    miscgroup.add_argument('--max_files', required=False, default=500,
            type=int, help='max files per search iteration')
    # ------------------------------------------------------------------

    context = parser.parse_args()
    context.vlevel = 99 if context.debug else (context.verbose-context.quiet)

    #if context.all_sites:
    #    context.want_sites = copy.deepcopy(nres_sites)

##--------------------------------------------------------------------------##
## Disable progress reporting in 'cronjob' mode:
if context.cronjob:
    #_fancy_downloading = False
    fdl.noprogress()


### Double-check valid site selection:
#if context.one_site and context.skipped_sites:
#    sys.stderr.write("Error: can't skip AND choose sits ...\n")
#    sys.exit(1)

### Some site specification is required:
#if not context.want_sites:
#    sys.stderr.write(BRED + "\nError: no site(s) specified!\n" + ENDC)
#    sys.stderr.write("\nTry `%s --help` for more information.\n\n"
#            % os.path.basename(__file__))
#    sys.exit(1)

## Halt in case of negative NDAYS lookback:
if (context.ndays < 0.0):
    sys.stderr.write("Error: NDAYS must be positive!\n")
    sys.exit(1)

## Halt in case of non-ISO start/end:
for datestr in (context.start, context.end):
    if datestr and not date_is_iso(datestr):
        sys.stderr.write("\nError: non-ISO date string: '%s'\n" % datestr)
        sys.exit(1)

## Time range precedence warning:
if (context.ndays > 0.0) and context.start:
    sys.stderr.write("WARNING: option --ndays supercedes --start\n")

##--------------------------------------------------------------------------##
##------------------ Sanity Checks and Option Consequences  ----------------##
##--------------------------------------------------------------------------##

## Slim down options in debug mode:
if context.debug:
    sys.stderr.write("Slimming down frames list for debug mode!\n")
    context.max_depth = 1
    context.max_files = 50
    context.do_download = False

## Disable downloads if requested (debugging):
fdl.allow_download(context.do_download)

## Root output folder must already exist:
if not os.path.isdir(context.save_root):
    sys.stderr.write("Error: output folder not found:\n"
                    + "--> %s\n\n" % context.save_root)

## Frame RLEVEL must be one of the allowed values:
rlevel_dirs = {0:'raw', 91:'specproc'}
if (context.data_rlevel != None):
    if not context.data_rlevel in rlevel_dirs.keys():
        sys.stderr.write("\n"
                        + "Unrecognized RLEVEL: %d\n\n" % context.data_rlevel
                        + "You should not see this ... please fix!!\n\n")
        sys.exit(1)

##--------------------------------------------------------------------------##



##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Frame look-up:
#headers = {'Authorization': 'Token ' + context.archive_token}
#base_url = 'https://archive-api.lco.global/'
#frame_url = base_url + 'frames/'
params = {'limit':context.max_files}
params['basename'] = 'nrs'              # only NRES files!!
#params['OBSTYPE'] = 'EXPOSE'
#params[ 'covers'] = lcofrm.wkt_from_coord((314.809246, +43.629033))

## Specific RLEVEL (if provided):
if (context.data_rlevel != None):
    params['RLEVEL'] = context.data_rlevel

## Optionally restrict search to single site:
params['SITEID'] = context.one_site

## Cheesy selection of calibration or science data:
if context.cal_only:
    params['PROPID'] = 'calibrate'
if context.sci_only:
    params['OBSTYPE'] = 'TARGET'

## Optional data time range (NOTE: --ndays supercedes --start):
if date_is_iso(context.end):
    params['end'  ] = promote_full_iso(context.end)
if date_is_iso(context.start):
    params['start'] = promote_full_iso(context.start)
if (context.ndays > 0.0):
    look_back = astt.TimeDelta(context.ndays, format='jd')
    params['start'] = (astt.Time.now() - look_back).iso

## Initial pass to get total frame count:
#get_cmd = {'url':frame_url, 'headers':headers, 'params':params}
try:
    #total = lcofrm.count_results(get_cmd)
    total = lcofrm.count_frames_by_params(params)
except:
    sys.stderr.write("Error: frame count failed!\n")
    sys.stderr.write("traceback:\n\n%s\n" % traceback.format_exc())
    sys.exit(1)
sys.stderr.write("Search identified %d frames in the LCO archive.\n" % total)

## Fetch frames (gets newest first):
results = []
try:
    depth, rcount = lcofrm.find_frames_by_params(params, results,
            maxdepth=context.max_depth, vlevel=0)
    #depth, rcount = lcofrm.recursive_request(get_cmd, results, 
    #        maxdepth=context.max_depth)
except:
    sys.stderr.write("Error: recursive frame retrieval failed!\n")
    sys.stderr.write("traceback:\n\n%s\n" % traceback.format_exc())
    sys.exit(1)
nhits = len(results)

sys.exit(0)
## Order by observation date:
results.sort(key=lambda x:x['DATE_OBS'])

## Optionally reverse list (to oldest first):
if not context.oldest_first:
    results.reverse()

## Ensure all sites are actually NRES sites:
for frame in results:
    if not (frame['SITEID'] in nres_sites):
        sys.stderr.write("Error: unsupported site: %s\n" % frame['SITEID'])
        sys.exit(1)

## Pre-download summary:
sys.stderr.write("Final frame list has %d entries.\n\n" % nhits)

##--------------------------------------------------------------------------##
## Multi-method download wrapper with built-in retry:
#def attempt_download(frame, filename, tmpfile):
#    if _fancy_downloading:
#        try:
#            fdl.fetch(frame['url'], tmpfile, resume=True, progress=True)
#        except:
#            return False
#    else:
#        try:
#            with open(tmpfile, 'wb') as f:
#                f.write(requests.get(frame['url']).content)
#        except (KeyboardInterrupt, SystemExit):
#            sys.stderr.write("Exit requested ...\n")
#            sys.exit(1)
#        except:
#            return False
#    sys.stderr.write("moving ... ")
#    shutil.move(tmpfile, filename)
#    sys.stderr.write("done.\n")
#    return True

#def download_with_retry(frame, filename, tmpfile, maxtry=5):
#    n_iters = 0
#    success = False
#    while not success:
#        success = attempt_download(frame, filename, tmpfile)
#        n_iters += 1
#        if (n_iters >= maxtry):
#            return False
#        pass
#    return True

##--------------------------------------------------------------------------##
## LCO path conventions:
def make_rel_daydir(frame):
    lsite = frame['SITEID']
    nrcam = nres_cam_map.get(lsite, None)
    if not nrcam:
        return (False, None)
    ibase = os.path.basename(frame['filename'])
    obs_day = ibase.split('-')[2]
    return os.path.join(lsite, nrcam, obs_day)

## Assemble local full path for specified frame:
def get_local_path_from_frame(frame, local_root):
    ibase = os.path.basename(frame['filename'])
    daydir_path = os.path.join(local_root, make_rel_daydir(frame))
    rlevel_path = os.path.join(daydir_path, rlevel_dirs[frame['RLEVEL']])
    return os.path.join(rlevel_path, ibase)

## Check if local frame copy exists:
def already_have_frame(frame, local_root):
    return os.path.isfile(get_local_path_from_frame(frame, local_root))

def frame_not_found(frame, local_root):
    return (not os.path.isfile(get_local_path_from_frame(frame, local_root)))

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Identify not-yet-fetched images:
sys.stderr.write("Identify remaining frames ... ")
still_need = [x for x in results if frame_not_found(x, context.save_root)]
sys.stderr.write("done.\n")
sys.stderr.write("Frames remaining: %d\n" % len(still_need))

## Stop if nothing to fetch:
if not still_need:
    #sys.stderr.write("No new frames to retrieve!\n")
    sys.stderr.write("Nothing to do!\n")
    sys.exit(0)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##


## Make spec list for download:
fetch_specs = []
for frame in still_need:
    full_save = get_local_path_from_frame(frame, context.save_root)
    temp_save = os.path.basename(full_save)
    fetch_specs.append((frame['url'], full_save, temp_save))

## Enable/skip downloads:
fdl.allow_download(context.do_download)

## Download files:
try:
    outcome = fdl.smart_fetch_bulk(fetch_specs)
except KeyboardInterrupt as exc:
    sys.stderr.write("Caught interrupt!\n")
    sys.exit(1)
except SystemExit as exc:
    sys.stderr.write("Caught system exit!\n")
    sys.exit(1)
except:
    sys.stderr.write("Unspecified failure ...\n")
    sys.exit(1)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

#sys.exit(0)
#
### Download files:
#total = len(still_need)
#ndownloaded = 0
#for i,frame in enumerate(still_need, 1):
#
#    # Get local path, check for file:
#    ibase = os.path.basename(frame['filename'])
#    isave = get_local_path_from_frame(frame, context.save_root)
#    sys.stderr.write("\rFetching %s (%d of %d) ... " % (ibase, i, total))
#    if os.path.isfile(isave):
#        sys.stderr.write("already retrieved!   ")
#        continue
#
#    # Ensure output folder exists:
#    save_folder = os.path.dirname(isave)
#    mkdir_p(save_folder)    # ensure existence
#
#    #isave = os.path.join(save_folder, ibase)
#
#    sys.stderr.write("not yet downloaded!  \nDownloading ... ")
#    itemp = 'dltemp_' + ibase
#    if context.do_download:
#        if not download_with_retry(frame, isave, itemp):
#            sys.stderr.write("Trouble downloading ... move on!\n")
#            continue
#    else:
#        sys.stderr.write("skipped (download disabled)!\n") 
#    ndownloaded += 1
#
#sys.stderr.write("\nAll downloads completed.\n")



######################################################################
# CHANGELOG (bin/fetch-nres-data.py):
#---------------------------------------------------------------------
#
#  2018-09-11:
#     -- Increased __version__ to 0.1.5.
#     -- Added site selection options -s, --site and --skip-site.
#
#  2018-09-05:
#     -- Increased __version__ to 0.1.0.
#     -- First created bin/fetch-nres-data.py.
#
