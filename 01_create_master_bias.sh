#!/bin/bash
#
#    Create a master bias and dark for the specified camera and day-obs.
#
# Rob Siverd
# Created:      2017-07-10
# Last updated: 2018-11-19
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Default options:
debug=0 ; clobber=0 ; force=0 ; timer=0 ; vlevel=0
script_version="0.64"
this_prog="${0##*/}"
#shopt -s nullglob
# Propagate errors through pipelines: set -o pipefail
# Exit if uninitialized variable used (set -u): set -o nounset
# Exit in case of nonzero status (set -e): set -o errexit

## Program options:
days_prior=0
days_after=0
keep_clean=0          # if true (=1), save temp files in daydir for later use
access_mode="symlink" # (symlink|copy) how to access existing 'clean' files

## Standard scratch files/dirs:
tmp_name="$(date +%Y%m%d.%H%M%S).$$.$(whoami)"
tmp_root="/tmp"
[ -d /dev/shm ] && [ -w /dev/shm ] && tmp_root="/dev/shm"
#tmp_dir="$tmp_root"
tmp_dir="$tmp_root/$tmp_name"
mkdir -p $tmp_dir
foo="$tmp_dir/foo_$$.fits"
bar="$tmp_dir/bar_$$.fits"
baz="$tmp_dir/baz_$$.fits"
qux="$tmp_dir/qux_$$.fits"
jnk="$foo $bar $baz $qux"  # working copy
def_jnk="$jnk"             # original set
dir_cleanup='(echo -e "\nAutomatic clean up ... " ; cmde "rm -vrf $tmp_dir")'
jnk_cleanup='for X in $jnk ; do [ -f $X ] && cmde "rm -vf $X" ; done'
trap "$jnk_cleanup" EXIT
##trap '[ -d $tmp_dir ] && cmde "rm -vrf $tmp_dir"' EXIT
#trap "[ -d $tmp_dir ] && $dir_cleanup" EXIT
trap "[ -d $tmp_dir ] && $dir_cleanup ; $jnk_cleanup" EXIT
#trap 'oops=$? ; echo ; exit $oops' HUP INT TERM

## Required programs:
declare -a need_exec
need_exec+=( awk cat fitsperc FuncDef getopt hdrtool imhget )
need_exec+=( kimstat medianize nres-cdp-trim-oscan sed tr )
for need in ${need_exec[*]}; do
   if ! ( /usr/bin/which $need >& /dev/null ); then
      echo "Error:  can't find '$need' in PATH !!" >&2
      exit 1
   fi
done

## Helper function definitions:
fd_args="--argchk --colors --cmde --echo"
#fd_args+=" --Critical"
fd_args+=" --rowwrite"
fd_args+=" --timers"
fd_args+=" --warnings"
FuncDef $fd_args >/dev/null || exit $?
eval "$(FuncDef $fd_args)"  || exit $?

##--------------------------------------------------------------------------##
## Syntax / help menu:
usage () {
   cat << EOH

Usage: $this_prog [options] camera YYYYMMDD
Create master bias for specified camera and day-obs.
Version: $script_version

Program options:
    -A, --after=N       include N days after fdate in stack [def: $days_after]
    -B, --before=N      include N days before fdate in stack [def: $days_prior]
    -c, --clobber       allow overwrite of output file
    -k, --keep_clean    preserve individual 'clean' files for later use
    -o, --output=FILE   send output to FILE

Other options:
        --deps          return list of required programs (for parent script)
    -h, --help          display this page
    -q, --quiet         less on-screen info and progress reporting
    -t, --timer         measure and report script execution time
    -v, --verbose       more on-screen info and progress reporting

EOH
#    -r, --random        randomize processing order (for parallel operation)
}

##--------------------------------------------------------------------------##
## Parse command line with getopt (reorders and stores CL args):
source aux/00_arg_parsing.sh

## Check for an appropriate number of arguments:
if [ -z "$2" ]; then
   usage >&2
   exit 1
fi
camid="$1"
fdate="$2"

[ $debug -eq 1 ] && vlevel=3
[ $vlevel -gt 1 ] && echo "Verbosity: $vlevel" >&2

##**************************************************************************##
##==========================================================================##
##--------------------------------------------------------------------------##

## Load shared functions:
func_files=( `ls func/??_*.sh 2>/dev/null` )
vyecho "Loading function files: ${func_files[*]}\n"
for item in ${func_files[*]}; do
   source $item
done

##--------------------------------------------------------------------------##
## Load configuration files:
conf_file="config.sh"
cams_list="cameras.txt"
good_cams="$tmp_dir/cleaned_camera_list.txt"
[ -f $conf_file ] || ErrorAbort "Can't find file: $conf_file"
[ -f $cams_list ] || ErrorAbort "Can't find file: $cams_list" 
cmde "source $conf_file"
declare -A cam_storage
declare -A cam_startdate
grep -v '^#' $cams_list > $good_cams   # remove comments
exec 10<$good_cams
while read cam subfolder startdate <&10; do
   cam_storage[$cam]="${arch_root}/${subfolder}"
   cam_startdate[$cam]="${startdate}"
done
exec 10>&-

## Check for output folder:
[ -d $save_root ] || ErrorAbort "Can't find directory: $save_root"

echo "Known cameras: ${!cam_storage[@]}"
#echo "cam_storage: ${cam_storage[*]}"

## Verify that clean/stack versions are available:
if [ ${#min_data_versions[*]} != 3 ]; then
   ErrorAbort "min_data_versions[] not defined!" 99
fi
if [ ${#min_code_versions[*]} != 3 ]; then
   ErrorAbort "min_code_versions[] not defined!" 99
fi

## Set up image version requirements:
need_data_versions=( `get_version_subset -b ${min_data_versions[*]}` )
need_code_versions=( `get_version_subset -b ${min_code_versions[*]}` )

##--------------------------------------------------------------------------##
## Validate camera:
use_arch="${cam_storage[$camid]}"
if [ -z "$use_arch" ]; then
   Recho "\nUnrecognized camera: '$camid'\n\n" >&2
   exit 1
fi

## Validate fdate:
earliest="${cam_startdate[$camid]}"
if [ $fdate -lt $earliest ]; then
   Recho "$fdate is outside the allowed date range for $camid.\n" >&2
   Recho "Earliest supported data are from: ${earliest}\n\n" >&2
   exit 99
fi

## Check for data folder and input files:
nite_dir="$use_arch/$fdate"
nite_raw="$nite_dir/raw"
echo "nite_raw: $nite_raw"
if [ ! -d $nite_raw ]; then
   recho "Folder not found: $nite_raw \n\n" >&2
   exit 1
fi

## List of nite_dirs:
recent_fdate_sec=$(fdate2unix $fdate)
#echo "recent_fdate_sec: $recent_fdate_sec"

fdate_list=()
for (( n = -$days_prior; n <= $days_after; n++ )); do
   unixtime_then=$(( recent_fdate_sec + 86400 * n ))
   fdate_list+=( `unix2fdate $unixtime_then` )
done
echo "fdate_list: ${fdate_list[*]}"

## Date range tag:
ndays=$(( days_prior + 1 + days_after ))
drtag=`printf '%02dd_%02db_%02da' $ndays $days_prior $days_after`
echo "drtag: $drtag"

##--------------------------------------------------------------------------##
## Count biases and darks:
bias_list=()
for nite in ${fdate_list[*]}; do
   bias_list+=( `ls $use_arch/$nite/raw/*b00.fits* 2>/dev/null` )
done

## Ensure uniqueness using DATE-OBS:
bias_list=( `pick_unique_subset ${bias_list[*]}` )
nbias=${#bias_list[*]}
echo "nbias: $nbias"

## For now, abort in case of missing files:
if [ $nbias -lt 2 ]; then
   recho "Too few biases found, nothing to process!\n"
   vcmde "rm -rf $tmp_dir"
   exit 0
fi

## Output files:
nite_folder="$save_root/$camid/$fdate"
nite_bias="$nite_folder/med_${camid}_bias_${fdate}_${drtag}.fits.fz"
cmde "mkdir -p $nite_folder" || exit $?

##--------------------------------------------------------------------------##
##                Existing Image Removal: Barrier Check                     ##
##--------------------------------------------------------------------------##

## Load helpers:
source func/01_barriers.sh

## Add fdates to list, one per line:
nlist="$tmp_dir/nite_list.$$.txt"
echo ${fdate_list[*]} | tr ' ' '\n' > $nlist

## Check fdate list against known barriers:
yecho "Barrier check ... "
if ( barrier_check_pass $camid bias $nlist ); then
   gecho "pass!\n"
   bcheck="PASS"
else
   recho "FAIL!\n"
   bcheck="FAIL"
fi

## Check for existing barrier violation:
if [ "$bcheck" = "FAIL" ]; then
   if [ -f $nite_bias ]; then
      yecho "Found existing (bogus) bias:\n"
      recho "--> $nite_bias \n\n"
      cmde "rm $nite_bias"
   fi

   rm -rf $tmp_dir
   exit 1
fi

##--------------------------------------------------------------------------##
##                Existing Image Removal: Version Check                     ##
##--------------------------------------------------------------------------##

if [ -f $nite_bias ]; then
   echo "Data version requirements: ${need_data_versions[*]}"
   echo "Code version requirements: ${need_code_versions[*]}"
   if ( header_gain_is_unity $nite_bias ) && \
      ( data_version_pass $nite_bias ${need_data_versions[*]} ) && \
      ( code_version_pass $nite_bias ${need_code_versions[*]} ); then
      Gecho "Existing $nite_bias passed version check!\n"
   else
      Recho "Existing $nite_bias FAILED version check!\n"
      cmde "rm $nite_bias"
   fi
else
   recho "No existing stacked bias found.\n"
fi

##--------------------------------------------------------------------------##
##          Existing Image Removal: Remove if wrong input count             ##
##--------------------------------------------------------------------------##

if [ -f $nite_bias ]; then
   numstack=$(imhget $nite_bias NUMSTACK)
   nlisted="${#bias_list[@]}"
   echo "Existing image inputs: $numstack"
   echo "Expected input images: $nlisted"
   if [ $numstack -eq $nlisted ]; then
      Gecho "Existing $nite_bias has matching inputs ($numstack images).\n"
   else
      Recho "Existing $nite_bias has input mismatch!  Rebuild needed ...\n"
      cmde "rm $nite_bias"
   fi
fi

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Create master bias (if not present):
if [ -f $nite_bias ]; then
   Gecho "$nite_bias already exists!\n"
else

   # Pretreat individual bias frames as needed:
   timer start
   yecho "Overscan correction etc... \n"
   for image in "${bias_list[@]}"; do
      ibase="${image##*/}"
      ifits="${ibase%.fz}"
      cbase="clean_${ifits}.fz"
      isave="$tmp_dir/$cbase"

      # Get 'clean' bias (process only if necessary):
      icheck="$(get_save_folder $image)/$cbase"
      if [ -f $icheck ]; then
         yecho "Checking ${icheck##*/} ... "
         if ( header_gain_is_unity $icheck ) && \
            ( data_version_pass $icheck ${need_data_versions[*]} ) && \
            ( code_version_pass $icheck ${need_code_versions[*]} ); then
            Gecho "version check PASSED!\n"
            gecho "Using existing temp-bias (${access_mode}): ${icheck}\n"
            case $access_mode in
               copy)    vcmde "cp -f $icheck $isave" || exit $?  ;;
               symlink) vcmde "ln -s $icheck $isave" || exit $?  ;;
               *) ErrorAbort "Unhandled access_mode: '$access_mode'" ;;
            esac
            continue
         else
            recho "version check FAILED, rebuild!\n"
            #[ $keep_clean -eq 1 ] && cmde "rm $icheck"
         fi
      fi

      # Existing 'clean' file unavailable, make it:
      echo
      cmde "nres-cdp-trim-oscan -q $image -o $foo"             || exit $?
      cmde "record_code_version $foo -b $script_version"       || exit $?
      cmde "record_data_version $foo -b $script_version"       || exit $?
      cmde "update_output_header $foo $camid BIAS 0.0 none"    || exit $?
      cmde "fpack -F -Y -qt 32 $foo"                           || exit $?
      cmde "mv -f $foo $isave"                                 || exit $?

      # Preserve files (if requested):
      if [ $keep_clean -eq 1 ]; then
         vcmde "mkdir -p $(dirname $icheck)"                   || exit $?
         cmde "cp -f $isave $icheck"                           || exit $?
      fi
      echo
   done
   timer

   min_bias_data_vers=$(find_min_cal_version -b $tmp_dir/clean*fits.fz)
   echo "min_bias_data_vers: $min_bias_data_vers"
   min_bias_code_vers=$(find_min_cal_version -B $tmp_dir/clean*fits.fz)
   echo "min_bias_code_vers: $min_bias_code_vers"

   # Combine biases with outlier rejection (stack-args in config.sh):
   mecho "\n`RowWrite 75 -`\n"
   opts="$bias_stack_args"
   cmde "medianize $opts $tmp_dir/clean*fits.fz -o '!$foo'" || exit $?
   timer

   # Add stats and identifiers to header:
   cmde "fitsperc -qS $foo --parallel3"                     || exit $?
   cmde "kimstat -qSC9 $foo"                                || exit $?
   cmde "record_data_version $foo -b $min_bias_data_vers"   || exit $?
   #cmde "record_code_version $foo -b $min_bias_code_vers"   || exit $?
   cmde "record_code_version $foo -b $script_version"       || exit $?
   cmde "update_output_header $foo $camid BIAS 0.0 $drtag"  || exit $?
   cmde "fpack -F -Y -qt 32 $foo"                           || exit $?
   cmde "mv -f $foo $nite_bias"                             || exit $?

   ## Preserve files (if requested):
   #if [ $keep_clean -eq 1 ]; then
   #   yecho "Preserving stack ...\n"
   #   for item in $tmp_dir/clean*fits; do
   #      dst_file="$(get_save_folder $item)/${item##*/}"
   #      if [ ! -f $dst_file ]; then
   #         vcmde "mkdir -p $(dirname $dst_file)"           || exit $?
   #         cmde "mv -f $item $dst_file"                    || exit $?
   #      fi
   #   done
   #fi
   timer
fi

##--------------------------------------------------------------------------##
## Clean up:
[ -d $tmp_dir ] && [ -O $tmp_dir ] && rm -rf $tmp_dir
[ -f $foo ] && rm -f $foo
[ -f $bar ] && rm -f $bar
[ -f $baz ] && rm -f $baz
[ -f $qux ] && rm -f $qux
exit 0

######################################################################
# CHANGELOG (01_create_master_bias.sh):
#---------------------------------------------------------------------
#
#  2018-08-07:
#     -- Increased script_version to 0.60.
#     -- Now check for GAIN=1.0 along with code/data version requirements in
#           both clean_ image inspection and for existing image removal.
#
#  2018-08-06:
#     -- Increased script_version to 0.59.
#     -- Added camid to median-stacked master bias file name.
#
#  2018-08-03:
#     -- Increased script_version to 0.58.
#     -- Gain is now corrected in nres-cdp-trim-overscan.py. Calib build
#           script versions are increased to force data rebuild.
#
#  2018-01-07:
#     -- Increased script_version to 0.57.
#     -- Added message indicating when no stacked bias exists yet.
#     -- Implemented separate for eff/scr versioning.
#     -- Increased script_version to 0.56.
#     -- Implemented separate version requirements for clean and stack data.
#     -- Added check for in-bounds fdate.
#     -- BIASVERS in master output file is now the lowest BIASVERS of its
#           input images.
#
#  2018-01-05:
#     -- Increased script_version to 0.55.
#     -- Now check existing master and 'clean' biases for version compliance.
#           Offending files are automatically removed and rebuilt.
#     -- Current script_version now recorded to BIASVERS keyword.
#     -- Now use new 'aux' and 'func' locations for common code and routines.
#
#  2017-08-07:
#     -- Increased script_version to 0.50.
#     -- Implemented stacking barriers.
#
#  2017-07-26:
#     -- Increased script_version to 0.45.
#     -- Now ensure existence of destination folder before attempting to save
#           'clean' temporary files.
#     -- Now ensure unique files are selected as input using DATE-OBS. This
#           protects against unexpected results if .fits and .fits.fz exist.
#
#  2017-07-24:
#     -- Increased script_version to 0.40.
#     -- To increase efficiency and provide compatibility with symlinks,
#           individual 'clean' files are only stored if not already present.
#     -- Added timers, compared access modes (copy vs. symlink). Symlink wins.
#     -- Stacking can now leverage existing 'clean' files when available.
#
#  2017-07-23:
#     -- Increased script_version to 0.35.
#     -- Moved argument parsing into separate 00_arg_parsing.sh code file
#           (shared with 01_create_master_dark.sh).
#     -- Removed no-longer-used dark stacking code (now in next script).
#
#  2017-07-20:
#     -- Increased script_version to 0.30.
#     -- days_prior and days_after are now selectable on the command line
#           with -B, --before=N and -A, --after=N options.
#     -- Finished script update, multi-day biases now build successfully.
#
#  2017-07-13:
#     -- Increased script_version to 0.20.
#     -- Renamed script to 01_create_master_bias.sh and moved dark creation
#           to a separate script, 01_create_master_dark.sh. The desire to use
#           multi-day bias stacks to correct the individual dark frames prior
#           to stacking requires that these scripts be separated.
#           
#
#  2017-07-12:
#     -- Increased script_version to 0.15.
#     -- Added code to select a range of FDATEs for stacking.
#
#  2017-07-10:
#     -- Increased script_version to 0.10.
#     -- First created 01_create_master_bias_dark.sh.
#
