#!/bin/bash
#
#    Create a master bias and dark for the specified camera and day-obs.
#
# Rob Siverd
# Created:      2017-07-10
# Last updated: 2017-07-23
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Default options:
debug=0 ; clobber=0 ; force=0 ; timer=0 ; vlevel=0
script_version="0.35"
this_prog="${0##*/}"
#shopt -s nullglob
# Propagate errors through pipelines: set -o pipefail
# Exit if uninitialized variable used (set -u): set -o nounset
# Exit in case of nonzero status (set -e): set -o errexit

## Program options:
#save_file=""
#shuffle=0
#confirmed=0
days_prior=0
days_after=0
keep_clean=0      # if true (>0), save temporary files in daydir for later use

## Standard scratch files/dirs:
tmp_name="$(date +%Y%m%d.%H%M%S).$$.$(whoami)"
tmp_root="/tmp"
[ -d /dev/shm ] && [ -w /dev/shm ] && tmp_root="/dev/shm"
tmp_dir="$tmp_root"
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
#fd_args+=" --timers"
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

###--------------------------------------------------------------------------##
### Parse command line with getopt (reorders and stores CL args):
#s_opts="A:B:co:rhqtv"
#l_opts="after:,before:,clobber,output:,random"
#l_opts+=",debug,deps,help,quiet,timer,verbose"
#args=`getopt -n $this_prog -o $s_opts -l $l_opts -- "$@"` ; failed=$?
#
### Check for parsing errors:
#if [ $failed -ne 0 ]; then 
#   echo "Try \`$this_prog --help' for more information." >&2
#   exit 1
#fi
#
### Change the arguments in the environment (e.g., $1, $2, ... $N):
#eval set -- "$args"

## Loop through arguments (shift removes them from environment):
source 00_arg_parsing.sh

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

### List of recent nights:
#nights=5
#declare -a nite_list
#for (( x = $nights; x >= -1; x-- )); do
#   nite_list+=( $(date -u +%Y%m%d --date="$x days ago") )
#done

##--------------------------------------------------------------------------##
## Load configuration files:
conf_file="config.sh"
cams_list="cameras.txt"
[ -f $conf_file ] || ErrorAbort "Can't find file: $conf_file"
[ -f $cams_list ] || ErrorAbort "Can't find file: $cams_list" 
cmde "source $conf_file"
declare -A cam_storage
exec 10<$cams_list
while read cam folder <&10; do
   cam_storage[$cam]="$folder"
done
exec 10>&-

## Check for output folder:
[ -d $save_root ] || ErrorAbort "Can't find directory: $save_root"

echo "Known cameras: ${!cam_storage[@]}"
#echo "cam_storage: ${cam_storage[*]}"

##--------------------------------------------------------------------------##
## Validate camera:
use_arch="${cam_storage[$camid]}"
if [ -z "$use_arch" ]; then
   Recho "\nUnrecognized camera: '$camid'\n\n" >&2
   exit 1
fi

## Check for data folder and input files:
nite_dir="$use_arch/$fdate"
nite_raw="$nite_dir/raw"
echo "nite_raw: $nite_raw"
if [ ! -d $nite_raw ]; then
   recho "Folder not found: $nite_raw \n\n" >&2
   exit 1
fi

## Convert between FDATE and UNIX seconds:
fdate2unix () { date +%s --date=$1; }
unix2fdate () { date +%Y%m%d --date="@$1"; }

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
#dark_list=()
for nite in ${fdate_list[*]}; do
   bias_list+=( `ls $use_arch/$nite/raw/*b00.fits* 2>/dev/null` )
   #dark_list+=( `ls $use_arch/$nite/raw/*d00.fits* 2>/dev/null` )
done

#bias_list=( `ls $nite_raw/*b00.fits* 2>/dev/null` )
#dark_list=( `ls $nite_raw/*d00.fits* 2>/dev/null` )
nbias=${#bias_list[*]}
#ndark=${#dark_list[*]}
echo "nbias: $nbias"
#echo "ndark: $ndark"
#exit 0

## For now, abort in case of missing files:
if [ $nbias -lt 2 ]; then
   recho "Too few biases found, nothing to process!\n"
   exit 0
fi

## Output files:
nite_folder="$save_root/$camid/$fdate"
nite_bias="$nite_folder/med_bias_${fdate}_${drtag}.fits"
#nite_dark="$nite_folder/med_dark_${fdate}_${drtag}.fits"
cmde "mkdir -p $nite_folder" || exit $?

##--------------------------------------------------------------------------##
## Create master bias (if not present):
if [ -f $nite_bias ]; then
   Gecho "$nite_bias already exists!\n"
else

   # Pretreat individual bias frames as needed:
   yecho "Overscan correction etc... \n"
   for image in "${bias_list[@]}"; do
      ibase="${image##*/}"
      ifits="${ibase%.fz}"
      isave="$tmp_dir/clean_${ifits}"
      #isave="$tmp_dir/clean_${drtag}_${ifits}"
      cmde "nres-cdp-trim-oscan -q $image -o $foo"       || exit $?
      #cmde "hdrtool $foo  -U CAMID --value='$camid'"     || exit $?
      #cmde "hdrtool $foo -CU DRTAG --value='$drtag'"     || exit $?
      update_output_header $foo $camid $drtag BIAS 0.0   || exit $?
      cmde "mv -f $foo $isave"                           || exit $?
   done

   # Combine biases with outlier rejection:
   opts="-c3 --itable"
   cmde "medianize $opts $tmp_dir/clean*fits -o '!$foo'" || exit $?

   # Add stats and identifiers to header:
   cmde "fitsperc -qS $foo"                              || exit $?
   cmde "kimstat -qSC9 $foo"                             || exit $?
   #cmde "hdrtool $foo  -U OBSTYPE --value='BIAS'"        || exit $?
   #cmde "hdrtool $foo  -U EXPTIME --value=0.0"           || exit $?
   #cmde "hdrtool $foo -CU DRTAG   --value='$drtag'"      || exit $?
   update_output_header $foo $camid $drtag BIAS 0.0      || exit $?
   cmde "mv -f $foo $nite_bias"                          || exit $?

   # Clean up:
   if [ $keep_clean -eq 1 ]; then
      yecho "Preserving stack ...\n"
      for item in $tmp_dir/clean*fits; do
         echo "item: $item"
         dst_dir=`get_save_folder $item`
         cmde "mv -f $item $dst_dir/." || exit $?
      done
   else
      yecho "Cleaning up temp files ...\n"
      cmde "rm $tmp_dir/clean*fits"
   fi
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
