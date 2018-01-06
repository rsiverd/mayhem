#!/bin/bash
#
#    Create a master bias and dark for the specified camera and day-obs.
#
# Rob Siverd
# Created:      2017-07-10
# Last updated: 2017-07-12
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Default options:
debug=0 ; clobber=0 ; force=0 ; timer=0 ; vlevel=0
script_version="0.15"
this_prog="${0##*/}"
#shopt -s nullglob
# Propagate errors through pipelines: set -o pipefail
# Exit if uninitialized variable used (set -u): set -o nounset
# Exit in case of nonzero status (set -e): set -o errexit

## Program options:
#save_file=""
#shuffle=0
#confirmed=0

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
#need_exec+=( shuf shuffle sort ) # for randomization
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
    -c, --clobber       allow overwrite of output file
    -o, --output=FILE   send output to FILE
    -r, --random        randomize processing order (for parallel operation)

Other options:
        --deps          return list of required programs (for parent script)
    -h, --help          display this page
    -q, --quiet         less on-screen info and progress reporting
    -t, --timer         measure and report script execution time
    -v, --verbose       more on-screen info and progress reporting

EOH
#        --debug         extra verbosity to assist bug hunting
#    -f, --force         force redo of whatever this script does
#    -f, --force         allow clobbering of output file
}

##--------------------------------------------------------------------------##
## Parse command line with getopt (reorders and stores CL args):
s_opts="co:rhqtv" # f
l_opts="START,clobber,output:,random"
l_opts+=",debug,deps,help,quiet,timer,verbose" # force
args=`getopt -n $this_prog -o $s_opts -l $l_opts -- "$@"` ; failed=$?

## Check for parsing errors:
if [ $failed -ne 0 ]; then 
   echo "Try \`$this_prog --help' for more information." >&2
   exit 1
fi

## Change the arguments in the environment (e.g., $1, $2, ... $N):
eval set -- "$args"

## Loop through arguments (shift removes them from environment):
while true ; do
   case $1 in
      #-------------------------------------------------------------------
     #--START)
     #   confirmed=1
     #   shift
     #   ;;
      #-------------------------------------------------------------------
      -c|--clobber)
         [ $vlevel -ge 0 ] && yecho "Enabled output clobber!\n" >&2
         clobber=1
         shift
         ;;
      #-------------------------------------------------------------------
     #-n|--number)
     #   case $2 in
     #      -*)
     #         msg="Option -n|--number requires an argument!"
     #         #msg="Option -n|--number needs a positive integer argument!"
     #         #msg="Option -n|--number needs a positive numerical argument!"
     #         Recho "\n${msg}\n" >&2
     #         usage >&2
     #         exit 1
     #         ;;
     #      *)
     ###       if !( num_check_pass $2 ); then
     ###       if !( num_check_pass $2 ) || (is_negative $2); then
     ###       if !( int_check_pass $2 ) || [ $2 -lt 0 ]; then
     #            Recho "Invalid value: " >&2 ; Yecho "$2 \n\n" >&2
     #            exit 1
     #         fi
     #         num_val=$2
     #         ;;
     #   esac
     #   [ $vlevel -ge 0 ] && yecho "Using value: ${num_val}\n" >&2
     #   shift 2
     #   ;;
      -o|--output)
         case $2 in
            -*)
               Recho "\nOption -o|--output requires an argument!\n" >&2
               usage >&2
               exit 1
               ;;
            *)
               save_file=$2
               # check value here ...
               if [ $clobber -eq 0 ] && [ -f $save_file ]; then
                  Recho "\nFile already exists: " >&2
                  Yecho "$save_file \n\n" >&2
                  exit 1
               fi
               ;;
         esac
         [ $vlevel -ge 0 ] && yecho "Output to: $save_file \n" >&2
         shift 2
         ;;
      -r|--random)
         [ $vlevel -ge 0 ] && yecho "Randomizing order!\n" >&2
         shuffle=1
         shift
         ;;
      #-------------------------------------------------------------------
      # Additional options (output control etc.):
      --debug)
         yecho "Debugging mode enabled!\n"
         debug=1
         shift
         ;;
      --deps)
         echo ${need_exec[*]}
         exit 0
         ;;
     #-f|--force)
     #   [ $vlevel -ge 0 ] && yecho "Output clobbering enabled!\n" >&2
     #   [ $vlevel -ge 0 ] && yecho "Forcing <WHAT THIS SCRIPT DOES>!\n" >&2
     #   clobber=1
     #   force=1
     #   shift
     #   ;;
      -h|--help)
         usage
         exit 0
         ;;
      -q|--quiet)
         (( vlevel-- ))
         shift
         ;;
      -t|--timer)
         [ $vlevel -ge 1 ] && yecho "Timing script execution!\n" >&2
         timer=1
         shift
         ;;
      -v|--verbose)
         (( vlevel++ ))
         shift
         ;;
      #-------------------------------------------------------------------
      --)
         shift
         break 
         ;;
      *)
         echo -e "\n$this_prog error: unhandled option: $1 \n" >&2
         exit 1 
         ;;
      #-------------------------------------------------------------------
   esac
done

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
conf_file="config.txt"
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
days_before=2
days_after=0
recent_fdate_sec=$(fdate2unix $fdate)
echo "recent_fdate_sec: $recent_fdate_sec"

fdate_list=()
for (( n = -$days_before; n <= $days_after; n++ )); do
   unixtime_then=$(( recent_fdate_sec + 86400 * n ))
   fdate_list+=( `unix2fdate $unixtime_then` )
done
echo "fdate_list: ${fdate_list[*]}"

## Date range tag:
ndays=$(( days_before + 1 + days_after ))
drtag=`printf '%02dd_%02db_%02da' $ndays $days_before $days_after`
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

   # Pretreat individual bias frames:
   yecho "Overscan correction etc... \n"
   for image in "${bias_list[@]}"; do
      ibase="${image##*/}"
      ifits="${ibase%.fz}"
      isave="$tmp_dir/clean_$ifits"
      cmde "nres-cdp-trim-oscan -q $image -o $isave"     || exit $?
   done

   # Combine biases with outlier rejection:
   cmde "medianize -c3 $tmp_dir/clean*fits -o '!$foo'"   || exit $?

   # Add stats and identifiers to header:
   cmde "fitsperc -qS $foo"                              || exit $?
   cmde "kimstat -qSC9 $foo"                             || exit $?
   cmde "hdrtool $foo -U OBSTYPE --value='BIAS'"         || exit $?
   cmde "hdrtool $foo -U EXPTIME --value=0.0"            || exit $?
   cmde "mv -f $foo $nite_bias"                          || exit $?

   # Clean up:
   cmde "rm $tmp_dir/clean*fits"
fi

exit 0

##--------------------------------------------------------------------------##
## Stop here if too few darks:
if [ $ndark -lt 2 ]; then
   recho "Too few darks found, nothing to process!\n"
   exit 0
fi

##--------------------------------------------------------------------------##
## Create master dark (if not present):
if [ -f $nite_dark ]; then
   Gecho "$nite_dark already exists!\n"
else

   # Pretreat individual dark frames:
   yecho "Overscan correction and bias subtraction ...\n"
   for image in "${dark_list[@]}"; do
      ibase="${image##*/}"
      ifits="${ibase%.fz}"
      isave="$tmp_dir/clean_$ifits"

      # Subtract overscan:
      cmde "nres-cdp-trim-oscan -q $image -o $foo"       || exit $?

      # Subtract bias and divide out exposure time:
      darkexp="$(imhget EXPTIME $foo)"
      echo "darkexp: $darkexp"
      cmde "fitsarith -qHi $foo -S $nite_bias -d $darkexp -o '!$bar'" || exit $?
      cmde "mv -f $bar $isave"
   done

   # Combine biases with outlier rejection:
   cmde "medianize -c3 $tmp_dir/clean*fits -o '!$foo'"   || exit $?

   # Add stats and identifiers to header:
   cmde "fitsperc -qS $foo"                              || exit $?
   cmde "kimstat -qSC9 $foo"                             || exit $?
   cmde "hdrtool $foo -U OBSTYPE --value='DARK'"         || exit $?
   cmde "hdrtool $foo -U EXPTIME --value=1.0"            || exit $?
   cmde "mv -f $foo $nite_dark"                          || exit $?

   # Clean up:
   cmde "rm $tmp_dir/clean*fits"
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
