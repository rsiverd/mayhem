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
#        --debug         extra verbosity to assist bug hunting
#    -f, --force         force redo of whatever this script does
#    -f, --force         allow clobbering of output file
#    -r, --random        randomize processing order (for parallel operation)
}

###--------------------------------------------------------------------------##
### Parse command line with getopt (reorders and stores CL args):
#s_opts="A:B:cko:rhqtv" # f
#l_opts="after:,before:,clobber,keep_clean,output:,random"
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
dark_list=()
for nite in ${fdate_list[*]}; do
   dark_list+=( `ls $use_arch/$nite/raw/*d00.fits* 2>/dev/null` )
done

ndark=${#dark_list[*]}
echo "ndark: $ndark"

## Output files:
nite_folder="$save_root/$camid/$fdate"
nite_bias="$nite_folder/med_bias_${fdate}_${drtag}.fits"
nite_dark="$nite_folder/med_dark_${fdate}_${drtag}.fits"
cmde "mkdir -p $nite_folder" || exit $?

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
      # Identify best current bias frame:
      #use_bias=`get_bias_path $image $camid`
      use_bias=$(pick_best_bdcal $image $camid bias) || exit $?
      bias_drtag=$(imhget DRTAG $use_bias)
      #echo "bias_drtag: $bias_drtag"
      if [ -z "$use_bias" ] || [ ! -f $use_bias ]; then
         Recho "Blank name or file missing: '$use_bias'\n\n" >&2
         exit 1
      fi

      # Temporary 'clean' file name (includes DRTAG of associated bias):
      ibase="${image##*/}"
      ifits="${ibase%.fz}"
      cbase="clean_${bias_drtag}_${ifits}"
      isave="$tmp_dir/$cbase"

      # Use existing cleaned dark if possible:
      icheck="$(get_save_folder $image)/$cbase"
      if [ -f $icheck ]; then
         gecho "Using existing temp-dark (symlink)!\n"
         cmde "ln -s $icheck $isave" || exit $?
         continue
      fi
      echo

      # --------------------------------------------
      # Otherwise, create cleaned dark for stacking:
      # --------------------------------------------

      # Subtract overscan:
      cmde "nres-cdp-trim-oscan -q $image -o $foo"       || exit $?

      # Subtract bias and divide out exposure time:
      darkexp="$(imhget EXPTIME $foo)"
      echo "darkexp: $darkexp"
      cmde "fitsarith -qHi $foo -S $use_bias -d $darkexp -o '!$bar'" || exit $?
      update_output_header $bar $camid $drtag DARK 1.0   || exit $?
      cmde "mv -f $bar $isave"                           || exit $?

      ##echo "ibase: $ibase"
      ##echo "fdate: $fdate"
      #echo "cbase: $cbase"
      #echo "isave: $isave"
      #exit 0

      ## Subtract bias and divide out exposure time:
      #darkexp="$(imhget EXPTIME $foo)"
      #echo "darkexp: $darkexp"
      #cmde "fitsarith -qHi $foo -S $use_bias -d $darkexp -o '!$bar'" || exit $?
      ##cmde "hdrtool $bar  -U OBSTYPE  --value='DARK'"    || exit $?
      ##cmde "hdrtool $bar  -U EXPTIME  --value=1.0"       || exit $?
      #update_output_header $bar $camid $drtag DARK 1.0   || exit $?
      #cmde "mv -f $bar $isave"                           || exit $?
   done

   # Combine darks with outlier rejection:
   opts="-r410 -c3 --itable"
   cmde "medianize $opts $tmp_dir/clean*fits -o '!$foo'" || exit $?

   # Add stats and identifiers to header:
   cmde "fitsperc -qS $foo"                              || exit $?
   cmde "kimstat -qSC9 $foo"                             || exit $?
   #cmde "hdrtool $foo  -U OBSTYPE  --value='DARK'"       || exit $?
   #cmde "hdrtool $foo  -U EXPTIME  --value=1.0"          || exit $?
   update_output_header $foo $camid $drtag DARK 1.0      || exit $?
   cmde "mv -f $foo $nite_dark"                          || exit $?

   # Clean up:
   if [ $keep_clean -eq 1 ]; then
      yecho "Preserving stack ...\n"
      for item in $tmp_dir/clean*fits; do
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
# CHANGELOG (01_create_master_dark.sh):
#---------------------------------------------------------------------
#
#  2017-07-23:
#     -- Increased script_version to 0.35.
#     -- Moved argument parsing into separate 00_arg_parsing.sh code file
#           (shared with 01_create_master_dark.sh).
#     -- pick_bias() function is now part of config.txt.
#     -- Now record DRTAG in FITS header for all potential output images.
#     -- Added -k, --keep_clean option to preserve individal stack files.
#     -- Added $drtag string to 'clean' filename for individual darks.
#
#  2017-07-20:
#     -- Increased script_version to 0.30.
#     -- days_prior and days_after are now selectable on the command line
#           with -B, --before=N and -A, --after=N options.
#     -- Finished script update, multi-day darks now build successfully.
#
#  2017-07-13:
#     -- Increased script_version to 0.20.
#     -- Renamed a copy of old script 01_create_master_bias_dark.sh to 
#           01_create_master_dark.sh. This is capable of producing multi-day
#           master darks and using multi-day master biases.
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
