#!/bin/bash
#
#    Create a master bias and dark for the specified camera and day-obs.
#
# Rob Siverd
# Created:      2017-07-24
# Last updated: 2017-07-26
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Default options:
debug=0 ; clobber=0 ; force=0 ; timer=0 ; vlevel=0
script_version="0.25"
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
calib_type=""     # name of lamp+calibration channel selected (REQUIRED)
access_mode="symlink" # (symlink|copy) how to access existing 'clean' files

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
need_exec+=( kimstat medianize nres-cdp-trim-oscan nres-labcal sed tr )
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
Create master lamp flat or arc for specified camera and day-obs.
Version: $script_version

Program options:
    -A, --after=N       include N days after fdate in stack [def: $days_after]
    -B, --before=N      include N days before fdate in stack [def: $days_prior]
    -c, --clobber       allow overwrite of output file
    -k, --keep_clean    preserve individual 'clean' files for later use
    -o, --output=FILE   send output to FILE

Lamp type (REQUIRED):
        --thar01        Thorium Argon DOUBLE, fiber0 + fiber1
        --thar12        Thorium Argon DOUBLE, fiber1 + fiber2
        --tung01        Tungsten Halogen FLAT, fiber0 + fiber1
        --tung12        Tungsten Halogen FLAT, fiber1 + fiber2

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

##--------------------------------------------------------------------------##
## Parse command line with getopt (reorders and stores CL args):
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

## Abort in case calib/channel was not specified:
if [ -z "$calib_type" ]; then
   Recho "Error: no lamp/channel specified!\n" >&2
   usage >&2
   exit 1
fi

## Load calibration config (file name specified in config.sh):
[ -f $calib_conf ] || ErrorAbort "Can't find file: $calib_conf"

declare -A objects_keyword obstype_keyword 
declare -A filename_suffix filename_prefix
exec 10<$calib_conf
while read ctype suffix objskey obstype prefix <&10; do
   filename_suffix[$ctype]=$suffix
   filename_prefix[$ctype]=$prefix
   objects_keyword[$ctype]=$objskey
   obstype_keyword[$ctype]=$obstype
done
exec 10>&-

## Convert calibration choice into OBJECTS keyword:
lampobj=${objects_keyword[$calib_type]}
fsuffix=${filename_suffix[$calib_type]}
obstype=${obstype_keyword[$calib_type]}
lamppre=${filename_prefix[$calib_type]}
#echo "lampobj: $lampobj"
#echo "obstype: $obstype"
#echo "fsuffix: $fsuffix"

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
## Create list of all lamp calibration files (both channels):
lamp_files="$tmp_dir/lamp_files.txt"
rm $lamp_files 2>/dev/null
for nite in ${fdate_list[*]}; do
   ls $use_arch/$nite/raw/*${fsuffix}.fits* 2>/dev/null >> $lamp_files
done
lamp_list=( `imhget -l $lamp_files OBJECTS \
   | awk -v want=$lampobj '{ if ($2 == want) print $1 }'` )

## Ensure uniqueness using DATE-OBS:
lamp_list=( `pick_unique_subset ${lamp_list[*]}` )
nlamp=${#lamp_list[*]}
echo "nlamp: $nlamp"

## Output files:
nite_folder="$save_root/$camid/$fdate"
nite_lampsave="$nite_folder/med_${calib_type}_${fdate}_${drtag}.fits"
cmde "mkdir -p $nite_folder" || exit $?

##--------------------------------------------------------------------------##
## Stop here if too few darks:
if [ $nlamp -lt 2 ]; then
   recho "Too few calibs found, nothing to process!\n"
   vcmde "rm -rf $tmp_dir"
   exit 0
fi

##--------------------------------------------------------------------------##
## Create master dark (if not present):
if [ -f $nite_lampsave ]; then
   Gecho "$nite_lampsave already exists!\n"
else

   # Pretreat individual dark frames:
   timer start
   yecho "Overscan correction and bias/dark subtraction ...\n"
   for image in "${lamp_list[@]}"; do
      # Identify best current bias frame:
      yecho "Selecting good bias ... "
      use_bias=$(pick_best_bdcal $image $camid bias --prev) || exit $?
      bdr_tag=$(imhget DRTAG $use_bias)
      if [ -z "$use_bias" ] || [ ! -f $use_bias ]; then
         Recho "Blank name or file missing: '$use_bias'\n\n" >&2
         exit 1
      fi
      gecho "done.\n"
      echo "use_bias: $use_bias"

      # Identify best current dark frame:
      yecho "Selecting good dark ... "
      use_dark=$(pick_best_bdcal $image $camid dark --prev) || exit $?
      ddr_tag=$(imhget DRTAG $use_dark)
      if [ -z "$use_dark" ] || [ ! -f $use_dark ]; then
         Recho "Blank name or file missing: '$use_dark'\n\n" >&2
         exit 1
      fi
      gecho "done.\n"
      echo "use_dark: $use_dark"

      # Temporary 'clean' file name (includes DRTAG of associated calibs):
      ibase="${image##*/}"
      ifits="${ibase%.fz}"
      cbase="clean_${drtag}_${ifits}"
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
      # Otherwise, create cleaned file for stacking:
      # --------------------------------------------

      # Subtract overscan:
      cmde "nres-cdp-trim-oscan -q $image -o $foo"       || exit $?
      lampexp="$(imhget EXPTIME $foo)"
      echo "lampexp: $lampexp"

      # Debias and subtract scaled dark:
      cmde "nres-labcal $foo -b $use_bias -d $use_dark -o $bar"   || exit $?
      cmde "hdrtool $bar --add_hist='use_bias ${use_bias##*/}'"   || exit $?
      cmde "hdrtool $bar --add_hist='use_dark ${use_dark##*/}'"   || exit $?
      update_output_header $bar $camid $obstype $lampexp $drtag   || exit $?
      cmde "mv -f $bar $isave"                                    || exit $?

   done
   timer

   # Combine darks with outlier rejection (stack_args in config.sh):
   mecho "\n`RowWrite 75 -`\n"
   opts="$dark_stack_args"
   cmde "medianize $opts $tmp_dir/clean*fits -o '!$foo'" || exit $?
   timer

   # Add stats and identifiers to header:
   cmde "fitsperc -qS $foo"                              || exit $?
   cmde "kimstat -qSC9 $foo"                             || exit $?
   update_output_header $foo $camid $drtag DARK 1.0      || exit $?
   cmde "mv -f $foo $nite_lampsave"                      || exit $?

   # Preserve stack files (if requested):
   if [ $keep_clean -eq 1 ]; then
      yecho "Preserving stack ...\n"
      for item in $tmp_dir/clean*fits; do
         dst_file="$(get_save_folder $item)/${item##*/}"
         if [ ! -f $dst_file ]; then
            cmde "mv -f $item $dst_file" || exit $?
         fi
      done
   fi
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
# CHANGELOG (03_create_master_lamp.sh):
#---------------------------------------------------------------------
#
#  2017-07-26:
#     -- Increased script_version to 0.25.
#     -- Now pare input list down to unique files using DATE-OBS keyword.
#
#  2017-07-24:
#     -- Increased script_version to 0.21.
#     -- Now record use_bias and use_dark in 'clean' image HISTORY cards.
#     -- Increased script_version to 0.20.
#     -- The $drtag in filenames is now exclusively the 'master' value set by the
#           script (values from biases and darks are ignored).
#     -- Increased script_version to 0.10.
#     -- First created 03_create_master_lamp.sh from 03_create_master_lamp.sh.
#
