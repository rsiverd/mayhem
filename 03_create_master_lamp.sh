#!/bin/bash
#
#    Create a master bias and dark for the specified camera and day-obs.
#
# Rob Siverd
# Created:      2017-07-24
# Last updated: 2018-11-19
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Default options:
debug=0 ; clobber=0 ; force=0 ; timer=0 ; vlevel=0
script_version="0.50"
this_prog="${0##*/}"
#shopt -s nullglob
# Propagate errors through pipelines: set -o pipefail
# Exit if uninitialized variable used (set -u): set -o nounset
# Exit in case of nonzero status (set -e): set -o errexit

## Program options:
days_prior=0
days_after=0
keep_clean=0          # if true (>0), save temp files in daydir for later use
calib_type=""         # name of lamp+calibration channel selected (REQUIRED)
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
need_data_versions=( `get_version_subset -l ${min_data_versions[*]}` )
need_code_versions=( `get_version_subset -l ${min_code_versions[*]}` )

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
## Create list of all lamp calibration files (both channels):
lamp_files="$tmp_dir/lamp_files.txt"
rm $lamp_files 2>/dev/null
for nite in ${fdate_list[*]}; do
   ls $use_arch/$nite/raw/*${fsuffix}.fits* 2>/dev/null >> $lamp_files
done
nfiles=$(cat $lamp_files | wc -l)
if [ $nfiles -lt 2 ]; then
   recho "Not enough calibration files to build master!\n"
   vcmde "rm -rf $tmp_dir"
   exit 0
fi

## Extract OBJECTS keyword:
lamp_list=( `imhget -l $lamp_files OBJECTS \
   | awk -v want=$lampobj '{ if ($2 == want) print $1 }'` )

## Ensure uniqueness using DATE-OBS:
lamp_list=( `pick_unique_subset ${lamp_list[*]}` )
nlamp=${#lamp_list[*]}
echo "nlamp: $nlamp"

## Output files:
nite_folder="$save_root/$camid/$fdate"
nite_lampsave="$nite_folder/med_${camid}_${calib_type}_${fdate}_${drtag}.fits.fz"
cmde "mkdir -p $nite_folder" || exit $?

##--------------------------------------------------------------------------##
## Stop here if too few lamp calibs:
if [ $nlamp -lt 2 ]; then
   recho "Too few calibs found, nothing to process!\n"
   vcmde "rm -rf $tmp_dir"
   exit 0
fi

##--------------------------------------------------------------------------##
##                FIXME: IMPLEMENT CCD TEMPERATURE CHECK                    ##
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
##                Existing Image Removal: Barrier Check                     ##
##--------------------------------------------------------------------------##

## Add fdates to list, one per line:
nlist="$tmp_dir/nite_list.$$.txt"
echo ${fdate_list[*]} | tr ' ' '\n' > $nlist

## Check fdate list against known barriers:
yecho "Barrier check ... "
if ( barrier_check_pass $camid lamp $nlist ); then
   gecho "pass!\n"
   bcheck="PASS"
else
   recho "FAIL!\n"
   bcheck="FAIL"
fi

## Check for existing barrier violation:
if [ "$bcheck" = "FAIL" ]; then
   if [ -f $nite_lampsave ]; then
      yecho "Found existing (bogus) lamp:\n"
      recho "--> $nite_lampsave \n\n"
      cmde "rm $nite_lampsave"
   fi

   rm -rf $tmp_dir
   exit 1
fi

##--------------------------------------------------------------------------##
##                Existing Image Removal: Version Check                     ##
##--------------------------------------------------------------------------##

if [ -f $nite_lampsave ]; then
   echo "Data version requirements: ${need_data_versions[*]}"
   echo "Code version requirements: ${need_code_versions[*]}"
   if ( header_gain_is_unity $nite_lampsave ) && \
      ( data_version_pass $nite_lampsave ${need_data_versions[*]} ) && \
      ( code_version_pass $nite_lampsave ${need_code_versions[*]} ); then
      # image is OK!
      Gecho "Existing $nite_lampsave passed version check!\n"
   else
      Recho "Existing $nite_lampsave FAILED version check!\n"
      cmde "rm $nite_lampsave"
   fi
fi

##--------------------------------------------------------------------------##
##             Existing Image Removal: Better Data Available                ##
##--------------------------------------------------------------------------##

# After multiple iterations calibration build procedures, newer versions of
# some files may appear. The purpose of this section is two-fold:
# 1) For each listed dark, remove non-current 'clean' intermediate files. This
#     ensures there is only one cleaned copy of each image.
# 2) If a nite_lampsave already exists, check its content (the 'quality' of cleaned
#     images that went into making it) against clean images available now. If
#     today's files are different and *sufficiently* better, the existing
#     master should be removed and rebuilt using these superior files.

echo 
echo
echo
echo "lamp_list:"
for item in "${lamp_list[@]}"; do echo "--> $item"; done
echo
echo
#exit

if [ -f $nite_lampsave ]; then
   yecho "\nChecking calibration inputs of $nite_lampsave ...\n"
   have_cals_list="$tmp_dir/have_cals.txt"
   want_cals_list="$tmp_dir/want_cals.txt"
   both_cals_list="$tmp_dir/both_sets.txt"
   rm $have_cals_list $want_cals_list 2>/dev/null

   # List of calibrations from nite_lampsave history:
   get_history_calibs $nite_lampsave > $have_cals_list
   #cmde "cat $have_cals_list"

   # Go file-by-file and build a list of expected calibration inputs:
   for image in "${lamp_list[@]}"; do
      ibase="${image##*/}"
      for ctype in bias dark; do
         imcal=$(pick_best_bdcal $image $camid $ctype)
         cbase="${imcal##*/}"
         echo "$cbase" >> $want_cals_list
      done
   done
   paste -d' ' $have_cals_list $want_cals_list > $both_cals_list
   #cmde "cat $both_cals_list"
   exec 10<$both_cals_list
   while read chave cwant <&10; do
      have_width=$(get_calib_width $chave)
      want_width=$(get_calib_width $cwant)
      if [ $want_width -gt $have_width ]; then
         #Recho "DISCREPANCY DETECTED!\n"
         Recho "Out-of-date inputs detected, rebuild needed!\n"
         echo "chave: $chave, have_width: $have_width"
         echo "cwant: $cwant, want_width: $want_width"
         echo
         cmde "rm $nite_lampsave"
         break
      fi
   done
   exec 10>&-
   rm $have_cals_list $want_cals_list $both_cals_list 2>/dev/null
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

      # FIRST, identify best current bias frame:
      yecho "Selecting good bias ... "
      use_bias=$(pick_best_bdcal $image $camid bias --prev) || exit $?
      if [ -z "$use_bias" ] || [ ! -f $use_bias ]; then
         Recho "Blank name or file missing: '$use_bias'\n\n" >&2
         exit 1
      fi
      bdr_tag=$(imhget DRTAG $use_bias)
      gecho "done.\n"
      echo "use_bias: $use_bias"

      # NEXT, identify best current dark frame:
      yecho "Selecting good dark ... "
      use_dark=$(pick_best_bdcal $image $camid dark --prev) || exit $?
      if [ -z "$use_dark" ] || [ ! -f $use_dark ]; then
         Recho "Blank name or file missing: '$use_dark'\n\n" >&2
         exit 1
      fi
      ddr_tag=$(imhget DRTAG $use_dark)
      gecho "done.\n"
      echo "use_dark: $use_dark"

      # Temporary 'clean' file name:
      ibase="${image##*/}"
      ifits="${ibase%.fz}"
      #cbase="clean_${drtag}_${ifits}.fz"
      cbase="clean_${ifits}.fz"
      isave="$tmp_dir/$cbase"

      # Use existing cleaned lamp if possible:
      icheck="$(get_save_folder $image)/$cbase"
      if [ -f $icheck ]; then
         yecho "\nChecking ${icheck##*/} ... "
         if ( header_gain_is_unity $icheck ) && \
            ( clean_width_check_passed $icheck $use_bias bias ) && \
            ( clean_width_check_passed $icheck $use_dark dark ) && \
            ( data_version_pass $icheck ${need_data_versions[*]} ) && \
            ( code_version_pass $icheck ${need_code_versions[*]} ); then
            Gecho "version/width check PASSED!\n"
            gecho "Using existing temp-lamp (${access_mode}): ${icheck}\n"
            case $access_mode in
               copy)    vcmde "cp -f $icheck $isave" || exit $?  ;;
               symlink) vcmde "ln -s $icheck $isave" || exit $?  ;;
               *) ErrorAbort "Unhandled access_mode: '$access_mode'" ;;
            esac
            continue
         else
            recho "version/width check FAILED, rebuild!\n\n"
            #[ $keep_clean -eq 1 ] && cmde "rm $icheck"
         fi
      fi


      # --------------------------------------------
      # Otherwise, create cleaned file for stacking:
      # --------------------------------------------

      # Versions from input images:
      min_bias_data_vers=$(find_min_cal_version -b $use_bias $use_dark)
      min_bias_code_vers=$(find_min_cal_version -B $use_bias $use_dark)
      min_dark_data_vers=$(find_min_cal_version -d $use_dark)
      min_dark_code_vers=$(find_min_cal_version -D $use_dark)
      echo "min_bias_data_vers: $min_bias_data_vers"
      echo "min_bias_code_vers: $min_bias_code_vers"
      echo "min_dark_data_vers: $min_dark_data_vers"
      echo "min_dark_code_vers: $min_dark_code_vers"

      # Subtract overscan:
      cmde "nres-cdp-trim-oscan -q $image -o $foo"       || exit $?
      lampexp="$(imhget EXPTIME $foo)"
      echo "lampexp: $lampexp"

      # Debias and subtract scaled dark:
      cmde "nres-labcal $foo -b $use_bias -d $use_dark -o $bar"   || exit $?
      cmde "hdrtool $bar --add_hist='use_bias ${use_bias##*/}'"   || exit $?
      cmde "hdrtool $bar --add_hist='use_dark ${use_dark##*/}'"   || exit $?
      cmde "hdrtool $bar -U OBJECTS --value='$lampobj'"           || exit $?

      # Store version information:
      cmde "record_code_version $bar -b $min_bias_code_vers"      || exit $?
      cmde "record_code_version $bar -d $min_dark_code_vers"      || exit $?
      cmde "record_code_version $bar -l ${script_version}"        || exit $?
      cmde "record_data_version $bar -b $min_bias_data_vers"      || exit $?
      cmde "record_data_version $bar -d $min_dark_data_vers"      || exit $?
      cmde "record_data_version $bar -l ${script_version}"        || exit $?

      #cmde "record_cal_version $foo -l $script_version"           || exit $?
      hargs=( $camid $obstype $lampexp "none" )
      cmde "update_output_header $foo ${hargs[*]}"                || exit $?
      cmde "fpack -F -Y -qt 32 $foo"                              || exit $?
      cmde "mv -f $bar $isave"                                    || exit $?

      # Preserve files (if requested):
      if [ $keep_clean -eq 1 ]; then
         vcmde "mkdir -p $(dirname $icheck)"                      || exit $?
         cmde "cp -f $isave $icheck"                              || exit $?
      fi
      echo
   done
   timer

   # Identify minimum code/data versions from input file collection:
   min_bias_data_vers=$(find_min_cal_version -b $tmp_dir/clean*fits.fz)
   min_bias_code_vers=$(find_min_cal_version -B $tmp_dir/clean*fits.fz)
   min_dark_data_vers=$(find_min_cal_version -d $tmp_dir/clean*fits.fz)
   min_dark_code_vers=$(find_min_cal_version -D $tmp_dir/clean*fits.fz)
   min_lamp_data_vers=$(find_min_cal_version -l $tmp_dir/clean*fits.fz)
   min_lamp_code_vers=$(find_min_cal_version -L $tmp_dir/clean*fits.fz)
   echo "min_bias_data_vers: $min_bias_data_vers"
   echo "min_bias_code_vers: $min_bias_code_vers"
   echo "min_dark_data_vers: $min_dark_data_vers"
   echo "min_dark_code_vers: $min_dark_code_vers"
   echo "min_lamp_data_vers: $min_lamp_data_vers"
   echo "min_lamp_code_vers: $min_lamp_code_vers"

   # Get lamp exposure time from input frames:
   #echo "STACKED EXPTIME ..."
   avg_expt=$(calc_average_exptime $tmp_dir/clean*fits.fz)  || exit $?
   echo "avg_expt: $avg_expt"
   #exit 1

   # Combine darks with outlier rejection (stack_args in config.sh):
   mecho "\n`RowWrite 75 -`\n"
   opts="$lamp_stack_args"
   cmde "medianize $opts $tmp_dir/clean*fits.fz -o '!$foo'" || exit $?
   append_input_histories $foo $tmp_dir/clean*fits.fz       || exit $?
   timer

   # Add stats and identifiers to header:
   cmde "fitsperc -qS $foo --parallel3"                     || exit $?
   cmde "kimstat -qSC9 $foo"                                || exit $?
   cmde "hdrtool $foo -U OBJECTS --value='$lampobj'"        || exit $?
   #cmde "record_cal_version $foo -l $script_version"        || exit $?

   cmde "record_data_version $foo -b $min_bias_data_vers"   || exit $?
   cmde "record_code_version $foo -b $min_bias_code_vers"   || exit $?
   cmde "record_data_version $foo -d $min_dark_data_vers"   || exit $?
   cmde "record_code_version $foo -d $min_dark_code_vers"   || exit $?
   cmde "record_data_version $foo -l $min_lamp_data_vers"   || exit $?
   cmde "record_code_version $foo -l $script_version"       || exit $?

   hargs=( $camid $obstype $avg_expt $drtag )
   cmde "update_output_header $foo ${hargs[*]}"             || exit $?
   cmde "fpack -F -Y -qt 32 $foo"                           || exit $?
   cmde "mv -f $foo $nite_lampsave"                         || exit $?

   ## Preserve stack files (if requested):
   #if [ $keep_clean -eq 1 ]; then
   #   yecho "Preserving stack ...\n"
   #   for item in $tmp_dir/clean*fits; do
   #      dst_file="$(get_save_folder $item)/${item##*/}"
   #      if [ ! -f $dst_file ]; then
   #         cmde "mv -f $item $dst_file" || exit $?
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
# CHANGELOG (03_create_master_lamp.sh):
#---------------------------------------------------------------------
#
#  2018-08-07:
#     -- Increased script_version to 0.44.
#     -- Now check for GAIN=1.0 along with code/data version requirements in
#           both clean_ image inspection and for existing image removal.
#
#  2018-08-06:
#     -- Increased script_version to 0.43.
#     -- Output master/stacked file names now include camid.
#
#  2018-08-03:
#     -- Increased script_version to 0.42.
#     -- Gain is now corrected in nres-cdp-trim-overscan.py. Calib build
#           script versions are increased to force data rebuild.
#
#  2018-02-11:
#     -- Increased script_version to 0.41.
#     -- Added additional check and sensible messages in case of no (or few)
#           calibrations found of *any* type (before OBJECTS selection).
#
#  2018-02-09:
#     -- Increased script_version to 0.40.
#     -- Implemented data/code version propagation scheme.
#     -- Implemented dual code/data check on existing stacked and clean images.
#
#  2018-01-07:
#     -- Increased script_version to 0.36.
#     -- Implemented separate version requirements for clean and stack data.
#     -- Added check for in-bounds fdate.
#
#  2018-01-05:
#     -- Increased script_version to 0.35.
#     -- Now check existing master and 'clean' images for version compliance.
#           Offending files are automatically removed and rebuilt.
#
#  2018-01-05:
#     -- Increased script_version to 0.33.
#     -- Current script_version now recorded to LAMPVERS keyword.
#     -- Now use new 'aux' and 'func' locations for common code and routines.
#
#  2017-08-07:
#     -- Increased script_version to 0.30.
#     -- Implemented stacking barriers.
#
#  2017-07-27:
#     -- Increased script_version to 0.26.
#     -- Now propagate OBJECTS into all output files.
#     -- Fixed mis-ordered arguments to update_output_header (yikes!).
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
