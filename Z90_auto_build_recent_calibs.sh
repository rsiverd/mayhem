#!/bin/bash
#
#    Rebuild calibrations (Z81 driver script) for the specified camera and
# lookback time (days). Inputs:
#  * camera_id (e.g., nres01, nres03)
#  # lookback_days (e.g., 5)
#
# Rob Siverd
# Created:      2017-08-14
# Last updated: 2019-02-19
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Default options:
debug=0 ; clobber=0 ; force=0 ; timer=0 ; vlevel=0
script_version="0.20"
this_prog="${0##*/}"
#shopt -s nullglob
# Propagate errors through pipelines: set -o pipefail
# Exit if uninitialized variable used (set -u): set -o nounset
# Exit in case of nonzero status (set -e): set -o errexit

## Program options:
#save_file=""
#shuffle=0
#confirmed=0
#days_back=5

## Standard scratch files/dirs:
tmp_name="$(date +%Y%m%d.%H%M%S).$$.$(whoami)"
tmp_root="/tmp"
[ -d /dev/shm ] && [ -w /dev/shm ] && tmp_root="/dev/shm"
tmp_dir="$tmp_root"
#tmp_dir="$tmp_root/$tmp_name"
#mkdir -p $tmp_dir
foo="$tmp_dir/foo_$$.txt"
bar="$tmp_dir/bar_$$.txt"
baz="$tmp_dir/baz_$$.fits"
qux="$tmp_dir/qux_$$.fits"
jnk="$foo $bar $baz $qux"  # working copy
def_jnk="$jnk"             # original set
dir_cleanup='(echo -e "\nAutomatic clean up ... " ; cmde "rm -vrf $tmp_dir")'
jnk_cleanup='for X in $jnk ; do [ -f $X ] && cmde "rm -vf $X" ; done'
trap "$jnk_cleanup" EXIT
##trap '[ -d $tmp_dir ] && cmde "rm -vrf $tmp_dir"' EXIT
#trap "[ -d $tmp_dir ] && $dir_cleanup" EXIT
#trap "[ -d $tmp_dir ] && $dir_cleanup ; $jnk_cleanup" EXIT
#trap 'oops=$? ; echo ; exit $oops' HUP INT TERM

## Required programs:
declare -a need_exec
need_exec+=( awk cat pipeline-dayobs-logic FuncDef sed tr )
#need_exec+=( shuf shuffle sort ) # for randomization
for need in ${need_exec[*]}; do
   if ! ( /usr/bin/which $need >& /dev/null ); then
      echo "Error:  can't find '$need' in PATH !!" >&2
      exit 1
   fi
done

## Helper function definitions:
fd_args="--argchk --nocolors --cmde --echo"
#fd_args+=" --Critical"
fd_args+=" --rowwrite"
#fd_args+=" --timers"
fd_args+=" --warnings"
FuncDef $fd_args >/dev/null || exit $?
eval "$(FuncDef $fd_args)"  || exit $?

## Check for arguments:
usage () { 
   #Recho "\nSyntax: $this_prog --START\n\n"
   Recho "\nSyntax: $this_prog camera_id days_back\n\n"
}
#if [ "$1" != "--START" ]; then
if [ -z "$2" ]; then
   usage >&2
   exit 1
fi
camid="$1"
ndays="$2"

##--------------------------------------------------------------------------##
## Date/time manipulation:
date_funcs="func/01_time_and_date.sh"
[ -f $date_funcs ] || ErrorAbort "Can't find file: $date_funcs"
vcmde "source $date_funcs"

## Post-restart time buffer (number of hours spent taking arcs/flats):
buff_hrs=6

## Disk/CPU kindness:
prefix="nice -19 ionice -c3"

##**************************************************************************##
##==========================================================================##
##--------------------------------------------------------------------------##

## Calib-builder script is required:
cal_script="Z81_build_all_calibs_range.sh"
[ -f $cal_script ] || ErrorAbort "Can't find file: $cal_script"

## Load site<-->camera mapplings:
site_cams="aux/site_camera_mappings.sh"
[ -f $site_cams ] || ErrorAbort "Can't find file: $site_cams" 
eval "$(cat $site_cams)"

#echo "camera_from_site:"
#echo "${camera_from_site[@]}"

lsite=${site_from_camera[$camid]}
echo "Selected camera:  $camid"
echo "LCO network site: $lsite"
echo "Cal buffer hours: $buff_hrs"

## Ingest results from choose-dayobs. The proc_fdate result should be the
## most recent fdate that meets 'done' criteria for the specified camera:
arrname="dchoice"
pdl_opts="-q --nres-cal-delay $buff_hrs -s $lsite --bash dchoice"
#eval "$(pipeline-dayobs-logic -q -b $buff_hrs -s $lsite --bash dchoice)"
eval "$(pipeline-dayobs-logic $pdl_opts)"
final_fdate="${dchoice[latest_nres_cal_fdate]}"
#echo "final_fdate: $final_fdate"
start_fdate=$(day_change $final_fdate -$ndays)
echo "Building date range: $start_fdate ... $final_fdate"
sleep 5

## Build all calibrations for the specified range:
echo
cmde "$prefix ./$cal_script $camid $start_fdate $final_fdate"

##--------------------------------------------------------------------------##
## Clean up:
#[ -d $tmp_dir ] && [ -O $tmp_dir ] && rm -rf $tmp_dir
#[ -f $foo ] && rm -f $foo
#[ -f $bar ] && rm -f $bar
#[ -f $baz ] && rm -f $baz
#[ -f $qux ] && rm -f $qux
exit 0

######################################################################
# CHANGELOG (Z90_auto_build_recent_calibs.sh):
#---------------------------------------------------------------------
#
#  2018-02-19:
#     -- Increased script_version to 0.20.
#     -- Script essentially works as intended. Hooray!
#
#  2018-02-18:
#     -- Increased script_version to 0.15.
#     -- Resumed implementation using new choose-dayobs Python script.
#
#  2017-08-14:
#     -- Increased script_version to 0.10.
#     -- First created Z90_auto_build_recent_calibs.sh.
#
