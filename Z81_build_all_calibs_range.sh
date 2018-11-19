#!/bin/bash
#
#    Build master biases, darks, and lamps for the past few days.
#
# Rob Siverd
# Created:      2017-07-26
# Last updated: 2017-07-26
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Default options:
debug=0 ; clobber=0 ; force=0 ; timer=0 ; vlevel=0
script_version="0.10"
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
need_exec+=( awk cat FuncDef sed tr )
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

## Check for arguments:
usage () { 
   #Recho "\nSyntax: $this_prog --START\n\n"
   Recho "\nSyntax: $this_prog camera day-obs1 day-obs2\n\n"
   yecho "Input date range should correspond to the day-obs of lamp/target \n"
   yecho "frames (i.e., after DAY-OBS rollover). The corresponding biases and\n"
   yecho "darks for lamp/target spectra are obtained on the previous DAY-OBS.\n"
}
#if [ "$1" != "--START" ]; then
if [ -z "$3" ]; then
   usage >&2
   exit 1
fi
camid="$1"
fdate1="$2"
fdate2="$3"

## FDATE sanity checks:
for item in $fdate1 $fdate2 ; do
   if [ ${#item} -ne 8 ]; then
      recho "Invalid FDATE: $item \n" >&2
      exit 1
   fi
done

##--------------------------------------------------------------------------##
## Date/time manipulation:
date_funcs="func/01_time_and_date.sh"
[ -f $date_funcs ] || ErrorAbort "Can't find file: $date_funcs"
vcmde "source $date_funcs"

##**************************************************************************##
##==========================================================================##
##--------------------------------------------------------------------------##

### List of recent nights:
#nights=5
#declare -a nite_list
#for (( x = $nights; x >= -1; x-- )); do
#   nite_list+=( $(date -u +%Y%m%d --date="$x days ago") )
#done
#
#echo "today: $today"

## Check for required script(s):
build_range="./Z01_rebuild_range.sh"
[ -f $build_range ] || ErrorAbort "Can't find file: $build_range"

##--------------------------------------------------------------------------##
## Versions to build (look-back days, keep in order):
#lookbacks=( 0 2 4 )
lookbacks=( 0 2 4 6 )
max_back=$(echo ${lookbacks[*]} | tr ' ' '\n' | sort -nk1 | tail -n1)
echo "max_back: $max_back"

## What types of lamp frames to build:
lamptypes=( tung01 tung12 thar01 thar12 )

##--------------------------------------------------------------------------##
## Bias/dark DAY-OBS is shifted by -1 relative to lamp/target:
bdday1=$(get_prev_day $fdate1)
bdday2=$(get_prev_day $fdate2)

## Flat/double DAY-OBS matches the TARGET date range:
lampday1=$fdate1
lampday2=$fdate2

## ALSO need master biases and darks:

##--------------------------------------------------------------------------##
## Build master biases (compensate for look-back time):
for nprev in ${lookbacks[*]}; do
   delta=$(( max_back - nprev ))
   day1=$(day_change $bdday1 -$delta)
   day2=$bdday2
   use_args="bias -k -B$nprev"
   cmde "$build_range $camid $day1 $day2 $use_args" || exit $?
done

## Build master darks:
for nprev in ${lookbacks[*]}; do
   delta=$(( max_back - nprev ))
   day1=$(day_change $bdday1 -$delta)
   day2=$bdday2
   use_args="dark -k -B$nprev"
   cmde "$build_range $camid $day1 $day2 $use_args" || exit $?
done

##--------------------------------------------------------------------------##
## Build master lamps (tung01):
#cmde "$build_range $camid $lampday1 $lampday2 lamp -k -B0 --tung01" || exit $?
#cmde "$build_range $camid $lampday1 $lampday2 lamp -k -B2 --tung01" || exit $?
#cmde "$build_range $camid $lampday1 $lampday2 lamp -k -B4 --tung01" || exit $?

## Build master lamps (tung12):
#cmde "$build_range $camid $lampday1 $lampday2 lamp -k -B0 --tung12" || exit $?
#cmde "$build_range $camid $lampday1 $lampday2 lamp -k -B2 --tung12" || exit $?
#cmde "$build_range $camid $lampday1 $lampday2 lamp -k -B4 --tung12" || exit $?

## Build master lamps (thar01):
#cmde "$build_range $camid $lampday1 $lampday2 lamp -k -B0 --thar01" || exit $?
#cmde "$build_range $camid $lampday1 $lampday2 lamp -k -B2 --thar01" || exit $?
#cmde "$build_range $camid $lampday1 $lampday2 lamp -k -B4 --thar01" || exit $?

## Build master lamps (thar12):
#cmde "$build_range $camid $lampday1 $lampday2 lamp -k -B0 --thar12" || exit $?
#cmde "$build_range $camid $lampday1 $lampday2 lamp -k -B2 --thar12" || exit $?
#cmde "$build_range $camid $lampday1 $lampday2 lamp -k -B4 --thar12" || exit $?

##--------------------------------------------------------------------------##
## Build master lamps (all types):
for ltype in ${lamptypes[*]}; do
   for nprev in ${lookbacks[*]}; do
      use_args="lamp -k -B$nprev --$ltype"
      cmde "$build_range $camid $lampday1 $lampday2 $use_args" || exit $?
   done
done

##--------------------------------------------------------------------------##
## Clean up:
#[ -d $tmp_dir ] && [ -O $tmp_dir ] && rm -rf $tmp_dir
#[ -f $foo ] && rm -f $foo
#[ -f $bar ] && rm -f $bar
#[ -f $baz ] && rm -f $baz
#[ -f $qux ] && rm -f $qux
exit 0

######################################################################
# CHANGELOG (Z81_build_all_calibs_range.sh):
#---------------------------------------------------------------------
#
#  2017-07-26:
#     -- Increased script_version to 0.10.
#     -- First created Z81_build_all_calibs_range.sh.
#
