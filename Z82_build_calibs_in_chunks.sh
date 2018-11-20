#!/bin/bash
#
# This script adds a layer between the user and Z81 calib building. It
# performs two main functions:
# 1) Break the input date range into chunks of some maximum size. This helps
#     the machine(s) involved better exploit disk caching by returning to
#     files before too much time has elapsed or other files are used.
# 2) After a chunk is fully processed, 'clean' images can be removed right
#     away. This significantly decreases the amount of disk space needed for
#     intermediate images while calibrations are building.
#
# Rob Siverd
# Created:      2018-11-19
# Last updated: 2018-11-19
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
trap "$jnk_cleanup >&2" EXIT
##trap '[ -d $tmp_dir ] && cmde "rm -vrf $tmp_dir"' EXIT
#trap "[ -d $tmp_dir ] && $dir_cleanup >&2" EXIT
#trap "[ -d $tmp_dir ] && $dir_cleanup >&2; $jnk_cleanup >&2" EXIT
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
   Recho "\nSyntax: $this_prog arg1\n\n"
}
#if [ "$1" != "--START" ]; then
if [ -z "$1" ]; then
   usage >&2
   exit 1
fi

## Check for arguments:
usage () {
   Recho "\nSyntax: $this_prog camera day-obs1 day-obs2 chunkdays\n\n"
   yecho "Input date range should correspond to the day-obs of lamp/target \n"
   yecho "frames (i.e., after DAY-OBS rollover). The corresponding biases and\n"
   yecho "darks for lamp/target spectra are obtained on the previous DAY-OBS.\n"
   yecho "\nOther arguments:\n"
   yecho " * chunkdays --> chunk size for date range breakdown. For decent\n"
   yecho "        performance, need chunkdays > maximum lookback\n"
   yecho "\n"
}
#if [ "$1" != "--START" ]; then
if [ -z "$4" ]; then
   usage >&2
   exit 1
fi
camid="$1"
fdate1="$2"
fdate2="$3"
chunkdays="$4"

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

## Child script and args:
child_script="./Z81_build_all_calibs_range.sh"
[ -f $child_script ] || PauseAbort "Can't find file: $child_script"

##--------------------------------------------------------------------------##
## Split date range into chunks. Each line gets start/end of a chunk range:
#chunksize=15
fdate_list=( `list_dates_between $fdate1 $fdate2` )
nites_to_chunks_firstlast $chunkdays ${fdate_list[*]} > $foo
#cmde "cat $foo"

## Execute child script using per-chunk date ranges:
## FIXME / TODO: add cleanup after processing each chunk.
exec 10<$foo
while read chunk_start chunk_end <&10; do
   #echo "cdate1,cdate2: $cdate1, $cdate2"
   echo "$child_script $camid $chunk_start $chunk_end"
done
exec 10>&-

## List of recent nights:
#nights=5
#declare -a nite_list
#for (( x = $nights; x >= -1; x-- )); do
#   nite_list+=( $(date -u +%Y%m%d --date="$x days ago") )
#done


##--------------------------------------------------------------------------##
## Clean up:
#[ -d $tmp_dir ] && [ -O $tmp_dir ] && rm -rf $tmp_dir
[ -f $foo ] && rm -f $foo
#[ -f $bar ] && rm -f $bar
#[ -f $baz ] && rm -f $baz
#[ -f $qux ] && rm -f $qux
exit 0

######################################################################
# CHANGELOG (Z82_build_calibs_in_chunks.sh):
#---------------------------------------------------------------------
#
#  2018-11-19:
#     -- Increased script_version to 0.10.
#     -- First created Z82_build_calibs_in_chunks.sh.
#
