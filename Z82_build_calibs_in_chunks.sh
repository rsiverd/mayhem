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
   yecho " * chunkdays --> [DAYS, REQUIRED] number of days to use per chunk\n"
   yecho "        when splitting the processing date range. To prevent bad\n"
   yecho "        performance, we ensure that (chunkdays > maximum lookback)\n"
   yecho " * erase_lag --> [DAYS, OPTIONAL] If specified, this integer is taken\n"
   yecho "        to be an offset (in days) between the chunk in process and a\n"
   yecho "        the nights that have been fully processed (for which 'clean'\n"
   yecho "        images can be safely removed without performance loss)\n"
   yecho "        erase_nites = chunk_nites - erase_lag\n"
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
## Note erase_lag (if provided):
delete_old=0
if [ -n "$5" ]; then
   delete_old=1
   erase_lag=$5

   # Erase lag should be larger than chunk size:
   if [ $erase_lag -lt $chunkdays ]; then
      recho "Error: need erase_lag >= chunk_size ...\n" >&2
      exit 1
   fi
   yecho "User-specified erase lag: ${erase_lag}\n"
fi

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Date/time manipulation:
date_funcs="func/01_time_and_date.sh"
[ -f $date_funcs ] || ErrorAbort "Can't find file: $date_funcs"
vcmde "source $date_funcs"

## Storage and pipeline configuration:
conf_file="config.sh"
yecho "Loading storage and pipeline config ...\n"
[ -f $conf_file ] || ErrorAbort "Can't find file: $conf_file"
cmde "source $conf_file"

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
   cmde "$child_script $camid $chunk_start $chunk_end"
   if [ $delete_old -eq 1 ]; then
      # Start/end of date range to wipe:
      erase_start=`day_change $chunk_start -$erase_lag`
      erase_end=`day_change $chunk_end -$erase_lag`
      #echo "chunk: $chunk_start -> $chunk_end"
      #echo "erase: $erase_start -> $erase_end"
 
      # Full list of FDATEs marked for clean-up:
      #echo "erase_fdates: ${erase_fdates[*]}"
      #echo; echo

      # Wipe clean files in each of the listed FDATEs:
      erase_fdates=( `list_dates_between $erase_start $erase_end` )
      for edate in ${erase_fdates[*]}; do
         proc_save="$save_root/$camid/$edate"
         echo "proc_save: $proc_save"
         if [ -d $proc_save ]; then
            wipe_list=( `ls $proc_save/clean_*fits* 2>/dev/null` )
            tmp_files="${#wipe_list[@]}"
            echo "Found $tmp_files images to remove ..."
            vcmde "rm ${wipe_list[@]}"
         fi
         echo
      done
      #echo "save_root: $save_root"
   fi
   exit 1
done
exec 10>&-

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
