#!/bin/bash
#
#    Rebuild biases going back to July 1.
#
# Rob Siverd
# Created:      2017-07-20
# Last updated: 2018-01-07
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
   Recho "\nSyntax: $this_prog camera day-obs1 day-obs2 cal_type\n\n"
   Yecho "cal_type must be one of: bias, dark\n\n"
}
#if [ "$1" != "--START" ]; then
if [ -z "$3" ]; then
   usage >&2
   exit 1
fi

##**************************************************************************##
##==========================================================================##
##--------------------------------------------------------------------------##

camid="$1"
fdate1=$2
fdate2=$3
ctype="$4"
shift 4
cal_args=( "$@" )

## Required external files:
bias_script="./01_create_master_bias.sh"
dark_script="./02_create_master_dark.sh"
lamp_script="./03_create_master_lamp.sh"
[ -f $bias_script ] || ErrorAbort "Can't find file: $bias_script"
[ -f $dark_script ] || ErrorAbort "Can't find file: $dark_script"
[ -f $lamp_script ] || ErrorAbort "Can't find file: $lamp_script"

## Validate calib type:
case $ctype in
   bias)
      run_script="$bias_script"
      ;;
   dark)
      run_script="$dark_script"
      ;;
   lamp)
      run_script="$lamp_script"
      ;;
   *)
      Recho "Unsupported calibration type: '$ctype'\n\n" >&2
      exit 1
      ;;
esac

##--------------------------------------------------------------------------##
## Convert between FDATE and UNIX seconds:
fdate2unix () { date +%s --date=$1; }
unix2fdate () { date +%Y%m%d --date="@$1"; }

## Make list:
next_unixtime=$(fdate2unix $fdate1)
fdate_list=( `unix2fdate $next_unixtime` )
while true; do
   next_unixtime=$(( next_unixtime + 86400 ))
   next_fdate=`unix2fdate $next_unixtime`
   if [ $next_fdate -le $fdate2 ]; then
      fdate_list+=( $next_fdate )
   else
      break
   fi
done

##--------------------------------------------------------------------------##
## Process each listed night:
for fdate in ${fdate_list[*]}; do
   wecho "`RowWrite 75 -`\n"
   cmde "$run_script $camid $fdate ${cal_args[*]}"; status=$?
   if [ $status -gt 0 ]; then
      recho "Got an error ....\n" >&2
   fi
   if [ $status -ge 10 ]; then
      recho "High status encountered ... time to stop.\n" >&2
      exit $status
   fi
done

##--------------------------------------------------------------------------##
## Clean up:
#[ -d $tmp_dir ] && [ -O $tmp_dir ] && rm -rf $tmp_dir
[ -f $foo ] && rm -f $foo
[ -f $bar ] && rm -f $bar
[ -f $baz ] && rm -f $baz
[ -f $qux ] && rm -f $qux
exit 0

######################################################################
# CHANGELOG (Z01_rebuild_range.sh):
#---------------------------------------------------------------------
#
#  2018-01-07:
#     -- Increased script_version to 0.15.
#     -- Now abort processing in case of exit status >= 10.
#
#  2017-07-23:
#     -- Increased script_version to 0.10.
#     -- Added better dividers between tasks to simplify debugging.
#     -- Renamed script Z01_rebuild_biases.sh to Z01_rebuild_range.sh
#           to reflect actual function of current script version.
#
#  2017-07-20:
#     -- Increased script_version to 0.01.
#     -- First created Z01_rebuild_biases.sh.
#
