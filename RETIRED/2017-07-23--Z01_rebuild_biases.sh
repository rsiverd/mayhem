#!/bin/bash
#
#    Rebuild biases going back to July 1.
#
# Rob Siverd
# Created:      2017-07-20
# Last updated: 2017-07-20
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Default options:
debug=0 ; clobber=0 ; force=0 ; timer=0 ; vlevel=0
script_version="0.01"
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

bias_script="./01_create_master_bias.sh"
dark_script="./02_create_master_dark.sh"

##--------------------------------------------------------------------------##
## Process each listed night:
for fdate in ${fdate_list[*]}; do
   cmde "$bias_script    -B2 $camid $fdate"  || exit $?
   cmde "$bias_script    -B0 $camid $fdate"  || exit $?
   cmde "$dark_script -k -B2 $camid $fdate"  || exit $?
   cmde "$dark_script -k -B0 $camid $fdate"  || exit $?
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
# CHANGELOG (Z01_rebuild_biases.sh):
#---------------------------------------------------------------------
#
#  2017-07-20:
#     -- Increased script_version to 0.01.
#     -- First created Z01_rebuild_biases.sh.
#
