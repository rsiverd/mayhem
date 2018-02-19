#!/bin/bash
#
#    Rebuild calibrations (Z81 driver script) for the specified camera and
# lookback time (days). Inputs:
#  * camera_id (e.g., nres01, nres03)
#  # lookback_days (e.g., 5)
#
# Rob Siverd
# Created:      2017-08-14
# Last updated: 2019-02-18
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
need_exec+=( awk cat choose-dayobs FuncDef sed tr )
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

##**************************************************************************##
##==========================================================================##
##--------------------------------------------------------------------------##

## Load site<-->camera mapplings:
site_cams="aux/site_camera_mappings.sh"
[ -f $site_cams ] || ErrorAbort "Can't find file: $site_cams" 
eval "$(cat $site_cams)"

#echo "camera_from_site:"
#echo "${camera_from_site[@]}"

lsite=${site_from_camera[$camid]}
echo "Selected camera:  $camid"
echo "LCO network site: $lsite"

### Run choose-dayobs and process dict() output:
#choose-dayobs -q -s $lsite | tr -d \{ | tr -d \} | tr ',' '\n' \
#   | tr -d ' ' | tr -d \' > $foo
##cmde "cat $foo"
#proc_fdate=$(grep process_fdate $foo | cut -d: -f2)
##echo "rstring: $rstring"
#echo $rstring | tr ',' '\n' | tr -d ' ' | tr -d "'"

## Run choose-dayobs and process dict() output:
#arrname="dchoice"
#eval "$(choose-dayobs -q -s lsc \
#   | sed -e "s/^{'/'/" -e "s/'}$/'/" -e "s/', '/'\n'/g" \
#   | awk -v arr="$arrname" -F\' '
#   BEGIN { printf "declare -A %s\n", arr }
#   { 
#      printf "%s+=( [\"%s\"]=\"%s\" )\n", arr, $2,$4
#   }')"
#echo "dchoice:"
#for i in "${!dchoice[@]}"; do
#   printf "%20s:%20s\n" $i ${dchoice[$i]}
#done

## Ingest results from choose-dayobs:
arrname="dchoice"
eval "$(choose-dayobs -q -s $lsite --bash dchoice)"

proc_fdate="${dchoice[process_fdate]}"
echo "proc_fdate: $proc_fdate"


## List of recent nights:
declare -a nite_list
for (( x = $ndays; x >= 1; x-- )); do
   nite_list+=( $(date -u +%Y%m%d --date="$x days ago") )
done

echo "nite_list: ${nite_list[*]}"

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
#  2018-02-18:
#     -- Increased script_version to 0.15.
#     -- Resumed implementation using new choose-dayobs Python script.
#
#  2017-08-14:
#     -- Increased script_version to 0.10.
#     -- First created Z90_auto_build_recent_calibs.sh.
#
