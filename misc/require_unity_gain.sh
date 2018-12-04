#!/bin/bash
#
# Inspect images in the specified folder and remove those with non-unity
# gain (reprocessing required).
#
# Rob Siverd
# Created:      2018-08-07
# Last updated: 2018-08-07
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
trap "$jnk_cleanup >&2" EXIT
##trap '[ -d $tmp_dir ] && cmde "rm -vrf $tmp_dir"' EXIT
#trap "[ -d $tmp_dir ] && $dir_cleanup >&2" EXIT
#trap "[ -d $tmp_dir ] && $dir_cleanup >&2; $jnk_cleanup >&2" EXIT
#trap 'oops=$? ; echo ; exit $oops' HUP INT TERM

## Required programs:
declare -a need_exec
need_exec+=( awk cat FuncDef imhget sed tr )
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
   Recho "\nSyntax: $this_prog image_folder\n\n"
}
#if [ "$1" != "--START" ]; then
if [ -z "$1" ]; then
   usage >&2
   exit 1
fi
img_dir="$1"
[ -d $img_dir ] || PauseAbort "Can't find directory: $img_dir"

##**************************************************************************##
##==========================================================================##
##--------------------------------------------------------------------------##

Yecho "\nChecking GAIN values in $img_dir ...\n"

## List images:
find $img_dir -type f -name "*fits*" > $foo
nfiles=$(cat $foo | wc -l)
if [ $nfiles -eq 0 ]; then
   gecho "No images to check!\n"
   rm $foo $bar 2>/dev/null
   exit 0
fi

## Mine headers:
cmde "imhget -n --progress GAIN -l $foo -o $bar" || exit $?
cmde "mv -f $bar $foo"                           || exit $?
#cat $bar

## List images with non-unity gain:
#awk '{ if ($2 != 1.0) { print }}' $foo
awk '/fits/ { if ($2 != 1.0) { print $1 }}' $foo > $bar
nbad=$(cat $bar | wc -l)
echo "nbad: $nbad"
#cat $bar
cmde "cat $bar | xargs -n 50 rm"

##--------------------------------------------------------------------------##
## Clean up:
#[ -d $tmp_dir ] && [ -O $tmp_dir ] && rm -rf $tmp_dir
[ -f $foo ] && rm -f $foo
[ -f $bar ] && rm -f $bar
#[ -f $baz ] && rm -f $baz
#[ -f $qux ] && rm -f $qux
exit 0

######################################################################
# CHANGELOG (require_unity_gain.sh):
#---------------------------------------------------------------------
#
#  2018-08-07:
#     -- Increased script_version to 0.01.
#     -- First created require_unity_gain.sh.
#
