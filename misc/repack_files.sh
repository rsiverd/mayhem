#!/bin/bash
#
# fpack images that got skipped.
#
# Rob Siverd
# Created:      2018-11-19
# Last updated: 2018-11-19
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
   Recho "\nSyntax: $this_prog --START\n\n"
   #Recho "\nSyntax: $this_prog arg1\n\n"
}
if [ "$1" != "--START" ]; then
#if [ -z "$1" ]; then
   usage >&2
   exit 1
fi

##**************************************************************************##
##==========================================================================##
##--------------------------------------------------------------------------##

fp_args="-D -Y -v -qt 32"

yecho "Finding uncompressed files ... "
img_list=( `find . -type f -name "clean*fits" | sort | tee $foo` )
dir_list=( `find . -type f -name "clean*fits" | cut -d/ -f2 | sort -u` )
gecho "done.\n"
total=${#img_list[*]}
ndirs=${#dir_list[*]}
echo "Found $total images to compress in $ndirs folders."
count=0
ntodo=0
#for image in ${img_list[*]}; do
#   (( count++ ))
#   echo "image: $image"
#   cmde "fpack -D -Y -v -qt 32 $image"
#   [ $ntodo -gt 0 ] && [ $count -ge $ntodo ] && break
#done
for folder in ${dir_list[*]}; do
   (( count++ ))
   echo "folder: $folder"
   #f_images=( `ls $folder/clean*fits 2>/dev/null` )
   f_images=( `grep $folder $foo` )
   nfound=${#f_images[*]}
   echo "Found $nfound images to compress ..."
   if [ $nfound -gt 0 ]; then
      cmde "fpack $fp_args ${f_images[*]}" || exit $?
   fi
   #cmde "ls $folder/clean*fits"
   [ $ntodo -gt 0 ] && [ $count -ge $ntodo ] && break
done

##--------------------------------------------------------------------------##
## Clean up:
#[ -d $tmp_dir ] && [ -O $tmp_dir ] && rm -rf $tmp_dir
[ -f $foo ] && rm -f $foo
#[ -f $bar ] && rm -f $bar
#[ -f $baz ] && rm -f $baz
#[ -f $qux ] && rm -f $qux
exit 0

######################################################################
# CHANGELOG (repack_files.sh):
#---------------------------------------------------------------------
#
#  2018-11-19:
#     -- Increased script_version to 0.01.
#     -- First created repack_files.sh.
#
