#!/bin/bash
#
# Look for empty files in the specified directories and remove any found.
#
# Rob Siverd
# Created:      2018-11-05
# Last updated: 2018-11-05
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
   Recho "\nSyntax: $this_prog folder1 [... folderN]\n\n"
}
#if [ "$1" != "--START" ]; then
if [ -z "$1" ]; then
   usage >&2
   exit 1
fi

##**************************************************************************##
##==========================================================================##
##--------------------------------------------------------------------------##

#ntodo=2
#count=0
#total=$#
#for imdir in "$@"; do
#   (( ++count ))
#   yecho "imdir: $imdir ($count of $total) ... "
#   images=( `ls $imdir/*fits* 2>/dev/null` )
#   nfits=${#images[*]}
#   if [ $nfits -eq 0 ]; then
#      echo "no images!"
#      continue
#   fi
#   vcmde "du -k $imdir/*fits* > $foo" || exit $?
#   empties=( `awk '$1 == 0 { print $2 }' $foo` )
#   #cmde "file $imdir/*fits* > $foo" || exit $?
#   #empties=( `grep ': empty' $foo | cut -d: -f1` )
#   nempty="${#empties[@]}"
#   echo "found $nempty empty file(s)."
#   if [ $nempty -gt 0 ]; then
#      cmde "rm ${empties[*]}"
#   fi
#done

ntodo=2
count=0
total=$#
for imdir in "$@"; do
   (( ++count ))
   yecho "imdir: $imdir ($count of $total) ... "
   #cleans=( `ls $imdir/clean*fits 2>/dev/null` )
   #nfits=${#cleans[*]}
   images=( `ls $imdir/*fits 2>/dev/null` )
   nfits=${#images[*]}
   if [ $nfits -eq 0 ]; then
      echo "no images!"
      continue
   fi

   ## quick check for already-packed images:
   ##token='HDU #2  Array:  NAXIS = 2,  BITPIX = -32'
   #token="^HDU #2  Array:  NAXIS = 2,  BITPIX = -32"
   #for item in ${images[*]}; do
   #   #cmde "liststruc $item | grep -- '$token'"
   #   liststruc $item | grep "$token"
   #done
      
   for item in ${images[*]}; do
      if ! ( vcmde "fitsverify -q -e $item" ); then
         recho "Image has a problem!!\n"
         exit
      fi
   done
   vcmde "du -k $imdir/*fits* > $foo" || exit $?
   empties=( `awk '$1 == 0 { print $2 }' $foo` )
   #cmde "file $imdir/*fits* > $foo" || exit $?
   #empties=( `grep ': empty' $foo | cut -d: -f1` )
   nempty="${#empties[@]}"
   echo "found $nempty empty file(s)."
   if [ $nempty -gt 0 ]; then
      cmde "rm ${empties[*]}"
   fi
done


##--------------------------------------------------------------------------##
## Clean up:
#[ -d $tmp_dir ] && [ -O $tmp_dir ] && rm -rf $tmp_dir
[ -f $foo ] && rm -f $foo
[ -f $bar ] && rm -f $bar
#[ -f $baz ] && rm -f $baz
#[ -f $qux ] && rm -f $qux
exit 0

######################################################################
# CHANGELOG (purge_empty_image_files.sh):
#---------------------------------------------------------------------
#
#  2018-11-05:
#     -- Increased script_version to 0.10.
#     -- First created purge_empty_image_files.sh.
#
