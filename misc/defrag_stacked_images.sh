#!/bin/bash
#
# Check for images with many extents in the specified folder. Attempt to
# defrag the worst offenders by moving them to /tmp (or /dev/shm) and back.
#
# Rob Siverd
# Created:      2018-12-03
# Last updated: 2018-12-03
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
need_exec+=( awk cat filefrag FuncDef sed tr )
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
   Recho "\nSyntax: $this_prog stack_folder [dfg_extents]\n\n"
}
#if [ "$1" != "--START" ]; then
if [ -z "$1" ]; then
   usage >&2
   exit 1
fi
stack_dir="$1"
max_frags=${2:-10}
[ -d $stack_dir ] || PauseAbort "Can't find directory: $stack_dir"

if ! int_check_pass $max_frags ; then
   Recho "Non-integer fragment count not allowed: '$max_frags'\n" >&2
   exit 1
fi
echo "Fragments limit: $max_frags"

##**************************************************************************##
##==========================================================================##
##--------------------------------------------------------------------------##

##--------------------------------------------------------------------------##
## Count stacked frames:
#Mecho "\n`RowWrite 75 =`\n"
yecho "Checking files ... "
vcmde "ls $stack_dir/med*fz 2>/dev/null > $foo"
nfiles=$(cat $foo | wc -l)
if [ $nfiles -eq 0 ]; then
   yecho "nothing to check!\n"
   rm $foo
   exit 0
fi
yecho "found $nfiles stacked images.\n"

##--------------------------------------------------------------------------##
## Check for fragmentation:
yecho "Checking fragmentation ... "
vcmde "filefrag `cat $foo` > $bar"  || exit $?
vcmde "mv -f $bar $foo"             || exit $?
gecho "done.\n"

## Select files for adjustment:
nbad=$(awk -v fmax=$max_frags '$2 >= fmax' $foo | tee $bar | wc -l)
echo "Over-fragged files: $nbad"
if [ $nbad -eq 0 ]; then
   echo "Nothing to do!"
   rm $foo $bar $baz 2>/dev/null
   exit 0
fi

##--------------------------------------------------------------------------##
## Attempt to defrag:
frag_list=( `cut -d: -f1 $bar` )
for image in "${frag_list[@]}"; do
   ibase="${image##*/}"
   xfile="${stack_dir}/${ibase}.xfer"
   echo "image: $image"
   echo "xfile: $xfile"
   cmde "cp $image $baz"
   cmde "mv -f $baz $xfile" || exit $?

   # Just to be sure ...
   #cmde "filefrag $image $xfile"
   if ( cmde "cmp $image $xfile" ); then
      cmde "mv -f $xfile $image" || exit $?
   else
      PauseAbort "YIKES!\n"
   fi
done

##--------------------------------------------------------------------------##
## Re-check after committment to disk:
cmde "sync" || exit $?
#cmde "filefrag `cat $foo`"

mecho "\nWorst 5 after defrag attempt:\n"
cmde "filefrag $stack_dir/med*fz | sort -rnk2 | head -5"

##--------------------------------------------------------------------------##
## Clean up:
#[ -d $tmp_dir ] && [ -O $tmp_dir ] && rm -rf $tmp_dir
[ -f $foo ] && rm -f $foo
[ -f $bar ] && rm -f $bar
[ -f $baz ] && rm -f $baz
[ -f $qux ] && rm -f $qux
exit 0

######################################################################
# CHANGELOG (defrag_stacked_images.sh):
#---------------------------------------------------------------------
#
#  2018-12-03:
#     -- Increased script_version to 0.10.
#     -- First created defrag_stacked_images.sh.
#
