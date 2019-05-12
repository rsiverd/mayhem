#!/bin/bash
#
# Inspect FITS (*.fits and *.fits.fz) files specified in input list
# with fitsverify tool. Remove all files that fail verification.
#
# Rob Siverd
# Created:      2019-05-11
# Last updated: 2019-05-11
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
baz="$tmp_dir/baz_$$.txt"
qux="$tmp_dir/qux_$$.txt"
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
need_exec+=( awk cat fitsverify FuncDef sed tr xargs )
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
   #Recho "\nSyntax: $this_prog folder_with_images\n\n"
   Recho "\nSyntax: $this_prog list_of_FITS_files.txt\n\n"
}
#if [ "$1" != "--START" ]; then
if [ -z "$1" ]; then
   usage >&2
   exit 1
fi
#imdir="$1"
#[ -d $imdir ] || PauseAbort "Can't find directory: $imdir"
imlist="$1"
[ -f $imlist ] || PauseAbort "Can't find file: $imlist"

##**************************************************************************##
##==========================================================================##
##--------------------------------------------------------------------------##

### List images:
#yecho "Listing files ... "
#find $imdir -type f -name "*.fits*" | sort > $foo
#total=$(cat $foo | wc -l)
#gecho "done. Found $total files.\n"
#if [ $total -eq 0 ]; then
#   echo "Nothing to do."
#   exit 0
#fi

## If image list given directly:
total=$(cat $imlist | wc -l)

##--------------------------------------------------------------------------##
## Verify files:
yecho "Verifying data ...\n"
unbuffer cat $imlist | xargs -n10 fitsverify -e -q | tee $bar | \
   awk -v total=$total '{
      printf "\rChecked %d of %d ...  ", NR, total
   } END { printf "done.\n" }'

## Count failures:
grep -v 'OK:' $bar | cut -d' ' -f3 | cut -d, -f1 > $baz
nfail=$(cat $baz | wc -l)
echo "Found $nfail bad images"
if [ $nfail -gt 0 ]; then
   yecho "Removing files ...\n"
   cat $baz | xargs -n10 rm -v
fi

##--------------------------------------------------------------------------##
## Clean up:
#[ -d $tmp_dir ] && [ -O $tmp_dir ] && rm -rf $tmp_dir
[ -f $foo ] && rm -f $foo
[ -f $bar ] && rm -f $bar
[ -f $baz ] && rm -f $baz
[ -f $qux ] && rm -f $qux
exit 0

######################################################################
# CHANGELOG (remove-invalid-fits-files.sh):
#---------------------------------------------------------------------
#
#  2019-05-11:
#     -- Increased script_version to 0.10.
#     -- First created remove-invalid-fits-files.sh.
#
