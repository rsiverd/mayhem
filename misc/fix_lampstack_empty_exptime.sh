#!/bin/bash
#
# Examine stacked lamp images (tung01, tung12, thar01, thar12) in the
# specified folder. Identify stacks with empty EXPTIME and insert sensible
# value (based on 01d variant).
#
# Rob Siverd
# Created:      2018-12-03
# Last updated: 2018-12-03
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
   Recho "\nSyntax: $this_prog stack_folder\n\n"
}
#if [ "$1" != "--START" ]; then
if [ -z "$1" ]; then
   usage >&2
   exit 1
fi
stack_dir="$1"
[ -d $stack_dir ] || PauseAbort "Can't find directory: $stack_dir"

##**************************************************************************##
##==========================================================================##
##--------------------------------------------------------------------------##

## Count stacked frames:
Mecho "\n`RowWrite 75 =`\n"
yecho "Checking files ... "
vcmde "ls $stack_dir/med*fz 2>/dev/null > $foo"
nfiles=$(cat $foo | wc -l)
if [ $nfiles -eq 0 ]; then
   yecho "nothing to check!\n"
   rm $foo
   exit 0
fi
yecho "found $nfiles stacked images.\n"

## Look for problem images:
vcmde "imhget EXPTIME DRTAG -l $foo -o $bar" || exit $?
nbad=$(awk 'NF != 3' $bar | wc -l)
#cmde "cat $bar"
if [ $nbad -eq 0 ]; then
   gecho "everything looks good!\n"
   rm $foo $bar 2>/dev/null
   exit 0
else
   recho "found $nbad images in need of help ...\n"
fi

##--------------------------------------------------------------------------##
## Process each flavor independently (different exposure times):
flavors=( tung01 tung12 thar01 thar12 )
for iflav in ${flavors[*]}; do
   yecho "iflav: $iflav \n"
   cmde "ls $stack_dir/med_*_${iflav}_*fits.fz 2>/dev/null > $foo"
   nimgs=$(cat $foo | wc -l)
   echo "nimgs: $nimgs"
   if [ $nimgs -eq 0 ]; then
      yecho "No $iflav images to process ...\n"
      continue
   fi
   cmde "imhget -l $foo EXPTIME DRTAG -o $bar"
   
   # get DRTAGs:
   grep -o '[0-9][0-9]d_[0-9][0-9]b_[0-9][0-9]a' $foo > $baz

   # count problem files:
   bogus=( `paste $bar $baz | awk 'NF != 4 || $3 != $4' | awk '{print $1}'` )
   nfail=${#bogus[*]}
   if [ $nfail -eq 0 ]; then
      gecho "No failed $iflav images found!\n"
      continue
   fi
   recho "Have $nfail bogus $iflav images ...\n"

   valid=( $(paste $bar $baz | awk 'NF == 4 && $3 == $4' | head -1) )
   nhits=${#valid[*]}
   echo "valid: ${valid[*]}"
   echo "nhits: $nhits"
   if [ $nhits -eq 0 ]; then
      #PauseAbort "NOTHING IS SAFE!"
      Recho "Nothing is safe ... no choice but deletion =(\n"
      #read pause
      cmde "rm ${bogus[*]}"
      #exit
      continue
   fi
   [ $nhits -ne 4 ] && PauseAbort "SOMETHING WENT AWRY ..."

   expsec=${valid[1]}
   echo "expsec: $expsec"

   # Fix files that need fixing:
   recho "BROKEN:\n"
   for image in ${bogus[*]}; do
      echo "image: $image"
      drtag=$(echo $image | grep -o '[0-9][0-9]d_[0-9][0-9]b_[0-9][0-9]a')
      echo "drtag: $drtag"
      #cmde "imhget DRTAG $image"
      cmde "hdrtool $image -U   DRTAG --str_data='$drtag'"
      cmde "hdrtool $image -U EXPTIME --hdr_data=${expsec}"
   done
   echo
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
# CHANGELOG (fix_lampstack_empty_exptime.sh):
#---------------------------------------------------------------------
#
#  2018-12-03:
#     -- Increased script_version to 0.01.
#     -- First created fix_lampstack_empty_exptime.sh.
#
