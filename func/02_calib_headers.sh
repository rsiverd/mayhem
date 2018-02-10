#!/bin/bash
#
#    This file contains helper routines related to setting and checking of
# key header content in calibration output files.
#
# Rob Siverd
# Created:      2018-01-05
# Last updated: 2018-02-09
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

##--------------------------------------------------------------------------##
## Output image header updates:
update_output_header () {
   local image="$1"
   local camid="$2"
   local obstype="$3"
   local exptime="$4"
   local drtag="$5"
   vcmde "hdrtool $image -d"                                || return $?
   vcmde "hdrtool $image -U OBSTYPE  --value='$obstype'"    || return $?
   vcmde "hdrtool $image -U EXPTIME  --value='$exptime'"    || return $?
   vcmde "hdrtool $image -U CAMID    --value='$camid'"      || return $?
   vcmde "hdrtool $image -U DRTAG    --value='$drtag'"      || return $?
   vcmde "hdrtool $image -U EXTNAME  --value='SPECTRUM'"    || return $?
   vcmde "hdrtool $image -U TELESCOP --value='$camid'"      || return $?
   vcmde "hdrtool $image -Cd"                               || return $?
}

##--------------------------------------------------------------------------##
## Set script version:
record_scr_version () {
   local image="$1"
   local version="$3"
   [ -z "$3" ] && ErrorAbort "Too few args for ${FUNCNAME[0]} ..." 99
   case $2 in
      -b) caltype="bias"; kword="BSCRIPT" ;;
      -d) caltype="dark"; kword="DSCRIPT" ;;
      -l) caltype="lamp"; kword="LSCRIPT" ;;
      *) ErrorAbort "Unhandled argument: '$2'" 99 ;;
   esac
   ctext="$caltype script version"
   hargs="-U $kword --hdr_data=$version --comment='$ctext'"
   vcmde "hdrtool $image $hargs" || return $?
}

## Get script_version:
get_scr_versions () {
   results=( `imhget -u $1 BSCRIPT DSCRIPT LSCRIPT` ) || exit $?
   echo ${results[*]} | sed 's/___/0.0/g'
}

##--------------------------------------------------------------------------##
## Set effective data version:
record_eff_version () {
   local image="$1"
   local version="$3"
   [ -z "$3" ] && ErrorAbort "Too few args for ${FUNCNAME[0]} ..." 99
   case $2 in
      -b) name="bias"; kword="BIASVERS" ;;
      -d) name="dark"; kword="DARKVERS" ;;
      -l) name="lamp"; kword="LAMPVERS" ;;
      *) ErrorAbort "Unhandled argument: '$2'" 99 ;;
   esac
   ctext="effective $name input version"
   hargs="-U $kword --hdr_data=$3 --comment='$ctext'"
   vcmde "hdrtool $image $hargs" || return $?
}

## Get effective data version:
get_eff_versions () {
   results=( `imhget -u $1 BIASVERS DARKVERS LAMPVERS` ) || exit $?
   echo ${results[*]} | sed 's/___/0.0/g'
}

##--------------------------------------------------------------------------##
## Select appropriate version dependencies for specified cal type:
get_version_subset () {
   [ $# -ne 4 ] && ErrorAbort "Wrong args for ${FUNCNAME[0]} ..." 99
   case $1 in
      -b) name="bias"; nkept=1 ;;
      -d) name="dark"; nkept=2 ;;
      -l) name="lamp"; nkept=3 ;;
      *) ErrorAbort "Unhandled argument: '$2'" 99 ;;
   esac
   shift
   keep=()
   for (( x = 0; x < $nkept; x++ )); do
      keep+=( $1 ); shift
   done
   echo ${keep[*]}
}

##--------------------------------------------------------------------------##
## Find minimum version for an ensemble of images:
find_min_cal_version () {
   #local image="$1"
   #local version="$3"
   [ -z "$tmp_dir" ] && ErrorAbort "Temp directory unknown!!" 99
   [ -z "$2" ] && ErrorAbort "Too few args for ${FUNCNAME[0]} ..." 99
   case $1 in
      -b) name="bias"; kword="BIASVERS" ;;
      -d) name="dark"; kword="DARKVERS" ;;
      -l) name="lamp"; kword="LAMPVERS" ;;
      -B) name="bias"; kword="BSCRIPT"  ;;
      -D) name="dark"; kword="DSCRIPT"  ;;
      -L) name="lamp"; kword="LSCRIPT"  ;;
      *) ErrorAbort "Unhandled argument: '$2'" 99 ;;
   esac
   shift
   imlist=( $* )
   #cmde "imhget -n $kword ${imlist[*]}"
   now_sec=$(date +%s.%N)
   tmp_cals="$tmp_dir/available_cal_versions_${now_sec}.txt"
   #cmde "imhget -n $kword ${imlist[*]} | awk '{print \$2}' | sed 's/___/0.0/g'"
   vcmde "imhget -n $kword ${imlist[*]} -o $tmp_cals"
   awk '{print $2}' $tmp_cals | sed 's/___/0.0/g' | sort -n | head -1
   vcmde "rm $tmp_cals"
   #cmde "imhget -n $kword ${imlist[*]} | awk '{print \$2}' | sed 's/___/0.0/g'"
   #imhget -n $kword ${imlist[*]} | awk '{print \$2}' | sed 's/___/0.0/g' |

}

##--------------------------------------------------------------------------##
## Check calib versions:
eff_version_pass () {
   [ -z "$2" ] && ErrorAbort "Too few args for ${FUNCNAME[0]} ..." 99
   [ $# -gt 4 ] && ErrorAbort "Too many args for ${FUNCNAME[0]} ..." 99
   local image=$1; shift
   needvers=( $* )
   versions=( `get_eff_versions $image` )
   for (( i = 0; i < $#; i++ )); do
      need="${needvers[i]}"
      have="${versions[i]}"
      okay=`bc <<< "$have >= $need"`
      #echo "i: $i, need: $need, have: $have, okay: $okay"
      [ $okay -eq 0 ] && return 1
   done
   return 0 # everything passed
}

##--------------------------------------------------------------------------##
## Collect unique input image history from a set of input images:
append_input_histories () {
   local dst_image="$1"; shift
   [ -z "$tmp_dir" ] && ErrorAbort "Temp directory unknown!!" 99
   tmp_hist="$tmp_dir/input_histories_${now_sec}.txt"
   for image in $*; do
      vcmde "listhead $image | grep '^HISTORY use_' | colrm 1 8 >> $tmp_hist"
   done
   exec 10<$tmp_hist
   while read item <&10; do
      vcmde "hdrtool $foo --add_hist='$item'" || return $?
   done
   exec 10>&-
   vcmde "hdrtool -d $foo" || return $?
   vcmde "rm $tmp_hist"
}


######################################################################
# CHANGELOG (02_calib_headers.sh):
#---------------------------------------------------------------------
#
#  2018-02-09:
#     -- find_min_cal_version() now works for script version, not just data.
#     -- Changes to improve clarity and distinguish script/data versions:
#
#  2018-01-07:
#     -- Added dividers to start/end of update_output_header().
#     -- Changes to improve clarity and distinguish script/data versions:
#           --> record_cal_version() now called record_eff_version()
#           --> get_cal_versions() now called get_eff_versions()
#           --> cal_version_pass() now called eff_version_pass()
#     -- Created record_scr_version() and get_scr_versions() routines.
#     -- Created find_min_cal_version() function. This retrieves the requested
#           version (bias/dark/lamp) from a set of FITS files and returns the
#           lowest number seen.
#
#  2018-01-05:
#     -- Added update_output_header() from config.sh.
#     -- First created 02_calib_headers.sh.
#
