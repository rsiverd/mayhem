#!/bin/bash
#
#    This file contains helper routines related to setting and checking of
# key header content in calibration output files.
#
# Rob Siverd
# Created:      2018-01-05
# Last updated: 2018-01-05
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
   #local mversion="$6"
   vcmde "hdrtool $image  -U OBSTYPE  --value='$obstype'"   || return $?
   vcmde "hdrtool $image  -U EXPTIME  --value='$exptime'"   || return $?
   vcmde "hdrtool $image  -U CAMID    --value='$camid'"     || return $?
   vcmde "hdrtool $image  -U DRTAG    --value='$drtag'"     || return $?
   #vcmde "hdrtool $image  -U MVERSION --value='$mversion'"  || return $?
   #vcmde "hdrtool $image  -U BDRTAG   --value='$bdr_tag'"   || return $?
   #vcmde "hdrtool $image  -U DDRTAG   --value='$ddr_tag'"   || return $?
   vcmde "hdrtool $image  -U EXTNAME  --value='SPECTRUM'"   || return $?
   vcmde "hdrtool $image -CU TELESCOP --value='$camid'"     || return $?
}

##--------------------------------------------------------------------------##
## Set versioning keywords:
record_cal_version () {
   local image="$1"
   local version="$3"
   [ -z "$3" ] && ErrorAbort "Too few args for record_cal_version() ..."
   case $2 in
      -b) name="bias"; kword="BIASVERS" ;;
      -d) name="dark"; kword="DARKVERS" ;;
      -l) name="lamp"; kword="LAMPVERS" ;;
      *) ErrorAbort "Unhandled argument: '$2'" ;;
      #*) Recho "lolwut: '$2'"; exit 1 ;;
   esac
   ctext="Master $name script version"
   hargs="-U $kword --hdr_data=$3 --comment='$ctext'"
   #echo "hargs: $hargs"
   vcmde "hdrtool $image $hargs" || return $?
}

##--------------------------------------------------------------------------##
## Get calib versions:
get_cal_versions () {
   results=( `imhget -u $1 BIASVERS DARKVERS LAMPVERS` ) || exit $?
   echo ${results[*]} | sed 's/___/0.0/g'
}

## Check calib versions:
cal_version_pass () {
   [ -z "$2" ] && ErrorAbort "Too few args for cal_version_pass() ..."
   [ $# -gt 4 ] && ErrorAbort "Too many args for cal_version_pass() ..."
   local image=$1; shift
   needvers=( $* )
   versions=( `get_cal_versions $image` )
   for (( i = 0; i < $#; i++ )); do
      need="${needvers[i]}"
      have="${versions[i]}"
      okay=`bc <<< "$have >= $need"`
      #echo "i: $i, need: $need, have: $have, okay: $okay"
      [ $okay -eq 0 ] && return 1
   done
   return 0 # everything passed
}

######################################################################
# CHANGELOG (02_calib_headers.sh):
#---------------------------------------------------------------------
#
#  2018-01-05:
#     -- Added update_output_header() from config.sh.
#     -- First created 02_calib_headers.sh.
#
