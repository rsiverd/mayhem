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
#record_cal_version () {
#   local image="$1"
#   local version="$3"
#   [ -z "$3" ] && return 1
#   case $2 in
#      -b) vkey="BIASVERS"
#}


######################################################################
# CHANGELOG (02_calib_headers.sh):
#---------------------------------------------------------------------
#
#  2018-01-05:
#     -- Added update_output_header() from config.sh.
#     -- First created 02_calib_headers.sh.
#
