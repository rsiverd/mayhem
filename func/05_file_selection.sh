#!/bin/bash
#
#    File selection routines for mayhem pipeline. Most of this content was
# formerly kept in config.sh. Keeping these in a separate source file is
# necessary due to dependency ordering constraints.
#
# Rob Siverd
# Created:      2018-02-19
# Last updated: 2018-08-06
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

##--------------------------------------------------------------------------##
## Get destination path for 'clean' files:
get_save_folder () {
   keys=( "TELESCOP" "DAY-OBS" )
   nhave=`imhget -c ${keys[*]} $1`
   if [ $nhave -ne ${#keys[*]} ]; then
      ErrorAbort "Image $1 missing one of ${keys[*]} keys ..." 99
   fi
   hdata=( `imhget ${keys[*]} $1` )
   local camid=${hdata[0]}
   local dayobs=${hdata[1]}
   echo "$save_root/$camid/$dayobs"
}

## Retrieve the path for saving/checking cleaned files:
get_clean_ipath () {
   local image="$1"
   local ibase="${image##*/}"
   clean_dir=`get_save_folder $image`
   echo "clean_dir: $clean_dir"
   clean_image="$clean_dir/$ibase"
}

##--------------------------------------------------------------------------##
## Get specific bias frame path:
get_bias_path () {
   local image="$1"
   local camid="$2"
   local drtag="$3"
   local dayobs=$(imhget DAY-OBS $image)
   echo "$save_root/$camid/$dayobs/med_bias_${dayobs}_${drtag}.fits"
}

##--------------------------------------------------------------------------##
## Select unique files based on DATE-OBS:
pick_unique_subset () {
   local tmpfile="$tmp_dir/hdrdata.$$.txt"
   imhget DATE-OBS $* -o $tmpfile
   unique_dates=( `awk '{ print $2 }' $tmpfile | sort -u` )
   for obsdate in ${unique_dates[*]}; do
      grep -m1 $obsdate $tmpfile | awk '{ print $1 }'
   done
   rm $tmpfile
}

##--------------------------------------------------------------------------##
## Pick the best available bias/dark calibration image (go back up to 30 days):
pick_best_bdcal () {
   local image="$1"
   local camid="$2"
   local ctype="$3"
   local delta="$4"
   local dayobs=$(imhget DAY-OBS $image)  || exit $?
   local maxiter=30  # maximum days back to search:

   # validate inputs:
   if [ -z "$camid" ]; then
      echo -e "\npick_best_bdcal error: camid cannot be blank!" >&2
      return 1
   fi
   case $ctype in
      bias|dark) ;;
      *)
         echo -e "\npick_best_bdcal error: unsupported ctype: '$ctype'" >&2
         return 1
         ;;
   esac

   # offset dayobs in 'yesterday' mode:
   if [ "$delta" = "--prev" ]; then
      dayobs=$(get_prev_day $dayobs)
   fi

   # identify latest-and-greatest master bias/dark:
   local cal_folder
   for (( x = 0; x < $maxiter; x++ )); do
      # select trial day-obs, check for daydir:
      cal_day=$(day_change $dayobs -$x)
      #echo "x: $x, $cal_day"
      cal_folder="$save_root/$camid/$cal_day"
      if [ ! -d $cal_folder ]; then
         echo -e "\npick_best_bdcal warning: folder not found:" >&2
         echo "--> $cal_folder" >&2
         continue
      fi

      # list available files, report latest-and-greatest:
      hits=( `ls $cal_folder/med_${camid}_${ctype}_${cal_day}_*.fits 2>/dev/null` )
      nfound=${#hits[*]}
      #echo "nfound: $nfound" >&2
      if [ $nfound -gt 0 ]; then
         echo "${hits[-1]}"
         return 0
      fi
   done

   # search failure:
   echo -e "\npick_best_bdcal error: no potential images found!" >&2
   return 1
}



##**************************************************************************##
##==========================================================================##
##--------------------------------------------------------------------------##

######################################################################
# CHANGELOG (05_file_selection.sh):
#---------------------------------------------------------------------
#
#  2018-08-06:
#     -- pick_best_bdcal now knows to check for camid in file name.
#
#  2018-02-19:
#     -- Imported most of the functions from config.sh.
#     -- First created 05_file_selection.sh.
#
