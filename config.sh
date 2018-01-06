
arch_root="/archive/engineering"
save_root="/scratch/rsiverd/mayhem"
calib_conf="./calib_types.txt"
barrier_list="./barriers.txt"

## Camera storage locations:
#declare -A cam_storage
#cam_storage["nres"]="/archive/engineering/bpl/nres"
#cam_storage["nres01"]="/archive/engineering/lsc/nres01"

bias_stack_args="-r410 --itable"
dark_stack_args="-r410 --itable -c5"

#declare -A objects_keyword
#objects_keyword["thar01"]="thar&thar&none"
#objects_keyword["thar12"]="none&thar&thar"
#objects_keyword["tung01"]="tung&tung&none"
#objects_keyword["tung12"]="none&tung&tung"

##--------------------------------------------------------------------------##
## Calibration version requirements:
min_biasvers="0.55"     # biases with version older than this are rebuilt
min_darkvers="0.60"     #  darks with version older than this are rebuilt
min_lampvers="0.35"     #  lamps with version older than this are rebuilt

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Date/time convenience functions:

## Convert between FDATE and UNIX seconds:
fdate2unix () { date +%s --date=$1; }
unix2fdate () { date +%Y%m%d --date="@$1"; }

day_change () {
   local old_fdate="$1"
   local day_delta="$2"
   unixtime_old=$(fdate2unix $old_fdate)
   unixtime_new=$(( unixtime_old + $day_delta * 86400 ))
   unix2fdate $unixtime_new
}

get_next_day () { day_change $1  1; }
get_prev_day () { day_change $1 -1; }

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Path and file selection routines:


##--------------------------------------------------------------------------##
## Get destination path for 'clean' files:
get_save_folder () {
   hdata=( `imhget TELESCOP DAY-OBS $1` )
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
      hits=( `ls $cal_folder/med_${ctype}_${cal_day}_*.fits 2>/dev/null` )
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

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
## Check for fdate in array:
date_is_listed () {
   local fdate="$1"
   shift
   for cmpdate in $*; do
      [ "$fdate" = "$cmpdate" ] && return 0
   done
   return 1
}

dates_are_listed () {
   local tdate1="$1"
   local tdate2="$2"
   shift 2
   others=( $* )
   for nite in $tdate1 $tdate2 ; do
      if !( date_is_listed $nite ${others[*]} ); then
         return 1
      fi
   done
   return 0
}

######################################################################
# CHANGELOG (config.sh):
#---------------------------------------------------------------------
#
#  2018-01-05:
#     -- Introduced min_biasvers=0.55, min_darkvers=0.60, min_lampvers=0.35.
#     -- Moved barrier_check_pass() into separate aux/01_barriers.sh file.
#     -- update_output_header now also sets MVERSION for revision tracking.
#     -- Added this change log.
#

