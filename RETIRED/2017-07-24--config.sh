arch_root="/archive/engineering"
save_root="/scratch/rsiverd/mayhem"

## Camera storage locations:
#declare -A cam_storage
#cam_storage["nres"]="/archive/engineering/bpl/nres"
#cam_storage["nres01"]="/archive/engineering/lsc/nres01"


# -----------------------------------------------------------------------
# Image picking routines:


# Get destination path for 'clean' files:
get_save_folder () {
   hdata=( `imhget CAMID DRTAG DAY-OBS $1` )
   local camid=${hdata[0]}
   local drtag=${hdata[1]}
   local dayobs=${hdata[2]}
   echo $save_root/$camid/$dayobs
}

## Get specific bias frame path:
get_bias_path () {
   local image="$1"
   local camid="$2"
   local drtag="$3"
   local dayobs=$(imhget DAY-OBS $image)
   echo $save_root/$camid/$dayobs/med_bias_${dayobs}_${drtag}.fits
}

## Pick the best available bias/dark calibration image:
pick_best_bdcal () {
   local image="$1"
   local camid="$2"
   local ctype="$3"
   local dayobs=$(imhget DAY-OBS $image)  || exit $?

   # validate inputs:
   if [ -z "$camid" ]; then
      echo "pick_best_bdcal error: camid cannot be blank!" >&2
      return 1
   fi
   case $ctype in
      bias|dark) ;;
      *)
         echo "pick_best_bdcal error: unsupported calib type: '$ctype'" >&2
         return 1
         ;;
   esac

   # check for data directory:
   local cal_folder="$save_root/$camid/$dayobs"
   if [ ! -d $cal_folder ]; then
      echo "pick_best_bdcal error: folder not found:" >&2
      echo "--> $cal_folder" >&2
      return 1
   fi
   hits=(`ls $save_root/$camid/$dayobs/med_bias_${dayobs}_*.fits 2>/dev/null`)
   nfound=${#hits[*]}
   #echo "nfound: $nfound" >&2
   if [ $nfound -eq 0 ]; then
      echo "pick_best_bdcal error: no potential images found!" >&2
      return 1
   fi
   echo "${hits[-1]}"
}

