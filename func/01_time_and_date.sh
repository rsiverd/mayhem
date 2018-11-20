#!/bin/bash
#
#    Time and date look-up, arithmetic, and modification routines.
#
# Rob Siverd
# Created:      2018-02-19
# Last updated: 2018-02-19
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

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

## List all dates between starting/ending FDATE:
list_dates_between () {
   local fdate1="$1"
   local fdate2="$2"
   next_unixtime=$(fdate2unix $fdate1)
   fdate_list=( `unix2fdate $next_unixtime` )
   while true; do
      next_unixtime=$(( next_unixtime + 86400 ))
      next_fdate=`unix2fdate $next_unixtime`
      if [ $next_fdate -le $fdate2 ]; then
         fdate_list+=( $next_fdate )
      else
         break
      fi
   done
   echo ${fdate_list[*]}
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
# CHANGELOG (01_time_and_date.sh):
#---------------------------------------------------------------------
#
#  2018-02-19:
#     -- First created 01_time_and_date.sh.
#
