#!/bin/bash
#
#    Barrier-checking routines for master calibration builder scripts.
#
# Rob Siverd
# Created:      2018-01-05
# Last updated: 2018-01-05
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Load barriers:
barrier_check_pass () {
   local camid=$1
   local otype=$2
   local nlist=$3
   matches=0
   if [ ! -f $barrier_list ]; then
      echo "No barrier list found!"
      return 0
   fi
   barriers=( `awk -v camid=$camid -v otype=$otype '{
      if (($1 == camid) && ($2 == otype)) print $3"_"$4
      }' $barrier_list` )

   #echo "barriers: ${barriers[*]}"
   for dpair in ${barriers[*]}; do
      nites=( `echo $dpair | tr '_' ' '` )
      #echo "nites: ${nites[*]}"
      #grep -e ${nites[0]} -e ${nites[1]} $nlist
      hits=$(grep -e ${nites[0]} -e ${nites[1]} $nlist | wc -l)
      #echo "hits: $hits"
      if [ $hits -eq 2 ]; then
         #echo "BARRIER CROSSED! ${nites[*]}" >&2
         return 1
      fi
   done
   return 0
}

######################################################################
# CHANGELOG (01_barriers.sh):
#---------------------------------------------------------------------
#
#  2018-01-05:
#     -- First created 01_barriers.sh.
#
