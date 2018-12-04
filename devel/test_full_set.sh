#!/bin/bash

tdir="scratch"
mkdir -p $tdir

tung_list=( `ls images/med*tung*fz` )
thar_list=( `ls images/med*thar*fz` )

baff_file="vtx-baffle-lsc.fits.fz"
extr_args="-Y handy/spec_yprofile_med.fits.fz -B $baff_file"
extr_args="-Y handy/spec_yprofile_med.fits.fz -B vtx-baffle-lsc.fits.fz"

for image in ${tung_list[*]}; do
   ibase="${image##*/}"
   ifits="${ibase%.fz}"
   trsave="$tdir/trace1_${ifits}"
   bgsave="$tdir/bkg_${ifits}"
   cln_save="$tdir/clean_${ifits}"
   cmde "./mayhem-find-spec-trace.py $extr_args $image -o $trsave"
   status=$?
   [ $status -ne 0 ] && break

   # Trace-assisted background subtraction:
   bkg_args="-T $trsave $extr_args"
   cmde "./mayhem-estimate-background.py $image $bkg_args -o $bgsave"
   status=$?
   [ $status -ne 0 ] && break

   # Remove background:
   cmde "fitsarith -qHi $image -S $bgsave -o $cln_save" || break
done


