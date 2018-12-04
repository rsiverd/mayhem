
fitsarith -i scratch/clean_med_nres01_tung01_20181125_07d_06b_00a.fits \
          -M scratch/clean_med_nres01_tung12_20181125_07d_06b_00a.fits \
          -D handy/spec_yprofile_med.fits.fz \
          -D handy/spec_yprofile_med.fits.fz \
          -o product_tung01_tung12.fits

fpeg -i product_tung01_tung12.fits -B 1.0 --sqrt \
     -o \!geomean_tung01_tung12.fits

