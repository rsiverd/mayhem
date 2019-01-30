import os, sys, time
import numpy as np
import astropy.io.fits as pf

## A ~solved LSC spectrum as example:
specfile = 'lscnrs01-fa09-20181125-0031-e91.fits.fz'
save_txt = 'stats_specproc_wlen.csv'
if os.path.isfile(specfile):
    wavespec = pf.getdata(specfile, extname='WAVESPEC')
    wavethar = pf.getdata(specfile, extname='WAVETHAR')
    tharflat = pf.getdata(specfile, extname='THARFLAT')
else:
    sys.stderr.write("File not found: %s\n" % specfile)
    sys.exit(1)

## -----------------------------------------------------------------------
## Try to visualize image coordinates of central wavelengths on CCD:
ctrwl_data = np.genfromtxt('./lam_ctr_data.csv',
        dtype=None, delimiter=',', names=True)


