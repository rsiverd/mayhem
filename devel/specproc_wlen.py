import os, sys, time
import numpy as np
import astropy.io.fits as pf

## A ~solved LSC spectrum as example:
specfile = 'lscnrs01-fa09-20181125-0031-e91.fits.fz'
wavespec = pf.getdata(specfile, extname='WAVESPEC')
wavethar = pf.getdata(specfile, extname='WAVETHAR')
tharflat = pf.getdata(specfile, extname='THARFLAT')

def mmm(a):
    return np.min(a), np.average(a), np.max(a)

## Typical wavelength step in each order:
for iord,ordwl in enumerate(wavethar):
    msgtxt  = "iord %3d:\n" % iord
    wstep   = np.diff(ordwl)
    msgtxt += "--> wlen min/avg/max: %10.5f, %10.5f, %10.5f\n" % mmm(ordwl)
    msgtxt += "--> dlam min/avg/max: %10.5f, %10.5f, %10.5f\n" % mmm(wstep)
    sys.stderr.write("%s\n" % msgtxt)

