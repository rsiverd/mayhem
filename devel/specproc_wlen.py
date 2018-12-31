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

def mmm(a):
    return (np.min(a), np.average(a), np.max(a))

## Typical wavelength step in each order:
results = []
for iord,ordwl in enumerate(wavethar):
    msgtxt  = "iord %3d:\n" % iord
    wstep   = np.diff(ordwl)
    msgtxt += "--> wlen min/avg/max: %10.5f, %10.5f, %10.5f\n" % mmm(ordwl)
    msgtxt += "--> dlam min/avg/max: %10.5f, %10.5f, %10.5f\n" % mmm(wstep)
    sys.stderr.write("%s\n" % msgtxt)
    results.append([iord] + list(mmm(ordwl)) + list(mmm(wstep)))

## Save a copy to file:
with open(save_txt, 'w') as f:
    f.write("order,wl_min,wl_avg,wl_max,dl_min,dl_avg,dl_max\n")
    for thing in results:
        f.write(','.join([str(x) for x in thing]) + '\n')


sys.exit(0)

