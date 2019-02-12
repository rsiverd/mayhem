import os, sys, time
import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
import wavelength_reference
reload(wavelength_reference)
wlr = wavelength_reference

nist_data = wlr.load_nist_argon_pd()
lope_data = wlr.load_lovis_pepe_thar()

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

### Save a copy to file:
#with open(save_txt, 'w') as f:
#    f.write("order,wl_min,wl_avg,wl_max,dl_min,dl_avg,dl_max\n")
#    for thing in results:
#        f.write(','.join([str(x) for x in thing]) + '\n')

## ----------------------------------------------------------------------- ##
## ----------------------------------------------------------------------- ##
## ----------------------------------------------------------------------- ##
## ----------------------------------------------------------------------- ##
## Line selection with some smarts:
def get_nist_lines(wl1, wl2, reltol=0.001,
        lamcol='lam_obs_nm', flxcol='rel_intensity'):
    twlen = nist_data[lamcol]
    which = (wl1 < twlen) & (twlen < wl2)
    neato = nist_data[which]
    nkept = which.sum()
    sys.stderr.write("NIST Argon: %d lines found.\n" % nkept)
    #for ww,ff in neato[[lamcol, flxcol]].values:
    #    sys.stderr.write("%10.5f --- %10.3f\n" % (ww, ff))

    # stick to brightest of lines found:
    thresh = neato[flxcol].max() * reltol
    #sys.stderr.write("thresh: %10.5f\n" % thresh)
    smart = (neato[flxcol] >= thresh)   # relative to high peak
    bright = neato[smart]
    sys.stderr.write("After peak-rel-cut, have %d lines.\n" % smart.sum())
    return bright[lamcol].values

def get_lope_lines(wl1, wl2, nmax=20):
    twlen = lope_data['lam_vac_nm']
    tflux = lope_data['flux']
    flcut = np.percentile(tflux, 75)
    which = (wl1 < twlen) & (twlen < wl2) & (tflux > flcut)
    neato = lope_data[which]
    sys.stderr.write("Lovis+Pepe: %d lines found.\n" % neato.size)
    top_few_idx = np.argsort(neato['flux'])[-nmax:]
    sys.stderr.write("Selecting top %d with highest flux ...\n" % nmax)
    return neato['lam_vac_nm'][top_few_idx]

def ocheck(oidx):
    #fig = plt.figure(1, figsize=(8,5))
    fig = plt.gcf()
    fig.clf()
    thar_flux = tharflat[oidx]
    thar_wlen = wavethar[oidx]
    wl1, wl2 = thar_wlen.min(), thar_wlen.max()
    sys.stderr.write("Order %2d: %.3f < lambda < %.3f\n" % (oidx, wl1, wl2))
    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    ax1.plot(thar_wlen, thar_flux, c='b', label='order %d'%oidx)

    ar_lines = get_nist_lines(wl1, wl2)
    for line in ar_lines[:-1]:
        ax1.axvline(line, ls=':', c='r')
    ax1.axvline(ar_lines[-1], ls=':', c='r', label='NIST Argon')

    thar_lines = get_lope_lines(wl1, wl2)
    if (thar_lines.size > 0):
        for line in thar_lines[:-1]:
            ax1.axvline(line, ls=':', c='g')
        ax1.axvline(thar_lines[-1], ls=':', c='g', label='Lovis & Pepe (2007)')

    ax1.legend()
    fig.tight_layout()
    plt.draw()


