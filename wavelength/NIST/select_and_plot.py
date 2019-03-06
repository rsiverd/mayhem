
import os, sys, time
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import scipy.signal as ssig

nist_top = pf.getdata('./NIST_spectrum.all.fits')
nist_all = pf.getdata('./NIST_spectrum.all.fits')

wlen_nm = nist_all['wlen'] / 10.0


def find_peak_wls(lamvec, flxvec, pkmin=200):
    #offset = np.array([-2, -1, 0, 1, 2])
    offset = np.arange(5) - 2
    kern = np.array([-1.0, 0.0, 1.0])
    high_idxs = (flxvec > pkmin).nonzero()[0]    # on a bright emission line
    high_wide = np.unique([offset+x for x in high_idxs])
    #sys.stderr.write("flxvec[high_wide]: %s\n" % str(flxvec[high_wide]))
    hits = ssig.argrelmax(flxvec[high_wide], order=3)[0]
    #sys.stderr.write("hits: %s\n" % str(hits))
    #peak_idx = high_vals.nonzero()[0][hits]     # indexes into original arrays
    peak_idx = high_wide[hits]                   # indexes into original arrays
    peak_ctr_wlen = []
    peak_tot_flux = []

    for ctr_idx in peak_idx:
        which = slice(ctr_idx - 1, ctr_idx + 2)
        pklam = lamvec[which]
        pkflx = flxvec[which]
        peak_ctr_wlen.append(np.sum(pklam * pkflx) / np.sum(pkflx))
        peak_tot_flux.append(np.sum(pkflx))
    peak_ctr_wlen = np.array(peak_ctr_wlen)
    peak_tot_flux = np.array(peak_tot_flux)
    peak_lam = lamvec[peak_idx]
    #sys.stderr.write("peaks: %s\n" % str(peak_lam))
    #sys.stderr.write("refined: %s\n" % str(peak_ctr_wlen))
    wlr_string = ', '.join(['%.6f'%x for x in peak_ctr_wlen])
    sys.stderr.write("use this:\n--> np.array([%s])\n" % wlr_string)
    return peak_lam

def get_wl_section(wl1_nm, wl2_nm):
    which = (wl1_nm * 10 < nist_all['wlen']) & (nist_all['wlen'] < wl2_nm * 10)
    subset = nist_all[which]
    tlam = subset['wlen'] / 10.0    # in nm
    tflx = subset['flux']           # in whatever
    fig = plt.gcf()
    fig.clf()
    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    ax1.plot(tlam, tflx)
    #ax1.set_yscale('symlog')
    #ax1.set_ylim(bottom=3)
    wl_peaks = find_peak_wls(tlam, tflx)
    [ax1.axvline(x, ls=':', c='r', lw=1.0) for x in wl_peaks]
    return tlam, tflx, wl_peaks

#tlam, tflx, wlpeaks = get_wl_section(745, 759)
