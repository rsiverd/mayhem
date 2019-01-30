import os, sys, time
import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt

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
## Load central wavelength data:
ctrwl_data = np.genfromtxt('./lam_ctr_data.csv',
        dtype=None, delimiter=',', names=True)

use_lamcen = ctrwl_data['lam_dumb'] * 1e3   # in nm
use_ordsep = ctrwl_data['sep_dumb']

use_lamcen = ctrwl_data['lam_iter'] * 1e3   # in nm
use_ordsep = ctrwl_data['sep_iter']

## -----------------------------------------------------------------------
## Compute X- and Y-positions of central wavelength for each order:
xpixel = 1.0 + np.arange(wavespec.shape[1])
wavekind = wavethar
#wavekind = wavespec
owlmin = np.min(wavekind, axis=1)
owlmax = np.max(wavekind, axis=1)

#def central_wavelength_xy(lam_cen_nm, ord_sep_pix, wave_sol_nm):
#    ctr_coo = []
#    for oi,wlc in enumerate(use_lamcen):
#        if (wlc < wavekind.min()) or (wlc > wavekind.max()):
#            continue
#        for tim_oid,ordwl in enumerate(wavekind):
#            if (ordwl.min() < wlc) and (wlc < ordwl.max()):
#                sys.stderr.write("tim_oid: %d\n" % tim_oid) 
#                wlc_x = np.interp(wlc, ordwl, xpixel)
#                wlc_y = use_ordsep[tim_oid+1]
#                ctr_coo.append((wlc, wlc_x, wlc_y))

ctr_coo = []
centerx = 2048.5
tim_idx = np.arange(wavekind.shape[0])
for oi,wlc in enumerate(use_lamcen):
    which = (owlmin <= wlc) & (wlc <= owlmax)
    nhits = np.sum(which)
    sys.stderr.write("nhits = %d for oi=%2d\n" % (nhits, oi))
    if nhits == 0:
        continue
    #possibles = [(ti, np.interp(wlc, ordwl, xpixel)) for ti,ordwl in \
    #        zip(tim_idx[which], wavekind[which])]
    xpixels = [np.interp(wlc, ordwl, xpixel) for ordwl in wavekind[which]]
    #tmp_ti, tmp_xx = zip(*possibles)
    idx_mid = np.argmin(np.abs(np.array(xpixels) - centerx))
    ctr_coo.append((oi, xpixels[idx_mid], use_ordsep[oi]))
    
def central_wavelength_xy(lam_cen_nm, ord_sep_pix, wave_sol_nm,
        midpct=(25, 75), xleft=1000, xright=3000):

    midpoints = []
    wls1, wls2 = np.percentile(wave_sol_nm, midpct, axis=1)
    for oi,wlc in enumerate(use_lamcen):
        which = (wls1 <= wlc) & (wlc <= wls2)
        nhits = np.sum(which)
        if nhits == 0:
            continue
        if nhits >= 2:
            sys.stderr.write("LOLWUT?!?!?!\n")
            return None
        ordwl = wave_sol_nm[which][0]
        thisx = np.interp(wlc, ordwl, xpixel)
        thisy = ord_sep_pix[oi]
        midpoints.append((wlc, thisx, thisy))
    return np.array(midpoints)

    #break
    #if nhits > 1:
    #    print(possibles)
    #    break
    #ctr_y = use_ordsep[oi]
    #break

sys.exit(0)
ctr_coo = np.array(ctr_coo)
_, cwl_xpix, cwl_ypix = ctr_coo.T

## -----------------------------------------------------------------------
## Plot something sensible ...

fig = plt.figure(1, figsize=(8,8))
fig.clf()
ax1 = fig.add_subplot(111, aspect='equal')
ax1.grid()

ax1.scatter(cwl_xpix, cwl_ypix, s=35, lw=0)
ax1.set_xlim(0, 4100)
ax1.set_ylim(0, 4100)
ax1.set_xlabel("X Pixel")
ax1.set_ylabel("Y Pixel")

fig.tight_layout()
plt.draw()

