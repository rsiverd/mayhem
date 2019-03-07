import os, sys, time
import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
import theil_sen as ts
import fov_rotation
reload(fov_rotation)
r3d = fov_rotation.Rotate3D()

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

#ctr_coo = []
#centerx = 2048.5
#tim_idx = np.arange(wavekind.shape[0])
#for oi,wlc in enumerate(use_lamcen):
#    which = (owlmin <= wlc) & (wlc <= owlmax)
#    nhits = np.sum(which)
#    sys.stderr.write("nhits = %d for oi=%2d\n" % (nhits, oi))
#    if nhits == 0:
#        continue
#    #possibles = [(ti, np.interp(wlc, ordwl, xpixel)) for ti,ordwl in \
#    #        zip(tim_idx[which], wavekind[which])]
#    xpixels = [np.interp(wlc, ordwl, xpixel) for ordwl in wavekind[which]]
#    #tmp_ti, tmp_xx = zip(*possibles)
#    idx_mid = np.argmin(np.abs(np.array(xpixels) - centerx))
#    ctr_coo.append((oi, xpixels[idx_mid], use_ordsep[oi]))
    
def central_wavelength_xy(lam_cen_nm, ord_sep_pix, wave_sol_nm,
        midpct=(25, 75)):

    midpoints = []
    wls1, wls2 = np.percentile(wave_sol_nm, midpct, axis=1)
    for oi,wlc in enumerate(lam_cen_nm):
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

#ctr_coo = np.array(ctr_coo)

wl_midpoints_dumb = central_wavelength_xy(1e3 * ctrwl_data['lam_dumb'],
        ctrwl_data['sep_dumb'], wavethar)

wl_midpoints_iter = central_wavelength_xy(1e3 * ctrwl_data['lam_iter'],
        ctrwl_data['sep_iter'], wavethar)

_, ccdx_dumb, ccdy_dumb = wl_midpoints_dumb.T
_, ccdx_iter, ccdy_iter = wl_midpoints_iter.T

## -----------------------------------------------------------------------
## Save midpoint coordinates to file:
save_file_iter = "lamcen_iter_coords.dat"
save_file_dumb = "lamcen_dumb_coords.dat"
with open(save_file_iter, 'w') as f:
    f.write("lamcen_nm  ccd_xpix  ccd_ypix\n")
    for row in wl_midpoints_iter:
        f.write("%9.4f %9.3f %9.3f\n" % tuple(row))
    pass
with open(save_file_dumb, 'w') as f:
    f.write("lamcen_nm  ccd_xpix  ccd_ypix\n")
    for row in wl_midpoints_dumb:
        f.write("%9.4f %9.3f %9.3f\n" % tuple(row))
    pass

## -----------------------------------------------------------------------
## Linear fit to central wavelengths:
icept, slope = ts.linefit(ccdx_iter, ccdy_iter)
fitted_x = np.array([0, 4000])
fitted_y = icept + slope * fitted_x

rotation_deg = 90.0 - np.degrees(np.arctan(slope))
sys.stderr.write("Tilt of central wavelength ridge: %.3f deg\n" % rotation_deg)

## Attempt to return to native spectrograph coordinates:
camera_xyz = np.vstack((ccdx_iter, ccdy_iter, np.zeros_like(ccdx_iter)))
native_xyz = r3d.zrot(np.radians(rotation_deg), camera_xyz) # nice!
spcx, spcy, _ = native_xyz

## -----------------------------------------------------------------------
## Transform CCD -> spectrograph coordinates:
def ccd2spec_xy(ccdx, ccdy, rot_deg):
    if ccdx.shape != ccdy.shape:
        raise ValueError("CCD coordinate arrays have mismatched shape\n:"
                + "%s != %s\n" % (str(ccdx.shape), str(ccdy.shape)))
    if len(ccdx.shape) != 1:
        raise ValueError("Expected 1-D input, have shape %s" % str(ccdx.shape))
    # MEMORY INTENSIVE!!
    #old_dim = ccdx.shape
    #ccd_xyz = np.vstack((ccdx.flatten(), 
    #                     ccdy.flatten(),
    #                     np.zeros(ccdx.size)))
    ccd_xyz = np.vstack((ccdx, ccdy, np.zeros(ccdx.size)))
    sx, sy, _ = r3d.zrot(np.radians(rot_deg), ccd_xyz)
    return sx.A1, sy.A1
    #return np.squeeze(np.asarray(sx)), np.squeeze(np.asarray(sy))
    #return np.array(sx), np.array(sy)
    #return sx.reshape(old_dim), sy.reshape(old_dim)

ny, nx = 4096, 4096
x_list = (0.5 + np.arange(nx)) / nx - 0.5            # relative (centered)
y_list = (0.5 + np.arange(ny)) / ny - 0.5            # relative (centered)
#xx, yy = np.meshgrid(x_list, y_list)                 # relative (centered)
xx, yy = np.meshgrid(nx*x_list, ny*y_list)           # absolute (centered)

spec_coo = []
for tx,ty in zip(xx, yy):
    spec_coo.append(ccd2spec_xy(tx, ty, rotation_deg))
sxx, syy = (np.vstack(coo) for coo in zip(*spec_coo))
del spec_coo


## -----------------------------------------------------------------------
## -----------------------------------------------------------------------
## -----------------------------------------------------------------------
## Plot something sensible ...

fig = plt.figure(1, figsize=(8,8))
fig.clf()
ax1 = fig.add_subplot(111, aspect='equal')
ax1.grid()

ax1.plot(fitted_x, fitted_y, lw=1, ls=':', c='b')
ax1.scatter(ccdx_dumb, ccdy_dumb, s=15, lw=0, c='r', label='dumb')
ax1.scatter(ccdx_iter, ccdy_iter, s=15, lw=0, c='g', label='iter')
ax1.set_xlim(0, 4100)
ax1.set_ylim(0, 4100)
ax1.set_xlabel("X Pixel")
ax1.set_ylabel("Y Pixel")
ax1.legend(loc='best')

fig.tight_layout()
plt.draw()

