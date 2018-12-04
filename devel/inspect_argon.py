
import numpy as np
import astropy.io.fits as pf


## Argon lines:
data = np.genfromtxt('clean_argon.txt', dtype=None, names=True)

lam_obs  = data['lam_obs']
lam_ritz = data['lam_ritz']
rel_flux = data['rel_flux']

lam_diff = np.roll(lam_obs, -1) - lam_obs
lam_diff[-1] = 9999.9999

brightest = rel_flux.max()
which_bri = (rel_flux >= brightest / 5.0)


## Find a pair of bright lines close in wavelength (same order):
lam_sep = 2.0 # nm
#which1 = np.where(np.diff(lam_obs) < lam_sep)[0]    # idx of first in pair
#which1 = np.where(np.diff(lam_obs) < lam_sep)[0]    # idx of first in pair
lam_sep = 0.5
which1 = (lam_diff <= lam_sep).nonzero()[0]

bri_cut = 3000.0
lbright = (rel_flux[which1+0] >= bri_cut)
rbright = (rel_flux[which1+1] >= bri_cut)
#bripair = lbright & rbright
#success = data[lbright & rbright]
hit_idx = which1[lbright & rbright]
argon1 = data[hit_idx + 0]
argon2 = data[hit_idx + 1]

print("Expect the two adjacent lines to be:")
print("%s" % str(argon1))
print("%s" % str(argon2))

