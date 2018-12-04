import os, sys, time
import numpy as np
import pandas as pd
import astropy.io.fits as pf

def matched_diff(data, endpt):
    diffs = np.roll(data, -1) - data
    diffs[-1] = endpt
    return diffs

nist_lines = pd.read_csv('../wavelength/linelist.csv')

lam_obs  = nist_lines['wlen_meas'].values / 10.0     # now in nm
rel_flx = nist_lines['rel_intensity'].values
lam_dif = matched_diff(lam_obs, 9999.9999)
lam_dif[-1] = 9999.9999

bri2dim  = np.argsort(rel_flx)[::-1]

ntop = 100
lam_top  = np.sort(lam_obs[bri2dim[:ntop]])

lam_top_diff = matched_diff(lam_top, 9999.9999)

# look for spacings of roughly 1nm:
reasonable = (0.5 < lam_top_diff) & (lam_top_diff < 1.5)

# brightest things on the red side:
ntop = 50
reddish = (lam_obs > 650)
red_idx = reddish.nonzero()[0]
#red_dat = np.column_stack((lam_obs, rel_flx))[reddish]
#red_lam = lam_obs[red_idx]
#red_flx = rel_flx[red_idx]
get_few = np.argsort(rel_flx[red_idx])[::-1][:ntop]
top_idx = red_idx[get_few]
#top_few = np.argsort(red_dat[:, 1])[::-1][:ntop]



