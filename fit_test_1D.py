#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os, sys, time

def ycalc(xvals, slope):
    return slope*xvals

def fit_yx_ratio(xvals, yvals, weights=None):
    """Fit Y = a*X"""
    #if (weights == None):

    return np.sum(xvals * yvals) / np.sum(xvals * xvals)

# adopt sqrt(y) as scatter:
def make_fake_yvals(xvals, slope, sigma):
    y_true = slope * xvals
    y_vars = 0.0 * y_true + sigma**2                # baseline variance
    pwhich = (y_true > 0.0).nonzero()[0]            # which values are positive
    y_vars[pwhich] += y_true[pwhich]                # include sqrt(N) component
    y_errs = np.random.randn(y_true.size) * np.sqrt(y_vars)
    return y_true + y_errs

use_slope = 63.23
err_sigma = 10.0
npoints   =  50
xval_rng  =  75.0
xval_min  =  -5.0

clean_xvals = xval_min + xval_rng * np.random.random(npoints)
noisy_yvals = make_fake_yvals(clean_xvals, use_slope, err_sigma)
#noisy_yvals = make_fake_yvals(clean_xvals, use_slope)

# fit attempt:
#result = np.linalg.lstsq(clean_xvals[:, None], noisy_yvals)
result = np.linalg.lstsq(clean_xvals[:, None], noisy_yvals, rcond=-1)
params = result[0]
sys.stderr.write("fitted slope: %s\n" % str(params))

easyslope = fit_yx_ratio(clean_xvals, noisy_yvals)

fig = plt.figure(1)
fig.clf()
ax1 = fig.add_subplot(111)
ax1.grid(True)
ax1.scatter(clean_xvals, noisy_yvals, lw=0)

