#Imports
import matplotlib.pyplot as plt         #python / matlab
import pylab
import random                           #random generator package
import pyfits
import os
import numpy as np
import matplotlib.cm as cm
import bisect as bis
import matplotlib.image as mpimg
import random

#leastsquare package
from scipy.optimize.minpack import leastsq
from scipy import interpolate
from math import cos, sin, acos, asin, pi, atan, degrees, sqrt

#Astro Libraries
from astLib import astSED

minLambda=0.5886        #min wavelength
maxLambda=0.59          #max wavelength
deltaLambda =0.0001     #step interval
maxLambda+=deltaLambda

#Can plot orders from 146 to 73 (about 390 to 795 nm)
minOrder = 146
maxOrder =  73
deltaOrder = -1
maxOrder += deltaOrder
booLog = 6
pixelSize= 5.4


def main_errors(p, mainArgs):

    x, y, waveList, xSig, ySig = readCalibrationData(mainArgs[2])

    hdulist = pyfits.open('../c_noFlat_sky_0deg_460_median.fits')
    imWidth = hdulist[0].header['NAXIS1']
    imHeight = hdulist[0].header['NAXIS2']

    x = x - imWidth / 2
    y = y - imHeight / 2

    x_model, y_model, Lambda = main(p, mainArgs)

    x_best = x.copy()
    y_best = y.copy()
    for k in range(0, len(waveList)):
        ix, = np.where(waveList[k] == Lambda)
        if (ix.size == 0):
            x_best[k] = 0
            y_best[k] = 0
        else:
            best = ix[np.argmin(np.abs(y_model[ix] - y[k]))]
            x_best[k] = x_model[best]
            y_best[k] = y_model[best]

    return np.hstack([(x_best - x) / xSig, (y_best - y) / ySig]), waveList

