#!/usr/bin/env python
# vim: set fileencoding=utf-8 ts=4 sts=4 sw=4 et tw=80 :
#
# Separate module for optical material types.
#
# Rob Siverd
# Created:       2019-03-30
# Last modified: 2019-03-30
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.1.0"

### Python version-agnostic module reloading:
#try:
#    reload                              # Python 2.7
#except NameError:
#    try:
#        from importlib import reload    # Python 3.4+
#    except ImportError:
#        from imp import reload          # Python 3.0 - 3.3

## Modules:
import os
import sys
import time
import numpy as np

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## Calculating index of refraction for several glass types:
class Glass(object):

    def __init__(self, glasstype):
        self._bcoeffs = {
            'SiO2':np.array([0.67071081e0, 0.433322857e0, 0.877379057e0]),
            'LLF1':np.array([1.21640125e0, 1.33664540e-1, 8.83399468e-1]),
            'PBM2':np.array([1.39446503e0, 1.59230985e-1, 2.45470216e-1]),
             'LF5':np.array([1.28035628e0, 1.6350597e-1,  8.93930112e-1]), }
        self._ccoeffs = {
            'SiO2':np.array([0.00449192312e0, 0.0132812976e0, 95.8899878e0]),
            'LLF1':np.array([8.57807248e-3,   4.20143003e-2,   1.07593060e+2]),
            'PBM2':np.array([1.10571872e-2,   5.07194882e-2,   3.14440142e1]),
             'LF5':np.array([9.29854416e-3,   4.49135769e-2,   1.10493685e2]), }
        if self._unknown_glass(glasstype):
            raise
        self._gtype = glasstype
        self._coeffs = zip(self._bcoeffs[self._gtype],
                           self._ccoeffs[self._gtype])
        return

    def _unknown_glass(self, glasstype):
        if not glasstype in self._bcoeffs.keys():
            sys.stderr.write("Unknown glass type: %s\n" % glasstype)
            return True
        else:
            return False

    # Squared index of refraction for specified wavelengths:
    def refraction_index_squared(self, wlen_um):
        lam_um_sq = wlen_um**2
        n_squared = np.ones_like(wlen_um, dtype='float')
        for bb,cc in self._coeffs:
            n_squared += (lam_um_sq * bb) / (lam_um_sq - cc)
        return n_squared

    def refraction_index(self, wlen_um):
        return np.sqrt(self.refraction_index_squared(wlen_um))

    def glass_dn_dlambda_easy(self, wlen_um, glasstype, epsfrac=1e-5):
        #if self._unknown_glass(glasstype):
        #    raise

        wlen_lower = wlen_um * (1.0 - epsfrac)
        wlen_upper = wlen_um * (1.0 + epsfrac)

        nn_lower = self.refraction_index(wlen_lower) #, glasstype)
        nn_upper = self.refraction_index(wlen_upper) #, glasstype)

        return (nn_upper - nn_lower) / (wlen_upper - wlen_lower)




######################################################################
# CHANGELOG (glass.py):
#---------------------------------------------------------------------
#
#  2019-03-30:
#     -- Increased __version__ to 0.0.1.
#     -- First created glass.py.
#
