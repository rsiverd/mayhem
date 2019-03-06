
import numpy as np
import pandas as pd
fdata = pd.read_csv('./nist_spectrum.csv', comment='#')
ldata = fdata[::-1]  # ordered by increasing wavelength

#cutoff = 5000 #10000
bright = ldata[(ldata['intensity'] > 5000)]
booming = ldata[(ldata['intensity'] > 10000)]
#booming = (ldata['intensity'] > cutoff)
#peaks = ldata[booming]

## -----------------------------------------------------------------------
## FITSify output:
import astropy.io.fits as pf

twlen = ldata['wlen'].values
tflux = ldata['intensity'].values

## Make full FITS table:
columns = [pf.Column(name='wlen', format='D', unit='Angstrom', array=twlen),
            pf.Column(name='flux', format='D', unit='counts', array=tflux)]
tbhdu = pf.BinTableHDU.from_columns(pf.ColDefs(columns))
tbhdu.writeto('NIST_spectrum.all.fits', overwrite=True)

## Make 'bright' FITS table:
which = (tflux >= np.percentile(tflux, 95))
columns = [pf.Column(name='wlen', format='D', unit='Angstrom', array=twlen[which]),
            pf.Column(name='flux', format='D', unit='counts', array=tflux[which])]
tbhdu = pf.BinTableHDU.from_columns(pf.ColDefs(columns))
tbhdu.writeto('NIST_spectrum.top.fits', overwrite=True)

