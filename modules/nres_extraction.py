#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
#    Module providing useful routines for NRES spectrum extraction.
#
# Rob Siverd
# Created:       2017-08-14
# Last modified: 2018-12-24
#--------------------------------------------------------------------------
#**************************************************************************
#--------------------------------------------------------------------------

## Current version:
__version__ = "0.4.2"

## Modules:
#import argparse
#import mimetypes
#import linecache
#import getopt
#import shutil
import signal
#import glob
import os
import sys
import time
#import ephem
import numpy as np
#import datetime as dt
#import scipy.signal as ssig
import scipy.ndimage as ndi
#import scipy.optimize as opti
#import scipy.interpolate as stp
#import scipy.spatial.distance as ssd
import statsmodels.api as sm
import theil_sen as ts

## FITS I/O:
try:
    import astropy.io.fits as pf
except ImportError:
    try:
       import pyfits as pf
    except ImportError:
        sys.stderr.write("\nError!  No FITS I/O module found!\n"
               "Install either astropy.io.fits or pyfits and try again!\n\n")
        sys.exit(1)

##--------------------------------------------------------------------------##
## Catch interruption cleanly:
#def signal_handler(signum, frame):
#    sys.stderr.write("\nInterrupted!\n\n")
#    sys.exit(1)
#
#signal.signal(signal.SIGINT, signal_handler)

##--------------------------------------------------------------------------##
## Save FITS image with clobber (astropy / pyfits):
#def qsave(iname, idata, header=None, padkeys=1000, **kwargs):
#    sys.stderr.write("Writing to '%s' ... " % iname)
#    if header:
#        while (len(header) < padkeys):
#            header.append() # pad header
#    if os.path.isfile(iname):
#        os.remove(iname)
#    pf.writeto(iname, idata, header=header, **kwargs)
#    sys.stderr.write("done.\n")

##--------------------------------------------------------------------------##
#def argnear(vec, val):
#    return (np.abs(vec - val)).argmin()

## Robust location/scale estimate using median/MAD:
#def calc_ls_med_MAD(a, axis=None):
#    """Return median and median absolute deviation of *a* (scaled to normal)."""
#    med_val = np.median(a, axis=axis)
#    sig_hat = (1.482602218 * np.median(np.abs(a - med_val), axis=axis))
#    return (med_val, sig_hat)

## Robust location/scale estimate using median/IQR:
#def calc_ls_med_IQR(a, axis=None):
#    """Return median and inter-quartile range of *a* (scaled to normal)."""
#    pctiles = np.percentile(a, [25, 50, 75], axis=axis)
#    med_val = pctiles[1]
#    sig_hat = (0.741301109 * (pctiles[2] - pctiles[0]))
#    return (med_val, sig_hat)

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##


##--------------------------------------------------------------------------##
## Flat-relative spectrum solver:
#def flat_rel_solver(schunk, fchunk, rcond=-1):
#    """First axis is Y (cross-dispersion), second is X (CCD column)."""
#    # Cross-dispersion size and number of CCD columns in block:
#    xdpixels, ncolumns = schunk.shape
#    rel_spec = np.zeros(ncolumns, dtype='float32')
#    for i in range(xdpixels):
#       fcol = fchunk[i, :].reshape(-1, 1)
#       srow = schunk[i, :]
#       rel_spec[i] = np.linalg.lstsq(fcol, srow, rcond=rcond)[0]
#    #for i,(fcol,scol) in enumerate(zip(fchunk, schunk)):
#    #    rel_spec[i] = np.linalg.lstsq(fcol.reshape(-1, 1), scol, rcond=None)[0]
#    return rel_spec

##--------------------------------------------------------------------------##
## Flat-relative extraction class:
class FlatRelativeExtraction(object):

    def __init__(self):
        return

    # Fit for Y/X using direct summation (single slope):
    @staticmethod
    def fit_yxratio_direct_sums(xvals, yvals):
        """Fit Y = a*X. Returns `a'."""
        return np.sum(xvals * yvals) / np.sum(xvals * xvals)


    # Bulk fitting of yxratio with direct method:
    @staticmethod
    def flat_rel_solver(lblob, sblob, weights=None):
        """
        Inputs:
            lblob   <-- lampflat pixels
            sblob   <-- spectrum pixels
            weights <-- per-pixel weights/errors
                        FIXME: need to double-check what to use here
        """
    
        # set uniform weights if none provided:
        if (weights == None):
            wvals = np.ones_like(lblob)
        else:
            wvals = weights
        if (lblob.shape != sblob.shape) or (lblob.shape != wvals.shape):
            sys.stderr.write("Mismatched shapes:\n")
            sys.stderr.write("lblob: %s\n" % str(lblob.shape))
            sys.stderr.write("sblob: %s\n" % str(sblob.shape))
            sys.stderr.write("wvals: %s\n" % str(wvals.shape))
            raise
        if (len(lblob.shape) != 2):
            sys.stderr.write("Unexpected shape: %s\n" % str(lblob.shape))
            sys.stderr.write("Expected 2-D input ....\n")
            raise
        narrower = np.argmin(lblob.shape)   # skinny axis should be summed over
        #sys.stderr.write("narrower: %d\n" % narrower)
        numer = np.sum(wvals * lblob * sblob, axis=narrower)
        denom = np.sum(wvals * lblob * lblob, axis=narrower)
        return numer / denom

##--------------------------------------------------------------------------##
## Fit for Y/X using numpy linalg.lstsq:
def fit_yxratio_numpy_lstsq(xvals, yvals, rcond=-1, full=False):
    result = np.linalg.lstsq(xvals[:, None], yvals, rcond=rcond)
    if full:
        return result       # everything we know
    else:
        return result[0][0] # parameters only

##--------------------------------------------------------------------------##
## Fit for Y/X using direct summation (single slope):
def fit_yxratio_direct_sums(xvals, yvals):
    """Fit Y = a*X. Returns `a'."""
    return np.sum(xvals * yvals) / np.sum(xvals * xvals)

## Bulk fitting of yxratio with direct method:
def flat_rel_solver(xblob, yblob, weights=None):
    """
    Inputs:
        xblob <-- lampflat pixels
        yblob <-- spectrum pixels
    """

    # set uniform weights if none provided:
    if (weights == None):
        wvals = np.ones_like(xblob)
    else:
        wvals = weights
    if (xblob.shape != yblob.shape) or (xblob.shape != wvals.shape):
        sys.stderr.write("Mismatched shapes:\n")
        sys.stderr.write("xblob: %s\n" % str(xblob.shape))
        sys.stderr.write("yblob: %s\n" % str(yblob.shape))
        sys.stderr.write("wvals: %s\n" % str(wvals.shape))
        raise
    if (len(xblob.shape) != 2):
        sys.stderr.write("Unexpected shape: %s\n" % str(xblob.shape))
        sys.stderr.write("Expected 2-D input ....\n")
        raise
    narrower = np.argmin(xblob.shape)   # skinny axis should be summed over
    #sys.stderr.write("narrower: %d\n" % narrower)
    numer = np.sum(wvals * xblob * yblob, axis=narrower)
    denom = np.sum(wvals * xblob * xblob, axis=narrower)
    return numer / denom

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##


##--------------------------------------------------------------------------##
## Overwrite image border:
def borderize(image, width, value):
    new_img = image.copy()
    new_img[:, 0:width] = value
    new_img[:, -width:] = value
    new_img[0:width, :] = value
    new_img[-width:, :] = value
    return new_img

##--------------------------------------------------------------------------##
## Basic extraction routine:
#mask_image = mask_pix
#y_offset = 10
def extract_flux(sci_image, mask_image, y_offset):
    ny, nx = sci_image.shape
    rows = np.arange(ny)
    cols = np.arange(nx)
    xx, yy = np.meshgrid(cols, rows)   # absolute

    raw_spectrum = []
    sys.stderr.write("Fudging mask ... ")
    use_mask = borderize(np.roll(mask_image, y_offset, axis=0), 5, 0)
    sys.stderr.write("finding objects ... ")
    img_labs, numL = ndi.label(use_mask.byteswap().newbyteorder())
    obj_areas = ndi.find_objects(img_labs)
    obj_npix = np.bincount(img_labs.ravel())
    sys.stderr.write("done.\n")
    sys.stderr.write("numL: %d\n" % numL)
    #qsave(img_labs, 'yomama.%d.fits' % y_offset)
    sys.stderr.write("\n")
    for i, aslice in enumerate(obj_areas):
        lvalue = i + 1
        #sys.stderr.write("lvalue: %d\n" % lvalue)
        img_pix = sci_image[aslice]
        lab_pix = img_labs[aslice]
        xcoords = xx[aslice]
        ycoords = yy[aslice]
        keepers = (lab_pix == lvalue)
        im_mask = np.zeros_like(lab_pix, dtype='int')
        #im_mask = np.zeros_like(lab_pix, dtype='float')
        im_mask[ keepers] = 1
        #im_mask[~keepers] = np.nan

        # order-specific pixels (still 2D):
        #use_pix = img_pix * im_mask
        all_pix = img_pix * im_mask
        #use_pix = np.ma.masked_array(all_pix, np.isnan(nmasked))
        use_pix = np.ma.masked_array(all_pix, ~keepers)
        x_pixel = xcoords[0]
        y_ridge = np.sum(use_pix * ycoords, axis=0) / np.sum(use_pix, axis=0)
        avg_flx = np.mean(use_pix, axis=0)
        sum_flx = np.sum(use_pix, axis=0)
        #ord_flx = np.average(use_pix, axis=0)
        #ord_flx = np.sum(use_pix, axis=0)
        #pix_flx = np.vstack((x_pixel, y_ridge, avg_flx, sum_flx))
        pix_flx = np.vstack((x_pixel, y_ridge, sum_flx))
        raw_spectrum.append(pix_flx)

    return raw_spectrum

##--------------------------------------------------------------------------##
## How to save spectra:
def fitsify_spectrum(data, filename):
    sys.stderr.write("Writing table '%s' ... " % filename)
    hdu_list = pf.HDUList()
    #for i, (xpix, ypix, aflux, sflux) in enumerate(data):
    for i, (xpix, ypix, sflux) in enumerate(data):
        #sys.stderr.write("i: %2d\n" % i)
        col_list = []
        col_list.append(pf.Column(name='X_PIXEL', format='J', array=xpix))
        col_list.append(pf.Column(name='Y_PIXEL', format='D', array=ypix))
        col_list.append(pf.Column(name='COUNTS', format='D', array=sflux))
        #col_list.append(pf.Column(name='AVG_COUNTS', format='D', array=aflux))
        #col_list.append(pf.Column(name='SUM_COUNTS', format='D', array=sflux))
        col_defs = pf.ColDefs(col_list)
        tab_hdu = pf.BinTableHDU.from_columns(col_defs)
        tab_hdu.name = 'ORDER_%d' % i
        hdu_list.append(tab_hdu)
    hdu_list.writeto(filename, clobber=True)
    sys.stderr.write("done.\n")



##--------------------------------------------------------------------------##
##                      Loading and storing of trace data                   ##
##--------------------------------------------------------------------------##

class TraceData(object):

    def __init__(self, trace_list, metadata):
        self._trace_list = trace_list
        self._metadata = metadata
        self._imshape = self._get_imshape(self._metadata)
        return

    # Look for image dimensions in metadata:
    @staticmethod
    def _get_imshape(mdata):
        if ('SRC_XPIX' in mdata.keys()) and ('SRC_YPIX' in mdata.keys()):
            return (mdata['SRC_YPIX'], mdata['SRC_XPIX'])
        else:
            return None

    # Return raw trace parameters:
    def get_trace_list(self):
        """Return raw trace fit parameters."""
        return self._trace_list

    # Return full primary data header:
    def get_metadata(self):
        return self._metadata

    # Return traces as pixel masks (requires appropriate metadata):
    def get_trace_masks(self, vlevel=0):
        """Returns traces as pixel masks."""
        if not self._imshape:
            sys.stderr.write("Image dimensions not available!\n")
            #return None
            raise
        return self._mask_from_traces(self._imshape, self._trace_list, vlevel)

    # Build pixel masks corresponding to listed traces:
    @staticmethod
    def _mask_from_traces(imshape, trace_list, vlevel=0):
        mask_image = np.zeros(imshape, dtype='bool')
        trace_coords = []
        n_traces = len(trace_list)
        for i,trace_fit in enumerate(trace_list, 1):
            if (vlevel >= 1):
                sys.stderr.write("\rAdding trace %d of %d ... " % (i, n_traces))
            xlist = np.arange(trace_fit['xmin'],
                            trace_fit['xmax']).astype('uint16')
            ordfit_ycoord = ridge_eval(trace_fit['params'], xlist)
            ylower = np.int_(np.floor(ordfit_ycoord))
            yc_list, xc_list = [], []
            apron_pix = trace_fit['apron']
            for offset in range(-apron_pix + 1, apron_pix + 1):
                xc_list.append(xlist)
                yc_list.append(ylower + offset)
                pass
            trace_coords.append((np.vstack(yc_list), np.vstack(xc_list)))
        return trace_coords

##--------------------------------------------------------------------------##
##--------------------------------------------------------------------------##

## In this (initial) format, each trace will be given its own HDU. That HDU has
## a single 'params' column with the polynomial fit coefficients. Each HDU also
## has a few keywords providing useful metadata.

## Metadata keyword/comment mapping ():
_trace_hkey_spec = [
        ( 'xmin',  'XMIN', '[pixel] trace lower X limit (left side)'),
        ( 'xmax',  'XMAX', '[pixel] trace upper X limit (right side)'),
        ('apron', 'APRON', '[pixel] apron size used for tracing'),
        ]

_metadata_order = ['EXTRVERS', 'TR_IMAGE', 'SRC_XPIX', 'SRC_YPIX',
            'TRMETHOD', 'HALFBPIX', 'BAFFMASK', 'YPROFILE', ]

class TraceIO(object):

    def __init__(self):
        self._divcmt = pf.Card("COMMENT", 65*'-')
        return

    def _header_from_dict(self, fit_data):
        c_list = [self._divcmt, self._divcmt]
        for dkey,fkey,cmnt in _trace_hkey_spec:
            c_list.append(pf.Card(fkey, fit_data[dkey], comment=cmnt))
        return pf.Header(c_list)

    def _trace_to_HDU(self, fit_data):
        header  = self._header_from_dict(fit_data)
        pcolumn = pf.Column(name="params", format='D', unit=None, 
                        array=fit_data['params'])
        return pf.BinTableHDU.from_columns([pcolumn,], header=header)

    def _trace_from_HDU(self, trace_HDU):
        fit_data = {'params':trace_HDU.data['params']}
        for dkey,fkey,cmnt in _trace_hkey_spec:
            fit_data[dkey] = trace_HDU.header[fkey]
        return fit_data

    # Save a list of traces to a FITS table:
    def store_traces(self, filename, traces_list, hdata=None):
        tables = []
        prihdr = pf.Header()
        prihdr['TRIOVERS'] = (__version__, 'TraceIO code version')
        if hdata:
            # Standard keys go in first:
            for kk in _metadata_order:
                if kk in hdata.keys():
                    prihdr[kk] = tuple(hdata.pop(kk))
            prihdr.append(self._divcmt)

            # Dump in anything else:
            if len(hdata):
                prihdr.update({k:tuple(v) for k,v in hdata.items()})
                prihdr.append(self._divcmt)
        prihdu = pf.PrimaryHDU(header=prihdr)

        tables.append(prihdu)
        for trace in traces_list:
            tables.append(self._trace_to_HDU(trace))

        hdu_list = pf.HDUList(tables)
        hdu_list.writeto(filename, overwrite=True)
        return

    # Load traces from the specified file:
    def load_traces(self, filename):
        traces_list = []
        with pf.open(filename) as hdu_list:
            all_pri_keys = hdu_list[0].header
            use_pri_keys = all_pri_keys.copy(strip=True)
            for hdu in hdu_list[1:]:
                traces_list.append(self._trace_from_HDU(hdu))
        return TraceData(traces_list, use_pri_keys)

##--------------------------------------------------------------------------##
##                      overplotting of traces onto image                   ##
##--------------------------------------------------------------------------##

## Trimmer for 'drawing' ridges:
#def trim_to_image_dims(xcoo, ycoo, imshape):
def trim_to_image_dims(xcoo, ycoo, image):
    ny, nx = image.shape
    useful = (0 <= xcoo) & (xcoo < nx) & (0 <= ycoo) & (ycoo < ny)
    return (ycoo[useful], xcoo[useful])

def overplot_traces(idata, trace_list, vlevel=0):
    n_traces  = len(trace_list)
    tmp_image = np.copy(idata)
    for i,trace_fit in enumerate(trace_list, 1):
        if (vlevel >= 0):
            sys.stderr.write("\rPainting trace %d of %d ... " % (i, n_traces))
        #ordfit_params = nrex.fit_polynomial(xlist, ylist, fit_degree)['params']
        xlist = np.arange(trace_fit['xmin'], trace_fit['xmax']).astype('uint16')
        ordfit_ycoord = ridge_eval(trace_fit['params'], xlist)
        ylower = np.int_(np.floor(ordfit_ycoord))
        #ylower_safe = trim_to_image_dims(xlist, ylower + 0, tmp_image)
        #yupper_safe = trim_to_image_dims(xlist, ylower + 1, tmp_image)
        tmp_image[trim_to_image_dims(xlist, ylower + 0, tmp_image)] = np.nan
        tmp_image[trim_to_image_dims(xlist, ylower + 1, tmp_image)] = np.nan
    sys.stderr.write("done.\n")
    return tmp_image

def get_trace_xya(trace_fit):
    """
    Return X position, ridge Y position, and pixel apron from the specified
    trace parameter dictionary. 
    
    NOTE: positions are in array coordinates (0-indexed)
    """
    ridge_x = np.arange(trace_fit['xmin'], trace_fit['xmax']).astype('uint16')
    ridge_y = ridge_eval(trace_fit['params'], ridge_x)
    return (ridge_x, ridge_y, trace_fit['apron'])


def mask_from_traces(imshape, trace_list, vlevel=0):
    mask_image = np.zeros(imshape, dtype='bool')
    trace_coords = []
    n_traces = len(trace_list)
    for i,trace_fit in enumerate(trace_list, 1):
        if (vlevel >= 1):
            sys.stderr.write("\rAdding trace %d of %d ... " % (i, n_traces))
        xlist = np.arange(trace_fit['xmin'], trace_fit['xmax']).astype('uint16')
        ordfit_ycoord = ridge_eval(trace_fit['params'], xlist)
        ylower = np.int_(np.floor(ordfit_ycoord))
        yc_list, xc_list = [], []
        apron_pix = trace_fit['apron']
        for offset in range(-apron_pix + 1, apron_pix + 1):
            xc_list.append(xlist)
            yc_list.append(ylower + offset)
            pass
        trace_coords.append((np.vstack(yc_list), np.vstack(xc_list)))
    return trace_coords


###--------------------------------------------------------------------------##
### Trace object for echelle orders:
#class Trace(object):
#
#    def __init__(self):
#        pass
#
#    @staticmethod
#    def _column_average(imvals, ypixels):
#        return np.sum(ypixels * imvals, axis=0) / np.sum(imvals, axis=0)

##--------------------------------------------------------------------------##
## Ridge fitting and evaluation:
def fit_polynomial(xvec, yvec, poly=2, rlm=True):
    results = {'xmin':xvec.min(), 'xmax':xvec.max()}
    design_matrix = np.column_stack([xvec**i for i in range(poly + 1)])
    if rlm:
        best_fit = sm.RLM(yvec, design_matrix).fit()
    else:
        best_fit = sm.OLS(yvec, design_matrix).fit()
    results['params'] = best_fit.params
    results['fitting'] = best_fit
    return results

def theil_sen_fit(xvec, yvec, poly=1, rlm=False):
    results = {'xmin':xvec.min(), 'xmax':xvec.max()}
    results['params'] = ts.linefit(xvec, yvec, weighted=False, joint=True)
    results['fitting'] = None
    return results

def fit_yridge(spectrum_order, poly=2, rlm=True):
    xcoo, ycoo, flux = spectrum_order
    return fit_polynomial(xcoo, ycoo, poly=poly, rlm=rlm)

## Evaluate a polynomial ridge fit:
def ridge_eval(model, xpos):
    return np.sum([cc*np.float_(xpos)**i for i,cc in enumerate(model)], axis=0)

## Return ridge x,y array coordinates:
def ridge_pos_2d(rspec):
    xcoords = np.arange(rspec['xmin'], rspec['xmax'] + 1, dtype='float32')
    ycoords = ridge_eval(rspec['params'], xcoords)
    return xcoords, ycoords

## Fit polynomials to all orders in a spectrum:
def fit_spectrum_ridges_fluxes(spectrum_orders, ypoly=2, fpoly=2, vlevel=0):
    ridges = []
    fluxes = []
    n_orders = len(spectrum_orders)
    for i, (xcoo, ycoo, flux) in enumerate(spectrum_orders):
        if vlevel >= 1:
            sys.stderr.write("\rFitting order %d of %d ... " % (i+1, n_orders))
        ridges.append(fit_polynomial(xcoo, ycoo, poly=ypoly))
        fluxes.append(fit_polynomial(xcoo, flux, poly=fpoly))
    if vlevel >= 1:
        sys.stderr.write("done.\n")
    return (ridges, fluxes)

## Evaluate all orders:
def splat_orders_onto_image(image, ridge_list, fluxes_list, dtype='float32'):
    orderpos = np.zeros_like(image, dtype=dtype)
    for rr,ff in zip(ridge_list, fluxes_list):
        xvalues, yvalues = ridge_pos_2d(rr)
        xvalues, fvalues = ridge_pos_2d(ff)
        xcoords = np.int_(xvalues)
        y_upper = np.int_(np.ceil(yvalues))
        y_lower = np.int_(np.floor(yvalues))
        orderpos[y_upper+1, xcoords] = fvalues
        orderpos[y_upper+0, xcoords] = fvalues
        orderpos[y_lower-0, xcoords] = fvalues
        orderpos[y_lower-1, xcoords] = fvalues
    return orderpos

def splat_orders_onto_image(image, ridge_list, fluxes_list, dtype='float32'):
    orderpos = np.zeros_like(image, dtype=dtype)
    for rr,ff in zip(ridge_list, fluxes_list):
        xvalues, yvalues = ridge_pos_2d(rr)
        xvalues, fvalues = ridge_pos_2d(ff)
        xcoords = np.int_(xvalues)
        y_upper = np.int_(np.ceil(yvalues))
        y_lower = np.int_(np.floor(yvalues))
        orderpos[y_upper+1, xcoords] = fvalues
        orderpos[y_upper+0, xcoords] = fvalues
        orderpos[y_lower-0, xcoords] = fvalues
        orderpos[y_lower-1, xcoords] = fvalues
    return orderpos




##--------------------------------------------------------------------------##
## Ridge object for tracing echelle orders:
class Ridge(object):

    def __init__(self, image, bmask=None):
        """
        Initialize ridge detector. Inputs:
        image       -- 2D image with spectrum to trace
        bounds      -- [optional] where to stop extracting (e.g., baffle mask)
        """
        self.idata = image
        self.bmask = bmask
        pass

    # ---------------------------------------------
    # Follow-the-ridge driver routine:
    def extract(self, yapprox, xcolumns, apron, nshift=40, 
            mincounts=None, maxdepth=0, vlevel=0):
        """
        Main extraction driver routine. Spreads outwards from initial
        guess, following flux 'ridge' until the signal is lost or an edge
        is reached. Returns X,Y (0-indexed) coordinates of identified ridge.

        yapprox     -- approximate Y-coord of ridge (array coords)
        xcolumns    -- slice with column range for initial guess/fit
        apron       -- half-size of re-centroiding box (pixels)
        maxdepth    -- [debugging] limit the number of extension iterations
        vlevel      -- verbosity control
        """

        # MAXDEPTH warning:
        if maxdepth > 0:
            sys.stderr.write("WARNING: maxdepth in use: %d\n" % maxdepth)

        # Get data, perform initial linear fit:
        ysection = self._make_yslice(yapprox, apron)
        r1x, r1y = self._guess_ridge_from_slices(self.idata, xcolumns, ysection)
        #initial_fit = nrex.fit_polynomial(r1x, r1y, poly=1)
        initial_fit = fit_polynomial(r1x, r1y, poly=1)
        #sys.stderr.write("r1x: %s\n" % str(r1x))
        #sys.stderr.write("r1y: %s\n" % str(r1y))
        #sys.stderr.write("initial_fit: %s\n" % str(initial_fit))

        ## DEBUGGING: evaluate initial_fit for inspection:
        #sys.stderr.write("-------------------------------------------\n")
        #dbg_fitted_y = ridge_eval(initial_fit['params'], r1x)
        #sys.stderr.write("initial_fit_debug:\n")
        #sys.stderr.write("%s\n" % str(np.vstack((r1x, dbg_fitted_y)).T))
        #sys.stderr.write("-------------------------------------------\n")


        # Refine ridge position using symmetric apron, refit:
        r2x, r2y, r2counts = self._recenter_ridge_ixmodel(self.idata,
                r1x, initial_fit['params'], apron)
        #starter_fit = nrex.fit_polynomial(r2x, r2y, poly=1)
        starter_fit = fit_polynomial(r2x, r2y, poly=1)
        #sys.stderr.write("starter_fit: %s\n" % str(starter_fit))
        #sys.stderr.write("r2x: %s\n" % str(r2x))
        #sys.stderr.write("r2y: %s\n" % str(r2y))
        ##asdf = raw_input()

        ## DEBUGGING: evaluate starter_fit for inspection:
        #sys.stderr.write("-------------------------------------------\n")
        #dbg_fitted_y = ridge_eval(starter_fit['params'], r1x)
        #sys.stderr.write("starter_fit_debug:\n")
        #sys.stderr.write("%s\n" % str(np.vstack((r1x, dbg_fitted_y)).T))
        #sys.stderr.write("-------------------------------------------\n")
        ##return (r2x, dbg_fitted_y)

        # Extend initial fit in both directions:
        extkw = {'apron':apron, 'mincounts':mincounts, 
                'maxdepth':maxdepth, 'vlevel':vlevel}
        rsegs = self._extend_ridge_to_edge(self.idata, starter_fit, 
                nudgepix=nshift, **extkw)
        lsegs = self._extend_ridge_to_edge(self.idata, starter_fit, 
                nudgepix=-nshift, **extkw)

        # Combine segments:
        segments = [(r2x, r2y)]
        segments.extend(rsegs)
        segments.extend(lsegs)

        # Separate coordinates, return sorted results:
        xlist, ylist = zip(*segments)
        xlist, ylist = np.hstack(xlist), np.hstack(ylist)
        order = np.argsort(xlist)
        return xlist[order], ylist[order]

    # ---------------------------------------------
    @staticmethod
    def _make_yslice(ycenter, apron):
        """Make slice centered on ycenter (array coords) with +/- apron."""
        ylower = ycenter - apron
        n_rows = 2*apron + 1
        yupper = ylower + n_rows
        ylower = ylower if ylower >= 0 else 0
        return slice(ylower, yupper)

    @staticmethod
    def _get_xyi_from_slices(image, xslice, yslice):
        """Extract arrays of X-position, Y-position, and image values for the
        region specified by the input slices."""

        xpixels = np.arange(image.shape[1])[xslice]
        ypixels = np.arange(image.shape[0])[yslice]
        return (xpixels, ypixels, image[yslice, xslice])

    # Basic y-centroid routine:
    @staticmethod
    def _ysum_ridge(imvals, ypixels):
        ytotals = np.sum(imvals, axis=0)
        ycenter = np.sum(ypixels * imvals, axis=0) / ytotals
        return (ycenter, ytotals)
    
    # ---------------------------------------------
    # Initial ridge guess (slices for input):
    def _guess_ridge_from_slices(self, image, xslice, yslice):
        xvec, yvec, vals = self._get_xyi_from_slices(image, xslice, yslice)
        yridge, ytotals = self._ysum_ridge(vals, yvec[:, None])
        return (xvec, yridge)

    # ---------------------------------------------
    # Intersection finder for collision check:
    def _edge_crossed(self, bdata):
        if np.any(bdata):
            crossings = np.diff(bdata)
            edge_index = crossings.nonzero()[0]    # array index of edge
            if edge_index.size != 1:
                raise   # insane!
            edge_index = edge_index[0]
            return ( True, edge_index)
        else:
            return (False, None)
    
    # Collision/interference prediction:
    def _trim_to_edge(self, bmask, xpixels, rmodel, apron):
        """This routine trims the input x-segment as needed to avoid collision
        with a boundary mask (i.e., NRES baffle). The trimmed x-segment and
        updated 'inbounds' criterion are returned.
        
        If a boundary is crossed, x-segment is truncated, inbounds=False
        If no boundary is crossed, x-segment is unmodified, inbounds=True
        """
        fitted_y = ridge_eval(rmodel, xpixels)
        ylower = np.int_(np.floor(fitted_y))
        lo_bound = bmask[ylower + 0, xpixels]
        hi_bound = bmask[ylower + 1, xpixels]

        # check if boundary is left or right:
        b_side = 'left' if (lo_bound[0] or hi_bound[0]) else 'right'
        result = [self._edge_crossed(x) for x in [lo_bound, hi_bound]] 
        edgepix = [pp for cc,pp in result if cc]
        ncrossed = len(edgepix)

        if ncrossed > 0:
            edge_index = np.sort(edgepix)[-1]
            edge_x_pix = xpixels[edge_index]
            if b_side == 'left':
                left_edge = edge_x_pix + 1
                truncated_xpixels = np.arange(left_edge, xpixels.max()+1)
            else:
                right_edge = edge_x_pix
                truncated_xpixels = np.arange(xpixels.min(), right_edge+1)
            return (truncated_xpixels, False)
        return (xpixels, True)

    # ---------------------------------------------
    # Model-assisted refinement of ridge position:
    def _recenter_ridge_ixmodel(self, image, xpixels, rmodel, apron):
        """
        Re-centroid ridge Y-coordinate using apron ~centered on best-fit model.
    
        image       --  2D image with extractable signal
        xpixels     --  numpy array of X positions [array coords, 0 ... N-1]
        rmodel      --  polynomial fit parameters
        apron       --  half-size of re-centroiding box (pixels)
        """
        fitted_y = ridge_eval(rmodel, xpixels)
        #fitted_y = nrex.ridge_eval(rmodel, xpixels)
        y_apron = np.arange(-apron, apron)      # biased low, use with np.ceil
        new_yarea = np.int_(np.ceil(fitted_y) + y_apron[:, None])
        new_ivals = image[new_yarea, xpixels]
        new_yridge, ytotals = self._ysum_ridge(new_ivals, new_yarea)
        return (xpixels, new_yridge, ytotals)

    # ---------------------------------------------
    # Workhorse: iteratively extend initial guess to edge of image:
    def _extend_ridge_to_edge(self, image, segfit, nudgepix, apron, 
            bmask=None, mincounts=None, maxdepth=0, vlevel=0):
        """
        Given a linear segment as an initial guess, this routine traces a ridge
        along a spectrum order and returns the X,Y coordinates it finds. 
        Starting from the initial segment, the best-fit is iteratively extended
        in both +/- X directions. Extension ends when:
        * ridge reaches edge of image
        * scatter in ridge Y position exceeds half apron size
        * other stopping point reached (e.g., baffle) [NOT YET IMPLEMENTED]
        * ridge disappears (no light)                 [NOT YET IMPLEMENTED]
    
        Inputs:
        ----------
        image      ---  2D image with spectrum orders to trace
        segfit     ---  fit to initial line segment (dict)
        nudgepix   ---  shift distance between segments
        apron      ---  size of +/- apron in cross-dispersion direction
                        (should be large enough to enclose entire trace)
        bmask      ---  [optional] array containing 2D extraction boundary
                            in image format (e.g., baffle mask)
        mincounts  ---  [optional] count level below which extraction stops
        maxdepth   ---  [optional] limit the number of extension iterations
        vlevel     ---  verbosity control
        """

        niy, nix = image.shape
        old_xseg = np.arange(segfit['xmin'], segfit['xmax'])
        fit_pars = segfit['params']
        segments = []
        inbounds = True
        lo_bound = 1        # fixed at pixel 1 for now
        hi_bound = niy - 2  # stop short of rightmost pixel
        ext_iter = 0

        while inbounds:
            ext_iter += 1
            new_xseg = old_xseg + nudgepix
            xmin, xmax = new_xseg.min(), new_xseg.max()

            # Stop if segment is fully out-of-bounds:
            if (xmin >= hi_bound) or (xmax <= lo_bound):
                if vlevel >= 1:
                    sys.stderr.write("out-of-bounds!\n")
                break

            # Trim left edge to boundary if hit:
            if (xmin <= lo_bound):
                xmin = lo_bound
                new_xseg = np.arange(xmin, xmax+1)
                inbounds = False    # stop after this iteration

            # Trim right edge to boundary if hit:
            if (xmax >= hi_bound):
                xmax = hi_bound
                new_xseg = np.arange(xmin, xmax+1)
                inbounds = False    # stop after this iteration

            ## Predict collision:
            if isinstance(self.bmask, np.ndarray):
                new_xseg, inbounds = self._trim_to_edge(self.bmask, 
                        new_xseg, fit_pars, apron)
                xmin, xmax = new_xseg.min(), new_xseg.max()

            if vlevel >= 0:
                msg = "\rFitting segment %d <= X <= %d ... " % (xmin, xmax)
                sys.stderr.write(msg)
                #sys.stderr.write("%s ... " % str(fit_pars))

            # Stop if edge is reached:
            if ((xmin < 0) or (nix < xmax)):
                if vlevel >= 1:
                    sys.stderr.write("edge!\n")
                break

            # Find new ridge segment:
            nrx, nry, ytotals = \
                self._recenter_ridge_ixmodel(image, new_xseg, fit_pars, apron)
            cts_per_pix = ytotals / (2.0 * apron + 1.0)

            if (vlevel >= 3):
                sys.stderr.write("\n")
                sys.stderr.write("nrx: %s\n" % str(nrx))
                sys.stderr.write("nry: %s\n" % str(nry))
                sys.stderr.write("tot: %s\n" % str(ytotals))
                sys.stderr.write("cts: %s\n" % str(cts_per_pix))
                sys.stderr.write("\n")

            # Look for scatter in ridge Y position:
            ypos_avg = np.average(nry)
            ry_scatter = np.std(nry)
            ypos_avg, ypos_std = np.average(nry), np.std(nry)
            #ypos_med, ypos_mad = calc_ls_med_MAD(nry)
            #ypos_med, ypos_mad = calc_ls_med_IQR(nry)
            if vlevel >= 1:
                sys.stderr.write(" ypos_avg=%.3f" % ypos_avg)
                sys.stderr.write(" ypos_std=%.3f" % ypos_std)
                #sys.stderr.write(" ypos_med=%.3f" % ypos_med)
                #sys.stderr.write(" ypos_mad=%.3f" % ypos_mad)
                sys.stderr.write(" ytotals=%.3f " % np.average(ytotals))
                sys.stderr.write("\n")
            #sys.stderr.write("ytotals: %s\n" % str(ytotals))
            ypos_scatter = ypos_std
            #ypos_scatter = ypos_mad
            if ((mincounts != None) and (np.average(ytotals) < mincounts)):
                if vlevel >= 1:
                    sys.stderr.write("Signal lost (flux cut)!\n")
                break
    
            if (ypos_scatter >= 0.5*apron):
                if vlevel >= 1:
                    sys.stderr.write("Signal lost!\n")
                break

            #fit_results = nrex.fit_polynomial(nrx, nry, poly=1)
            #fit_results = nrex.theil_sen_fit(nrx, nry)
            fit_results = theil_sen_fit(nrx, nry)
            fit_pars = fit_results['params']
            fit_stat = fit_results['fitting']
            segments.append((nrx, nry))
            old_xseg = new_xseg

            # Stop early if maximum iterations reached:
            if (maxdepth > 0) and (ext_iter >= maxdepth):
                if (vlevel >= 1):
                    sys.stderr.write("depth break!!\n")
                break
            pass

        if vlevel == 0:
            sys.stderr.write("done.\n")
        return segments

##--------------------------------------------------------------------------##

######################################################################
# CHANGELOG (nres_extraction.py):
#---------------------------------------------------------------------
#
#  2018-12-02:
#     -- Increased __version__ to 0.4.1.
#     -- Now include a comment divider before parameters in trace HDUs to
#           simplify visual inspection.
#
#  2018-11-30:
#     -- Increased __version__ to 0.4.0.
#     -- Added flat_rel_solver() method (FOX extraction method).
#
#  2018-08-03:
#     -- Increased __version__ to 0.3.7.
#     -- Commented out unused ylower_safe and yupper_safe calculation in
#           overplot_traces() method. These were doubling the calculation
#           time for pixel identification (may speed up overall process).
#
#  2018-05-09:
#     -- Increased __version__ to 0.3.6.
#     -- Added trim_to_image_dims() and overplot_traces() routines.
#
#  2018-05-08:
#     -- Increased __version__ to 0.3.5.
#     -- Added new TraceIO class.
#
#  2017-10-30:
#     -- Increased __version__ to 0.3.0.
#     -- Automatic extraction seems to largely work!
#     -- Implemented left-boundary adjustment in extend_ridge_to_edge.
#     -- Added mincounts parameter to extend_to_edge routine.
#     -- Added Ridge() object tracing assistant.
#
#  2017-08-15:
#     -- Increased __version__ to 0.2.1.
#     -- Now import theil_sen estimation module (testing).
#
#  2017-08-15:
#     -- Increased __version__ to 0.2.0.
#     -- Added Trace object.
#
#  2017-08-14:
#     -- Increased __version__ to 0.1.0.
#     -- First created nres_extraction.py.
#
