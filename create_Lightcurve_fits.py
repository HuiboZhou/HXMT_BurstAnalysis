#!/usr/bin/env python
from astropy.io import fits 
import numpy as np
from astropy.time import Time

def create_lightcurve_file(time, rate, error, outfile, **kwargs):
    
    # Table
    c1 = fits.Column(name='TIME', array=time, format='1D')
    c2 = fits.Column(name='RATE', array=rate, format='1E')
    c3 = fits.Column(name='ERROR', array=error, format='1E')

    tb = fits.BinTableHDU.from_columns([c1,c2,c3])
    # Prmary Header
    header = fits.Header()
    primary_hdr = fits.Header()
    primary_hdr['comments'] = 'FITS(Flexible Image Transport System)'
    primary_hdu = fits.PrimaryHDU(header=primary_hdr)
    
    hdul = fits.HDUList([primary_hdu, tb])
    hdul.writeto(outfile,overwrite=True)
    
    # write Keywords
    TRUE = np.bool(True)
    FALSE = np.bool(False)
    hdulist = fits.open(outfile)
    hdulist[1].header['TIMESYS'] = 'TT'
    hdulist[1].header['MJDREFI'] = 55927
    hdulist[1].header['MJDREFF'] = 7.6601852e-4
    hdulist[1].header['TIMEREF'] = 'LOCAL'
    hdulist[1].header['TASSIGN'] = 'SATELLITE'
    hdulist[1].header['TIMEUNIT']= 's'
    hdulist[1].header['TIMEDEL'] = kwargs['binsize']
    hdulist[1].header['TIERRELA']= 1e-6
    hdulist[1].header['TIERABSO']= 5e-5
    hdulist[1].header['CLOCKAPP']= FALSE
    hdulist[1].header['TSTART']  = kwargs['starttime']
    hdulist[1].header['TSTOP']   = kwargs['stoptime']

    hdulist.writeto(outfile,overwrite=True)

if __name__ == "__main__":
    data = np.loadtxt("testdata/LE_lixb_ch106-1170_saturation_corrected.dat")
    time = np.array([x[0] for x in data])
    rate = np.array([x[1] for x in data])
    error= np.array([x[2] for x in data])
    create_lightcurve_file(time, rate, error, "testdata/test_LE_LC.FITS", starttime=min(time), stoptime=max(time), binsize=0.005)


