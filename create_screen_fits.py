#!/usr/bin/env python
from astropy.io import fits 
import numpy as np
from astropy.time import Time

def create_screen_file(time, outfile, **kwargs):
    
    # Table
    c1 = fits.Column(name='TIME', array=time, format='1D')

    tb = fits.BinTableHDU.from_columns([c1])
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
    hdulist[1].header['HDUCLASS'] = 'OGIP'
    hdulist[1].header['HDUCLAS1'] = 'EVENTS'
    hdulist[1].header['HDUCLAS2'] = 'ALL'
    hdulist[1].header['LONGSTRN'] = 'OGIP 1.0'
    hdulist[1].header['TIMESYS'] = 'TT'

    hdulist[1].header['MJDREFI'] = 55927
    hdulist[1].header['MJDREFF'] = 7.6601852e-4
    hdulist[1].header['TIMEREF'] = 'LOCAL'
    hdulist[1].header['TASSIGN'] = 'SATELLITE'
    hdulist[1].header['TIMEUNIT']= 's'

    hdulist.writeto(outfile,overwrite=True)

if __name__ == "__main__":
    data = np.loadtxt("/Users/tuoyouli/Downloads/tt.txt")
    time = data
    create_screen_file(time, "/Users/tuoyouli/Downloads/tt.fits")


