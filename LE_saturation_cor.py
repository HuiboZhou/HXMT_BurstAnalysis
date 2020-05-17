#!/usr/bin/env python
from astropy.io import fits 
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
import sys

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
    hdulist[1].writeto(outfile,overwrite=True)

def get_args(argv):
    keywords = []
    values  = []
    for i in range(len(argv)):
        if i == 0: continue
        keywords.append(argv[i].split('=')[0])
        values.append(argv[i].split('=')[1])
    arguments = dict(zip(keywords, values))
    return arguments

def genlc(data,binsize=1,fig=False,rate=True,pannel=True):
    N = (max(data)-min(data))/binsize
    N = int(N)
    lc = np.histogram(data,N)[0]
    if rate:
        lc = lc/binsize # calculate counts rate instead of counts
    lc_time = np.histogram(data,N)[1][0:-1]
    #null = np.where(lc == 0)
    #lc = np.delete(lc,null)
    #lc_time = np.delete(lc_time,null)
    #print null
    return lc_time,lc

def genlc_bin(data, bins, rate=True):
    time_resolution = bins[1]-bins[0]
    lc = np.histogram(data,bins=bins)[0]
    if rate:
        binsize=bins[1]-bins[0]
        lc = lc/binsize # calculate counts rate instead of counts
    lc_time = np.histogram(data,bins=bins)[1][0:-1] + time_resolution/2
    #null = np.where(lc == 0)
    #lc = np.delete(lc,null)
    #lc_time = np.delete(lc_time,null)
    #print null
    return lc_time,lc


def select_oneboxtime(time, boxnum=0, **kwargs):
    le_small_detid = np.array([0,2,3,4,6,7,8,9,10,12,
        14,20,22,23,24,25,26,28,30,32,34,35,36,38,39,40,
        41,42,44,46,52,54,55,56,57,58,60,61,62,64,66,67,
        68,70,71,72,73,74,76,78,84,86,88,89,90,92,93,94])
    #data:
    detid = kwargs['detid']
    detbox= kwargs['detbox']
    channel=kwargs['pi']
    #parameters
    channelmin = kwargs['minPI']
    channelmax = kwargs['maxPI']

    new_time = time[(np.isin(detid,le_small_detid))&
            (channel>=channelmin)&
            (channel<=channelmax)&
            (detbox==boxnum)]
    return new_time

def genlc_forcedtrigger(time_evt, evttype, detbox, **kwargs):
    #read parameter
    tstart = kwargs['starttime']
    tstop  = kwargs['stoptime']
    time_resolution = kwargs['binsize']
    boxnum = kwargs['boxnum']
    #--
    detbox = detbox[(time_evt>=tstart)&(time_evt<=tstop)]
    evttype = evttype[(time_evt>=tstart)&(time_evt<=tstop)]
    time_evt = time_evt[(time_evt>=tstart)&(time_evt<=tstop)]
    forcedtrigger_time = time_evt[(evttype==1)&(detbox==boxnum)] #select forced trigger evt for one detBox 
    forcedtrigger_lc_x, forcedtrigger_lc_y = genlc_bin(forcedtrigger_time, 
            bins=np.arange(tstart, tstop+time_resolution, time_resolution))
    return forcedtrigger_lc_x, forcedtrigger_lc_y

def coeff_forcedtrigger_rate(time, rate, detboxid=0):
    coeff = np.array([1000./x if x>0 else 0 for x in rate])
    return coeff

def coeff_detnumnorm(time, fte_rate, detboxid=0):
    if detboxid == 0:
        detnum = 19
    if detboxid == 1:
        detnum = 20
    if detboxid == 2:
        detnum = 19
    coeff = np.array([1/detnum if x>0 else 0 for x in fte_rate])
    return coeff

def detnum_cor(time, *args):

    for i,rate in enumerate(args):
        if i == 0:
            coeff0 = np.array([1 if x>0 else 0 for x in rate])
        elif i == 1:
            coeff1 = np.array([1 if x>0 else 0 for x in rate])
        elif i == 2:
            coeff2 = np.array([1 if x>0 else 0 for x in rate])
    detnum_eff = coeff0 + coeff1 + coeff2
    return detnum_eff


def saturation_cor(screenfile, evtfile, **kwargs):
    if kwargs['binsize'] <0.005:
        raise ValueError("For the saturation correction of LE, the time resolution has to be greater than or equal to 5 milliseconds\n")
    hdulist = fits.open(screenfile)
    time_screen = hdulist[1].data.field("TIME")
    pi          = hdulist[1].data.field("PI")
    detid       = hdulist[1].data.field("DET_ID")
    detbox_screen = hdulist[1].data.field("DetBox")
    hdulist = fits.open(evtfile)
    time_evt    = hdulist[1].data.field("TIME")
    evttype     = hdulist[1].data.field("Event_Type")
    detbox_evt  = hdulist[1].data.field("DetBox")
    #FINISH read data
    #read parameters
    tstart = kwargs['starttime'] #fix LC time bins at tstart to tstop 
    tstop  = kwargs['stoptime']

    #good time selection for screen
    time_detbox0 = select_oneboxtime(time_screen, boxnum=0, detid=detid, detbox=detbox_screen, pi=pi, minPI=kwargs['minPI'], maxPI=kwargs['maxPI'])
    time_detbox1 = select_oneboxtime(time_screen, boxnum=1, detid=detid, detbox=detbox_screen, pi=pi, minPI=kwargs['minPI'], maxPI=kwargs['maxPI'])
    time_detbox2 = select_oneboxtime(time_screen, boxnum=2, detid=detid, detbox=detbox_screen, pi=pi, minPI=kwargs['minPI'], maxPI=kwargs['maxPI'])
    print("Entries for detBox0: ", len(time_detbox0))
    print("Entries for detBox1: ", len(time_detbox1))
    print("Entries for detBox2: ", len(time_detbox2))

    #generate LC for each detbox
    time_resolution = kwargs['binsize']
    lc_x_detbox0, lc_y_detbox0 = genlc_bin(time_detbox0, bins=np.arange(tstart, tstop+time_resolution, time_resolution))
    lc_x_detbox1, lc_y_detbox1 = genlc_bin(time_detbox1, bins=np.arange(tstart, tstop+time_resolution, time_resolution))
    lc_x_detbox2, lc_y_detbox2 = genlc_bin(time_detbox2, bins=np.arange(tstart, tstop+time_resolution, time_resolution))

#    fig, (ax1,ax2) = plt.subplots(2,1, sharex=True)
#    ax1.errorbar(lc_x_detbox0, lc_y_detbox0, drawstyle="steps-mid")
#    ax1.errorbar(lc_x_detbox1, lc_y_detbox1, drawstyle="steps-mid")
#    ax1.errorbar(lc_x_detbox2, lc_y_detbox2, drawstyle="steps-mid")

    #generate LC for forced trigger event based on bins of each detbox LC
    fte_lc_x_detbox0, fte_lc_y_detbox0 = genlc_forcedtrigger(time_evt, evttype, detbox=detbox_evt, starttime=tstart, stoptime=tstop,
            binsize=time_resolution, boxnum=0)
    fte_lc_x_detbox1, fte_lc_y_detbox1 = genlc_forcedtrigger(time_evt, evttype, detbox=detbox_evt, starttime=tstart, stoptime=tstop,
            binsize=time_resolution, boxnum=1)
    fte_lc_x_detbox2, fte_lc_y_detbox2 = genlc_forcedtrigger(time_evt, evttype, detbox=detbox_evt, starttime=tstart, stoptime=tstop,
            binsize=time_resolution, boxnum=2)


    coeff0 = coeff_forcedtrigger_rate(fte_lc_x_detbox0, fte_lc_y_detbox0, detboxid=0)
    coeff1 = coeff_forcedtrigger_rate(fte_lc_x_detbox1, fte_lc_y_detbox1, detboxid=1)
    coeff2 = coeff_forcedtrigger_rate(fte_lc_x_detbox2, fte_lc_y_detbox2, detboxid=2)
    coeff_detnorm0 = coeff_detnumnorm(fte_lc_x_detbox0, fte_lc_y_detbox0, detboxid=0)
    coeff_detnorm1 = coeff_detnumnorm(fte_lc_x_detbox1, fte_lc_y_detbox1, detboxid=1)
    coeff_detnorm2 = coeff_detnumnorm(fte_lc_x_detbox2, fte_lc_y_detbox2, detboxid=2)

    detnum_cor_coeff = detnum_cor(fte_lc_x_detbox0, fte_lc_y_detbox0, fte_lc_y_detbox1, fte_lc_y_detbox2)

    time_corrected = lc_x_detbox0
    rate_corrected = (lc_y_detbox0 * coeff0 * coeff_detnorm0 + 
        lc_y_detbox1 * coeff1 * coeff_detnorm1+ 
        lc_y_detbox2 * coeff2 * coeff_detnorm2)*58/detnum_cor_coeff

    #ERROR
    error_detbox0 = np.sqrt(lc_y_detbox0/time_resolution)
    error_detbox1 = np.sqrt(lc_y_detbox1/time_resolution)
    error_detbox2 = np.sqrt(lc_y_detbox2/time_resolution)
    error_corrected = 58*np.sqrt((error_detbox0*coeff0*coeff_detnorm0/detnum_cor_coeff)**2 + 
            (error_detbox1*coeff1*coeff_detnorm1/detnum_cor_coeff)**2 + 
            (error_detbox2*coeff2*coeff_detnorm2/detnum_cor_coeff)**2 )
    #NOTE:use propogating error instead of counts error.
    #error_corrected = np.sqrt(rate_corrected/time_resolution)

    ##
    if 'outfile' in kwargs:
        outfile = kwargs['outfile']
        create_lightcurve_file(time_corrected, rate_corrected, error_corrected, outfile, binsize=binsize, starttime=tstart, stoptime=tstop)
        print("-------------------------")
        print("\n Lightcurve file {} is generated Successfully \n".format(outfile))

    return time_corrected, rate_corrected


if __name__ == "__main__":
    arguments = get_args(sys.argv)
    evtfile = arguments['evtfile']
    screenfile= arguments['screenfile']
    binsize   = np.float(arguments['binsize'])
    starttime = np.float(arguments['starttime'])
    stoptime  = np.float(arguments['stoptime'])
    minPI = np.int(arguments['minPI'])
    maxPI = np.int(arguments['maxPI'])
    if 'outfile' in arguments:
        outfile = arguments['outfile']

    saturation_cor(screenfile, evtfile, minPI=minPI, maxPI=maxPI, binsize=binsize, starttime=starttime, stoptime=stoptime, outfile=outfile)

