from __future__ import division
from astropy.io import fits 
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time


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
    lc = np.histogram(data,bins=bins)[0]
    if rate:
        binsize=bins[1]-bins[0]
        lc = lc/binsize # calculate counts rate instead of counts
    lc_time = np.histogram(data,bins=bins)[1][0:-1]
    #null = np.where(lc == 0)
    #lc = np.delete(lc,null)
    #lc_time = np.delete(lc_time,null)
    #print null
    return lc_time,lc

def read_one_box(filename, evtfilename):
    data = np.loadtxt(filename)
    evttype = np.array([x[0] for x in data])
    detid   = np.array([x[1] for x in data])
    channel = np.array([x[2] for x in data])
    data = np.loadtxt(evtfilename)
    time = data
    return time, detid, channel, evttype

def select(time, detid, channel, evttype):
    le_small_detid = np.array([0,2,3,4,6,7,8,9,10,12,
        14,20,22,23,24,25,26,28,30,32,34,35,36,38,39,40,41,42,44,46,52,54,55,56,57,58,60,61,62,64,66,67,68,70,71,72,73,74,76,78,84,86,88,89,90,92,93,94])
    new_time =       time[(evttype==0)&(np.isin(detid,le_small_detid))&(channel>=400)&(channel<=4000)]
    new_evttype = evttype[(evttype==1)&(np.isin(detid,le_small_detid))&(channel>=400)&(channel<=4000)]
    new_detid =     detid[(evttype==1)&(np.isin(detid,le_small_detid))&(channel>=400)&(channel<=4000)]
    new_channel=  channel[(evttype==1)&(np.isin(detid,le_small_detid))&(channel>=400)&(channel<=4000)]
    return new_time, new_evttype, new_detid, new_channel

def czcf_lc(time, evttype, resolution):
    print(len(time), len(evttype))
    time0, czcf_rate = genlc1(time[evttype==1], bins=np.arange(262708466.999120, 262708469.876226,resolution))
    return time0, czcf_rate

def coeff(time, czcf_rate0, czcf_rate1, czcf_rate2):
    x0 = 19
    x1 = 20
    x2 = 19
    print("LLLLLLLL", czcf_rate0[czcf_rate0<1000])
    coeff0 = np.array([x0*1000./x if x>0 else 0 for x in czcf_rate0])
    coeff1 = np.array([x1*1000./x if x>0 else 0 for x in czcf_rate1])
    coeff2 = np.array([x2*1000./x if x>0 else 0 for x in czcf_rate2])
    print("LLLLLLLLL", coeff0, coeff1, coeff2)

    return coeff0, coeff1, coeff2

def cal_coeff1(time, czcf_rate0, czcf_rate1, czcf_rate2):
    coeff0 = np.array([19 if x>0 else 0 for x in czcf_rate0])
    coeff1 = np.array([20 if x>0 else 0 for x in czcf_rate1])
    coeff2 = np.array([19 if x>0 else 0 for x in czcf_rate2])

    return coeff0, coeff1, coeff2

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
            bins=np.arange(tstart, tstop, time_resolution))
    return forcedtrigger_lc_x, forcedtrigger_lc_y

def coeff_forcedtrigger_rate(time, rate, detboxid=0):
    if detboxid == 0:
        detnum = 19
    if detboxid == 1:
        detnum = 20
    if detboxid == 2:
        detnum = 19
    coeff = np.array([x*1000./x if x>0 else 0 for x in rate])
    return coeff

def detnum_cor(time, *args):

    for i,rate in enumerate(args):
        if i == 0:
            coeff0 = np.array([19 if x>0 else 0 for x in rate])
        elif i == 1:
            coeff1 = np.array([20 if x>0 else 0 for x in rate])
        elif i == 2:
            coeff2 = np.array([19 if x>0 else 0 for x in rate])
    detnum_eff = coeff0 + coeff1 + coeff2
    return detnum_eff


def saturation_cor(screenfile, evtfile, **kwargs):
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

    #generate LC for each detbox
    time_resolution = kwargs['binsize']
    lc_x_detbox0, lc_y_detbox0 = genlc_bin(time_detbox0, bins=np.arange(tstart, tstop, time_resolution))
    lc_x_detbox1, lc_y_detbox1 = genlc_bin(time_detbox1, bins=np.arange(tstart, tstop, time_resolution))
    lc_x_detbox2, lc_y_detbox2 = genlc_bin(time_detbox2, bins=np.arange(tstart, tstop, time_resolution))

    #generate LC for forced trigger event based on bins of each detbox LC
    fte_lc_x_detbox0, fte_lc_y_detbox0 = genlc_forcedtrigger(time_evt, evttype, detbox=detbox_evt, starttime=min(time_screen), stoptime=max(time_screen),
            binsize=time_resolution, boxnum=0)
    fte_lc_x_detbox1, fte_lc_y_detbox1 = genlc_forcedtrigger(time_evt, evttype, detbox=detbox_evt, starttime=min(time_screen), stoptime=max(time_screen),
            binsize=time_resolution, boxnum=1)
    fte_lc_x_detbox2, fte_lc_y_detbox2 = genlc_forcedtrigger(time_evt, evttype, detbox=detbox_evt, starttime=min(time_screen), stoptime=max(time_screen),
            binsize=time_resolution, boxnum=2)


    coeff0 = coeff_forcedtrigger_rate(fte_lc_x_detbox0, fte_lc_y_detbox0)
    coeff1 = coeff_forcedtrigger_rate(fte_lc_x_detbox1, fte_lc_y_detbox1)
    coeff2 = coeff_forcedtrigger_rate(fte_lc_x_detbox2, fte_lc_y_detbox2)

    detnum_cor_coeff = detnum_cor(fte_lc_x_detbox0, fte_lc_y_detbox0, fte_lc_y_detbox1, fte_lc_y_detbox2)


    time_corrected = lc_x_detbox0
    rate_corrected = (lc_y_detbox0 * coeff0/19 + 
        lc_y_detbox1 * coeff1/20 + 
        lc_y_detbox2 * coeff2/19)*58/(detnum_cor_coeff)  #NORMALIZED for detector number
    return time_corrected, rate_corrected


if __name__ == "__main__":
    time_resolution = 0.005
    x, y = saturation_cor("../sgr_262708467-262708468/P020402500803_LE_screen.fits",
            "../sgr_262708457-262708477_NEWLE/HXMT_P020402500803_LE-Evt_FFFFFF_V1_L1P_NEW.FITS",
            minPI=106, maxPI=1170, binsize=time_resolution, starttime=262708467, stoptime=262708468)
    with open("LE_lixb_ch106-1170_saturation_corrected.dat", 'w')as fout:
        for i in range(len(x)):
            fout.write("%.8f %f %f\n"%(x[i], y[i], np.sqrt(y[i])/time_resolution))
    plt.figure()
    plt.errorbar(x, y , np.sqrt(y)/0.005, ds='steps-mid')
    plt.show()

