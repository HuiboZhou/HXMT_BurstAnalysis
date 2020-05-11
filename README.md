# A burst analysis toolkit for HXMT data

## LE staturation correction

Correct the light curve for each DetBox of LE telescope based on the forced trigger event rate and of the three LE DetBox. The light curves of the three DetBoxs are weighted and averaged according to the number of detectors in each DetBox.

```LE_saturation_cor.py```

USAGE: 

```python LE_saturation_cor.py evtfile="/path/to/raw_Evt_file" screenfile="/path/to/LE_screen_file" binsize=0.005 minPI=106 maxPI=1170 starttime=262708467 stoptime=262708468 outfile=/path/to/output/Lightcurve_file.fits ```
        
