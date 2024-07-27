# from obspy.core import read, Stream, Trace, AttribDict, UTCDateTime
# from obspy.clients.fdsn import Client
# import os
# import argparse
# import obspy
# import obstools as obs
# from obspy import taup
# from obstools.atacr import DayNoise, TFNoise, EventStream, StaNoise, utils
# import obstools.atacr.plotting as atplot
# from obstools.scripts import comply_calculate, atacr_clean_spectra, atacr_correct_event, atacr_daily_spectra, atacr_download_data, atacr_download_event, atacr_transfer_functions
# from stdb.scripts import query_fdsn_stdb
# import matplotlib.pyplot as plt
# import gc
# import fnmatch
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, '/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/COMPS')
sys.path.insert(0, '/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python/OBStools')
sys.path.insert(0, '/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python')
import ObsQA
from comp_tools import *
dirs = ObsQA.TOOLS.io.dir_libraries('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR')[1]
eventsfolder = dirs['Py_CorrectedTraces']
catalog_full = pd.read_excel('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python/utilities/Janiszewski_etal_2023_StationList.xlsx')
catalog = pd.read_pickle(eventsfolder + '/sta_catalog_proxima_test.pkl')
ATaCR_Parent = dirs['Py_DataParentFolder']
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## ===============================================================================
## STEPS = [1,2,3,4,5,6,7] ##Absolutely every step - Downloading adds hour(s) or more to the process
## STEPS = [2,3] ##Everything but the download steps - About 4min for six stations.
## STEPS = [4,5,6,7] ##Everything but the download steps - About 4min for six stations
## -------------------------------------------------------
## Step-1: Station Metadata. Step a0 in ML-ATaCR. Always run this.
## Step-2: Download event data. Step a3 in ML-ATaCR.
## Step-3: Download day data. Step a2 in ML-ATaCR.
## Step-4: Daily Spectra. Step b1 in ML-ATaCR.
## Step-5: Clean and Average Daily Spectra. Step b2 in ML-ATaCR.
## Step-6: Calculate transfer functions. Step b3 ML-ATaCR.
## Step-7: Correct events. Step b4 in ML-ATaCR.
## ===============================================================================
cat = catalog.copy()
print(cat)
event_mode = True
Minmag,Maxmag=7.1,8.0
fork = True
STEPS = [6,7]
## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##
for STEP in STEPS:
    if STEP==-3:
        NoiseFolder = dirs['Py_RawDayData']
        print('Day Noise While Loop Mode')
        ObsQA.TOOLS.io.DayNoiseWhileLoop(cat,NoiseFolder,ATaCR_Parent,days=10,attempts=100)
    for ii,Station in enumerate(cat.iloc):
        ## StaFolder = Path(dirs['Py_RawDayData']) / Station.StaName
        ## Files = list(StaFolder.glob('*.SAC'))
        staname = Station.StaName
        subfolder = staname + '/'
        print('[//////////////////////////]'*2)
        print('----Station: ' + staname +  ' (' + str(ii+1) + ' of ' + str(len(cat)) + ')')
        icatalog = Station.to_frame().T
        print('[//////////////////////////]'*2)
        ObsQA.TOOLS.io.Run_ATaCR(icatalog,fork=fork,event_mode=event_mode, ATaCR_Parent = dirs['Py_DataParentFolder'],STEPS=[STEP],log_prefix=Station.StaName,Minmag=Minmag,Maxmag=Maxmag)
        if event_mode & (STEP==3):
            Origins = Station.Origin
            Starts = Origins
            if isinstance(Origins[0],obspy.core.event.origin.Origin):
                    Starts = [e.time for e in Origins]
            dateformat = '%Y.%j.%H.%M'
            hps_data_folder = (Path(dirs['Py_RawDayData']) / Station.StaName / 'HPS_Data')
            hps_data_folder.mkdir(exist_ok=True)
            days = list(np.unique([s.strftime('%Y.%j') for s in Starts]))
            [[shutil.move(fi,hps_data_folder / fi.name) for fi in list((Path(dirs['Py_RawDayData']) / Station.StaName).glob(d + '*.SAC'))] for d in days]
            # shutil.
        print('....done')
## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##