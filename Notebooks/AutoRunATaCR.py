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
# ------------------------------------------------------------------------------------------------------------------------
from pathlib import Path
project_path = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp')
import shutil
import numpy as np
import pandas as pd
import sys
sys.path.append(str(project_path / 'Packages'))
sys.path.insert(0, str(project_path / 'Packages' / 'ATaCR'))
sys.path.insert(0, str(project_path / 'Packages' / 'CompCode'))
sys.path.insert(0, str(project_path / 'Packages' / 'ATaCR'/ 'OBStools'))
import ObsQA
from comp_tools import *


# --------------------------------------------------------------------------------------------------
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# --------------------------------------------------------------------------------------------------
# Test-1,Test-2,Test-3
# 
# 
# 
# 
# 
# 
# --------------------------------------------------------------------------------------------------
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# --------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------
# ============================================ FOLDERS ===========================================
# ---------------------------------------------------------------------------------------------------
project_path = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/ATACR_HPS_Comp')
ATaCR_DataFolder = str(project_path / '_DataArchive' / 'ATaCR_Data')

dirs = OBS.TOOLS.io.dir_libraries(ATaCR_DataFolder)[1]
datafolder = dirs['Py_DataParentFolder']
eventsfolder = dirs['Py_CorrectedTraces']
eventsfolder = dirs['Py_CorrectedTraces']
ATaCR_Parent = dirs['Py_DataParentFolder']
catalog_full = pd.read_excel(str(project_path / '_DataArchive' / 'utilities' / 'Janiszewski_etal_2023_StationList.xlsx'))
catalog = pd.read_pickle(eventsfolder + '/sta_catalog_proxima_test.pkl')
# ---------------------------------------------------------------------------------------------------
# ============================================ LOAD DATA ===========================================
# ---------------------------------------------------------------------------------------------------
catalog = pd.read_pickle(eventsfolder + '/event_catalog_updated.pkl')
Station,evi = catalog.iloc[22],3
Event = Station.Events[evi]
# catalog = pd.read_pickle('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python/Metrics/EVENTS/EventMetrics_using_STA_avgTFs.pkl')
# catalog = pd.read_pickle('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python/EVENTS/event_catalog_updated.pkl')
# catalog = pd.read_pickle(eventsfolder + '/sta_catalog_evrecord_set_goodchans_updated.pkl')
# catalog = catalog.drop(index=29)
catalog = pd.read_pickle(eventsfolder + '/sta_catalog_proxima_test.pkl')
# evaudit = ObsQA.io.audit_events(eventsfolder)
evaudit = pd.read_pickle(Path(eventsfolder) / 'event_record_audit.pkl')
# [XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX]
# [XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX][XX]


# ________________________________________________________________________________________________________________________________________________________________________
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
catalog = catalog[catalog.Station=='M07A']
ev_ind = np.where(np.array(catalog.iloc[0].Events)=='2012.181.21.07')[0][0]
catalog.Origin = [[catalog.Origin.iloc[0][ev_ind]]]
catalog.Metadata = [[catalog.Metadata.iloc[0][ev_ind]]]
catalog.Magnitude_mw = [[catalog.Magnitude_mw.iloc[0][ev_ind]]]
catalog.Events = [[catalog.Events.iloc[0][ev_ind]]]
catalog.Files = [[catalog.Files.iloc[0][ev_ind]]]
catalog.Depth_KM = [[catalog.Depth_KM.iloc[0][ev_ind]]]


days = ['2011.310',
 '2011.355',
 '2011.357',
 '2012.011',
 '2012.037',
 '2012.043',
 '2012.104',
 '2012.121',
 '2012.153',
 '2012.164']

# days = ['2012.153',
#  '2012.175',
#  '2012.130',
#  '2011.318',
#  '2012.006',
#  '2012.133',
#  '2012.177',
#  '2012.042',
#  '2011.348',
#  '2011.333',
#  '2012.061',
#  '2012.012',
#  '2012.111',
#  '2012.054',
#  '2012.185']

days = [UTCDateTime.strptime(d,'%Y.%j') for d in days]
# days = []
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ________________________________________________________________________________________________________________________________________________________________________


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
## ===============================================================================
    # def __init__(self, daylist=None):

    #     def _load_dn(day):
    #         exmpl_path = Path(resource_filename('obstools', 'examples'))
    #         fn = '2012.'+day+'*.SAC'
    #         fn = exmpl_path / 'data' / fn
    #         st = read(str(fn))
    #         tr1 = st.select(component='1')[0]
    #         tr2 = st.select(component='2')[0]
    #         trZ = st.select(component='Z')[0]
    #         trP = st.select(component='H')[0]
    #         window = 7200.
    #         overlap = 0.3
    #         key = '7D.M08A'
    #         return DayNoise(tr1, tr2, trZ, trP, window, overlap, key)

    #     self.daylist = []
    #     self.initialized = False
    #     self.QC = False
    #     self.av = False
    #     self.direc = None

    #     if isinstance(daylist, DayNoise):
    #         daylist = [daylist]
    #     elif daylist == 'demo' or daylist == 'Demo':
    #         print("Uploading demo data - March 01 to 04, 2012, station " +
    #               "7D.M08A")
    #         self.daylist = [_load_dn('061'), _load_dn(
    #             '062'), _load_dn('063'), _load_dn('064')]
    #     if not daylist == 'demo' and daylist:
    #         self.daylist.extend(daylist)
    #         self.station_depth = self.daylist[0].station_depth


cat = catalog.copy()
# display(cat)
event_mode = False
Minmag,Maxmag=6.0,8.0
fork = False
STEPS = [4,5]
## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##
## =============================================================================== ## =============================================================================== ##
for STEP in STEPS:
    if STEP==-3:
        NoiseFolder = dirs['Py_RawDayData']
        print('Day Noise While Loop Mode')
        ObsQA.TOOLS.io.DayNoiseWhileLoop(cat,NoiseFolder,ATaCR_Parent,days=15,attempts=100)
    else:
        for ii,Station in enumerate(cat.iloc):
            ## StaFolder = Path(dirs['Py_RawDayData']) / Station.StaName
            ## Files = list(StaFolder.glob('*.SAC'))
            staname = Station.StaName
            subfolder = staname + '/'
            print('[//////////////////////////]'*2)
            print('----Station: ' + staname +  ' (' + str(ii+1) + ' of ' + str(len(cat)) + ')')
            icatalog = Station.to_frame().T
            print('[//////////////////////////]'*2)
            ObsQA.TOOLS.io.Run_ATaCR(icatalog,days=days,fork=fork,event_mode=event_mode, ATaCR_Parent = ATaCR_Parent,STEPS=[STEP],log_prefix=Station.StaName,Minmag=Minmag,Maxmag=Maxmag)
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