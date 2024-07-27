from ObsQA.imports import *
from obspy.core import UTCDateTime as _UTCDateTime
import glob as g
import pandas as pd
import numpy as np
import pickle as pkl
from obspy.clients.fdsn import Client as _Client
import datetime as _datetime
import os as os
from pathlib import Path
import concurrent
from concurrent.futures import wait
import ObsQA
from ObsQA import classes
import time
import logging
# from classes import OBSMetrics
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def Run_ATaCR(catalog, ATaCR_Parent = None, STEPS=[1,2,3,4,5,6,7], netsta_names=None, chan='H', Minmag=6.3, Maxmag=6.7, limit=1000, pre_event_min_aperture=1, pre_event_day_aperture=30, dailyspectra_flags='-O --figQC --figAverage --figCoh --save-fig', cleanspectra_flags='-O --figQC --figAverage --figCoh --figCross --save-fig', tf_flags='-O --figTF --save-fig', correctevents_flags='--figRaw --figClean --save-fig',logoutput_subfolder=None,log_prefix = '',staquery_output = './sta_query.pkl',fork=True,days=10,seed='MESSI_22FIFA_WORLD_CUP!',max_workers=1,event_mode=False,event_dt=None,taper_mode=0):
        # dailyspectra_flags='-O --figQC --figAverage --figCoh --save-fig'
        # STEPS = [1,2,3,4,5,6,7] #Absolutely every step - Downloading adds hour(s) or more to the process
        # STEPS = [2,3] #Everything but the download steps - About 4min for six stations.
        # STEPS = [4,5,6,7] #Everything but the download steps - About 4min for six stations
        #### -------------------------------------------------------
        # Step-1: Station Metadata. Step a0 in ML-ATaCR. Always run this.
        # Step-2: Download event data. Step a3 in ML-ATaCR.
        # Step-3: Download day data. Step a2 in ML-ATaCR.
        # Step-4: Daily Spectra. Step b1 in ML-ATaCR.
        # Step-5: Clean and Average Daily Spectra. Step b2 in ML-ATaCR.
        # Step-6: Calculate transfer functions. Step b3 ML-ATaCR.
        # Step-7: Correct events. Step b4 in ML-ATaCR.
        dirs = ObsQA.io.dir_libraries('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR')[1]
        if 1 in STEPS:
                print('Step 1/7 - BEGIN: Station Metadata')
                # C='?H?' #channels
                # !query_fdsn_stdb -N {','.join(N)} -C '{C}' -S {','.join(S)} ./Data/sta_query> ./Data/Step_1_7_StationMeta_logfile.log
                ObsQA.io.build_staquery(d=catalog,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
                print('Step 1/7 - COMPLETE: Station Metadata')
        if 2 in STEPS:
                print('Step 2/7 - BEGIN: Download Event Data')
                if logoutput_subfolder is None:
                        logoutput_subfolder = dirs['Py_Logs'] + '/2_7'
                ObsQA.io.DownloadEvents(catalog,ATaCR_Parent=ATaCR_Parent,netsta_names=netsta_names,Minmag=Minmag,Maxmag=Maxmag,limit=limit,pre_event_min_aperture=pre_event_min_aperture,logoutput_subfolder=logoutput_subfolder,staquery_output=staquery_output,chan=chan,log_prefix=log_prefix)
                print('Step 2/7 - COMPLETE: Download Event Data')
        if 3 in STEPS:
                print('Step 3/7 - BEGIN: Download Day Data')
                if logoutput_subfolder is None:
                        logoutput_subfolder = dirs['Py_Logs'] + '/3_7'
                # ObsQA.io.DownloadNoise(catalog,ATaCR_Parent=ATaCR_Parent,netsta_names=netsta_names,pre_event_day_aperture=pre_event_day_aperture,logoutput_subfolder=logoutput_subfolder,staquery_output=staquery_output,chan=chan,log_prefix=log_prefix)
                # ObsQA.io.DayNoiseWhileLoop(catalog,NoiseFolder,ATaCR_Parent,days=10,attempts=50)
                ObsQA.io.DownloadDayNoise(catalog,days=days,seed=seed,ATaCR_Parent=ATaCR_Parent,netsta_names=netsta_names,logoutput_subfolder=logoutput_subfolder,log_prefix = log_prefix,staquery_output=staquery_output,chan=chan,event_mode=event_mode)
                print('Step 3/7 - COMPLETE: Download Day Data')
        if 4 in STEPS:
                print('Step 4/7 - BEGIN: Quality Control Noise Data')
                if logoutput_subfolder is None:
                        logoutput_subfolder = dirs['Py_Logs'] + '/4_7'
                ObsQA.io.DailySpectra(catalog,ATaCR_Parent=ATaCR_Parent,netsta_names=netsta_names,extra_flags=dailyspectra_flags,logoutput_subfolder=logoutput_subfolder,staquery_output=staquery_output,chan=chan,log_prefix=log_prefix,fork=fork,max_workers=max_workers)
                print('Step 4/7 - COMPLETE: Quality Control Noise Data')
        if 5 in STEPS:
                print('Step 5/7 - BEGIN: Spectral Average of Noise Data')
                # !atacr_clean_spectra -O --figQC --figAverage --figCoh --figCross --save-fig --start='{SpecStart}' --end='{SpecEnd}' ./Data/sta_query.pkl> ./Data/Step_5_7_CleanSpectra_logfile.log
                if logoutput_subfolder is None:
                        logoutput_subfolder = dirs['Py_Logs'] + '/5_7'
                ObsQA.io.CleanSpectra(catalog,ATaCR_Parent=ATaCR_Parent,netsta_names=netsta_names,extra_flags=cleanspectra_flags,logoutput_subfolder=logoutput_subfolder,staquery_output=staquery_output,chan=chan,log_prefix=log_prefix,fork=fork,max_workers=max_workers)
                print('Step 5/7 - COMPLETE: Spectral Average of Noise Data')
        if 6 in STEPS:
                print('Step 6/7 - BEGIN: Calculate Transfer Functions')
                # !atacr_transfer_functions -O --figTF --save-fig ./Data/sta_query.pkl> ./Data/Step_6_7_CalcTFs_logfile.log
                if logoutput_subfolder is None:
                        logoutput_subfolder = dirs['Py_Logs'] + '/6_7'
                ObsQA.io.TransferFunctions(catalog,ATaCR_Parent=ATaCR_Parent,netsta_names=netsta_names,taper_mode=taper_mode,extra_flags=tf_flags,logoutput_subfolder=logoutput_subfolder,staquery_output=staquery_output,chan=chan,log_prefix=log_prefix,fork=fork,max_workers=max_workers)
                print('Step 6/7 - COMPLETE: Calculate Transfer Functions')
        if 7 in STEPS:
                print('Step 7/7 - BEGIN: Correct Event Data')
                # !atacr_correct_event --figRaw --figClean --save-fig ./Data/sta_query.pkl> ./Data/Step_7_7_CorrectEvents_logfile.log
                if logoutput_subfolder is None:
                        logoutput_subfolder = dirs['Py_Logs'] + '/7_7'
                ObsQA.io.CorrectEvents(catalog,ATaCR_Parent=ATaCR_Parent,netsta_names=netsta_names,extra_flags=correctevents_flags,logoutput_subfolder=logoutput_subfolder,staquery_output=staquery_output,chan=chan,log_prefix=log_prefix,fork=fork,max_workers=max_workers)
                print('Step 7/7 - COMPLETE: Correct Event Data')
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# eof