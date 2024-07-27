# ===================================================================================================
# ============================================ IMPORTS ==============================================
# ===================================================================================================
import sys
sys.path.insert(1,'/Users/charlesh/Documents/Codes/')
sys.path.append('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS')
sys.path.insert(0, '/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/COMPS')
sys.path.insert(0, '/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python/OBStools')
sys.path.insert(0, '/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python')
import requests 
from PIL import Image 
import math
import scipy
import numpy as np
import librosa
import os
import shutil
from scipy.signal import stft, detrend
os.environ['PYDEVD_WARN_SLOW_RESOLVE_TIMEOUT'] = '2'
from obspy import Trace
import librosa.display
import matplotlib.pyplot as plt
import matplotlib.image as img
import matplotlib.gridspec as gridspec
import matplotlib
import sys
import obspy
from obspy.signal import PPSD

from scipy import interpolate

import pickle as pkl
import glob as g
from obspy.clients.fdsn import Client
import datetime
import re
import math
from numpy import linalg as eigen
# import noisecut
import matplotlib.colors as mcolors
import matplotlib.cm as cm2
from scipy.stats import norm
import scipy.stats as stats
from scipy import fft
# import ntk
# from cmcrameri import cm

from scipy.interpolate import RBFInterpolator, InterpolatedUnivariateSpline #<----Experimental

import ObsQA as OBS
# from ObsQA.classes import OBSMetrics as OBSM
# from ObsQA.plots import qtp
from ObsQA import *
# OBSM = OBSMetrics
# import obstools as obs
import cmath
from comp_tools import *
from pathlib import Path
from scipy.signal import csd as _csd
import obspy.imaging.cm as cm

# ===================================================================================================
# ============================================  LOAD DATA ===========================================
# ===================================================================================================


NoiseFolder = '/Users/charlesh/Documents/Codes/OBS_Methods/NOISE'
CompFolder = NoiseFolder + '/COMPS/ATaCR_NC'
MethodsFolder = NoiseFolder + '/METHODS'
plotfolder = NoiseFolder + '/COMPS/FigureArchive/_GEN5'
compfolder = CompFolder
ATaCR_Py_DataFolder = OBS.TOOLS.io.dir_libraries(MethodsFolder + '/ATaCR')[1]
dirs = ATaCR_Py_DataFolder
datafolder = ATaCR_Py_DataFolder['Py_DataParentFolder']
eventsfolder = ATaCR_Py_DataFolder['Py_CorrectedTraces']
catalog = pd.read_pickle(eventsfolder + '/event_catalog_updated.pkl')
Station,evi = catalog.iloc[22],3
Event = Station.Events[evi]

catalog = pd.read_pickle(eventsfolder + '/sta_catalog_proxima_test.pkl')

Folder = Path(plotfolder) / 'MeetingFigs'
Folder.mkdir(exist_ok=True)
def smooth(d,k=10):
        return np.convolve(d, np.ones(k) / k, mode='same')

def metric_fractional_different(r,A,B):
    a_f,a_coh = A.Coherence(r)
    b_f,b_coh = B.Coherence(r)
    b_coh = interpolate.interp1d(b_f,b_coh)(a_f)
    frac = a_coh / b_coh
    return a_f,frac

# ===================================================================================================
# ===================================================================================================
# ===================================================================================================

datafolder = dirs['Py_DataParentFolder']
nsta_bar = 20
# evaudit = ObsQA.io.audit_events(eventsfolder)
evaudit = pd.read_pickle(Path(eventsfolder) / 'event_record_audit.pkl')
evaudit = evaudit[evaudit.nsta>=nsta_bar]
# evaudit = evaudit.iloc[100:]
# evsilike = ['2010.247.08.52','2013.143.17.19']
# evsilike = ['2010.247.08.52']
# evsilike = ['2011.242.06.57','2010.246.16.35','2011.191.00.57','2011.236.17.46','2010.163.19.26','2010.096.22.15']
# evsilike = ['2010.073.08.08','2010.216.07.15'] # '2010.151.19.51','2010.204.22.51'

# evsi = [np.where(evaudit.Event==d)[0][0] for d in evsilike]
# evaudit = evaudit.iloc[evsi]
# evaudit = evaudit.iloc[2:]
# evaudit = evaudit[evaudit.MW>=7.1]
# display(evaudit)
# folder = 'grouped_bands'
folder = 'separated_bands'
alumni_poster = Path('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/COMPS/FigureArchive/_GEN5/EventRecords/AlumniPoster/Events')
ysep_scl = 1.3
figsize = (14,13)
bands = [(1/10,1),(1/30,1/10),(1/100,1/30)]
trim = (10,7200)
tapers = [1]
# methods = ['PostATACR','PostHPS']
methods = ['PostATACR']
# methods = ['PostHPS']
channels = ['Z','1','2']
for correction_method in methods:
  current_method = correction_method.replace('Post','')
  coh_comp = correction_method.replace('PostHPS','HPS').replace('PostATACR','ATaCR')
  if correction_method=='PostHPS':
    return_hps = True
  else:
    return_hps = False
  if return_hps:
    return_atacr = False
  else:
    return_atacr = True
  for taper_mode in tapers:
    OutFolder = Path(plotfolder)
    SubFolders = Path('EventRecords') / ('Taper_' + str(taper_mode)) / correction_method / folder
    OutFolder = OutFolder / SubFolders
    # OutFolder = alumni_poster
    OutFolder.mkdir(parents=True,exist_ok=True)
    fig = plt.figure(figsize=(15,5))
    for evi,ev in enumerate(evaudit.iloc):
        print('=='*40)
        mirror_test = mirror_audit(ev.to_frame().T,datafolder=datafolder)[0]
        mirror_test = np.array(mirror_test[0]) & np.array(mirror_test[1])
        mirror_test = np.where(mirror_test)[0]
        if len(mirror_test)==0:
          print('|| No data mirrors between sets')
          continue
        else:
          print('|| >>> DATASETS MIRROR POPULATION: ' + str(len(mirror_test)))
        if len(mirror_test)<nsta_bar:
          print('Insufficient mirrored data (min:' + str(nsta_bar) + ')')
          continue
        event = ev.Event
        stations = ev.Stations
        networks = ev.Networks.tolist()
        networks = [networks[i] for i in mirror_test]
        stations = [stations[i] for i in mirror_test]
        evdepth = ev.depth
        post_record = Stream()
        pre_record = Stream()
        print('|| [' + str(evi) + '/' + str(len(evaudit)) + '] ' + event + ' | ')
        print('||---Begin load')
        # plt.figure(figsize=(15,5))
        for i,(net,sta) in enumerate(zip(networks,stations)):
            # if ev.Event=='2010.151.19.51':
            #     if sta=='C08W':
            #         _,_= stations.pop(i),networks.pop(i)
            #         continue
            #     if sta=='S01W':
            #         _,_ = stations.pop(i),networks.pop(i)
            #         continue
            Metrics,Comp = get_metrics_comp(net,sta,MethodsFolder,event,return_hps=return_hps,return_atacr=return_atacr,events_folder='EVENTS_Taper_' + str(taper_mode))

            fq,frac = metric_fractional_different('ZP',Metrics['ATaCR'],Metrics['ATaCR'].Noise)
            if np.sum(frac>=1)>(0.95*len(frac)):
              #  Metrics,Comp = get_metrics_comp(net,sta,MethodsFolder,event,return_hps=return_hps,return_atacr=return_atacr,events_folder='EVENTS_Taper_' + str(taper_mode))
               continue
            s = 2
            plt.scatter(fq[frac>=1],frac[frac>=1],c='b',s=s,marker='.')
            plt.scatter(fq[frac<1],frac[frac<1],c='r',s=s,marker='.')
            plt.axhline(1,linestyle=':',linewidth=0.5,c='k')
            # plt.xlim([fq.min(),fq.max()])
            # plt.ylim([frac.min(),frac.max()])
            if len(Metrics)==0:
                _ = stations.pop(i)
                _ = networks.pop(i)
                continue
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim([fq.min(),fq.max()])
        plt.ylim([frac.min(),frac.max()])
        plt.ylabel('Compliance Correction Misfit')
        # fig
plt.show()