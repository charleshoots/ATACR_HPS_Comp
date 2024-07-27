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
import obstools as obs
import cmath
from comp_tools import *
from pathlib import Path
from scipy.signal import csd as _csd
import obspy.imaging.cm as cm
from ObsQA.TOOLS.io import *
# ====================================================================================================================
# ====================================================================================================================
# ====================================================================================================================
NoiseFolder = '/Users/charlesh/Documents/Codes/OBS_Methods/NOISE'
CompFolder = NoiseFolder + '/COMPS/ATaCR_NC'
MethodsFolder = NoiseFolder + '/METHODS'
# ===================================================================================================
# ============================================  LOAD DATA ===========================================
# ===================================================================================================

plotfolder = NoiseFolder + '/COMPS/FigureArchive/_GEN5'
compfolder = CompFolder
ATaCR_Py_DataFolder = OBS.TOOLS.io.dir_libraries(MethodsFolder + '/ATaCR')[1]
dirs = ATaCR_Py_DataFolder
datafolder = ATaCR_Py_DataFolder['Py_DataParentFolder']
eventsfolder = ATaCR_Py_DataFolder['Py_CorrectedTraces']
catalog = pd.read_pickle(eventsfolder + '/event_catalog_updated.pkl')
Station,evi = catalog.iloc[22],3
Event = Station.Events[evi]
# catalog = pd.read_pickle('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python/Metrics/EVENTS/EventMetrics_using_STA_avgTFs.pkl')
# catalog = pd.read_pickle('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python/EVENTS/event_catalog_updated.pkl')
# catalog = pd.read_pickle(eventsfolder + '/sta_catalog_evrecord_set_goodchans_updated.pkl')
# catalog = catalog.drop(index=29)
catalog = pd.read_pickle(eventsfolder + '/sta_catalog_proxima_test.pkl')

Folder = Path(plotfolder) / 'MeetingFigs'
Folder.mkdir(exist_ok=True)


# evaudit = ObsQA.io.audit_events(eventsfolder)
evaudit = pd.read_pickle(Path(eventsfolder) / 'event_record_audit.pkl')


def smooth(d,k=10):
        return np.convolve(d, np.ones(k) / k, mode='same')
# ###=====================================================================================================================================================
# ###=====================================================================================================================================================
# ========================================================================================================================================================
# ================================================================CODE SNIPPETS===========================================================================
# ========================================================================================================================================================
# for a in cm._cmap_names_categorical:
#     display(cm.__dict__[a].resampled(4))
# # ======================================================================================================================================================
# # ======================================================================================================================================================
# for (Event,Station,Metrics,Comp) in OBS_Generator(catalog,ATaCR_Py_DataFolder['Py_DataParentFolder']):
#     print(Event)
# for i,(Event,Station,Metrics,Comp) in zip(range(1),OBS_Generator(catalog,ATaCR_Py_DataFolder['Py_DataParentFolder'])):
#     print(Event)
# # ======================================================================================================================================================
# # ======================================================================================================================================================
# display(catalog)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

nsta_bar = 10
# |||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||| Nx2 . COHERENCE RECORDS . PLOT CODE ||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||
# evaudit = ObsQA.io.audit_events(eventsfolder)
evaudit = pd.read_pickle(Path(eventsfolder) / 'event_record_audit.pkl')
evaudit = evaudit[evaudit.nsta>=nsta_bar]
# evsilike = ['2011.242.06.57','2010.246.16.35','2011.191.00.57','2011.236.17.46','2010.163.19.26','2010.096.22.15']
# evsilike = ['2011.191.00.57','2010.096.22.15']
# evsi = [np.where(evaudit.Event==d)[0][0] for d in evsilike]
# evaudit = evaudit.iloc[evsi]

# evaudit = evaudit[evaudit.MW>=7.1]
# display(evaudit)
tapers = [0]
# methods = ['PostATACR','PostHPS']
# methods = ['PostATACR']
methods = ['PostHPS','PostATACR']
ysep_scl = 1.3
figsize = (10,10)
return_noise = True
for correction_method in methods:
  coh_comp = correction_method.replace('PostHPS','HPS').replace('PostATACR','ATaCR')
  if correction_method=='PostHPS':
    return_hps = True
    return_atacr = False
  else:
    return_hps = False
    return_atacr = True
  
  OutFolder = Path(plotfolder)
  SubFolders = Path('EventRecords') / correction_method / 'coherence'
  OutFolder = OutFolder / SubFolders
  OutFolder.mkdir(parents=True,exist_ok=True)
  for evi,ev in enumerate(evaudit.iloc):
      event = ev.Event
      File = 'm' + str(ev.MW) + '_z' + str(ev.depth) + 'km' + '_' + event + '_' + correction_method.replace('Post','') + '_COH'
      title = File.replace('_',' | ').replace('z','z: ').replace('m','mag: m')

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
      Metrics = []
      print('=='*25)
      print('[' + str(evi) + '/' + str(len(evaudit)) + '] ' + File)
      print('|| [' + str(evi) + '/' + str(len(evaudit)) + '] ' + event + ' | ')
      print('||---Begin load')
      for i,(net,sta) in enumerate(zip(networks,stations)):
        try:
          M,Comp = get_metrics_comp(net,sta,MethodsFolder,event,return_hps=return_hps,return_atacr=return_atacr,events_folder='EVENTS')
          M['Noise'] = get_Noise(dirs['Py_DataParentFolder'],net,sta,'sta')['Noise']
          Metrics.append(M.copy())
        except:
          _ = stations.pop(i)
          _ = networks.pop(i)
          continue
        post_record += Comp[correction_method].copy()
        pre_record += Comp['RawZ'].copy()
        del Comp
      if len(Metrics)==0:
        continue
      evstream = post_record.copy()
      evstream_back = pre_record.copy()
      title = event
      sortindex = list(np.argsort([np.abs(st.stats.sac.stel*1000) for st in evstream]))
      sortindex.reverse()
      Metrics = [Metrics[s] for s in sortindex]
      fig, axes = plt.subplots(nrows=2, ncols=2,figsize=figsize,layout='constrained',squeeze=True,sharey='all',sharex='all')
      axes = axes.flatten()
      for ci,cmp in enumerate(['ZP','ZZ','1Z','2Z']):
        ax = axes[ci]
  # -------------
        if cmp=='ZP':
          [ax.scatter(M['Noise'].Coherence(cmp)[0],M['Noise'].Coherence(cmp)[1] + ysep*ysep_scl,c='gray',s=1) for ysep,M in enumerate(Metrics)]
          [ax.plot(M['Raw'].Coherence(cmp)[0],M['Raw'].Coherence(cmp)[1] + ysep*ysep_scl,c='b',linewidth=0.7) for ysep,M in enumerate(Metrics)]
        if cmp=='ZP':
          [ax.plot(M[coh_comp].Coherence(cmp)[0],M[coh_comp].Coherence(cmp)[1] + ysep*ysep_scl,c='r',linewidth=0.1) for ysep,M in enumerate(Metrics)]
        else:
          [ax.plot(M[coh_comp].Coherence(cmp)[0],M[coh_comp].Coherence(cmp)[1] + ysep*ysep_scl,c='r',linewidth=0.8) for ysep,M in enumerate(Metrics)]
  # -------------
        [ax.plot(M[coh_comp].Coherence(cmp)[0],M[coh_comp].Coherence(cmp)[1]*0 + ysep*ysep_scl,c='k',linewidth=0.1) for ysep,M in enumerate(Metrics)]
        fn = [fnotch(np.round(1000*abs(M['Raw'].traces.select(channel='*Z')[0].stats.sac.stel))) for M in Metrics]
        yticks = [ysep*ysep_scl for ysep,_ in enumerate(Metrics)]
        [ax.plot([f,f],[y,y+1],c='k',linewidth=0.4) for f,y in zip(fn,yticks)]
        [ax.text(f,y+1,'fn',horizontalalignment='center',fontsize=7,fontweight='bold') for f,y in zip(fn,yticks)]
        yticklabels = [str(round(1000*abs(M['Raw'].traces.select(channel='*Z')[0].stats.sac.stel))) + 'm [' + str( M['Raw'].traces.select(channel='*Z')[0].stats.network) + '] ' + str( M['Raw'].traces.select(channel='*Z')[0].stats.station) for M in Metrics]
        ax.set_xscale('log')
        f = Metrics[0][coh_comp].Coherence(cmp)[0]
        ax.set_xlim(f[1],f[-1])
        ax_title = cmp + ' Coherence'
        ax.set_title(ax_title,fontweight='bold')
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        ax.set_ylim(yticks[0],yticks[-1]+1.3)
      fig.suptitle(File.replace('_',' | '),fontweight='bold',fontsize=15)
      save_format = 'png'
      # fig.suptitle(fig_title.replace('_',' | ').replace('to',' to '),fontweight='bold',fontsize=15)
      # File = File
      outfile = OutFolder / (File + '.' + save_format)
      save_tight(outfile,format=save_format,dpi=500)
      # (OutFolder / 'PNG').mkdir(parents=True,exist_ok=True)
      # fig.savefig(OutFolder / 'PNG' / (File + '.png'),format='png',dpi=900)
      plt.close('all')
  print('000'*30)
  print(correction_method + ' ||---COHERENCE RECORDS COMPLETE---||')
  print('000'*30)
        # || 26min to complete