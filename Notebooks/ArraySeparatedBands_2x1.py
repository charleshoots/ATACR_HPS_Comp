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

# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||| 1x3 . SEPARATED BANDS . EVENT RECORDS . PLOT CODE |||||
# ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
datafolder = dirs['Py_DataParentFolder']
nsta_bar = 10
evaudit_focus = evaudit.copy()
evaudit_focus = evaudit_focus[evaudit_focus.nsta>=nsta_bar].copy()
# evaudit_focus = evaudit_focus.iloc[100:]
# evsilike = ['2010.247.08.52','2013.143.17.19']
# evsilike = ['2011.242.06.57','2010.246.16.35','2011.191.00.57','2011.236.17.46','2010.163.19.26','2010.096.22.15']
# evsilike = ['2010.073.08.08','2010.216.07.15'] # '2010.151.19.51','2010.204.22.51'
# evsi = [np.where(evaudit_focus.Event==d)[0][0] for d in evsilike]
# evaudit_focus = evaudit_focus.iloc[evsi]
# evaudit_focus = evaudit_focus.iloc[2:]
# evaudit_focus = evaudit_focus[evaudit_focus.MW>=7.1]
# display(evaudit_focus)
# folder = 'grouped_bands'
folder = 'separated_bands'
ysep_scl = 1.3
figsize = (14,13)
bands = [(1/10,1),(1/30,1/10),(1/100,1/30)]
trim = (10,7200)
# methods = ['PostATACR','PostHPS']
# methods = ['PostATACR']
methods = ['PostHPS','PostATACR']
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
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  OutFolder = Path(plotfolder)
  SubFolders = Path('EventRecords') / correction_method / folder
  OutFolder = OutFolder / SubFolders
  OutFolder.mkdir(parents=True,exist_ok=True)
  for evi,ev in enumerate(evaudit_focus.iloc):
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
      print('|| [' + str(evi) + '/' + str(len(evaudit_focus)) + '] ' + event + ' | ')
      print('||---Begin load')
      for i,(net,sta) in enumerate(zip(networks,stations)):
        if ev.Event=='2010.151.19.51':
          if sta=='C08W':
            _ = stations.pop(i)
            _ = networks.pop(i)
            continue
          if sta=='S01W':
            _ = stations.pop(i)
            _ = networks.pop(i)
            continue
        # try:
        Metrics,Comp = get_metrics_comp(net,sta,MethodsFolder,event,return_hps=return_hps,return_atacr=return_atacr,events_folder='EVENTS')
        if len(Metrics)==0:
            _ = stations.pop(i)
            _ = networks.pop(i)
            continue
        # except:
        #   _ = stations.pop(i)
        #   _ = networks.pop(i)
        #   continue
        post_record += Comp[current_method + '_Stream'].copy()
        pre_record += Comp['Raw'+ '_Stream'].copy()
        del Metrics
        del Comp
        if len(post_record)==0:
          continue
      print('||---Load complete')
      phases = ('P','S','SKS','PKiKP','SKiKS','SKSSKS',)
      # phases = ('P','S',)
      # phases=('ttall',)
      evstream = post_record.copy()
      evstream_back = pre_record.copy()
      facecolor=('b','r')
      title = event
      sortindex = None
      normscale = 0.7
      residual_fraction = 0.5
      for chan in channels:
        pre_record_chan = pre_record.select(channel='*'+chan).copy()
        post_record_chan = post_record.select(channel='*'+chan).copy()
        for bandi,(band,s) in enumerate(zip(bands,[[1,1],[1,1],[1,1]])):
          band_sec = np.sort([1/b for b in band])
          print('|| Plotting band: ' + str(band_sec[0]) + ' to ' + str(band_sec[-1]) + 's' + ' | channel:H' + chan)
          fig, axes = plt.subplots(nrows=1, ncols=2,figsize=figsize,layout='constrained',squeeze=False,sharey='all',sharex='all')
          # fig, axes = plt.subplots(nrows=1, ncols=3,figsize=figsize,layout='constrained',squeeze=False,sharey='all',sharex='all')
          # norms='postset'
          norms = 'trace'
          # norms = 'col'
          # norms = np.array([[abs(c.data).max() for c in s] for s in [pre_record,post_record]]).T.max(axis=1).tolist() #global max norm
          # norms = np.array([[abs(c.data).max() for c in s] for s in [pre_record]]).T.max(axis=1).tolist() #pre max norm
          # norms = np.array([[abs(c.data).max() for c in s] for s in [post_record]]).T.max(axis=1).tolist() #post max norm. >This works well at >m6.5
          for record_i,evstream in enumerate([pre_record_chan,post_record_chan]):
            # File = 'm' + str(ev.MW) + '_z' + str(ev.depth) + 'km' + '_' + event + '_' + str(band_sec[0]) + 'to' + str(band_sec[-1]) + 's_' + correction_method.replace('Post','') + '_T'
            File = 'm' + str(ev.MW) + '_z' + str(ev.depth) + 'km' + '_' + event + '_' + str(band_sec[0]) + 'to' + str(band_sec[-1]) + 's_' + correction_method.replace('Post','')
            fig_title = 'm' + str(ev.MW) + '_z' + str(ev.depth) + 'km' + '_' + '-'.join(event.split('.')[:2]) + ' ' + ':'.join(event.split('.')[2:]) + '_' + str(band_sec[0]) + 'to' + str(band_sec[-1]) + 's_Corrected using ' + correction_method.replace('Post','')
            ax = axes[0,record_i]
            # linewidth = [0.05,0.05,0.05][bandi]
            linewidth = [0.2,0.2,0.2][bandi]
            ax_title = ['{} Pre-Correction | Channel:{} \n '.format(current_method,chan),'{} Post-Correction | Channel:{} \n '.format(current_method,chan)][record_i]
            title = File.replace('_',' | ').replace('z','z: ').replace('m','mag: m')
            # ---------------
            # More reasonable black and white plots:
            ax = event_record_plot(evstream=evstream,evstream_back=None,norm=norms,scales=s,linewidth=linewidth,figsize=figsize,band=band,trim=trim,facecolor=facecolor,evdepth=evdepth,phases=phases,title=title,sortindex = sortindex,ax=ax,normscale=normscale,residual_fraction=residual_fraction)
            # Use this for those ugly fill mode plots:
            # ax = event_record_plot(evstream=evstream,evstream_back=None,norm=norm,scales=s,linewidth=linewidth,figsize=figsize,band=band,trim=trim,facecolor=facecolor,evdepth=evdepth,phases=phases,title=title,sortindex = sortindex,ax=ax,normscale=normscale,residual_fraction=residual_fraction)
            # ---------------

            ax.set_title(ax_title,fontweight='bold')          
# --------------------------------------------- Switch between grouped-band and separated bands plots
          # ##### Separated bands
          save_format = 'png'
          fig.suptitle(fig_title.replace('_',' | ').replace('to',' to '),fontweight='bold',fontsize=15)
          s_mag,s_depth,s_chan,s_bands,s_method = ['m' + str(ev.MW), 'z' + str(int(ev.depth)) + 'km' , 'H' + chan , 'to'.join([str(int(b)) for b in band_sec]) + 's',correction_method.replace('Post','')]
          tags = [s_chan,s_bands,s_method,s_mag+s_depth,event]
          File = '_'.join(tags) + '.{}'.format(save_format)
          outfile = OutFolder / File
          # print('{} saved to {}'.format(File,OutFolder))
          save_tight(outfile,format=save_format,dpi=200)
          # fig.savefig(outfile,format=save_format,dpi=400)
      # ##### Grouped bands
      # fig.suptitle(File.replace('_',' | ').replace('to',' to '),fontweight='bold',fontsize=15)
      # save_tight(OutFolder / (File + '.eps' ),dpi=600)
# --------------------------------------------
          plt.close('all')
  print('000'*30)
  print(correction_method + ' ||---EV RECORDS COMPLETE---||')
  print('000'*30)