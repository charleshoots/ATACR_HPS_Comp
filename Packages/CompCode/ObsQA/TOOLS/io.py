import ObsQA
# from ObsQA import *
# from ObsQA.imports import *
# from ObsQA import classes
# from ObsQA.noisecut import *
import glob as g
import numpy as np
import pandas as pd
import pickle as pkl
import obspy
from obspy import *
from obspy.core import UTCDateTime as UTCDateTime
from obspy.clients.fdsn import Client as _Client
import os
from pathlib import Path
import concurrent
from concurrent.futures import wait
import time
import logging
import matplotlib.pyplot as plt
from obstools.scripts import comply_calculate, atacr_clean_spectra, atacr_correct_event, atacr_daily_spectra, atacr_download_data, atacr_download_event, atacr_transfer_functions
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# from classes import OBSMetrics
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def audit_events(eventsfolder,Minmag=6.0,Maxmag=7.0):
    client = Client()
    timedelta = 60
    stations = ObsQA.TOOLS.io.get_event_catalog(eventsfolder)
    evlist = []
    [evlist.extend(stations.events[j]) for j in range(len(stations))]
    evlist = np.unique(evlist)
    evlist
    evaudit = dict()
    for ev in evlist:
        inds = [np.array([k for k in [np.where(np.array(stations.events[j])==ev)[0]]])[0] for j in range(len(stations))]
        for i in range(len(inds)):
            if len(inds[i])==0:
                inds[i] = None
            else:
                inds[i] = inds[i][0]
        idx = np.where([(j is not None) for j in inds])[0]
        stas = stations.Station[idx].tolist()
        evaudit[ev] = [stas]
        evaudit[ev] = [idx]
        ## evaudit['nsta'].append(len(stas))
    evaudit = pd.DataFrame.from_dict(evaudit,orient='index',columns=['Stations'])
    evaudit['Event'] = evaudit.index
    evaudit = evaudit.reset_index(drop=True)[['Event','Stations']]
    evaudit['nsta'] = [0 for _ in range(len(evaudit))]
    for j in range(len(evaudit)):
        evaudit.at[j,'Networks'] = np.array(list(stations.Network[evaudit.Stations[j]]))
        net = np.unique(np.array(list(stations.Network[evaudit.Stations[j]])))
        evaudit.at[j,'Stations'] = list(stations.Station[evaudit.Stations[j]])
        evaudit.at[j,'nsta'] = len(evaudit.iloc[j].Stations)
        
        print('[' + str(j) + '/' + str(len(evaudit)) + '] event: ' + evaudit.iloc[j].Event + ' || nsta: ' + str(evaudit.iloc[j].nsta))
        ev = evaudit.Event[j]
        starttime = UTCDateTime.strptime(str(ev),'%Y.%j.%H.%M')
        endtime = starttime + datetime.timedelta(seconds=timedelta)
        print('OBSPY Events: ' + str(starttime) + ' --> ' + str(endtime) + ' | Mags: m' + str(Minmag) +  '-' + str(Maxmag))
        cat = client.get_events(starttime=starttime, endtime=endtime,minmagnitude=Minmag, maxmagnitude=Maxmag)
        km = cat[0].origins[0].depth/1000
        mw = cat[0].magnitudes[0].mag
        lla = cat[0].origins[0].latitude,cat[0].origins[0].longitude
        evaudit.at[j,'MW'] = mw
        evaudit.at[j,'depth'] = km
        evaudit.at[j,'ev_lla'] = lla[0]
        evaudit.at[j,'ev_lon'] = lla[1]
    return evaudit
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def get_ATACR_HPS_Comp(d,atacr_TFs_used='ZP-21', win_length=200):
        RawP = d.Raw[0]['trP'].copy()
        RawZ = d.Raw[0]['trZ'].copy()
        PostATACR = d.Corrected[0][atacr_TFs_used].copy()
        PostATACR.stats.location = 'ATaCR (' + PostATACR.stats.location + ')'
        PostHPS, spectrograms = noisecut(RawZ.copy(), ret_spectrograms=True,win_length=win_length)
        PostBoth, spectrograms = noisecut(PostATACR.copy(), ret_spectrograms=True,win_length=win_length)
        return RawP,RawZ,PostATACR,PostHPS,PostBoth
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def get_arrivals(sta_llaz,ev_llaz,model = 'iasp91',phases=('ttall',)):
        ''''
        Simple function pulls arrival times and phase names for a given event observed at a given station
        sta_llaz = List object containing [Lat,Lon] of the station
        ev_llaz = List object containing [Lat,Lon] of the event
        phases = Tuple object containing a list of all desired phases. 'ttall'
        (default) will give every single phase available.
        -Charles Hoots,2022
        '''
        degdist = obspy.taup.taup_geo.calc_dist(ev_llaz[0],ev_llaz[1],sta_llaz[0],sta_llaz[1],6371,0)
        arrivals = obspy.taup.tau.TauPyModel(model=model).get_travel_times(ev_llaz[2], degdist,phase_list=phases)
        times = [[a.name,a.time] for a in arrivals]
        return times
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def obspyfft(d):
        signal = d.data
        fourier = np.abs(np.fft.rfft(signal))
        freq = np.abs(np.fft.rfftfreq(signal.size, d=1./d.stats.sampling_rate))
        return freq,fourier
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def AuditEventFolder(eventsfolder,subfolder='CORRECTED',parseby='pkl',Minmag=6.0,Maxmag=7.0):
        parse = {'SAC':'*Z.SAC','pkl':'*.pkl'}
        catalog = ObsQA.TOOLS.io.get_event_catalog(eventsfolder)
        stations = ObsQA.TOOLS.io.getstalist()
        # stations['StaName'] = stations.Network + '.' + stations.Station
        cols = stations.columns.tolist()
        client = Client()
        prefix = (catalog['Network'] + '.' + catalog['Station']).tolist()
        for ista in range(len(prefix)):
                sta = prefix[ista]
                folder = Path(eventsfolder) / sta
                folder = folder / subfolder
                fls = list(folder.glob(parse[parseby]))
                files = [str(fi).split('/')[-1] for fi in fls]
                evna = ['.'.join(f.split('.')[2:6]) for f in files if f.split('.')[-2]=='sta']
                mww = []
                depth_km = []
                origin_t = []
                event_meta = []
                averaging = [f.split('.')[-2] for f in files if f.split('.')[-2]=='sta']
                print(str(ista+1) + ' of ' + str(len(prefix)) + ' Sta: ' + sta + ', ' + str(len(evna)) + ' events found. Collecting metadata from IRIS..')
                for i,ev in enumerate(evna):
                        timedelta = 60
                        start = UTCDateTime.strptime(str(ev),'%Y.%j.%H.%M')
                        end = start + timedelta
                        cat = client.get_events(starttime=start, endtime=end,minmagnitude=Minmag, maxmagnitude=Maxmag)
                        mww.append(cat[0].magnitudes[0].mag)
                        depth_km.append(cat[0].origins[0].depth/1000)
                        origin_t.append(cat[0].origins[0].time)
                        event_meta.append(cat)
                stacat_id = np.where(stations.StaName==sta)[0][0]
                stations.iat[stacat_id,np.where(stations.columns=='Magnitude_mw')[0][0]] = mww
                stations.iat[stacat_id,np.where(stations.columns=='Depth_KM')[0][0]] = depth_km
                stations.iat[stacat_id,np.where(stations.columns=='Origin')[0][0]] = origin_t
                stations.iat[stacat_id,np.where(stations.columns=='Metadata')[0][0]] = event_meta
                stations.iat[stacat_id,np.where(stations.columns=='Averaging')[0][0]] = averaging
                stations.iat[stacat_id,np.where(stations.columns=='Events')[0][0]] = evna
                stations.iat[stacat_id,np.where(stations.columns=='Files')[0][0]] = files
                stations.iat[stacat_id,np.where(stations.columns=='n_events')[0][0]] = len(evna)
        catalog = stations
        catalog = catalog[catalog.n_events>0]
        if parseby=='SAC':
                catalog = catalog[cols]
        elif parseby=='pkl':
                cols.append('Averaging')
                catalog = catalog[cols]
        return catalog
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def get_event_catalog(eventsfolder,subfolder='CORRECTED',fmt = 'pkl'):
        evdict = dict()
        evdict['folder'] = [fld.split('/')[-2] for fld in g.glob(eventsfolder + '/*/')]
        evdict['Network'] = [fld.split('.')[0] for fld in evdict['folder']]
        evdict['Station'] = [fld.split('.')[1] for fld in evdict['folder']]
        evdict['StaName'] = [str(a) + '.' + str(b) for a,b in zip(evdict['Network'],evdict['Station'])]
        # evdict['folder'][0]
        evdict['n_events'] = list()
        evdict['events'] = list()
        for i in range(len(evdict['folder'])):
                staname =  evdict['StaName'][i]
                path = Path(eventsfolder)
                path = path / evdict['folder'][i] / subfolder
                files = list(path.glob('*.' + fmt))
                events = [str(f.name).replace(staname + '.','').replace('.sta','').replace('.day','').replace('.' + fmt,'') for f in files]
                # events = list(np.unique(['.'.join(files[g].split('.' + fmt)[0].split('.')[0:-1]) for g in range(len(files))]))
                evdict['events'].append(events)
                evdict['n_events'].append(len(events))
        catalog = pd.DataFrame(evdict)
        catalog = catalog.sort_values(by=['Network','Station'])
        catalog = catalog.reset_index(drop=True)
        return catalog
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def getstalist():
        current_path = os.path.dirname(__file__)
        excelfile = current_path + '/Janiszewski_etal_2023_StationList.xlsx'

        stations = pd.read_excel(excelfile)
        d = dict()
        d.update({a:b for a,b in zip(stations.columns,[c.replace('(','').replace(')','').replace(' ','_').replace('/','') for c in stations.columns])})
        stations = stations.rename(columns=d)
        stations.Station = [str(s) for s in stations.Station]
        stations.Network = [str(s) for s in stations.Network]
        stations['StaName'] = stations.Network.astype(str) + '.' + stations.Station.astype(str)
        allgood = np.in1d(stations[['Z_Is_Good','H1_Is_Good','H2_Is_Good','P_Is_Good']].sum(axis=1).tolist(),4)
        stations['Good_Channels'] = allgood
        stations = stations.assign(n_events=pd.Series())
        stations = stations.assign(Magnitude_mw=pd.Series())
        stations = stations.assign(Depth_KM=pd.Series())
        stations = stations.assign(Origin=pd.Series())
        stations = stations.assign(Metadata=pd.Series())
        stations = stations.assign(Averaging=pd.Series())
        stations = stations.assign(Events=pd.Series())
        stations = stations.assign(Files=pd.Series())
        stations = stations.sort_values(by=['Network','Station'])
        stations = stations.reset_index(drop=True)
        return stations
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def dir_libraries(CompFolder):
        ATaCR_ML_DataFolder = dict()
        ATaCR_ML_DataFolder['ML_ATaCR_Parent'] = CompFolder + '/ATaCR'
        ATaCR_ML_DataFolder['ML_DataParentFolder'] = CompFolder + '/ATaCR/DATA'
        ATaCR_ML_DataFolder['ML_RawDayData'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/datacache_day'
        ATaCR_ML_DataFolder['ML_PreProcDayData'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/datacache_day_preproc'
        ATaCR_ML_DataFolder['ML_RawEventData'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/datacache_event'
        ATaCR_ML_DataFolder['ML_PreProcEventData'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/datacache_event_preproc'
        ATaCR_ML_DataFolder['ML_StaSpecAvg'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/noisetc/AVG_STA'
        ATaCR_ML_DataFolder['ML_CorrectedTraces'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/noisetc/CORRSEIS'
        ATaCR_ML_DataFolder['ML_b1b2_StaSpectra'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/noisetc/SPECTRA'
        ATaCR_ML_DataFolder['ML_TransferFunctions'] = ATaCR_ML_DataFolder['ML_DataParentFolder'] + '/noisetc/TRANSFUN'

        ATaCR_Py_DataFolder = dict()
        ATaCR_Py_DataFolder['Py_DataParentFolder'] = CompFolder + '/ATaCR_Python'
        ATaCR_Py_DataFolder['Py_RawDayData'] = ATaCR_Py_DataFolder['Py_DataParentFolder'] + '/Data'
        # ATaCR_Py_DataFolder['Py_PreProcDayData']
        # ATaCR_Py_DataFolder['Py_RawEventData']
        # ATaCR_Py_DataFolder['Py_PreProcEventData']
        ATaCR_Py_DataFolder['Py_StaSpecAvg'] = ATaCR_Py_DataFolder['Py_DataParentFolder'] + '/AVG_STA'
        ATaCR_Py_DataFolder['Py_CorrectedTraces'] = ATaCR_Py_DataFolder['Py_DataParentFolder'] + '/EVENTS'
        ATaCR_Py_DataFolder['Py_b1b2_StaSpectra'] = ATaCR_Py_DataFolder['Py_DataParentFolder'] + '/SPECTRA'
        ATaCR_Py_DataFolder['Py_TransferFunctions'] = ATaCR_Py_DataFolder['Py_DataParentFolder'] + '/TF_STA'
        ATaCR_Py_DataFolder['Py_Logs'] = ATaCR_Py_DataFolder['Py_DataParentFolder'] + '/Logs'
        return ATaCR_ML_DataFolder,ATaCR_Py_DataFolder
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def mat2df(files):
        '''
        Converts all loaded matlab variables into a single dataframe
        '''
        if not isinstance(files,list):
                files = [files]
        files = sorted(files)
        out = []
        df = pd.DataFrame()
        FNUM = 0
        for f in files:
                df2 = pd.DataFrame()
                cfold = f[0:maxstrfind(f,'/')+1]
                cf = f[maxstrfind(f,'/')+1:len(f)]
                mat = spio.loadmat(f, simplify_cells=True)
                matvars = list(mat.keys())[3:len(list(mat.keys()))]
                oldkey = list(mat.keys())[-1] #Replaces Last Key with a single hardcoded name. Assumes a single variable (ie one struct) was saved to the matfile
                for k in range(3): #pop out the first three keys. they are artifacts from the mat import
                        mat.pop(list(mat.keys())[0])
                for k in matvars:
                        if isinstance(mat[k],list):
                                for dct in mat[k]:
                                        tmp = pd.DataFrame.from_dict(dct,orient='index').T
                                        tmp['FNUM'] = FNUM
                                        df2 = pd.concat([df2,tmp])
                                mat[k] = df2
                        else:
                                mat[k] = pd.DataFrame.from_dict(mat[k]).T
                                mat[k]['FNUM'] = FNUM
                if len(matvars)>1:
                        for i in range(len(matvars)-1):
                                mat[matvars[0]] = pd.merge(mat[matvars[0]],mat[matvars[i+1]],on='FNUM')
                fdf = mat[matvars[0]].set_index('FNUM')
                fdf['Folder'] = cfold
                fdf['File'] = cf
                out.append(fdf)
                FNUM += 1
        for d in out:
                df = pd.concat([df,d])
        df = df.groupby('FNUM',as_index=False,dropna=False).ffill().groupby('FNUM',as_index=False,dropna=False).bfill() #<---Fills data gaps using adjacent rows from only the same file.
        df = df.sort_values('File')
        return df
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def maxstrfind(s,p):
        '''
        Wild this doesnt exist in Python libraries. Finds LAST occurence of expression in a string.
        '''
        i = 0
        while i>-1:
                io = i
                i = s.find(p,i+1,len(s))
        return io
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def datenum_to_datetime64(dnum):
        '''
        Just some Matlab DateNum nonsense
        '''
        days = np.asarray(dnum) - 719529  # shift to unix epoch (1970-01-01)
        return np.round((days * 86400000)).astype("datetime64[ms]")
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def get_Noise(atacr_parent,net,sta,avg='Day',update=True):
        atacr_parent = Path(atacr_parent)
        if avg.lower()=='sta':
                specfold = 'AVG_STA'
                filetag = 'avg_sta'
        elif avg.lower()=='day':
                specfold = 'SPECTRA'
                filetag = 'spectra'
        subfolder = atacr_parent / specfold / '.'.join([net,sta])
        if update==False:
                out = pd.read_pickle(subfolder / 'Metrics' /  ('.'.join([net,sta])+ '.' + avg + 'Metrics.pkl'))
        else:
                # print(subfolder)
                s = loadpickles(subfolder)
                out = []
                for o,f in zip(s.Output,s.File):
                        keys = list(o.__dict__.keys())
                        keys.remove('cross')
                        keys.remove('power')
                        if avg.lower()=='day':
                                keys.remove('ft1')
                                keys.remove('ft2')
                                keys.remove('ftZ')
                                keys.remove('ftP')
                        vals = [[o.__dict__[k]] for k in keys]
                        noise = dict()
                        noise['StaNoise'] = [o]
                        noise['File'] = f
                        for k,v in zip(keys,vals):
                                noise[k] = v
                        if avg.lower()=='sta':
                                year = [int(s.split('.')[0]) for s in f.split('.avg_sta.pkl')[0].split('-')]
                                jday = [int(s.split('.')[1]) for s in f.split('.avg_sta.pkl')[0].split('-')]
                                noise['year'] = [year]
                                noise['julday'] = [jday]
                        if avg.lower()=='day':
                                ft = dict()
                                ft['1'] = o.__dict__['ft1']
                                ft['2'] = o.__dict__['ft2']
                                ft['Z'] = o.__dict__['ftZ']
                                ft['P'] = o.__dict__['ftP']
                                tr = [np.mean(np.real(np.fft.ifft(ft[k][o.__dict__['goodwins']])),axis=0) for k in list(ft.keys())]
                                GoodNoise = dict()
                                for k,t in zip(list(ft.keys()),tr):
                                        GoodNoise[k] = t
                                noise['ft'] = [ft]
                                noise['GoodNoise'] = [GoodNoise]
                        noise_out = pd.DataFrame.from_dict(noise)
                        if avg.lower()=='sta':
                                spec = o.power.__dict__
                                spec.update(o.cross.__dict__)
                                spec['f'] = o.__dict__['f']
                                noise['CSD'] = [spec]
                                noise_trace = dict()
                                noise_trace_avg = dict()
                                gd = noise_out.gooddays[0]
                                gw = ObsQA.TOOLS.io.get_Noise(atacr_parent,net,sta,avg='Day',update=True).goodwins[gd]
                                ft = ObsQA.TOOLS.io.get_Noise(atacr_parent,net,sta,avg='Day',update=True).ft[gd]

                                for c in ['1','2','Z','P']:
                                        noise_trace[c] = [np.mean(spec[c][g],axis=0) for spec,g in zip(ft.iloc,gw.iloc)]
                                        noise_trace_avg[c] = np.real(np.fft.ifft(np.mean(np.array(noise_trace[c]),axis=0)))
                                noise['Noise_Average'] = [noise_trace_avg]
                                noise['Noise'] = [noise_trace]
                                noise_out = pd.DataFrame.from_dict(noise)
                        if len(out)==0:
                                out = noise_out
                        else:
                                out = pd.concat([out,noise_out])
                if avg.lower()=='day':
                        out = out[['StaNoise','GoodNoise','File','key','tkey','julday','year','QC','av',
                                        'window','overlap','dt','npts','fs','ncomp','tf_list',
                                        'goodwins','rotation','f','ft']]
                        out = out.sort_values(by=['year','julday'],ascending=False)
                else:
                        out = out[['key','StaNoise','Noise','Noise_Average','File', 'year', 'julday','f', 'nwins']]
        out.reset_index(drop=True,inplace=True)
        if avg.lower()=='sta':
                Metrics = dict()
                Metrics_Noise = ObsQA.OBSM.classes.OBSMetrics(csd=out.StaNoise[0],f=out.StaNoise[0].f)
                # Metrics['Noise'] = out
                Metrics['Noise'] = Metrics_Noise
                out = Metrics
        return out
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def ClosestMLPreEventTF(tfdir,event_time):
        '''
        Replicates the ML ATaCR method for defining the 'nearest' TF to an event. It is a slightly different method than the Python version.
        '''
        files = g.glob(tfdir + '/*.mat')
        files = [ele for ele in files if '_AVERAGE_' not in ele] #<--Remove the station average from file list
        eventids = np.array(list(map(int,[f.split('/')[-1].split('_')[1] for f in files]))) #<Split file strings down to their event ids. Assumes format dir/stuff_EVENTID_stuff
        f = files[np.where(((np.array(event_time,dtype=int) - eventids)>0) & ((np.array(event_time,dtype=int) - eventids)==np.min((np.array(event_time,dtype=int) - eventids))))[0][0]]
        out = mat2df(f)
        return out
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def organize_evdata(evdata):
        '''
        Builds ObsPy trace objects inside the dataframes from the imported Matlab data
        '''
        itr1 = (evdata['channel'].squeeze().str.find('1')>0)
        itr2 = (evdata['channel'].squeeze().str.find('2')>0)
        itrZ = (evdata['channel'].squeeze().str.find('Z')>0)
        itrP = ~((evdata['channel'].squeeze().str.find('1')>0) + (evdata['channel'].squeeze().str.find('2')>0) + (evdata['channel'].squeeze().str.find('Z')>0))
        tmp = evdata[itr1]
        tr1 = Trace(data=tmp['data'].squeeze())
        tr1.stats.sampling_rate = tmp.sampleRate.squeeze()
        tr1.stats.network = tmp.network.squeeze()
        tr1.stats.station = tmp.station.squeeze()
        tr1.stats.channel = tmp.channel.squeeze()
        tr1.stats.location = 'ML-ATaCR-PreProc'
        tr1.stats.starttime = UTCDateTime( datenum_to_datetime64(tmp.startTime.squeeze()).tolist().strftime('%Y-%m-%dT%H:%M:%S.%f') )
        tmp = evdata[itr2]
        tr2 = Trace(data=tmp['data'].squeeze())
        tr2.stats.sampling_rate = tmp.sampleRate.squeeze()
        tr2.stats.network = tmp.network.squeeze()
        tr2.stats.station = tmp.station.squeeze()
        tr2.stats.channel = tmp.channel.squeeze()
        tr2.stats.location = 'ML-ATaCR-PreProc'
        tr2.stats.starttime = UTCDateTime( datenum_to_datetime64(tmp.startTime.squeeze()).tolist().strftime('%Y-%m-%dT%H:%M:%S.%f') )
        tmp = evdata[itrZ]
        trZ = Trace(data=tmp['data'].squeeze())
        trZ.stats.sampling_rate = tmp.sampleRate.squeeze()
        trZ.stats.network = tmp.network.squeeze()
        trZ.stats.station = tmp.station.squeeze()
        trZ.stats.channel = tmp.channel.squeeze()
        trZ.stats.location = 'ML-ATaCR-PreProc'
        trZ.stats.starttime = UTCDateTime( datenum_to_datetime64(tmp.startTime.squeeze()).tolist().strftime('%Y-%m-%dT%H:%M:%S.%f') )
        tmp = evdata[itrP]
        trP = Trace(data=tmp['data'].squeeze())
        trP.stats.sampling_rate = tmp.sampleRate.squeeze()
        trP.stats.network = tmp.network.squeeze()
        trP.stats.station = tmp.station.squeeze()
        trP.stats.channel = tmp.channel.squeeze()
        trP.stats.location = 'ML-ATaCR-PreProc'
        trP.stats.starttime = UTCDateTime( datenum_to_datetime64(tmp.startTime.squeeze()).tolist().strftime('%Y-%m-%dT%H:%M:%S.%f') )
        evdata = evdata.assign(Trace=0)
        evdata.iat[np.squeeze(np.where(itr1)).tolist(), -1] = tr1
        evdata.iat[np.squeeze(np.where(itr2)).tolist(), -1] = tr2
        evdata.iat[np.squeeze(np.where(itrZ)).tolist(), -1] = trZ
        evdata.iat[np.squeeze(np.where(itrP)).tolist(), -1] = trP
        evdata.tr1 = tr1
        evdata.tr2 = tr2
        evdata.trZ = trZ
        evdata.trP = trP
        return evdata
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def demean_detrend(stream):
        stream.detrend()
        for i in range(len(stream)):
            stream[i].data = stream[i].data - np.mean(stream[i].data)  
        stream.detrend()
        return stream
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def get_noisecut_output(net,sta,event,noisecut_datafolder):
  noisecut_datafolder = Path(noisecut_datafolder)
  file = '.'.join([net,sta,event,'pkl'])
  stafolder = '.'.join([net,sta])
  if len(list((noisecut_datafolder / stafolder).glob(file)))==0:
    print('{} not found'.format(file))
    return []
  else:
        st = pd.read_pickle(noisecut_datafolder / stafolder / file)
        for i,_ in enumerate(st.Raw[0]):
                st.Raw[0][i].stats.location = 'Raw'
        return st
  
def get_noisecut_event(parentfolder,staname,event,channel=['*1','*2','*Z','*H'],len_hrs=2,pre_trim=False,post_trim=True,win_length=200,width=None):
        # event_dt = 22  #<--This isn't necessary
        dateformat = '%Y.%j.%H.%M'
        event = UTCDateTime.strptime(str(event),dateformat)
        origin = event
        channel = [cc.replace('**','*') for cc in ['*' + c for c in np.array([channel]).flatten().tolist()]]
        folder = Path(str(parentfolder)) / 'Data' / staname / 'HPS_Data'  
        files = []
        # [files.extend(list(folder.glob((event-(3600*event_dt*i)).strftime('%Y.%j') + '*.SAC'))) for i in range(2)]
        [files.extend(list(folder.glob(event.strftime('%Y.%j') + '*.SAC'))) for i in range(1)]
        if len(files)==0:
                print('||--24hr event data not found')
                return []
        raw = Stream([read(f)[0] for f in files])
        raw = Stream([raw.select(channel=c)[0].copy() for c in channel])
        if pre_trim:
                raw.trim(origin,origin+3600*len_hrs)
        raw = demean_detrend(raw.copy())
        hps_out = [noisecut(s.copy(),ret_spectrograms=True,win_length=win_length,width=width) for s in raw]
        corrected = Stream([h[0].copy() for h in hps_out])

        if post_trim:
                corrected.trim(origin,origin+3600*len_hrs)
                raw.trim(origin,origin+3600*len_hrs)

        hps_spectrograms = [h[1] for h in hps_out]
        raw_collect = {'trZ':raw.select(channel='*Z')[0].copy(),'tr1':raw.select(channel='*1')[0].copy(),'tr2':raw.select(channel='*2')[0].copy(),'trP':raw.select(channel='*H')[0].copy()}
        corrected_collect = {'trZ':corrected.select(channel='*Z')[0].copy(),'tr1':corrected.select(channel='*1')[0].copy(),'tr2':corrected.select(channel='*2')[0].copy(),'trP':corrected.select(channel='*H')[0].copy()}
        out = pd.DataFrame({'Event':event,'Network':staname.split('.')[0],'Station':staname.split('.')[1],'Raw':[raw_collect],'Corrected':[corrected_collect],'Spectrograms':[hps_spectrograms]})
        return out

def run_noisecut(tr,win_length=200,width=None):
        hps, (S_full, S_background, S_hps, frequencies, times) = noisecut(tr, ret_spectrograms=True,win_length=win_length,width=width)
        return hps, (S_full, S_background, S_hps, frequencies, times)

def Get_ATaCR_CorrectedEvents(eventfolder,eventnames,net,sta,tfavg='sta',tf=''):
        if not isinstance(eventnames,list):
            eventnames = [eventnames]
        if not isinstance(net,list):
            net = [net]
        if not isinstance(sta,list):
            sta = [sta]
        dobspy,raw_collect,corrected = [],[],[]
        for i in range(len(eventnames)):
            neti = net[i]
            stai = sta[i]
            evi = eventnames[i]
            prefix = neti + '.' + stai
            folder = eventfolder / prefix / 'CORRECTED'
            f = list(folder.glob(prefix + '.' + evi + '*.pkl'))
            if len(f)==0:
                return []
            f = f[0]
            ri = pkl.load(open(f,'rb'))
            raw = {'tr1':ri.tr1.copy(),'tr2':ri.tr2.copy(),'trZ':ri.trZ.copy(),'trP':ri.trP.copy()}
            for k in list(raw.keys()):
                    raw[k].stats.location = 'Raw' #For better book-keeping. Label the raw traces.
            trcorr = {}
            raw_collect.append(raw)
            # rawz.append(ri.trZ.copy())
            dobspy.append(ri)
            for k in list(ri.correct.keys()): #Shape corrected traces into a list of ObsPy trace objects
                    tr = ri.trZ.copy()
                    tr.data = ri.correct[k]
                    tr.stats.location = k
                    trcorr[k] = tr
            corrected.append(trcorr)
        out = pd.DataFrame({'Event':eventnames,'Network':net,'Station':sta,'Raw':raw_collect,'Corrected':corrected})
        return out
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def GetML_EventData_and_TransferFunctions(event_time,network,sta,preprocfolder=None,correctedfolder=None,tffolder=None):
        evdata = None
        corrected_evdata = None
        ml_tf = None
        if preprocfolder is not None:
                folder = preprocfolder
                path = folder + '/' + event_time + '/*.mat'
                ml_files = g.glob(path)
                evdata = ObsQA.TOOLS.io.mat2df(ml_files)
                evdata = evdata[(evdata['Network']==network)&(evdata['Station']==sta)]
                evdata = ObsQA.TOOLS.io.organize_evdata(evdata)
        if correctedfolder is not None:
                folder = correctedfolder
                path = folder + '/' + network + '/' + sta + '/' + network + sta + '_' + str(event_time) + '_corrseis' + '.mat'
                ml_files = g.glob(path)
                corrected_evdata = ObsQA.TOOLS.io.mat2df(ml_files)
                corrected_evdata = corrected_evdata[(corrected_evdata['Network']==network)&(corrected_evdata['Station']==sta)]
                corrected_evdata = corrected_evdata.assign(Trace=0)
                for i in range(len(corrected_evdata)):
                        tmp = corrected_evdata.iloc[i]
                        tr = Trace(data=tmp['timeseries'].squeeze())
                        tr.stats.sampling_rate = 1/tmp['dt']
                        tr.stats.network = tmp['Network']
                        tr.stats.station = tmp['Station']
                        tr.stats.channel = tmp['label']
                        tr.stats.location = 'ML-ATaCR-Corrected'
                        if evdata is not None:
                                tr.stats.starttime = evdata.trZ.stats.starttime
                        corrected_evdata.iat[i, -1] = tr
        if tffolder is not None:
                ftf = tffolder + '/' + network + '/' + sta
                ml_tf = ObsQA.TOOLS.io.ClosestMLPreEventTF(ftf,event_time)
        ml_tf = ml_tf.reset_index()
        corrected_evdata = corrected_evdata.reset_index()
        evdata = evdata.reset_index()
        TF_list = {i : j for i, j in zip(corrected_evdata.label.to_list(), np.ones(len(corrected_evdata.label.to_list()),dtype=bool))}
        corrected_evdata['TF_list'] = TF_list
        for i in range(len(corrected_evdata)):
                corrected_evdata.at[i,'TF_list'] = TF_list
        corrected_evdata = corrected_evdata[['label','Network','Station','eventid','TFfilename', 'dt','NFFT', 'f', 'filop',
        'taxis', 'tf_op', 'spectrum', 'timeseries','isgood', 'Folder', 'File', 'TF_list','Trace']]
        return evdata,corrected_evdata,ml_tf
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def condensedata(input):
        collect = pd.DataFrame()
        for data in input:
                if isinstance(data,pd.DataFrame):
                        key = np.append(np.array(data.keys()[np.array(np.where(data.keys()=='channel')).flatten().tolist()]),np.array(data.keys()[np.array(np.where(data.keys()=='label')).flatten().tolist()]))[0]
                        labels = list(data[key])
                        if len(np.flatnonzero(np.core.defchararray.find(labels,'Z')!=-1))==1:
                                # isPreProc
                                data = data.rename(columns={'channel': 'label'})[['Trace','Network','Station','label','File']].copy().iloc[np.where(data[key]=='BHZ')]
                                data['State'] = 'PreProc'
                                data['CodeBase'] = 'Matlab'
                        else:
                                # isCorrected
                                data = data[['Trace','Network','Station','label','File']].copy()
                                data['State'] = 'Corrected'
                                data['CodeBase'] = 'Matlab'
                elif isinstance(data,obstools.atacr.classes.EventStream):
                        d1 = pd.DataFrame.from_dict({'Trace':data.trZ,'Network':data.trZ.stats.network,'Station':data.trZ.stats.station,'label':data.trZ.stats.channel,'File':data.prefix,'State':'PreProc','CodeBase':'Python'},orient='index').T
                        d2 = pd.DataFrame.from_dict({'Trace':list(evstream.correct.values()),'Network':data.trZ.stats.network,'Station':data.trZ.stats.station,'label':list(data.correct.keys()),'File':data.prefix,'State':'Corrected','CodeBase':'Python'})
                        data = pd.concat([d1,d2])
                collect = pd.concat([collect,data]).sort_values(by=['CodeBase'])
        return collect
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def MLtoDayNoise(folder,netlist,stalist):
        Data = pd.DataFrame({'Network':netlist,'Station':stalist,'DayNoise':[ [] for _ in range(len(stalist)) ],'NDays':np.zeros(len(stalist)),'Files':[ [] for _ in range(len(stalist)) ]})
        lls = []
        llf = []
        for ista in range(len(stalist)):
                csta=stalist[ista]
                cnet=netlist[ista]
                cfold = folder + '/' + cnet + '/' + csta
                df = mat2df(g.glob(cfold + '/*.mat'))
                files = list(df.File.unique())
                DayList = []
        for i in range(len(files)):
                cdf = organize_evdata(df.iloc[list(df.File==files[i])])
                window = cdf.iloc[0].Trace.stats.endtime - cdf.iloc[0].Trace.stats.starttime
                overlap = 0.3 #Default set in ML's ATaCR code (setup_parameter)
                key = cdf.iloc[0].Trace.stats.network + '.' + cdf.iloc[0].Trace.stats.station
                tr1 = list(cdf[list(cdf.channel=='BH1')].Trace)[0]
                tr2 = list(cdf[list(cdf.channel=='BH2')].Trace)[0]
                trZ = list(cdf[list(cdf.channel=='BHZ')].Trace)[0]
                trP = list(cdf[list(cdf.channel=='BDH')].Trace)[0]
                DN = DayNoise(tr1=tr1, tr2=tr2, trZ=trZ, trP=trP, window=window, key=key)
                DN.QC = True
                DN.av = True
                DayList.append(DN)
        Data.loc[ista,'NDays'] = len(DayList)
        lls.append(DayList)
        llf.append(files)
        Data['DayNoise'] = lls
        Data['Files'] = llf
        return Data
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def MLtoStaNoise(cnet,csta,preprocfolder,specfolder,level='b1'):
        DN = ObsQA.TOOLS.io.MLtoDayNoise(preprocfolder,[cnet],[csta])
        cfold = specfolder + '/' + cnet + '/' + csta + '/' + level.lower()
        files = g.glob(cfold + '/*.mat')
        df = ObsQA.TOOLS.io.mat2df(files)
        files = df.File.unique()
        udf = df.copy().iloc[0].to_frame().T.drop(0)
        for i in range(len(files)):
                udf = pd.concat([udf,df.iloc[list(df.File==files[i])].iloc[0].to_frame().T])
        df = udf.copy()
        for ifile in range(len(files)):
                cfile = files[ifile]
                cdf = df.iloc[list(df.File==cfile)]
                c11 = cdf.iloc[0].c11_stack.T
                c22 = cdf.iloc[0].c22_stack.T
                cZZ = cdf.iloc[0].czz_stack.T
                cPP = cdf.iloc[0].cpp_stack.T
                c12 = cdf.iloc[0].c12_stack.T
                c1Z = cdf.iloc[0].c1z_stack.T
                c1P = cdf.iloc[0].c1p_stack.T
                c2Z = cdf.iloc[0].c2z_stack.T
                c2P = cdf.iloc[0].c2p_stack.T
                cZP = cdf.iloc[0].cpz_stack.T
                cHH = cdf.iloc[0].chh_stack.T
                cHZ = cdf.iloc[0].chz_stack.T
                cHP = cdf.iloc[0].chp_stack.T
                tilt = cdf.iloc[0].rotor
                rotcoh = cdf.iloc[0].rotcoh
                f = cdf.iloc[0].f
                nwins = len(cdf.iloc[0].goodwins)
                direc = np.arange(0,360+10,10) #default in ML ATaCR code (b1)
                ncomp = sum(cdf.iloc[0].comp_exist)
                key = cnet + '.' + csta
                tf_list = {'ZP': True, 'Z1': True, 'Z2-1': True, 'ZP-21': True, 'ZH': True, 'ZP-H': True}
                power = obstools.atacr.classes.Power(c11 = c11, c22 = c22, cZZ = cZZ, cPP = cPP)
                cross = obstools.atacr.classes.Cross(c12=c12, c1Z=c1Z, c1P=c1P, c2Z=c2Z, c2P=c2P, cZP=cZP)
                rotation = obstools.atacr.classes.Rotation(cHH=cHH, cHZ=cHZ, cHP=cHP, coh=None, ph=None, tilt=tilt, coh_value=rotcoh, phase_value=None, direc=direc)
                goodwins = np.array(cdf.iloc[0].goodwins,dtype=bool)
                ind = np.where(np.char.replace(DN.iloc[0].Files,'.mat','_' + level + '_spectra.mat')==cfile)[0][0]
                DN.iloc[0].DayNoise[ind].power = power
                DN.iloc[0].DayNoise[ind].cross = cross
                DN.iloc[0].DayNoise[ind].rotation = rotation
                DN.iloc[0].DayNoise[ind].f = cdf.iloc[0].f
                DN.iloc[0].DayNoise[ind].goodwins = goodwins
        SN = StaNoise(daylist=DN.iloc[0].DayNoise)
        SN.init()
        SN.key = 'ML-ATaCR: ' + SN.key
        return SN
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def loadpickles(path):
        if not isinstance(path,Path):
                path = Path(path)
        py_files = list(path.glob('*.pkl'))
        if len(py_files)==0:
                raise Exception('Folder contains no .pkl files')
        out = pd.DataFrame({'Output':[ [] for _ in range(len(py_files)) ],'File':[ [] for _ in range(len(py_files)) ]})
        for i,f in enumerate(py_files):
                file = open(f, 'rb')
                pydata = pkl.load(file)
                out.iloc[i]['Output'] = pydata
                out.iloc[i]['File'] = f.name
                file.close()
        return out
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def randomdays(TStart,TEnd,seed='MESSI_22FIFA_WORLD_CUP',days=10,dateformat = '%Y-%m-%d %H:%M:%S'):
        TStart = datetime.datetime(year=TStart.year,month=TStart.month,day=TStart.day)
        TEnd = datetime.datetime(year=TEnd.year,month=TEnd.month,day=TEnd.day)
        Starts = [TStart + datetime.timedelta(days=j) for j in range((TEnd - TStart).days)]
        Ends = [TStart + datetime.timedelta(days=j+1) for j in range((TEnd - TStart).days)]
        if isinstance(seed,str):
                seed = int(''.join([str(ord(a)) for a in seed]))
        else:
                seed = int(str(seed).replace('.',''))
        rng = np.random.default_rng(seed=seed)
        idx = rng.choice([i for i in range(len(Starts))], size=days, replace=True)
        Starts = np.array(Starts)[idx].tolist()
        Ends = np.array(Ends)[idx].tolist()
        Starts = [str(UTCDateTime.strptime(str(a),dateformat)) for a in Starts]
        Ends = [str(UTCDateTime.strptime(str(a),dateformat)) for a in Ends]
        return Starts,Ends
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def build_staquery(d,chan='H',ATaCR_Parent=None,staquery_output='./sta_query.pkl'):
        out = dict()
        if not isinstance(staquery_output,Path):
                if ATaCR_Parent is None:
                        staquery_output = os.getcwd() + staquery_output.replace('./','/')
                else:
                        os.chdir(ATaCR_Parent)
                        staquery_output = ATaCR_Parent + '/' + staquery_output.replace('./','')
        for csta in d.iloc:
                net = csta.Network
                sta = csta.Station
                key = net + '.' + sta
                csta_dict = {'station':sta, 'network':net, 'altnet':[],'channel':chan, 'location':['--'], 'latitude':csta['Latitude'], 'longitude':csta['Longitude'], 'elevation':-csta['Water_Depth_m']/1000, 'startdate':UTCDateTime(csta.Start), 'enddate':UTCDateTime(csta.End), 'polarity':1.0, 'azcorr':0.0, 'status':'open'}
                out[key] = csta_dict
        output = pd.DataFrame.from_dict(out)
        if staquery_output is not None:
                output.to_pickle(staquery_output)
                print('Station Query File Written to: ' + staquery_output)
        else:
                return output
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def DownloadEvents(catalog,ATaCR_Parent=None,netsta_names=None,Minmag=6.3,Maxmag=6.7,limit=1000,pre_event_min_aperture=1,logoutput_subfolder='',log_prefix = '',staquery_output = './sta_query.pkl',chan='H'):
        # sys.stdout.flush()
        logfilename = '_Step_2_7_EventDownload_logfile.log'
        if logoutput_subfolder is not None:
                logoutput = logoutput_subfolder + '/' + log_prefix
        else:
                logoutput = '_Step_2_7_EventDownload_logfile.log'
        datafolder = './EVENTS/'
        dateformat = '%Y.%j.%H.%M'
        print('----Begin Event Download----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        for i in range(len(catalog)):
                csta = catalog.iloc[i]
                S = csta['Station']#station
                N = csta['Network']
                staname = str(N) + '.' + str(S)
                if os.path.isdir(datafolder + staname)==False:
                        os.system('mkdir ' + datafolder + '/' + staname)
                ObsQA.TOOLS.io.build_staquery(d=catalog[(catalog.Network==N) & (catalog.Station==S)],staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
                log_fout = logoutput  + logfilename
                # log_fout
                logoutput_subfolder = Path(logoutput_subfolder)
                logoutput_subfolder.mkdir(parents=True,exist_ok=True)
                # original = sys.stdout
                # sys.stdout = open(log_fout,'w+')
                logging.basicConfig(stream=log_fout)
                print('--' + staname + '--',flush=True)
                for j in range(len(csta.Events)):
                        ev = csta.Events[j]
                        print('---'*20)
                        print(staname + '| S:[' +str(i+1) + '/' + str(len(catalog)) + '] | E:[' + str(j+1) + '/' + str(len(csta.Events)) + '] | ' + str(ev),flush=True)
                        print('---'*20)
                        EventStart = UTCDateTime.strptime(ev,dateformat)
                        EventEnd = UTCDateTime.strptime(ev,dateformat) + datetime.timedelta(minutes=1)
                        print('Event ' + str(j+1) + '/' + str(len(csta.Events)) + ' | Start -> End : ' + str(EventStart) + ' -> ' + str(EventEnd))
                        if ATaCR_Parent is not None:
                                os.chdir(ATaCR_Parent)
                        # args = [staquery_output,'--start={}'.format(EventStart), '--end={}'.format(EventEnd),'--min-mag={}'.format(Minmag),'--max-mag={}'.format(Maxmag),'--limit={}'.format(limit)]
                        args = [staquery_output,'--start={}'.format(EventStart), '--end={}'.format(EventEnd),'--min-mag={}'.format(Minmag),'--max-mag={}'.format(Maxmag)]
                        # [print(ii) for ii in [staquery_output,'--start={}'.format(EventStart), '--end={}'.format(EventEnd),'--min-mag={}'.format(Minmag),'--max-mag={}'.format(Maxmag),'--limit={}'.format(limit)]]
                        # with open(log_fout, 'w') as sys.stdout:
                        atacr_download_event.main(atacr_download_event.get_event_arguments(args))
        print(' ')
        print('----Event Download Complete----')
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def DayNoiseWhileLoop(catalog,NoiseFolder,ATaCR_Parent,days=10,attempts=50,seed='MESSI_22FIFA_WORLD_CUP!'):
      NoiseFolder = Path(NoiseFolder)
      cat = catalog.copy()
      stafolders = [f for f in NoiseFolder.glob('*/') if f.is_dir()]
      stafolders = [stafolders[np.where(c==np.array([g.name for g in stafolders]))[0][0]] for c in cat.StaName.iloc]
      STEP = 3
      if isinstance(seed,str):
            seed = int(''.join([str(ord(a)) for a in seed]))
      else:
            seed = int(str(seed).replace('.',''))
      for ista,cstaf in enumerate(stafolders):
            print('--<>--'*20)
            nfiles = len([f for f in cstaf.glob('*Z.SAC')])
            nstart = nfiles
            print('Folder: [' + cstaf.name + '] (' + str(ista) + '/' + str(len(stafolders)) + ')'  + '  ||  Current day count: ' + str(nfiles))
            queries = 0
            if nfiles<days:
                  while (queries<attempts) and (nfiles<days):
                        queries+=1
                        seed*=2
                        query_day_len = days - nfiles
                        icatalog = cat[cat.Station==str(cstaf).split('/')[-1].split('.')[-1]]
                        print('Attempt: ' + str(queries) + ', Days needed: ' + str(query_day_len))
                        ObsQA.TOOLS.io.Run_ATaCR(icatalog, ATaCR_Parent = ATaCR_Parent,STEPS=[STEP],log_prefix=icatalog.StaName.iloc[0],days=query_day_len,seed=seed)
                        nfiles = len([f for f in cstaf.glob('*Z.SAC')])
                        if nfiles<days:
                                print(' || Day requirements not satisfied. Attempting new seed. ||')
                        else:
                                print('[' + str(ista) + '/' + str(len(stafolders)) + '][DAY REQUIREMENTS SATISFIED] | Days gained: ' + str(nfiles-nstart))
            else:
                    print('[' + str(ista) + '/' + str(len(stafolders)) + '][DAY REQUIREMENTS ALREADY SATISFIED] | Skipping.')
            if (nfiles<days) & (queries>=attempts):
                    print('Attempt ' + str(attempts) + '/' + str(attempts) + '. Maximum number of attempts exceeded. Skipping station.')
      failedattempts = [(cstaf.parts[-1], len([f for f in cstaf.glob('*Z.SAC')])) for cstaf in stafolders if len([f for f in cstaf.glob('*Z.SAC')])<days]        
      print('Complete')
      if len(failedattempts)>0:
            print('Stations that still do not meet requirements:')
            print(failedattempts)
      else:
            print('All station folders now meet day requirements')
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def DownloadDayNoise(catalog,days=[],end_delta=24,event_mode=False,seed='MESSI_22FIFA_WORLD_CUP!',ATaCR_Parent=None,netsta_names=None,logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H'):
        logfilename = log_prefix + '_Step_3_7_NoiseDownload_logfile.log'
        dateformat = '%Y.%j.%H.%M'
        event_dt = 22
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        datafolder = Path('./Data/')
        if logoutput_subfolder is not None:
                logfolder = Path(logoutput_subfolder)
                logfolder.mkdir(exist_ok=True,parents=True)
        else:
                logfolder = datafolder
        print('----Begin Noise Download----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        for i,Station in enumerate(catalog.iloc):
                if event_mode:
                        Origins = Station.Origin
                        Starts = Origins
                        if isinstance(Origins[0],obspy.core.event.origin.Origin):
                                Starts = [e.time for e in Origins]
                        Ends = [e+(3600*24) for e in Starts]
                elif isinstance(days,int):
                        Starts,Ends = ObsQA.TOOLS.io.randomdays(Station.Start,Station.End,seed=seed,days=days)
                else:
                        Starts = days
                        Ends = [datetime.datetime(year=s.year,month=s.month,day=s.day) + datetime.timedelta(hours=end_delta) for s in Starts]
                # --Log file stdout stuff--
                staname = Station.StaName
                sta_folder = Path(datafolder / staname)
                sta_folder.mkdir(exist_ok=True,parents=True)
                if logoutput_subfolder is not None:
                        log_fout = logfolder / logfilename
                else:
                        log_fout = sta_folder / logfilename
                # original = sys.stdout
                # sys.stdout = open(log_fout,'w+')
                logging.basicConfig(stream=log_fout)
                print('--' + staname + '--',flush=True)

                ObsQA.TOOLS.io.build_staquery(d=Station.to_frame().T,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
                for j,(NoiseStart,NoiseEnd) in enumerate(zip(Starts,Ends)):
                        print('<||>'*30)
                        print(staname + ' Station ' +str(i+1) + '/' + str(len(catalog)) + ' - Day ' + str(j+1) + '/' + str(len(Starts)),flush=True)
                        args = [staquery_output,'--start={}'.format(NoiseStart), '--end={}'.format(NoiseEnd)]
                        atacr_download_data.main(atacr_download_data.get_daylong_arguments(args))
        print(' ')
        print('----Noise Download Complete----')
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def DownloadEVNoise(catalog,ATaCR_Parent=None,netsta_names=None,pre_event_day_aperture=10,logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H'):
        # sys.stdout.flush()
        logfilename = log_prefix + '_Step_3_7_NoiseDownload_logfile.log'
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        datafolder = Path('./Data/')
        if logoutput_subfolder is not None:
                logfolder = Path(logoutput_subfolder)
                logfolder.mkdir(exist_ok=True,parents=True)
        else:
                logfolder = datafolder
        dateformat = '%Y.%j.%H.%M'

        print('----Begin Noise Download----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        for i in range(len(catalog)):
                csta = catalog.iloc[i]
                S = csta['Station']#station
                N = csta['Network']
                staname = str(N) + '.' + str(S)
                sta_folder = Path(datafolder / staname)
                sta_folder.mkdir(exist_ok=True,parents=True)
                ObsQA.TOOLS.io.build_staquery(d=catalog[(catalog.Network==N) & (catalog.Station==S)],staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
                if logoutput_subfolder is not None:
                        log_fout = logfolder / logfilename
                else:
                        log_fout = sta_folder / logfilename
                # original = sys.stdout
                # sys.stdout = open(log_fout,'w+')
                logging.basicConfig(stream=log_fout)
                print('--' + staname + '--',flush=True)
                for j in range(len(csta.Events)):
                        print(staname + ' Station ' +str(i+1) + '/' + str(len(catalog)) + ' - Event ' + str(j+1) + '/' + str(len(csta.Events)),flush=True)
                        ev = csta.Events[j]
                        NoiseStart = UTCDateTime.strptime(ev,dateformat) - datetime.timedelta(days=pre_event_day_aperture)
                        NoiseStart = NoiseStart - datetime.timedelta(hours = NoiseStart.hour, minutes = NoiseStart.minute, seconds=NoiseStart.second) #rounds down to the nearest day
                        NoiseEnd = UTCDateTime.strptime(ev,dateformat)
                        NoiseEnd = NoiseEnd - datetime.timedelta(hours = NoiseEnd.hour, minutes = NoiseEnd.minute, seconds=NoiseEnd.second) #rounds down to the nearest day
                        args = [staquery_output,'--start={}'.format(NoiseStart), '--end={}'.format(NoiseEnd)]
                        atacr_download_data.main(atacr_download_data.get_daylong_arguments(args))
        print(' ')
        print('----Noise Download Complete----')
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def DailySpectra(catalog,ATaCR_Parent=None,netsta_names=None,extra_flags = '-O --figQC --figAverage --figCoh --save-fig',logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H',fork=True,max_workers=1):
        # sys.stdout.flush()
        datafolder = './Data/'
        if logoutput_subfolder is not None:
                logoutput = logoutput_subfolder + '/' + log_prefix + '_Step_4_7_QCSpectra_logfile.log'
                log_fout = logoutput
        else:
                logoutput = '_Step_4_7_QCSpectra_logfile.log'
                log_fout = datafolder + logoutput
        SpecStart = catalog.Start.min().strftime("%Y-%m-%d, %H:%M:%S")
        SpecEnd = catalog.End.max().strftime("%Y-%m-%d, %H:%M:%S")
        args = [staquery_output]
        [args.append(flg) for flg in extra_flags.split()]
        [args.append(flg) for flg in ['--start={}'.format(SpecStart),'--end={}'.format(SpecEnd)]]
        # with open(log_fout, 'w') as sys.stdout:
        # original = sys.stdout
        # sys.stdout = open(log_fout,'w+')
        logging.basicConfig(stream=log_fout)
        # print('----Begin Daily Spectra----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.TOOLS.io.build_staquery(d=catalog,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        args = atacr_daily_spectra.get_dailyspec_arguments(args)
        if not fork:
                atacr_daily_spectra.main(args)
        else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as process_executor:
                        print('----Begin Daily Spectra----')
                        future = process_executor.submit(atacr_daily_spectra.main,args)
                        print('Waiting for tasks to complete')
                        wait([future])
                future.result()
                print('----Daily Spectra Complete----')
                # sys.stdout = original
        # atacr_daily_spectra.main(args)
        # print(' ')
        # print('----Daily Spectra Complete----')
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def Fork(fn, *args, **kwargs):
  """Submits a job to the concurrent.futures ThreadPoolExecutor.

  Args:
    fn: The function to be executed.
    *args: The arguments to be passed to the function.
    **kwargs: The keyword arguments to be passed to the function.

  Returns:
    A Future object representing the job.
  """
  max_workers = kwargs['max_workers']
  with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
    future = executor.submit(fn, *args, **kwargs)
    print('Waiting for tasks to complete')
    wait([future])
    return future


def CleanSpectra(catalog,ATaCR_Parent=None,netsta_names=None,extra_flags = '-O --figQC --figAverage --save-fig',logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H',fork=True,max_workers=1):
        # sys.stdout.flush()
        #  --figQC --figAverage --save-fig
        datafolder = './Data/'
        if logoutput_subfolder is not None:
                logoutput = logoutput_subfolder + '/' + log_prefix + '_Step_5_7_CleanSpectra_logfile.log'
                log_fout = logoutput
        else:
                logoutput = '_Step_5_7_CleanSpectra_logfile.log'
                log_fout = datafolder + logoutput
        SpecStart = catalog.Start.min().strftime("%Y-%m-%d, %H:%M:%S")
        SpecEnd = catalog.End.max().strftime("%Y-%m-%d, %H:%M:%S")
        args = [staquery_output]
        [args.append(flg) for flg in extra_flags.split()]
        [args.append(flg) for flg in ['--start={}'.format(SpecStart),'--end={}'.format(SpecEnd)]]
        # with open(log_fout, 'w') as sys.stdout:
        # original = sys.stdout
        log_fout = Path(log_fout)
        log_fout.parent.mkdir(parents=True,exist_ok=True)
        # sys.stdout = open(str(log_fout),'w+')
        logging.basicConfig(stream=log_fout)
        # print('----Begin Clean Spectra----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.TOOLS.io.build_staquery(d=catalog,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        args = atacr_clean_spectra.get_cleanspec_arguments(args)
        # atacr_clean_spectra.main(args)
        if not fork:
                atacr_clean_spectra.main(args)
        else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as process_executor:
                        # Start the load operations and mark each future with its URL
                        print('----Begin Clean Spectra----')
                        future = process_executor.submit(atacr_clean_spectra.main,args)
                        print('Waiting for tasks to complete')
                        wait([future])
                future.result()
                print('----Clean Spectra Complete----')
                # sys.stdout = original
        # print('----Clean Spectra Complete----')
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def TransferFunctions(catalog,ATaCR_Parent=None,netsta_names=None,extra_flags = '-O --figTF --save-fig'
                #       ,taper_mode=1
                      ,logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H',fork=True,max_workers=1):
        # sys.stdout.flush()
        datafolder = './Data/'
        if logoutput_subfolder is not None:
                logoutput = logoutput_subfolder + '/' + log_prefix + '_Step_6_7_CalcTFs_logfile.log'
                log_fout = logoutput
        else:
                logoutput = '_Step_6_7_CalcTFs_logfile.log'
                log_fout = datafolder + logoutput
        args = [staquery_output]
        extra_flags = extra_flags # + ' --taper ' + str(taper_mode)
        [args.append(flg) for flg in extra_flags.split()]
        # with open(log_fout, 'w') as sys.stdout:
        # original = sys.stdout
        # sys.stdout = open(log_fout,'w+')
        logging.basicConfig(stream=log_fout)
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.TOOLS.io.build_staquery(d=catalog,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        args = atacr_transfer_functions.get_transfer_arguments(args)
        # atacr_transfer_functions.main(args)
        if not fork:
                atacr_transfer_functions.main(args)
        else:
                print('----Begin Transfer Functions----')
                with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as process_executor:
                        future = process_executor.submit(atacr_transfer_functions.main,args)
                        print('Waiting for tasks to complete')
                        wait([future])
                future.result()
                print('----Transfer Functions Complete----')
                # sys.stdout = original
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def CorrectEvents(catalog,ATaCR_Parent=None,netsta_names=None,extra_flags = '--figRaw --figClean --save-fig',logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H',fork=True,max_workers=1):
        datafolder = './Data/'
        staquery_output = str(staquery_output)
        if logoutput_subfolder is not None:
                logoutput = logoutput_subfolder + '/' + log_prefix + '_Step_7_7_CorrectEvents_logfile.log'
                log_fout = logoutput
        else:
                logoutput = '_Step_7_7_CorrectEvents_logfile.log'
                log_fout = datafolder + logoutput
        args = [staquery_output]
        [args.append(flg) for flg in extra_flags.split()]
        # with open(log_fout, 'w') as sys.stdout:
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.TOOLS.io.build_staquery(d=catalog,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
        # original = sys.stdout
        # sys.stdout = open(log_fout,'w+')
        logging.basicConfig(stream=log_fout)
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        args = atacr_correct_event.get_correct_arguments(args)
        if not fork:
                atacr_correct_event.main(args)
        else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as process_executor:
                        future = process_executor.submit(atacr_correct_event.main,args)
                        print('Waiting for tasks to complete')
                        wait([future])
                future.result()
                print('----Correct Events Complete----')
                # sys.stdout = original
        # atacr_correct_event.main(args)
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def Run_ATaCR(catalog, ATaCR_Parent = None, STEPS=[1,2,3,4,5,6,7], netsta_names=None, chan='H', Minmag=6.3, Maxmag=6.7, limit=1000, pre_event_min_aperture=1, pre_event_day_aperture=30, dailyspectra_flags='-O --figQC --figAverage --figCoh --save-fig', cleanspectra_flags='-O --figQC --figAverage --figCoh --figCross --save-fig', tf_flags='-O --figTF --save-fig', correctevents_flags='--figRaw --figClean --save-fig',logoutput_subfolder=None,log_prefix = '',staquery_output = './sta_query.pkl',fork=True,days=10,seed='MESSI_22FIFA_WORLD_CUP!',max_workers=1,event_mode=False,event_dt=None,
        #       taper_mode=0
              ):
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
        dirs = ObsQA.TOOLS.io.dir_libraries('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR')[1]
        if 1 in STEPS:
                print('Step 1/7 - BEGIN: Station Metadata')
                # C='?H?' #channels
                # !query_fdsn_stdb -N {','.join(N)} -C '{C}' -S {','.join(S)} ./Data/sta_query> ./Data/Step_1_7_StationMeta_logfile.log
                ObsQA.TOOLS.io.build_staquery(d=catalog,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
                print('Step 1/7 - COMPLETE: Station Metadata')
        if 2 in STEPS:
                print('Step 2/7 - BEGIN: Download Event Data')
                if logoutput_subfolder is None:
                        logoutput_subfolder = dirs['Py_Logs'] + '/2_7'
                ObsQA.TOOLS.io.DownloadEvents(catalog,ATaCR_Parent=ATaCR_Parent,netsta_names=netsta_names,Minmag=Minmag,Maxmag=Maxmag,limit=limit,pre_event_min_aperture=pre_event_min_aperture,logoutput_subfolder=logoutput_subfolder,staquery_output=staquery_output,chan=chan,log_prefix=log_prefix)
                print('Step 2/7 - COMPLETE: Download Event Data')
        if 3 in STEPS:
                print('Step 3/7 - BEGIN: Download Day Data')
                if logoutput_subfolder is None:
                        logoutput_subfolder = dirs['Py_Logs'] + '/3_7'
                # ObsQA.TOOLS.io.DownloadNoise(catalog,ATaCR_Parent=ATaCR_Parent,netsta_names=netsta_names,pre_event_day_aperture=pre_event_day_aperture,logoutput_subfolder=logoutput_subfolder,staquery_output=staquery_output,chan=chan,log_prefix=log_prefix)
                # ObsQA.TOOLS.io.DayNoiseWhileLoop(catalog,NoiseFolder,ATaCR_Parent,days=10,attempts=50)
                ObsQA.TOOLS.io.DownloadDayNoise(catalog,days=days,seed=seed,ATaCR_Parent=ATaCR_Parent,netsta_names=netsta_names,logoutput_subfolder=logoutput_subfolder,log_prefix = log_prefix,staquery_output=staquery_output,chan=chan,event_mode=event_mode)
                print('Step 3/7 - COMPLETE: Download Day Data')
        if 4 in STEPS:
                print('Step 4/7 - BEGIN: Quality Control Noise Data')
                if logoutput_subfolder is None:
                        logoutput_subfolder = dirs['Py_Logs'] + '/4_7'
                ObsQA.TOOLS.io.DailySpectra(catalog,ATaCR_Parent=ATaCR_Parent,netsta_names=netsta_names,extra_flags=dailyspectra_flags,logoutput_subfolder=logoutput_subfolder,staquery_output=staquery_output,chan=chan,log_prefix=log_prefix,fork=fork,max_workers=max_workers)
                print('Step 4/7 - COMPLETE: Quality Control Noise Data')
        if 5 in STEPS:
                print('Step 5/7 - BEGIN: Spectral Average of Noise Data')
                # !atacr_clean_spectra -O --figQC --figAverage --figCoh --figCross --save-fig --start='{SpecStart}' --end='{SpecEnd}' ./Data/sta_query.pkl> ./Data/Step_5_7_CleanSpectra_logfile.log
                if logoutput_subfolder is None:
                        logoutput_subfolder = dirs['Py_Logs'] + '/5_7'
                ObsQA.TOOLS.io.CleanSpectra(catalog,ATaCR_Parent=ATaCR_Parent,netsta_names=netsta_names,extra_flags=cleanspectra_flags,logoutput_subfolder=logoutput_subfolder,staquery_output=staquery_output,chan=chan,log_prefix=log_prefix,fork=fork,max_workers=max_workers)
                print('Step 5/7 - COMPLETE: Spectral Average of Noise Data')
        if 6 in STEPS:
                print('Step 6/7 - BEGIN: Calculate Transfer Functions')
                # !atacr_transfer_functions -O --figTF --save-fig ./Data/sta_query.pkl> ./Data/Step_6_7_CalcTFs_logfile.log
                if logoutput_subfolder is None:
                        logoutput_subfolder = dirs['Py_Logs'] + '/6_7'
                ObsQA.TOOLS.io.TransferFunctions(catalog,ATaCR_Parent=ATaCR_Parent,netsta_names=netsta_names
                                                #  ,taper_mode=taper_mode
                                                 ,extra_flags=tf_flags,logoutput_subfolder=logoutput_subfolder,staquery_output=staquery_output,chan=chan,log_prefix=log_prefix,fork=fork,max_workers=max_workers)
                print('Step 6/7 - COMPLETE: Calculate Transfer Functions')
        if 7 in STEPS:
                print('Step 7/7 - BEGIN: Correct Event Data')
                # !atacr_correct_event --figRaw --figClean --save-fig ./Data/sta_query.pkl> ./Data/Step_7_7_CorrectEvents_logfile.log
                if logoutput_subfolder is None:
                        logoutput_subfolder = dirs['Py_Logs'] + '/7_7'
                ObsQA.TOOLS.io.CorrectEvents(catalog,ATaCR_Parent=ATaCR_Parent,netsta_names=netsta_names,extra_flags=correctevents_flags,logoutput_subfolder=logoutput_subfolder,staquery_output=staquery_output,chan=chan,log_prefix=log_prefix,fork=fork,max_workers=max_workers)
                print('Step 7/7 - COMPLETE: Correct Event Data')
        plt.close('all')
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# eof