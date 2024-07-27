import sys
sys.path.insert(1,'/Users/charlesh/Documents/Codes/')
sys.path.append('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS')
sys.path.insert(0, '/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/COMPS')
sys.path.insert(0, '/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python/OBStools')
sys.path.insert(0, '/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python')
import ObsQA

import ObsQA as OBS
# from ObsQA.classes import OBS.Metrics
from scipy.stats import norm
from cmcrameri import cm
# from ObsQA import *
import matplotlib.colors as mcolors
from obspy.clients.fdsn import Client


from ObsQA.imports import *
from obspy.core import UTCDateTime as _UTCDateTime
import glob as g
import pandas as pd
import numpy as np
import pickle as pkl
import datetime as _datetime
import os as os
from pathlib import Path
import concurrent
from concurrent.futures import wait
import ObsQA
# from ObsQA import classes
import time
import logging
import matplotlib.pyplot as plt
import requests 
from PIL import Image 

def OBS_Generator(catalog,parent,events_folder = 'EVENTS',avg='STA',tf='ZP-21',evireq=None,sta=None, net=None,hps_win_length=200,width=None):
        # datafolder = ATaCR_Py_DataFolder['Py_CorrectedTraces']
        nsta = len(catalog)
        cat = catalog.copy()
        staind = [i for i in range(len(cat))]
        if sta is not None:
            fi = np.atleast_1d(np.in1d(cat.Station,np.array(sta).squeeze()))
            cat = cat[fi].copy()
            staind = staind[np.squeeze(np.where(fi))]
        if net is not None:
            fi = np.in1d(cat.Network,np.array(net))
            cat = cat[fi].copy()
            staind = staind[np.squeeze(np.where(fi))]
        staind = np.atleast_1d(np.array(staind))
        for stai,Station in zip(staind,cat.iloc):
                StaName = Station.StaName
                evts = list(Station.Events)
                if evireq is not None:
                    evts = evts[evireq]
                evind = [i for i in range(len(evts))]
                nevents = len(evts)
                try:
                    Metrics = OBS.TOOLS.io.get_Noise(parent,Station.Network,Station.Station,avg='STA',update=True)
                except:
                    Metrics = dict()
                    print(StaName + ' || Missing Noise Data')
                for evi,event_request in zip(evind,evts):
                        event_request = '.'.join([s.zfill(2) for s in [event_request.split('.')][0]])
                        # Metrics,Comps = get_metrics_comp(Station.Network,Station.Station,parent,event=event_request,events_folder = events_folder,avg=avg,tf=tf,hps_win_length=hps_win_length,width=width,Metrics=Metrics)
                        try:
                            Metrics,Comps = get_metrics_comp(Station.Network,Station.Station,parent,event=event_request,events_folder = events_folder,avg=avg,tf=tf,hps_win_length=hps_win_length,width=width,Metrics=Metrics)
                        except:
                              print(StaName + '  || ' + event_request + ' || Missing Event Data')
                              continue
                        Notify = StaName + ' | ' + str(stai+1) + '/' + str(nsta) + ' | ' + str(evi+1) + '/' + str(nevents) + ' : ' + event_request
                        Event = event_request + '  ' + 'm' + str(Station.Magnitude_mw[evi]) + '  ' + str(Station.Depth_KM[evi]) + 'km'
                        print(Notify)
                        try:
                            Station.CurrentEvent = Station.Metadata[evi]
                        except:
                            Station.CurrentEvent = Station.Metadata[evi][0]
                        Station.Event_LLAZ = [Station.CurrentEvent.origins[0].latitude,Station.CurrentEvent.origins[0].longitude,Station.CurrentEvent.origins[0].depth/1000]
                        Station.STA_LLAZ = [Station.Latitude,Station.Longitude,Station.Water_Depth_m/1000]
                        yield (Event,Station,Metrics,Comps)

def event_generator(catalog,evireq=None,net=None,sta=None):
        nsta = len(catalog)
        cat = catalog.copy()
        staind = np.array([i for i in range(len(cat))])
        if sta is not None:
            fi = np.in1d(cat.Station,np.array(sta))
            # fi = np.atleast_1d(np.in1d(cat.Station,np.array(sta).squeeze()))
            cat = cat[fi].copy()
            staind = staind[fi]
        if net is not None:
            fi = np.in1d(cat.Network,np.array(net))
            cat = cat[fi].copy()
            staind = staind[fi]
        staind = np.atleast_1d(np.array(staind))
        for stai,Station in zip(staind,cat.iloc):
                StaName = Station.StaName
                evts = list(Station.Events)
                evind = [i for i in range(len(evts))]
                if evireq is not None:
                    evts = evts[evireq]
                    evind = evind[evireq]
                nevents = len(evts)
                net = Station.Network
                sta = Station.Station
                for evi,event_request in zip(evind,evts):
                        Notify = StaName + ' | ' + str(stai+1) + '/' + str(nsta) + ' | ' + str(evi+1) + '/' + str(nevents) + ' : ' + Station.Events[evi]
                        Station.CurrentEvent = Station.Metadata[evi][0]
                        Station.Event_LLAZ = [Station.CurrentEvent.origins[0].latitude,Station.CurrentEvent.origins[0].longitude,Station.CurrentEvent.origins[0].depth/1000]
                        Station.STA_LLAZ = [Station.Latitude,Station.Longitude,Station.Water_Depth_m/1000]
                        print(Notify)
                        yield net,sta,event_request,Station

def get_metrics_comp(net,sta,methodfolder,event,return_atacr=True,return_hps=True,return_raw=True,return_noise=True,events_folder = 'EVENTS',avg='STA',tf='ZP-21',hps_win_length=200,width=None,Metrics=dict()):
        Comps = dict()
        atacr_parent = Path(methodfolder) / 'ATaCR' / 'ATaCR_Python'
        atacr_data_folder = atacr_parent / events_folder
        hps_data_folder = Path(methodfolder) / 'NoiseCut' / 'Data' / 'Output'
        noise_data_folder = atacr_parent
        if return_atacr:
            if return_raw:
                d = OBS.TOOLS.io.Get_ATaCR_CorrectedEvents(atacr_data_folder,[event],net,sta,tfavg=avg.lower())
                Raw = d.Raw[0]
                Comps['RawP'] = Raw['trP'].copy()
                Comps['RawZ'] = Raw['trZ'].copy()
                Metrics['Raw'] = OBS.Metrics(tr1=Raw['tr1'].copy(),tr2=Raw['tr2'].copy(),trZ=Raw['trZ'].copy(),trP=Raw['trP'].copy())
        if return_hps:
            d_hps = OBS.TOOLS.io.get_noisecut_output(net,sta,event,hps_data_folder)
            if len(d_hps)==0:
                print('||--[HPS] Confirming No Data')
                return dict(),[]
        if return_atacr:
            d = OBS.TOOLS.io.Get_ATaCR_CorrectedEvents(atacr_data_folder,[event],net,sta,tfavg=avg.lower())
            Comps['d'] = d
            if len(d)==0:
                print('||--[ATaCR] Confirming No Data')
                return dict(),[]
            Comps['PostATACR'] = get_atacr(d,tf=tf)
            Metrics['ATaCR'] = OBS.Metrics(tr1=Raw['tr1'].copy(),tr2=Raw['tr2'].copy(),trZ=Comps['PostATACR'].copy(),trP=Raw['trP'].copy())
            if return_raw:
                Metrics['ATaCR'] = Metrics['ATaCR'] / Metrics['Raw'].copy()
        if return_hps:
            Comps['HPS_Stream'] = d_hps.HPS[0]
            Comps['Raw_Stream'] = d_hps.Raw[0]
            Comps['PostHPS'] = d_hps.HPS[0].select(channel='*Z')[0]
            Comps['PostHPS_H1'] = d_hps.HPS[0].select(channel='*1')[0]
            Comps['PostHPS_H2'] = d_hps.HPS[0].select(channel='*2')[0]
            Comps['RawP'] = d_hps.Raw[0].select(channel='*H')[0]
            Comps['RawZ'] = d_hps.Raw[0].select(channel='*Z')[0]
            Comps['Raw1'] = d_hps.Raw[0].select(channel='*1')[0]
            Comps['Raw2'] = d_hps.Raw[0].select(channel='*2')[0]
            # Comps['PostHPS'],Comps['PostHPS_H1'],Comps['PostHPS_H2'] = get_hps(Raw=Raw,win_length=hps_win_length,width=width) ## The old 2hr noisecut method

            Metrics['HPS'] = OBS.Metrics(tr1=Comps['PostHPS_H1'].copy(),tr2=Comps['PostHPS_H2'].copy(),trZ=Comps['PostHPS'].copy(),trP=Comps['RawP'].copy())
            if return_raw:
                Metrics['Raw'] = OBS.Metrics(tr1=Comps['Raw1'].copy(),tr2=Comps['Raw2'].copy(),trZ=Comps['RawZ'].copy(),trP=Comps['RawP'].copy())
                Metrics['HPS'] = Metrics['HPS'] / Metrics['Raw'].copy()
            if return_atacr:
                PostATACR = Comps['PostATACR']
                Comps['PostBoth'] = get_hps(PostATACR=PostATACR,win_length=hps_win_length,width=width)
                Metrics['ATaCR_HPS'] = OBS.Metrics(tr1=Comps['PostHPS_H1'].copy(),tr2=Comps['PostHPS_H2'].copy(),trZ=Comps['PostBoth'].copy(),trP=Comps['RawP'].copy())
                if return_raw:
                    Metrics['ATaCR_HPS'] = Metrics['ATaCR_HPS'] / Metrics['Raw'].copy()
                Metrics['ATaCR w/ HPS'] = Metrics['ATaCR'].copy() / Metrics['HPS'].copy()
        if return_noise:
            Metrics['Noise'] = OBS.TOOLS.io.get_Noise(noise_data_folder,net,sta,avg=avg,update=True)['Noise']
            for k in list(Metrics.keys()):
                 Metrics[k].Noise = Metrics['Noise']
        return Metrics,Comps

def get_atacr(d,tf='ZP-21'):
        PostATACR = d.Corrected[0][tf].copy()
        PostATACR.stats.location = 'ATaCR (' + PostATACR.stats.location + ')'
        return PostATACR

def get_hps(Raw=None,PostATACR=None,win_length=200,width=None):
        if Raw:
            PostHPS, spectrograms = noisecut(Raw['trZ'].copy(), ret_spectrograms=True,win_length=win_length,width=width)
            PostHPS.spectrograms = spectrograms #I just care about the vertical stft
            PostHPS_H1,spectrograms = noisecut(Raw['tr1'].copy(), ret_spectrograms=True,win_length=win_length,width=width)
            PostHPS_H2,spectrograms = noisecut(Raw['tr2'].copy(), ret_spectrograms=True,win_length=win_length,width=width)
            # PostHPS_H2.spectrograms = spectrograms
            # PostHPS_H1.spectrograms = spectrograms
        if PostATACR:
            PostBoth, spectrograms = noisecut(PostATACR.copy(), ret_spectrograms=True,win_length=win_length,width=width) #<move out of coditions?
            PostBoth.spectrograms = spectrograms

        if (PostATACR is not None) and (Raw is not None):
            return PostHPS,PostHPS_H1,PostHPS_H2,PostBoth
        elif Raw is not None:
            return PostHPS,PostHPS_H1,PostHPS_H2
        elif PostATACR is not None:
            return PostBoth


def get_obs_cmaps(light_min=0.4,light_max=0.9,extrema_n=4):
        cmaps = dict()
        raw_palette = 'grayC_categorical'
        atacr_palette = 'bilbao_categorical'
        nc_palette = 'devon_categorical'
        both_palette = 'acton_categorical'
        comp_palette = 'bamako_categorical'
        for k,p in zip(['raw','atacr','nc','both','comp'],[raw_palette,atacr_palette,nc_palette,both_palette,comp_palette]):
                colors = cm.cmaps[p].colors
                colors = cm.cmaps[p].colors[((cm.cmaps[p].colors**2).sum(axis=1)**(0.5))<light_max]
                colors = colors[(colors**2).sum(axis=1)**(0.5)>light_min]
                if k=='atacr':
                        colors[0][0]  = colors[0][0]**(0.5**(extrema_n+1))
                        for i,c in enumerate(colors):
                            c[0] = c[0]**(0.5**extrema_n)
                            colors[i] = c
                if k=='nc':
                        colors[0][2]  = colors[0][2]**(0.5**(extrema_n+1))
                        for i,c in enumerate(colors):
                            c[2] = c[2]**(0.5**extrema_n)
                            colors[i] = c
                if k=='raw':
                        # pass
                        colors[0] = [0,0,0]
                        # colors[0] = colors[0]**(2**extrema_n)
                        # for i,c in enumerate(colors):
                        #     c = c**(2**extrema_n)
                        #     colors[i] = c
                if k=='both':
                        colors[0][0]  = colors[0][0]**(0.5**(extrema_n+1))
                        colors[0][2]  = colors[0][2]**(0.5**(extrema_n+1))
                        for i,c in enumerate(colors):
                            c[0] = c[0]**(0.5**extrema_n)
                            c[2] = c[2]**(0.5**extrema_n)
                            colors[i] = c
                if k=='comp':
                        colors[0][1]  = colors[0][1]**(0.5**(extrema_n+1))
                        for i,c in enumerate(colors):
                            c[1] = c[1]**(0.5**extrema_n)
                            colors[i] = c
                colors = [mcolors.to_hex(c) for c in colors]
                cmaps[k] = colors
        return cmaps

def demean_detrend(stream):
        stream.detrend()
        for i in range(len(stream)):
            stream[i].data = stream[i].data - np.mean(stream[i].data)  
        stream.detrend()
        return stream

def preparetraces(stream,trim=(0,7200),band=None,sortindex=None,max_percentage=0.01,max_length=500,type='hann'):
        # Ideal order of operations: demean -> detrend -> taper -> FILTER -> taper -> demean -> detrend
        if sortindex is None:
            stream = [stream[i] for i in np.argsort([ev.stats.sac.dist for ev in stream])]
        else:
            stream = [stream[i] for i in sortindex]
        stream = Stream(traces=stream)
        for _ in range(2):
            stream = demean_detrend(stream.copy())
        stream.taper(max_percentage=max_percentage, type=type, max_length=max_length)
        if band:
            fmin,fmax = np.sort(band)
            [ev.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=4, zerophase=True) for ev in stream]
        stream.trim(stream[0].stats.starttime + trim[0],stream[0].stats.starttime + trim[1])
        stream.taper(max_percentage=max_percentage, type=type, max_length=max_length)
        for _ in range(2):
            stream = demean_detrend(stream.copy())
        return stream

def build_event_catalog(new_catalog,N_permag_perwin=7,mags=[[6.6,7.0],[6.4,6.59],[6.0,6.39]],windows=None,hardcap=None):
        # wins = [[catalog[catalog.Network==e].Start.min() ,catalog[catalog.Network==e].End.max()] for e in catalog.Network.unique()]
        cat = new_catalog.copy()
        client = Client()
        # Just a good all-around default set of windows for my dataset
        if windows is None:
            wins = [[UTCDateTime(month=12,year=2009,day=1),UTCDateTime(month=1,year=2011,day=1)],
            [UTCDateTime(month=1,year=2011,day=1),UTCDateTime(month=9,year=2011,day=1)],
            [UTCDateTime(month=9,year=2011,day=1),UTCDateTime(month=9,year=2012,day=1)],
            [UTCDateTime(month=9,year=2012,day=1),UTCDateTime(month=9,year=2013,day=1)],
            [UTCDateTime(month=9,year=2013,day=1),UTCDateTime(month=6,year=2014,day=1)],
            [UTCDateTime(month=6,year=2014,day=1),UTCDateTime(month=10,year=2015,day=1)],
            [UTCDateTime(month=5,year=2018,day=1),UTCDateTime(month=9,year=2019,day=1)]]
        elif windows=='sta':
            pass
        for stai,station in enumerate(cat.iloc):
            if windows=='sta':
                wins = [[station.Start,station.End]]
            N_permag_perwin = (20 - station.n_events)*2
            print('---'*20)
            # str(station.Start) + ' -> ' + str(station.End)
            print('[' + str(stai+1) + '/' + str(len(cat)) + ']     | ' + str(station.StaName))
            query = []
            skips = 0
            for wi,(Start,End) in enumerate(wins):
                for minmagnitude,maxmagnitude in (mags):
                    if len(query)<(N_permag_perwin*len(mags)):
                        events = client.get_events(starttime=Start, endtime=End,minmagnitude=minmagnitude,maxmagnitude=maxmagnitude,orderby='magnitude',limit=N_permag_perwin)
                    events = events.filter("time > " + str(station.Start),"time < " + str(station.End))
                    if query==[]:
                        query = events
                    else:
                        query +=events
                    if len(events)==0:
                        skips+=1
                print('Windows skipped: ' + str(skips) + '/' + str(len(mags)*len(wins)))
            if len(query)<(N_permag_perwin*len(mags)):
                print('Widening arc, too few events: ' + str(len(query)))
            events = client.get_events(starttime=station.Start, endtime=station.End,minmagnitude=min([min(m) for m in mags]),maxmagnitude=max([max(m) for m in mags]),limit=(N_permag_perwin*len(mags))-len(query))
            query+=events
            if len(query)==0:
                Exception('no events')
            query = obspy.Catalog([query[i] for i in np.unique([j.origins[0].time for j in query.events],return_index=True)[1]]) #<Remove duplicate events
            if hardcap:    
                if len(query)>hardcap:
                    query = obspy.Catalog(query[:(hardcap)])
            n_events = len(query)
            print('N-Events: ' + str(n_events))
            Magnitude_mw = [a.magnitudes[0].mag for a in query]
            Origin = [a.origins[0] for a in query]
            Metadata = query
            Averaging = []
            Events = [a.origins[0].time.strftime('%Y.%j.%H.%M') for a in query]
            Files = []
            Depth_KM = [a.origins[0].depth/1000 for a in query]
            cat.at[stai,'n_events'] = n_events
            cat.at[stai,'Magnitude_mw'] = Magnitude_mw
            cat.at[stai,'Origin'] = Origin
            cat.at[stai,'Metadata'] = Metadata
            cat.at[stai,'Averaging'] = Averaging
            cat.at[stai,'Events'] = Events
            cat.at[stai,'Files'] = Files
            cat.at[stai,'Depth_KM'] = Depth_KM
            cat = cat.reset_index(drop=True)
        print('N-Events Minimum: ' + str(cat.n_events.min()))
        return cat

def pull_from_catalog(add_stations,catalog_orig):
        catalog_orig['Station'] = [str(a) for a in catalog_orig.Station]
        qi = [np.where(catalog_orig.Station==a)[0][0] for a in add_stations]
        cat = catalog_orig.iloc[qi]
        W = dict()
        [W.update({x:y}) for x,y in zip(list(cat.columns),[a.replace(' ','_').replace('/','').replace('(','').replace(')','').replace('_deg','') for a in list(cat.columns)])]
        cat = cat.rename(columns=W)
        cat['Good_Channels'] = list(cat[['Z_Is_Good', 'H1_Is_Good','H2_Is_Good', 'P_Is_Good']].sum(axis=1)==4)
        cat = cat[['Station', 'Network', 'Latitude', 'Longitude', 'Experiment','Instrument_Design', 'Seismometer',
        'Environment', 'Pressure_Gauge','Water_Depth_m', 'Distance_from_Land_km','Distance_to_Plate_Boundary_km', 'Sediment_Thickness_m',
        'Surface_Current_ms', 'Crustal_Age_Myr', 'Start', 'End','Deployment_Length_days', 'Good_Channels']]
        cat['Network_Experiment'] = ['[' + a + '] ' + b for a,b in zip(cat.Network,cat.Experiment)]
        cat['StaName'] = [a + '.' + b for a,b in zip(cat.Network,cat.Station)]
        cat['n_events'] = 0
        cat['Magnitude_mw'] = None
        cat['Origin'] = None
        cat['Metadata'] = None
        cat['Averaging'] = None
        cat['Events'] = None
        cat['Files'] = None
        cat['Depth_KM'] = None
        if np.any(cat['Good_Channels']==False):
            print('You have selected stations with dead channels')
        else:
            return cat
        return cat

def get_mustang_noise_spec_url(net,sta,cha,start_t,end_t,quality = 'M',loc='--',title=None,db=None,palette='rainbow'):
    if title is None:
        title = '.'.join([net,sta,cha])
    req_parameters = dict()
    req_parameters['base_url'] = 'https://service.iris.edu/mustang/noise-spectrogram/1/query?'
    req_parameters['net'] = 'net=' + net
    req_parameters['sta'] = 'sta=' + sta
    req_parameters['loc'] = 'loc=' + loc
    req_parameters['cha'] = 'cha=' + cha
    req_parameters['quality'] = 'quality=' + quality
    req_parameters['starttime'] = 'starttime=' + start_t
    req_parameters['endtime'] = 'endtime=' + end_t
    req_parameters['output'] = 'output=power'
    req_parameters['out_format'] = 'format=plot'
    # req_parameters['plot_width'] = 'plot.width=2048'
    # req_parameters['plot_height'] = 'plot.height=1024'
    req_parameters['color_palette'] = 'plot.color.palette=' + palette
    req_parameters['plot_xticks'] = 'plot.time.tickunit=month'
    if db is not None:
        req_parameters['powerscale'] = 'plot.powerscale.range=' + str(min(db)) + ',' + str(max(db))
    req_parameters['title'] = 'plot.title=' + title
    req_url = req_parameters['base_url'] + '&'.join([req_parameters[k] for k in list(req_parameters.keys())[1:]])
    return req_url

def mustang_noise_spec_reports(cat,mustang_folder,channel='Z',loc='--',db=None,quality='M'):
    problem_stations = []
    for stai,Station in enumerate(cat.iloc):
        cha = 'BH'+channel
        print('--'*20)
        notify = str(stai+1) + '/' + str(len(cat)) + ' | Problem stations:' + str(len(problem_stations))
        net = Station.Network
        sta = Station.Station
        title = '.'.join([net,sta,cha]) + ' | Z:' + str(int(Station.Water_Depth_m)) + 'm'
        start_t = Station.Start.strftime('%Y-%m-%d')
        end_t = Station.End.strftime('%Y-%m-%d')
        file = '.'.join([str(net),str(sta),str(cha)])
        fmt = '.png'
        outfile = mustang_folder / (file + '.day_noise' + fmt)
        palette = 'rainbow'
        url = get_mustang_noise_spec_url(net=net,sta=sta,cha=cha,start_t = start_t,end_t= end_t,quality = quality,loc=loc,title=title,db=db,palette=palette)
        print(notify)
        data = requests.get(url)
        if not data.ok:
            print('|| [' + cha + ']-- Request: ' + file + ' | 404 Error')
            cha = '*'+channel
            title = '.'.join([net,sta,cha]) + ' | Z:' + str(int(Station.Water_Depth_m)) + 'm'
            url = get_mustang_noise_spec_url(net=net,sta=sta,cha=cha,start_t = start_t,end_t= end_t,quality = quality,loc=loc,title=title,db=db)
            data = requests.get(url)
        if not data.ok:
            print('|| [' + cha + ']-- Request: ' + file + ' | 404 Error')
            cha = '*H'+channel
            title = '.'.join([net,sta,cha]) + ' | Z:' + str(int(Station.Water_Depth_m)) + 'm'
            url = get_mustang_noise_spec_url(net=net,sta=sta,cha=cha,start_t = start_t,end_t= end_t,quality = quality,loc=loc,title=title,db=db)
            data = requests.get(url)
        if not data.ok:
            print('|| [' + cha + ']-- Request: ' + file + ' | 404 Error')
            cha = '*H*'+channel
            title = '.'.join([net,sta,cha]) + ' | Z:' + str(int(Station.Water_Depth_m)) + 'm'
            url = get_mustang_noise_spec_url(net=net,sta=sta,cha=cha,start_t = start_t,end_t= end_t,quality = quality,loc=loc,title=title,db=db)
            data = requests.get(url)
        if not data.ok:
            print('|| [' + cha + ']-- Request: ' + file + ' | 404 Error')
            cha = 'LH'+channel
            title = '.'.join([net,sta,cha]) + ' | Z:' + str(int(Station.Water_Depth_m)) + 'm'
            url = get_mustang_noise_spec_url(net=net,sta=sta,cha=cha,start_t = start_t,end_t= end_t,quality = quality,loc=loc,title=title,db=db)
            data = requests.get(url)
        if not data.ok:
            print('|| [' + cha + ']-- Request: ' + file + ' | 404 Error')
            cha = 'LD'+channel
            title = '.'.join([net,sta,cha]) + ' | Z:' + str(int(Station.Water_Depth_m)) + 'm'
            url = get_mustang_noise_spec_url(net=net,sta=sta,cha=cha,start_t = start_t,end_t= end_t,quality = quality,loc=loc,title=title,db=db)
            data = requests.get(url)
        if not data.ok:
            print('|| [' + cha + ']-- Request: ' + file + ' | 404 Error')
            cha = 'BD'+channel
            title = '.'.join([net,sta,cha]) + ' | Z:' + str(int(Station.Water_Depth_m)) + 'm'
            url = get_mustang_noise_spec_url(net=net,sta=sta,cha=cha,start_t = start_t,end_t= end_t,quality = quality,loc=loc,title=title,db=db)
            data = requests.get(url)
        if not data.ok:
            problem_stations.append(stai)
            _=[print('|| ' + p.split('?')[-1]) for p in url.split('&')]
            print('|| [' + cha + ']-- Request: ' + file + ' | 404 Error')
            continue
        _=[print('|| ' + p.split('?')[-1]) for p in url.split('&')]
        print('[' + cha + ']-- Request: ' + file + ' | Good')
        f = open(outfile,'wb') 
        f.write(data.content) 
        f.close()
    print('---'*30)
    print('---'*30)
    print('Problem stations::: ')
    display(cat.iloc[problem_stations])
    return problem_stations

def save_tight(filename,fig=None,dpi=200,format=None):
        # Saves figure to PDF with no margins. Do not modify
        # plt.gca().set_axis_off()
        # plt.subplots_adjust(top = 1, bottom = 0.0, right = 1, left = 0,hspace = 0.07, wspace = 0.03)
        plt.margins(0.1,0.1)
        # plt.gca().xaxis.set_major_locator(plt.NullLocator())
        # plt.gca().yaxis.set_major_locator(plt.NullLocator())
        if fig is None:
                plt.savefig(filename, bbox_inches = 'tight',pad_inches = 0.05,dpi=dpi,format=format)
        else:
                fig.savefig(filename, bbox_inches = 'tight',pad_inches = 0.05,dpi=dpi,format=format)
        return 'Complete'
def fnotch(d):
        '''The frequency knotch root function described in Crawford et al., 1998.
        depth (d) is in meters. Returned (f) is in Hz.'''
        g = 9.80665
        f = (g/(2*np.pi*d))**0.5
        return f

# _________________________________________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||||||||||||||| HPS Spectrogram Plots (Replicates Zali) |||||||||||||||||||||||||
# _________________________________________________________________________________________________________________
def hps_spectrograms(S_full, S_background, S_hps, frequencies, times,figsize=(11,12)):
        t = times
        f = frequencies
        s = S_full
        cmap = 'magma'
        xlabel = True
        yscale = 'log'
        ax = None
        vmin,vmax = None,None
        fig, axes = plt.subplots(nrows=3, ncols=1,figsize=figsize,height_ratios=[1,1,1],width_ratios=[1],layout='constrained',squeeze=False,sharey='col',sharex='row')
        titles = ['Raw','Noise','Noise Removed']
        for r,s in enumerate([S_full, S_background, S_hps]):
                gax = axes[r,0]
                pc = gax.pcolormesh(t,f, 10*np.log10(s), cmap = cmap, shading= 'auto')
                if (vmin is not None) and (vmax is not None):
                        pc.set_clim(vmin, vmax)
                else:
                        vmin,vmax = pc.get_clim()
                if ylabel:
                        gax.set_ylabel('Frequency (Hz)',fontweight='bold')
                if xlabel:
                        gax.set_xlabel('Time (s)',fontweight='bold')
                gax.set_yscale(yscale)
                gax.set_ylim(f[1],f[-1])
                gax.set_xlim(t[0],t[-1])
                gax.set_title(titles[r],fontweight='bold')
                if fig is not None:
                        fig.colorbar(pc, ax=gax, pad=0.01, label='dB')

# _________________________________________________________________________________________________________________
# ||||||||||||||||||||||||||||||||||||||||||||||| EVENT RECORD FUNCTION |||||||||||||||||||||||||||||||||||||||||||
# _________________________________________________________________________________________________________________
def event_record_plot(evstream,evstream_back=None,linewidth=0.2,trim = (None,None),scales = [1,1],band = None,facecolor=('b','r'),norm = 'trace',figsize = (20,13),phases = ('P','S','SKS','PKiKP','SKiKS','SKSSKS',),evdepth=None,title='',sortindex=None,ax=None,normscale=1.0,residual_fraction=1.0):
        # linewidth=0.2
        # prepare traces
        if evstream_back:
            sets = [evstream_back,evstream] # [UNCORRECT_SET , CORRECT_SET] 
        else:
            sets = [evstream]
        sets = [preparetraces(stream,trim=trim,band=band,sortindex=sortindex) for stream in sets]
        if len(sets)==2:
            residuals = [ev0.data - ev.data for ev0,ev in zip(sets[0],sets[1])]
            residuals = [normscale * residual_fraction * np.array(res) / np.max(np.abs(np.array(res))) for res in residuals]
        fig = None
        if ax is None:
            fig, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize,layout='constrained',squeeze=False,sharey='all',sharex='all')
            ax = axes[0,0]
        # norming for plots
        postsetmax = [np.max(ev.data) for ev in  sets[-1]]
        for si in range(len(sets)):
            for i in range(len(sets[0])): #norm the uncorrected set to the max in the corrected set
                # -----------
                if isinstance(norm,list):
                        sets[si][i].data = sets[si][i].data/norm[i]
                elif norm.lower()=='postset':
                        sets[si][i].data = sets[si][i].data/postsetmax[i]
                elif norm.lower()=='trace':
                        sets[si][i].data = sets[si][i].data/np.max(abs(sets[si][i].data))
                elif norm.lower()=='col':
                        # print('Trace norm scaling by r')
                        dist = [s.stats.sac.dist for s in sets[si]]
                        norms = [d/np.max(dist) for d in dist]
                        norms = [d**(1) for d in norms]
                        sets[si][i].data = sets[si][i].data/np.max(abs(sets[si][i].data))
                        sets[si][i].data = sets[si][i].data / norms[i]
                        colmax = np.max([abs(d.data.max()) for d in sets[si]])
                        sets[si][i].data = sets[si][i].data / colmax
                # -----------
        correctset = sets[-1]
        [ax.plot(tr.times(),tr.data + ysep,linewidth=linewidth,color='k') for ysep,tr in enumerate(correctset)]
        if len(sets)>1:
            for si,s in enumerate(sets):
                [ax.fill_between(tr.times(),tr.data*scales[si] + ysep,tr.data*0 + ysep, where=np.abs(tr.data)>=0, facecolor=facecolor[si]) for ysep,tr in enumerate(s)]
        [ax.plot(tr.times(),tr.data*0 + ysep,linewidth=0.4,color='k') for ysep,tr in enumerate(correctset)]
        if evdepth is not None:
            arrivals = [ObsQA.TOOLS.io.get_arrivals(sta_llaz=(sta.stats.sac.stla,sta.stats.sac.stlo,sta.stats.sac.stel),ev_llaz=(sta.stats.sac.evla,sta.stats.sac.evlo,evdepth),phases=phases) for sta in correctset]
            ardict = dict()
            corephase_dy = 0.2
            direcphase_dy = 0.03
            [[ardict.update({ph[0]:[]}) for ph in a] for a in arrivals]
            [[ardict[ph[0]].extend([ph[1]]) for ph in a] for a in arrivals]
            [ardict.update({k:np.max(ardict[k])}) for k in list(ardict.keys())]
            [[ax.vlines(ph[1], ysep, ysep+1, colors='b',linewidth=0.2) for ph in a] for ysep,a in enumerate(arrivals)]
            [[ax.vlines(ph[1], ysep, ysep+1, colors='r',linewidth=0.5,alpha=0.5) for ph in a if ph[0]=='P'] for ysep,a in enumerate(arrivals)]
            [[ax.vlines(ph[1], ysep, ysep+1, colors='r',linewidth=0.5,alpha=0.5) for ph in a if ph[0]=='S'] for ysep,a in enumerate(arrivals)]
            [ax.text(ardict[k], len(correctset) - corephase_dy, k, color='b',horizontalalignment='center') for k in list(ardict.keys()) if (k!='P') and (k!='S')]
            [ax.text(ardict[k], len(correctset) + direcphase_dy, k, color='r',horizontalalignment='center') for k in list(ardict.keys()) if ((k=='P') or (k=='S'))]
        yl = (-1,len(correctset))
        ax.set_ylim(yl)
        ax.set_xlim(correctset[0].times()[0],correctset[0].times()[-1])
        ax.set_yticks([i for i in range(len(correctset))])
        labels = [str(int(ev.stats.sac.dist)) +'km' + ' [' + ev.stats.network + '] ' + ev.stats.station  + '\n depth:' + str(int(abs(ev.stats.sac.stel*1000))) + 'm' for ev in correctset]
        ax.set_yticklabels(labels)
        ax.set_xlabel('Time(s)')
        if fig is not None:
            fig.suptitle(title,fontweight='bold',fontsize=15)
        return ax
def mirror_audit(evaudit,datafolder):
  audit = []
  for ev in evaudit.iloc:
    # pass
    # print(ev.Event)
    inside_atacr = []
    inside_hps = []
    for n,s in zip(ev.Networks,ev.Stations):
      stafolder = Path(datafolder) / 'Data' / (n + '.' + s)
      f = '.'.join(ev.Event.split('.')[:2]) + '*.SAC'
      N = len(list((stafolder / 'HPS_Data').glob(f)))
      if N!=4:
        inside_hps.append(False)
      else:
        inside_hps.append(True)
      f = n + '.' + s + '.' + ev.Event + '*.SAC'
      stafolder = Path(datafolder) / 'EVENTS' / (n + '.' + s)
      N = len(list((stafolder / 'CORRECTED').glob( '*' + f)))
      if N!=4:
        inside_atacr.append(False)
      else:
        inside_atacr.append(True)
    audit.append([np.array(inside_hps),np.array(inside_atacr)])
    return audit

def update_catalog(catalog,eventsfolder):
  for ii,s in enumerate(catalog.iloc):
    ids = np.where(np.array([len(np.where(np.unique(['.'.join((f.name.split('.')[:4])) for f in list((Path(eventsfolder) / s.StaName).glob('*.SAC'))])==e)[0]) for e in s.Events])>0)[0]
    Events = list(np.array(s.Events)[ids])
    catalog.at[ii,'Magnitude_mw'] = list(np.array(catalog.at[ii,'Magnitude_mw'])[ids])
    catalog.at[ii,'Origin'] = [catalog.at[ii,'Origin'][j] for j in ids]
    catalog.at[ii,'Metadata'] = obspy.Catalog([catalog.at[ii,'Metadata'][h] for h in ids])
    catalog.at[ii,'Events'] = Events
    catalog.at[ii,'n_events'] = len(Events)
    catalog.at[ii,'Depth_KM'] = [catalog.at[ii,'Depth_KM'][j] for j in ids]
  return catalog

def log_smoothing(x,y):
        # Peak Log-Smoothed trendline.
        # Math is inefficient but it works. Needs to be replaced with np.ceil(np.log10(Nyq/NFFT))
        dw = x[np.log10(x)<=np.ceil(np.log10(x[1]))].size
        dx = x[np.log10(x)<=np.ceil(np.log10(x[1]))].max()
        kernel_size = int(np.round((dx/x[-1])*x.size))
        kernel = np.ones(kernel_size) / kernel_size
        y = np.convolve(y, kernel, mode='valid')
        x = x[np.arange(0,y.size,dw)]
        y = y[np.arange(0,y.size,dw)]
        return x,y

# def save_tight(filename,fig=None,dpi=200):
#         # Saves figure to PDF with no margins. Do not modify
#         # plt.gca().set_axis_off()
#         plt.subplots_adjust(top = 1, bottom = 0.0, right = 1, left = 0,hspace = 0.07, wspace = 0.03)
#         plt.margins(0.1,0.1)
#         # plt.gca().xaxis.set_major_locator(plt.NullLocator())
#         # plt.gca().yaxis.set_major_locator(plt.NullLocator())
#         if fig is None:
#             plt.savefig(filename, bbox_inches = 'tight',pad_inches = 0.05,dpi=dpi)
#         else:
#             fig.savefig(filename, bbox_inches = 'tight',pad_inches = 0.05,dpi=dpi)
#         return 'Complete'