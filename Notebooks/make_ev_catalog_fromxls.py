# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------
# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------
datemin = datetime.datetime(year=2010,day=1,month=9)
datemax = datetime.datetime(year=2010,day=30,month=12)
n,nn=5,5
client = Client()
evts = client.get_events(starttime=datemin, endtime=datemax,minmagnitude=6.0,maxmagnitude=7.0,orderby='time-asc')
evts = [e for e in evts]
e = list(evts[0:n])
e.extend(evts[-nn:])
evts = e
Events = [str(e.origins[0].time.year) + '.' + str(e.origins[0].time.julday) + '.' + str(e.origins[0].time.hour) + '.' + str(e.origins[0].time.minute) for e in evts]
Magnitudes_mw = [e.magnitudes[0].mag for e in evts]
Origin = [e.origins[0] for e in evts]
Times = [e.origins[0].time for e in evts]
Metadata = evts
datemin = datetime.datetime(year=2012,day=1,month=2)
datemax = datetime.datetime(year=2012,day=1,month=4)
n,nn=1,5
client = Client()
evts = client.get_events(starttime=datemin, endtime=datemax,minmagnitude=6.0,maxmagnitude=7.0,orderby='time-asc')
e = list(evts[0:n])
e.extend(evts[-nn:])
evts = e
Events.extend([str(e.origins[0].time.year) + '.' + str(e.origins[0].time.julday) + '.' + str(e.origins[0].time.hour) + '.' + str(e.origins[0].time.minute) for e in evts])
Magnitudes_mw.extend([e.magnitudes[0].mag for e in evts])
Origin.extend([e.origins[0] for e in evts])
Times.extend([e.origins[0].time for e in evts])
Metadata.extend(evts)
datemin = datetime.datetime(year=2012,day=1,month=9)
datemax = datetime.datetime(year=2013,day=1,month=1)
n,nn=5,1
client = Client()
evts = client.get_events(starttime=datemin, endtime=datemax,minmagnitude=6.0,maxmagnitude=7.0,orderby='time-asc')
e = list(evts[0:n])
e.extend(evts[-nn:])
evts = e
Events.extend([str(e.origins[0].time.year) + '.' + str(e.origins[0].time.julday) + '.' + str(e.origins[0].time.hour) + '.' + str(e.origins[0].time.minute) for e in evts])
Magnitudes_mw.extend([e.magnitudes[0].mag for e in evts])
Origin.extend([e.origins[0] for e in evts])
Times.extend([e.origins[0].time for e in evts])
Metadata.extend(evts)
datemin = datetime.datetime(year=2013,day=1,month=9)
datemax = datetime.datetime(year=2013,day=1,month=10)
n,nn=5,5
client = Client()
evts = client.get_events(starttime=datemin, endtime=datemax,minmagnitude=6.0,maxmagnitude=7.0,orderby='time-asc')
e = list(evts[0:n])
e.extend(evts[-nn:])
evts = e
Events.extend([str(e.origins[0].time.year) + '.' + str(e.origins[0].time.julday) + '.' + str(e.origins[0].time.hour) + '.' + str(e.origins[0].time.minute) for e in evts])
Magnitudes_mw.extend([e.magnitudes[0].mag for e in evts])
Origin.extend([e.origins[0] for e in evts])
Times.extend([e.origins[0].time for e in evts])
Metadata.extend(evts)
event_catalog = dict()
event_catalog['Times'] = Times
event_catalog['Events'] = Events
event_catalog['Magnitudes_mw'] = Magnitudes_mw
event_catalog['Origin'] = Origin
event_catalog['Metadata'] = Metadata
event_catalog['n_sta'] = [0 for i in range(len(Metadata))]
event_catalog = pd.DataFrame.from_dict(event_catalog)
event_catalog = event_catalog.iloc[np.argsort(event_catalog.Times)]

catalog = pd.concat([catalog,pd.DataFrame(columns = list(event_catalog.keys()))])
catalog = pd.concat([catalog,pd.DataFrame(columns = ['n_events'])])
for jj,Station in enumerate(catalog.iloc):
    ind = (Station.End>event_catalog['Times']) & (Station.Start<event_catalog['Times'])
    idx = np.where(ind)[0]
    N = len(event_catalog[ind].Events.tolist())
    for k in event_catalog.keys():
        catalog.iat[jj,np.where(np.array(list(catalog.keys()))==k)[0][0]] = event_catalog[ind][k].tolist()
    catalog.iat[jj,np.where(np.array(list(catalog.keys()))=='n_events')[0][0]] = len(catalog.iloc[jj].Events)
    for i in idx:
        event_catalog.at[i,'n_sta']+=1
# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------
# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------# ###########---------
if len(catalog[catalog.n_events==0])>0:
    display('These stations are not used')
    catalog[catalog.n_events==0]
else:
    print('|| All ' +  str(len(catalog)) + ' stations are used with atleast ' + str(catalog.n_events.min()) + ' events ||')
    print('|| All ' + str(len(event_catalog)) + ' events are used with atleast ' + str(event_catalog.n_sta.min()) + ' stations ||')
display(event_catalog)
# event_catalog.to_pickle('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python/EVENTS/event_catalog_evrecord_set.pkl')
# catalog.to_pickle('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python/EVENTS/sta_catalog_evrecord_set.pkl')