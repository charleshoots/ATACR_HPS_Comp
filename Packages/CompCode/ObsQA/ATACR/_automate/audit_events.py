def audit_events(eventsfolder):
    client = Client()
    Minmag = 6.0
    Maxmag = 7.0
    timedelta = 60
    stations = ObsQA.io.get_event_catalog(eventsfolder)
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
        start = UTCDateTime.strptime(str(ev),'%Y.%j.%H.%M')
        cat = client.get_events(starttime=start, endtime=start + datetime.timedelta(seconds=timedelta),minmagnitude=Minmag, maxmagnitude=Maxmag)
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
