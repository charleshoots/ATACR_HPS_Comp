def AuditEventFolder(eventsfolder,subfolder='CORRECTED',parseby='pkl',Minmag=6.0,Maxmag=7.0):
        parse = {'SAC':'*Z.SAC','pkl':'*.pkl'}
        catalog = ObsQA.io.get_event_catalog(eventsfolder)
        stations = ObsQA.io.getstalist()
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
