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
