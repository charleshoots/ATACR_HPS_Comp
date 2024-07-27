def GetML_EventData_and_TransferFunctions(event_time,network,sta,preprocfolder=None,correctedfolder=None,tffolder=None):
        evdata = None
        corrected_evdata = None
        ml_tf = None
        if preprocfolder is not None:
                folder = preprocfolder
                path = folder + '/' + event_time + '/*.mat'
                ml_files = g.glob(path)
                evdata = ObsQA.io.mat2df(ml_files)
                evdata = evdata[(evdata['Network']==network)&(evdata['Station']==sta)]
                evdata = ObsQA.io.organize_evdata(evdata)
        if correctedfolder is not None:
                folder = correctedfolder
                path = folder + '/' + network + '/' + sta + '/' + network + sta + '_' + str(event_time) + '_corrseis' + '.mat'
                ml_files = g.glob(path)
                corrected_evdata = ObsQA.io.mat2df(ml_files)
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
                ml_tf = ObsQA.io.ClosestMLPreEventTF(ftf,event_time)
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
