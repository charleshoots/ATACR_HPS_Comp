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
