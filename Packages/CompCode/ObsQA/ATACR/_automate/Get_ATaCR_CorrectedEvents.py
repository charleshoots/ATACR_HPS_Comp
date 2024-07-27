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
            folder = Path(eventfolder + '/' + prefix + '/CORRECTED/')
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
