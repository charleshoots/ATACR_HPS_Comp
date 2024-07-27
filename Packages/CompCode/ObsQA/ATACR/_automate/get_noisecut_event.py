def get_noisecut_event(parentfolder,staname,event,channel=['*1','*2','*Z'],len_hrs=2,pre_trim=False,post_trim=True,win_length=200,width=None):
        event_dt = 22  #<--This isn't necessary
        origin = event
        channel = [cc.replace('**','*') for cc in ['*' + c for c in np.array([channel]).flatten().tolist()]]
        folder = Path(str(parentfolder)) / 'Data' / staname / 'HPS_Data'  
        files = []
        # [files.extend(list(folder.glob((event-(3600*event_dt*i)).strftime('%Y.%j') + '*.SAC'))) for i in range(2)]
        [files.extend(list(folder.glob((event-(3600*event_dt*i)).strftime('%Y.%j') + '*.SAC'))) for i in range(1)]
        if len(files)==0:
                print('24hr event data not found')
                return [],()
        raw = Stream([read(f)[0] for f in files])
        raw = Stream([raw.select(channel=c)[0].copy() for c in channel])
        # raw.trim(origin-(3600*event_dt),origin+(3600*2))
        if pre_trim:
                raw.trim(origin,origin+3600*len_hrs)

        hps_out = [run_noisecut(s.copy(),win_length=win_length,width=width) for s in raw]
        corrected = Stream([h[0].copy() for h in hps_out])

        if post_trim:
                corrected.trim(origin,origin+3600*len_hrs)

        hps_spectrograms = [h[1] for h in hps_out]
        raw_collect = {'trZ':raw.select(channel='*Z')[0].copy(),'tr1':raw.select(channel='*1')[0].copy(),'tr2':raw.select(channel='*2')[0].copy()}
        corrected_collect = {'trZ':corrected.select(channel='*Z')[0].copy(),'tr1':corrected.select(channel='*1')[0].copy(),'tr2':corrected.select(channel='*2')[0].copy()}
        out = pd.DataFrame({'Event':event,'Network':staname.split('.')[0],'Station':staname.split('.')[1],'Raw':raw_collect,'Corrected':corrected_collect,'Spectrograms':hps_spectrograms})
        return out