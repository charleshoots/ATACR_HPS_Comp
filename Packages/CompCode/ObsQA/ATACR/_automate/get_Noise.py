def get_Noise(atacr_parent,net,sta,avg='Day',update=True):
        if avg.lower()=='sta':
                specfold = 'AVG_STA'
                filetag = 'avg_sta'
        elif avg.lower()=='day':
                specfold = 'SPECTRA'
                filetag = 'spectra'
        subfolder = atacr_parent + '/' + specfold + '/' + net + '.' + sta
        if update==False:
                out = pd.read_pickle(atacr_parent + '/' + specfold + '/' + net + '.' + sta + '/Metrics/' + net + '.' + sta + '.' + avg + 'Metrics.pkl')
        else:
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
                        noise['Noise'] = [o]
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
                                gd = noise_out.gooddays[0]
                                gw = ObsQA.io.get_Noise(atacr_parent,net,sta,avg='Day',update=True).goodwins[gd]
                                ft = ObsQA.io.get_Noise(atacr_parent,net,sta,avg='Day',update=True).ft[gd]
                                noise_trace['1'] = np.real(np.fft.ifft(np.mean(np.array([np.mean(spec['1'][g],axis=0) for spec,g in zip(ft.iloc,gw.iloc)]),axis=0)))
                                noise_trace['2'] = np.real(np.fft.ifft(np.mean(np.array([np.mean(spec['2'][g],axis=0) for spec,g in zip(ft.iloc,gw.iloc)]),axis=0)))
                                noise_trace['Z'] = np.real(np.fft.ifft(np.mean(np.array([np.mean(spec['Z'][g],axis=0) for spec,g in zip(ft.iloc,gw.iloc)]),axis=0)))
                                noise_trace['P'] = np.real(np.fft.ifft(np.mean(np.array([np.mean(spec['P'][g],axis=0) for spec,g in zip(ft.iloc,gw.iloc)]),axis=0)))
                                noise['AVG_Noise_Trace'] = [noise_trace]
                                noise_out = pd.DataFrame.from_dict(noise)
                        if len(out)==0:
                                out = noise_out
                        else:
                                out = pd.concat([out,noise_out])
                if avg.lower()=='day':
                        out = out[['Noise','GoodNoise','File','key','tkey','julday','year','QC','av',
                                        'window','overlap','dt','npts','fs','ncomp','tf_list',
                                        'goodwins','rotation','f','ft']]
                        out = out.sort_values(by=['year','julday'],ascending=False)
                else:
                        out = out[['key','Noise','AVG_Noise_Trace','File', 'year', 'julday','f', 'nwins']]
        out.reset_index(drop=True,inplace=True)
        if avg.lower()=='sta':
                Metrics = dict()
                Metrics_Noise = ObsQA.classes.OBSMetrics(csd=out.Noise[0],f=out.Noise[0].f)
                # Metrics['Noise'] = out
                Metrics['Noise'] = Metrics_Noise
                out = Metrics
        return out
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
