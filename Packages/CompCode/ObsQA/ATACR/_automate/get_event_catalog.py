def get_event_catalog(eventsfolder,subfolder='CORRECTED',fmt = 'pkl'):
        evdict = dict()
        evdict['folder'] = [fld.split('/')[-2] for fld in g.glob(eventsfolder + '/*/')]
        evdict['Network'] = [fld.split('.')[0] for fld in evdict['folder']]
        evdict['Station'] = [fld.split('.')[1] for fld in evdict['folder']]
        evdict['StaName'] = [str(a) + '.' + str(b) for a,b in zip(evdict['Network'],evdict['Station'])]
        # evdict['folder'][0]
        evdict['n_events'] = list()
        evdict['events'] = list()
        for i in range(len(evdict['folder'])):
                staname =  evdict['StaName'][i]
                path = Path(eventsfolder)
                path = path / evdict['folder'][i] / subfolder
                files = list(path.glob('*.' + fmt))
                events = [str(f.name).replace(staname + '.','').replace('.sta','').replace('.day','').replace('.' + fmt,'') for f in files]
                # events = list(np.unique(['.'.join(files[g].split('.' + fmt)[0].split('.')[0:-1]) for g in range(len(files))]))
                evdict['events'].append(events)
                evdict['n_events'].append(len(events))
        catalog = pd.DataFrame(evdict)
        catalog = catalog.sort_values(by=['Network','Station'])
        catalog = catalog.reset_index(drop=True)
        return catalog
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
