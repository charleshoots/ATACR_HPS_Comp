def DayNoiseWhileLoop(catalog,NoiseFolder,ATaCR_Parent,days=10,attempts=50,seed='MESSI_22FIFA_WORLD_CUP!'):
      NoiseFolder = Path(NoiseFolder)
      cat = catalog.copy()
      stafolders = [f for f in NoiseFolder.glob('*/') if f.is_dir()]
      stafolders = [stafolders[np.where(c==np.array([g.name for g in stafolders]))[0][0]] for c in cat.StaName.iloc]
      STEP = 3
      if isinstance(seed,str):
            seed = int(''.join([str(ord(a)) for a in seed]))
      else:
            seed = int(str(seed).replace('.',''))
      for ista,cstaf in enumerate(stafolders):
            print('--<>--'*20)
            nfiles = len([f for f in cstaf.glob('*Z.SAC')])
            nstart = nfiles
            print('Folder: [' + cstaf.name + '] (' + str(ista) + '/' + str(len(stafolders)) + ')'  + '  ||  Current day count: ' + str(nfiles))
            queries = 0
            if nfiles<days:
                  while (queries<attempts) and (nfiles<days):
                        queries+=1
                        seed*=2
                        query_day_len = days - nfiles
                        icatalog = cat[cat.Station==str(cstaf).split('/')[-1].split('.')[-1]]
                        print('Attempt: ' + str(queries) + ', Days needed: ' + str(query_day_len))
                        ObsQA.io.Run_ATaCR(icatalog, ATaCR_Parent = ATaCR_Parent,STEPS=[STEP],log_prefix=icatalog.StaName.iloc[0],days=query_day_len,seed=seed)
                        nfiles = len([f for f in cstaf.glob('*Z.SAC')])
                        if nfiles<days:
                                print(' || Day requirements not satisfied. Attempting new seed. ||')
                        else:
                                print('[' + str(ista) + '/' + str(len(stafolders)) + '][DAY REQUIREMENTS SATISFIED] | Days gained: ' + str(nfiles-nstart))
            else:
                    print('[' + str(ista) + '/' + str(len(stafolders)) + '][DAY REQUIREMENTS ALREADY SATISFIED] | Skipping.')
            if (nfiles<days) & (queries>=attempts):
                    print('Attempt ' + str(attempts) + '/' + str(attempts) + '. Maximum number of attempts exceeded. Skipping station.')
      failedattempts = [(cstaf.parts[-1], len([f for f in cstaf.glob('*Z.SAC')])) for cstaf in stafolders if len([f for f in cstaf.glob('*Z.SAC')])<days]        
      print('Complete')
      if len(failedattempts)>0:
            print('Stations that still do not meet requirements:')
            print(failedattempts)
      else:
            print('All station folders now meet day requirements')
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
