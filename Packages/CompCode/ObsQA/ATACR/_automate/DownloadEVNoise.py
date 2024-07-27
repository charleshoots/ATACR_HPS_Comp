def DownloadEVNoise(catalog,ATaCR_Parent=None,netsta_names=None,pre_event_day_aperture=10,logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H'):
        # sys.stdout.flush()
        logfilename = log_prefix + '_Step_3_7_NoiseDownload_logfile.log'
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        datafolder = Path('./Data/')
        if logoutput_subfolder is not None:
                logfolder = Path(logoutput_subfolder)
                logfolder.mkdir(exist_ok=True,parents=True)
        else:
                logfolder = datafolder
        dateformat = '%Y.%j.%H.%M'

        print('----Begin Noise Download----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        for i in range(len(catalog)):
                csta = catalog.iloc[i]
                S = csta['Station']#station
                N = csta['Network']
                staname = str(N) + '.' + str(S)
                sta_folder = Path(datafolder / staname)
                sta_folder.mkdir(exist_ok=True,parents=True)
                ObsQA.io.build_staquery(d=catalog[(catalog.Network==N) & (catalog.Station==S)],staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
                if logoutput_subfolder is not None:
                        log_fout = logfolder / logfilename
                else:
                        log_fout = sta_folder / logfilename
                # original = sys.stdout
                # sys.stdout = open(log_fout,'w+')
                logging.basicConfig(stream=log_fout)
                print('--' + staname + '--',flush=True)
                for j in range(len(csta.Events)):
                        print(staname + ' Station ' +str(i+1) + '/' + str(len(catalog)) + ' - Event ' + str(j+1) + '/' + str(len(csta.Events)),flush=True)
                        ev = csta.Events[j]
                        NoiseStart = UTCDateTime.strptime(ev,dateformat) - datetime.timedelta(days=pre_event_day_aperture)
                        NoiseStart = NoiseStart - datetime.timedelta(hours = NoiseStart.hour, minutes = NoiseStart.minute, seconds=NoiseStart.second) #rounds down to the nearest day
                        NoiseEnd = UTCDateTime.strptime(ev,dateformat)
                        NoiseEnd = NoiseEnd - datetime.timedelta(hours = NoiseEnd.hour, minutes = NoiseEnd.minute, seconds=NoiseEnd.second) #rounds down to the nearest day
                        args = [staquery_output,'--start={}'.format(NoiseStart), '--end={}'.format(NoiseEnd)]
                        atacr_download_data.main(atacr_download_data.get_daylong_arguments(args))
        print(' ')
        print('----Noise Download Complete----')
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
