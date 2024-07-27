def DownloadEvents(catalog,ATaCR_Parent=None,netsta_names=None,Minmag=6.3,Maxmag=6.7,limit=1000,pre_event_min_aperture=1,logoutput_subfolder='',log_prefix = '',staquery_output = './sta_query.pkl',chan='H'):
        # sys.stdout.flush()
        logfilename = '_Step_2_7_EventDownload_logfile.log'
        if logoutput_subfolder is not None:
                logoutput = logoutput_subfolder + '/' + log_prefix
        else:
                logoutput = '_Step_2_7_EventDownload_logfile.log'
        datafolder = './EVENTS/'
        dateformat = '%Y.%j.%H.%M'
        print('----Begin Event Download----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        for i in range(len(catalog)):
                csta = catalog.iloc[i]
                S = csta['Station']#station
                N = csta['Network']
                staname = str(N) + '.' + str(S)
                if os.path.isdir(datafolder + staname)==False:
                        os.system('mkdir ' + datafolder + '/' + staname)
                ObsQA.io.build_staquery(d=catalog[(catalog.Network==N) & (catalog.Station==S)],staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
                log_fout = logoutput  + logfilename
                # log_fout
                logoutput_subfolder = Path(logoutput_subfolder)
                logoutput_subfolder.mkdir(parents=True,exist_ok=True)
                # original = sys.stdout
                # sys.stdout = open(log_fout,'w+')
                logging.basicConfig(stream=log_fout)
                print('--' + staname + '--',flush=True)
                for j in range(len(csta.Events)):
                        ev = csta.Events[j]
                        print('---'*20)
                        print(staname + '| S:[' +str(i+1) + '/' + str(len(catalog)) + '] | E:[' + str(j+1) + '/' + str(len(csta.Events)) + '] | ' + str(ev),flush=True)
                        print('---'*20)
                        EventStart = UTCDateTime.strptime(ev,dateformat)
                        EventEnd = UTCDateTime.strptime(ev,dateformat) + datetime.timedelta(minutes=1)
                        print('Event ' + str(j+1) + '/' + str(len(csta.Events)) + ' | Start -> End : ' + str(EventStart) + ' -> ' + str(EventEnd))
                        if ATaCR_Parent is not None:
                                os.chdir(ATaCR_Parent)
                        # args = [staquery_output,'--start={}'.format(EventStart), '--end={}'.format(EventEnd),'--min-mag={}'.format(Minmag),'--max-mag={}'.format(Maxmag),'--limit={}'.format(limit)]
                        args = [staquery_output,'--start={}'.format(EventStart), '--end={}'.format(EventEnd),'--min-mag={}'.format(Minmag),'--max-mag={}'.format(Maxmag)]
                        # [print(ii) for ii in [staquery_output,'--start={}'.format(EventStart), '--end={}'.format(EventEnd),'--min-mag={}'.format(Minmag),'--max-mag={}'.format(Maxmag),'--limit={}'.format(limit)]]
                        # with open(log_fout, 'w') as sys.stdout:
                        atacr_download_event.main(atacr_download_event.get_event_arguments(args))
        print(' ')
        print('----Event Download Complete----')
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
