def DownloadDayNoise(catalog,days=[],end_delta=24,event_mode=False,seed='MESSI_22FIFA_WORLD_CUP!',ATaCR_Parent=None,netsta_names=None,logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H'):
        logfilename = log_prefix + '_Step_3_7_NoiseDownload_logfile.log'
        dateformat = '%Y.%j.%H.%M'
        event_dt = 22
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        datafolder = Path('./Data/')
        if logoutput_subfolder is not None:
                logfolder = Path(logoutput_subfolder)
                logfolder.mkdir(exist_ok=True,parents=True)
        else:
                logfolder = datafolder
        print('----Begin Noise Download----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        for i,Station in enumerate(catalog.iloc):
                if event_mode:
                        Origins = Station.Origin
                        Starts = Origins
                        if isinstance(Origins[0],obspy.core.event.origin.Origin):
                                Starts = [e.time for e in Origins]
                        Ends = [e+(3600*24) for e in Starts]
                elif isinstance(days,int):
                        Starts,Ends = ObsQA.io.randomdays(Station.Start,Station.End,seed=seed,days=days)
                else:
                        Starts = days
                        Ends = [datetime.datetime(year=s.year,month=s.month,day=s.day) + datetime.timedelta(hours=end_delta) for s in Starts]
                # --Log file stdout stuff--
                staname = Station.StaName
                sta_folder = Path(datafolder / staname)
                sta_folder.mkdir(exist_ok=True,parents=True)
                if logoutput_subfolder is not None:
                        log_fout = logfolder / logfilename
                else:
                        log_fout = sta_folder / logfilename
                # original = sys.stdout
                # sys.stdout = open(log_fout,'w+')
                logging.basicConfig(stream=log_fout)
                print('--' + staname + '--',flush=True)

                ObsQA.io.build_staquery(d=Station.to_frame().T,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
                for j,(NoiseStart,NoiseEnd) in enumerate(zip(Starts,Ends)):
                        print('<||>'*30)
                        print(staname + ' Station ' +str(i+1) + '/' + str(len(catalog)) + ' - Day ' + str(j+1) + '/' + str(len(Starts)),flush=True)
                        args = [staquery_output,'--start={}'.format(NoiseStart), '--end={}'.format(NoiseEnd)]
                        atacr_download_data.main(atacr_download_data.get_daylong_arguments(args))
        print(' ')
        print('----Noise Download Complete----')
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
