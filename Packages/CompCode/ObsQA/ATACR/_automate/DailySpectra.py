def DailySpectra(catalog,ATaCR_Parent=None,netsta_names=None,extra_flags = '-O --figQC --figAverage --figCoh --save-fig',logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H',fork=True,max_workers=1):
        # sys.stdout.flush()
        datafolder = './Data/'
        if logoutput_subfolder is not None:
                logoutput = logoutput_subfolder + '/' + log_prefix + '_Step_4_7_QCSpectra_logfile.log'
                log_fout = logoutput
        else:
                logoutput = '_Step_4_7_QCSpectra_logfile.log'
                log_fout = datafolder + logoutput
        SpecStart = catalog.Start.min().strftime("%Y-%m-%d, %H:%M:%S")
        SpecEnd = catalog.End.max().strftime("%Y-%m-%d, %H:%M:%S")
        args = [staquery_output]
        [args.append(flg) for flg in extra_flags.split()]
        [args.append(flg) for flg in ['--start={}'.format(SpecStart),'--end={}'.format(SpecEnd)]]
        # with open(log_fout, 'w') as sys.stdout:
        # original = sys.stdout
        # sys.stdout = open(log_fout,'w+')
        logging.basicConfig(stream=log_fout)
        # print('----Begin Daily Spectra----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.io.build_staquery(d=catalog,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        args = atacr_daily_spectra.get_dailyspec_arguments(args)
        if not fork:
                atacr_daily_spectra.main(args)
        else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as process_executor:
                        print('----Begin Daily Spectra----')
                        future = process_executor.submit(atacr_daily_spectra.main,args)
                        print('Waiting for tasks to complete')
                        wait([future])
                future.result()
                print('----Daily Spectra Complete----')
                # sys.stdout = original
        # atacr_daily_spectra.main(args)
        # print(' ')
        # print('----Daily Spectra Complete----')
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
