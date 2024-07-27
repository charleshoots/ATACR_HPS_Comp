def CleanSpectra(catalog,ATaCR_Parent=None,netsta_names=None,extra_flags = '-O --figQC --figAverage --save-fig',logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H',fork=True,max_workers=1):
        # sys.stdout.flush()
        #  --figQC --figAverage --save-fig
        datafolder = './Data/'
        if logoutput_subfolder is not None:
                logoutput = logoutput_subfolder + '/' + log_prefix + '_Step_5_7_CleanSpectra_logfile.log'
                log_fout = logoutput
        else:
                logoutput = '_Step_5_7_CleanSpectra_logfile.log'
                log_fout = datafolder + logoutput
        SpecStart = catalog.Start.min().strftime("%Y-%m-%d, %H:%M:%S")
        SpecEnd = catalog.End.max().strftime("%Y-%m-%d, %H:%M:%S")
        args = [staquery_output]
        [args.append(flg) for flg in extra_flags.split()]
        [args.append(flg) for flg in ['--start={}'.format(SpecStart),'--end={}'.format(SpecEnd)]]
        # with open(log_fout, 'w') as sys.stdout:
        # original = sys.stdout
        log_fout = Path(log_fout)
        log_fout.parent.mkdir(parents=True,exist_ok=True)
        # sys.stdout = open(str(log_fout),'w+')
        logging.basicConfig(stream=log_fout)
        # print('----Begin Clean Spectra----')
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.io.build_staquery(d=catalog,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        args = atacr_clean_spectra.get_cleanspec_arguments(args)
        # atacr_clean_spectra.main(args)
        if not fork:
                atacr_clean_spectra.main(args)
        else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as process_executor:
                        # Start the load operations and mark each future with its URL
                        print('----Begin Clean Spectra----')
                        future = process_executor.submit(atacr_clean_spectra.main,args)
                        print('Waiting for tasks to complete')
                        wait([future])
                future.result()
                print('----Clean Spectra Complete----')
                # sys.stdout = original
        # print('----Clean Spectra Complete----')
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////