def TransferFunctions(catalog,ATaCR_Parent=None,netsta_names=None,extra_flags = '-O --figTF --save-fig',taper_mode=0,logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H',fork=True,max_workers=1):
        # sys.stdout.flush()
        datafolder = './Data/'
        if logoutput_subfolder is not None:
                logoutput = logoutput_subfolder + '/' + log_prefix + '_Step_6_7_CalcTFs_logfile.log'
                log_fout = logoutput
        else:
                logoutput = '_Step_6_7_CalcTFs_logfile.log'
                log_fout = datafolder + logoutput
        args = [staquery_output]
        extra_flags = extra_flags + ' --taper ' + str(taper_mode)
        [args.append(flg) for flg in extra_flags.split()]
        # with open(log_fout, 'w') as sys.stdout:
        # original = sys.stdout
        # sys.stdout = open(log_fout,'w+')
        logging.basicConfig(stream=log_fout)
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.io.build_staquery(d=catalog,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        args = atacr_transfer_functions.get_transfer_arguments(args)
        # atacr_transfer_functions.main(args)
        if not fork:
                atacr_transfer_functions.main(args)
        else:
                print('----Begin Transfer Functions----')
                with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as process_executor:
                        future = process_executor.submit(atacr_transfer_functions.main,args)
                        print('Waiting for tasks to complete')
                        wait([future])
                future.result()
                print('----Transfer Functions Complete----')
                # sys.stdout = original
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////