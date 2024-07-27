#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def CorrectEvents(catalog,ATaCR_Parent=None,netsta_names=None,extra_flags = '--figRaw --figClean --save-fig',logoutput_subfolder=None,log_prefix = '',staquery_output='./sta_query.pkl',chan='H',fork=True,max_workers=1):
        datafolder = './Data/'
        staquery_output = str(staquery_output)
        if logoutput_subfolder is not None:
                logoutput = logoutput_subfolder + '/' + log_prefix + '_Step_7_7_CorrectEvents_logfile.log'
                log_fout = logoutput
        else:
                logoutput = '_Step_7_7_CorrectEvents_logfile.log'
                log_fout = datafolder + logoutput
        args = [staquery_output]
        [args.append(flg) for flg in extra_flags.split()]
        # with open(log_fout, 'w') as sys.stdout:
        if netsta_names is not None:
                catalog = catalog[np.in1d((catalog.Network + '.' + catalog.Station),netsta_names)]
        ObsQA.io.build_staquery(d=catalog,staquery_output = staquery_output,chan=chan,ATaCR_Parent = ATaCR_Parent)
        # original = sys.stdout
        # sys.stdout = open(log_fout,'w+')
        logging.basicConfig(stream=log_fout)
        if ATaCR_Parent is not None:
                os.chdir(ATaCR_Parent)
        args = atacr_correct_event.get_correct_arguments(args)
        if not fork:
                atacr_correct_event.main(args)
        else:
                with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as process_executor:
                        future = process_executor.submit(atacr_correct_event.main,args)
                        print('Waiting for tasks to complete')
                        wait([future])
                future.result()
                print('----Correct Events Complete----')
                # sys.stdout = original
        # atacr_correct_event.main(args)
        # sys.stdout = original
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////