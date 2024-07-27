def MLtoStaNoise(cnet,csta,preprocfolder,specfolder,level='b1'):
        DN = ObsQA.io.MLtoDayNoise(preprocfolder,[cnet],[csta])
        cfold = specfolder + '/' + cnet + '/' + csta + '/' + level.lower()
        files = g.glob(cfold + '/*.mat')
        df = ObsQA.io.mat2df(files)
        files = df.File.unique()
        udf = df.copy().iloc[0].to_frame().T.drop(0)
        for i in range(len(files)):
                udf = pd.concat([udf,df.iloc[list(df.File==files[i])].iloc[0].to_frame().T])
        df = udf.copy()
        for ifile in range(len(files)):
                cfile = files[ifile]
                cdf = df.iloc[list(df.File==cfile)]
                c11 = cdf.iloc[0].c11_stack.T
                c22 = cdf.iloc[0].c22_stack.T
                cZZ = cdf.iloc[0].czz_stack.T
                cPP = cdf.iloc[0].cpp_stack.T
                c12 = cdf.iloc[0].c12_stack.T
                c1Z = cdf.iloc[0].c1z_stack.T
                c1P = cdf.iloc[0].c1p_stack.T
                c2Z = cdf.iloc[0].c2z_stack.T
                c2P = cdf.iloc[0].c2p_stack.T
                cZP = cdf.iloc[0].cpz_stack.T
                cHH = cdf.iloc[0].chh_stack.T
                cHZ = cdf.iloc[0].chz_stack.T
                cHP = cdf.iloc[0].chp_stack.T
                tilt = cdf.iloc[0].rotor
                rotcoh = cdf.iloc[0].rotcoh
                f = cdf.iloc[0].f
                nwins = len(cdf.iloc[0].goodwins)
                direc = np.arange(0,360+10,10) #default in ML ATaCR code (b1)
                ncomp = sum(cdf.iloc[0].comp_exist)
                key = cnet + '.' + csta
                tf_list = {'ZP': True, 'Z1': True, 'Z2-1': True, 'ZP-21': True, 'ZH': True, 'ZP-H': True}
                power = obstools.atacr.classes.Power(c11 = c11, c22 = c22, cZZ = cZZ, cPP = cPP)
                cross = obstools.atacr.classes.Cross(c12=c12, c1Z=c1Z, c1P=c1P, c2Z=c2Z, c2P=c2P, cZP=cZP)
                rotation = obstools.atacr.classes.Rotation(cHH=cHH, cHZ=cHZ, cHP=cHP, coh=None, ph=None, tilt=tilt, coh_value=rotcoh, phase_value=None, direc=direc)
                goodwins = np.array(cdf.iloc[0].goodwins,dtype=bool)
                ind = np.where(np.char.replace(DN.iloc[0].Files,'.mat','_' + level + '_spectra.mat')==cfile)[0][0]
                DN.iloc[0].DayNoise[ind].power = power
                DN.iloc[0].DayNoise[ind].cross = cross
                DN.iloc[0].DayNoise[ind].rotation = rotation
                DN.iloc[0].DayNoise[ind].f = cdf.iloc[0].f
                DN.iloc[0].DayNoise[ind].goodwins = goodwins
        SN = StaNoise(daylist=DN.iloc[0].DayNoise)
        SN.init()
        SN.key = 'ML-ATaCR: ' + SN.key
        return SN
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
