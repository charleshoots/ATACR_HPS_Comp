def mat2df(files):
        '''
        Converts all loaded matlab variables into a single dataframe
        '''
        if not isinstance(files,list):
                files = [files]
        files = sorted(files)
        out = []
        df = pd.DataFrame()
        FNUM = 0
        for f in files:
                df2 = pd.DataFrame()
                cfold = f[0:maxstrfind(f,'/')+1]
                cf = f[maxstrfind(f,'/')+1:len(f)]
                mat = spio.loadmat(f, simplify_cells=True)
                matvars = list(mat.keys())[3:len(list(mat.keys()))]
                oldkey = list(mat.keys())[-1] #Replaces Last Key with a single hardcoded name. Assumes a single variable (ie one struct) was saved to the matfile
                for k in range(3): #pop out the first three keys. they are artifacts from the mat import
                        mat.pop(list(mat.keys())[0])
                for k in matvars:
                        if isinstance(mat[k],list):
                                for dct in mat[k]:
                                        tmp = pd.DataFrame.from_dict(dct,orient='index').T
                                        tmp['FNUM'] = FNUM
                                        df2 = pd.concat([df2,tmp])
                                mat[k] = df2
                        else:
                                mat[k] = pd.DataFrame.from_dict(mat[k]).T
                                mat[k]['FNUM'] = FNUM
                if len(matvars)>1:
                        for i in range(len(matvars)-1):
                                mat[matvars[0]] = pd.merge(mat[matvars[0]],mat[matvars[i+1]],on='FNUM')
                fdf = mat[matvars[0]].set_index('FNUM')
                fdf['Folder'] = cfold
                fdf['File'] = cf
                out.append(fdf)
                FNUM += 1
        for d in out:
                df = pd.concat([df,d])
        df = df.groupby('FNUM',as_index=False,dropna=False).ffill().groupby('FNUM',as_index=False,dropna=False).bfill() #<---Fills data gaps using adjacent rows from only the same file.
        df = df.sort_values('File')
        return df
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////