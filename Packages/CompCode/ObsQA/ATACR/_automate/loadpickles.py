def loadpickles(path):
        if not isinstance(path,Path):
                path = Path(path)
        py_files = list(path.glob('*.pkl'))
        if len(py_files)==0:
                raise Exception('Folder contains no .pkl files')
        out = pd.DataFrame({'Output':[ [] for _ in range(len(py_files)) ],'File':[ [] for _ in range(len(py_files)) ]})
        for i,f in enumerate(py_files):
                file = open(f, 'rb')
                pydata = pkl.load(file)
                out.iloc[i]['Output'] = pydata
                out.iloc[i]['File'] = f.name
                file.close()
        return out
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
