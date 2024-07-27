def randomdays(TStart,TEnd,seed='MESSI_22FIFA_WORLD_CUP',days=10,dateformat = '%Y-%m-%d %H:%M:%S'):
        TStart = datetime.datetime(year=TStart.year,month=TStart.month,day=TStart.day)
        TEnd = datetime.datetime(year=TEnd.year,month=TEnd.month,day=TEnd.day)
        Starts = [TStart + datetime.timedelta(days=j) for j in range((TEnd - TStart).days)]
        Ends = [TStart + datetime.timedelta(days=j+1) for j in range((TEnd - TStart).days)]
        if isinstance(seed,str):
                seed = int(''.join([str(ord(a)) for a in seed]))
        else:
                seed = int(str(seed).replace('.',''))
        rng = np.random.default_rng(seed=seed)
        idx = rng.choice([i for i in range(len(Starts))], size=days, replace=True)
        Starts = np.array(Starts)[idx].tolist()
        Ends = np.array(Ends)[idx].tolist()
        Starts = [str(UTCDateTime.strptime(str(a),dateformat)) for a in Starts]
        Ends = [str(UTCDateTime.strptime(str(a),dateformat)) for a in Ends]
        return Starts,Ends
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
