def getstalist():
        current_path = os.path.dirname(__file__)
        excelfile = current_path + '/Janiszewski_etal_2023_StationList.xlsx'

        stations = pd.read_excel(excelfile)
        d = dict()
        d.update({a:b for a,b in zip(stations.columns,[c.replace('(','').replace(')','').replace(' ','_').replace('/','') for c in stations.columns])})
        stations = stations.rename(columns=d)
        stations.Station = [str(s) for s in stations.Station]
        stations.Network = [str(s) for s in stations.Network]
        stations['StaName'] = stations.Network.astype(str) + '.' + stations.Station.astype(str)
        allgood = np.in1d(stations[['Z_Is_Good','H1_Is_Good','H2_Is_Good','P_Is_Good']].sum(axis=1).tolist(),4)
        stations['Good_Channels'] = allgood
        stations = stations.assign(n_events=pd.Series())
        stations = stations.assign(Magnitude_mw=pd.Series())
        stations = stations.assign(Depth_KM=pd.Series())
        stations = stations.assign(Origin=pd.Series())
        stations = stations.assign(Metadata=pd.Series())
        stations = stations.assign(Averaging=pd.Series())
        stations = stations.assign(Events=pd.Series())
        stations = stations.assign(Files=pd.Series())
        stations = stations.sort_values(by=['Network','Station'])
        stations = stations.reset_index(drop=True)
        return stations
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
