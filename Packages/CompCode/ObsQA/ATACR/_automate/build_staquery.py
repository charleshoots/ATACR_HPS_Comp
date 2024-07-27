def build_staquery(d,chan='H',ATaCR_Parent=None,staquery_output='./sta_query.pkl'):
        out = dict()
        if not isinstance(staquery_output,Path):
                if ATaCR_Parent is None:
                        staquery_output = os.getcwd() + staquery_output.replace('./','/')
                else:
                        os.chdir(ATaCR_Parent)
                        staquery_output = ATaCR_Parent + '/' + staquery_output.replace('./','')
        for csta in d.iloc:
                net = csta.Network
                sta = csta.Station
                key = net + '.' + sta
                csta_dict = {'station':sta, 'network':net, 'altnet':[],'channel':chan, 'location':['--'], 'latitude':csta['Latitude'], 'longitude':csta['Longitude'], 'elevation':-csta['Water_Depth_m']/1000, 'startdate':UTCDateTime(csta.Start), 'enddate':UTCDateTime(csta.End), 'polarity':1.0, 'azcorr':0.0, 'status':'open'}
                out[key] = csta_dict
        output = pd.DataFrame.from_dict(out)
        if staquery_output is not None:
                output.to_pickle(staquery_output)
                print('Station Query File Written to: ' + staquery_output)
        else:
                return output
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
