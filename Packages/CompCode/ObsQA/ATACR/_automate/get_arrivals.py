def get_arrivals(sta_llaz,ev_llaz,model = 'iasp91',phases=('ttall',)):
        ''''
        Simple function pulls arrival times and phase names for a given event observed at a given station
        sta_llaz = List object containing [Lat,Lon] of the station
        ev_llaz = List object containing [Lat,Lon] of the event
        phases = Tuple object containing a list of all desired phases. 'ttall'
        (default) will give every single phase available.
        -Charles Hoots,2022
        '''
        degdist = obspy.taup.taup_geo.calc_dist(ev_llaz[0],ev_llaz[1],sta_llaz[0],sta_llaz[1],6371,0)
        arrivals = obspy.taup.tau.TauPyModel(model=model).get_travel_times(ev_llaz[2], degdist,phase_list=phases)
        times = [[a.name,a.time] for a in arrivals]
        return times
#### \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#### ---------------------------------------------------------------------------------------------------------------------------------------------------
#### ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////