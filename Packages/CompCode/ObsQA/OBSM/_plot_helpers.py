# ---------------------------------------------------------------------------------------------------------
# Plot Helpers --------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
def preparetraces(self,stream,trim=(None,None),band=None,max_percentage=0.01,max_length=500,type='hann',norm=True):
        stream = Stream(traces=stream)
        stream.detrend()
        for i in range(len(stream)):
                stream[i].data = stream[i].data - _np.mean(stream[i].data)  
        stream.detrend()
        stream.taper(max_percentage=max_percentage, type=type, max_length=max_length)
        if band:
                fmin,fmax = _np.sort(band)
                [ev.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=4, zerophase=True) for ev in stream]
        for i in range(len(stream)):
                stream[i].data = stream[i].data - _np.mean(stream[i].data)
        if _np.all(_np.array(trim)==None):  
                trim = 7200 #Might be wise to trim first and last few seconds from trace before normalizing. 
                # If a filter artifact is present at zero-time, for whatever reason, it will break the norm.
                stream.trim(stream[0].stats.starttime,stream[0].stats.starttime + trim)
        else:
                if isinstance(trim[0],UTCDateTime):
                        stream.trim(trim[0],trim[1])
                else:
                        stream.trim(stream[0].stats.starttime + trim[0],stream[0].stats.starttime + trim[1])
        stream.detrend()
        stream.taper(max_percentage=max_percentage, type=type, max_length=max_length)
        if norm:
                div = _np.abs(_np.max(stream[0].data))   
                stream[0].data = stream[0].data / div
        stream.detrend()
        return stream[0]
def get_arrivals(self,sta_llaz,ev_llaz,model = 'iasp91',phases=('ttall',)):
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