from ObsQA.imports import *

def get_arrivals(sta_llaz,ev_llaz,model = 'iasp91',phases=('ttall',)):
        degdist = obspy.taup.taup_geo.calc_dist(ev_llaz[0],ev_llaz[1],sta_llaz[0],sta_llaz[1],6371,0)
        arrivals = obspy.taup.tau.TauPyModel(model=model).get_travel_times(ev_llaz[2], degdist,phase_list=phases)
        times = [[a.name,a.time] for a in arrivals]
        return times