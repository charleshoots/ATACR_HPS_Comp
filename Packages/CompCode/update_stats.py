from obspy.clients.fdsn import Client
from obspy.core import read, Stream, Trace, AttribDict, UTCDateTime
client = Client('IRIS')

def update_stats(trace, stla=None, stlo=None, stel=None, cha=None, evla=None, evlo=None):
    """
    Function to include SAC metadata to :class:`~obspy.core.Trace` objects
    Parameters
    ----------
    trace : :class:`~obspy.core.Trace` object
        Trace object to update
    stla : float
        Latitude of station
    stlo : float
        Longitude of station
    stel : float
        Station elevation (m)
    cha : str
        Channel for component
    evla : float, optional
        Latitude of event
    evlo : float, optional
        Longitute of event
    Returns
    -------
    trace : :class:`~obspy.core.Trace` object
        Updated trace object
    """
    trace.stats.sac = AttribDict()
    if stla is None:
        inventory = client.get_stations(network=trace.stats.network, station=trace.stats.station,starttime=trace.stats.starttime,endtime=trace.stats.endtime)
        stla,stlo,stel = [inventory[0][0].__dict__[k] for k in ['_latitude','_longitude','_elevation']]

    trace.stats.sac.stla = stla
    trace.stats.sac.stlo = stlo
    trace.stats.sac.stel = stel

    if cha is not None:
        trace.stats.sac.kcmpnm = cha
        trace.stats.channel = cha
    if evla is not None and evlo is not None:
        trace.stats.sac.evla = evla
        trace.stats.sac.evlo = evlo
    return trace