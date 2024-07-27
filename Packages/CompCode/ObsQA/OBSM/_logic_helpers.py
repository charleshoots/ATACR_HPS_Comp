# ---------------------------------------------------------------------------------------------------------
# Logic Helpers -------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
import ObsQA as OBS
import numpy as _np
import obspy
import copy as cp
from obspy import Stream

def copy(self):
        return cp.deepcopy(self)
def append(self,other):
        self.__add__(other)
def __add__(self,other):
        if isinstance(other,OBS.OBSM.Metrics):
                self._add_traces(tr1=other.traces['1'][0],tr2=other.traces['2'][0],trZ=other.traces['Z'][0],trP=other.traces['P'][0])
        elif isinstance(other,Stream):
                self._add_traces(tr1=other.select(component='1'),tr2=other.select(component='2'),trZ=other.select(component='Z'),trP=other.select(component='H'))
        self._meta()
        self._updatespec()
        return self
def _add_traces(self,tr1=None,tr2=None,trZ=None,trP=None):
        self.traces = Stream([trZ.copy(),tr1.copy(),tr2.copy(),trP.copy()])
        self.A,self.B = self.traces.copy(),self.traces.copy()
        self.ntraces+=1
def __truediv__(self, other):
        self = self.copy()
        if isinstance(other,OBS.OBSM.Metrics):
                self.B = other.traces.copy()
        elif isinstance(other,Stream):
                self.B = other
        self._updatespec()
        return self.copy()
def _meta(self):
        self.dt = _np.unique([s.stats.delta for s in self.traces if isinstance(s,obspy.core.trace.Trace)])
        self.fs = 1/self.dt
        self.tlen = _np.size(self.traces.select(channel='*Z')[0].data)/self.fs
        if self.tlen>7200:
                self.window = 7200
        else:
                self.window = 500