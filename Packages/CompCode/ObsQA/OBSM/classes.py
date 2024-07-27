import numpy as _np
import obspy
from obspy.core import read, Stream, Trace, AttribDict, UTCDateTime
from scipy.signal import stft, detrend
import copy
import obstools
import matplotlib.pyplot as plt
import librosa

from obspy.imaging.cm import pqlx
from obspy.signal import PPSD
import obspy.imaging.cm as cm


from . import _logic_helpers
from . import _signal_qa
from . import _plot_helpers
from . import _plotters
from . import _signal_helpers

class OBSMetrics(object):

        def __init__(self,tr1=None, tr2=None, trZ=None, trP=None,overlap=0.3,csd=None,f=None):
                self.overlap = overlap
                self.traces=dict()
                self.traces_ppsd = dict()
                self.csd_pairs = ['1P','2P','1Z','2Z','ZP','PP','11','22','ZZ']
                self.traces['1'],self.traces['2'],self.traces['Z'],self.traces['P'] = [],[],[],[]
                self.csd = dict()
                self.csd['A'] = dict()
                self.csd['B'] = dict()
                self.csd['AB'] = dict()
                self.Noise = []
                for p in self.csd_pairs:
                        self.csd['A'][p] = []
                        self.csd['B'][p] = []
                        self.csd['AB'][p] = []
                self.ntraces = 0
                self.StaNoise = None
                self.f = f
                if tr1 is not None:
                        self._add_traces(tr1=tr1,tr2=tr2,trZ=trZ,trP=trP)
                        self.A = self.traces.copy()
                        self.B = self.traces.copy()
                        self._meta()
                        self._updatespec()
                if csd is not None:
                        if isinstance(csd,obstools.atacr.StaNoise):
                                dct = csd.power.__dict__.copy()
                                dct.update(csd.cross.__dict__.copy())
                                self.StaNoise = csd
                                csd = dct.copy()
                        for a in self.csd_pairs:
                                try:
                                        self.csd['A'][a] = csd['c'+a]
                                        self.csd['B'][a] = csd['c'+a]
                                        self.csd['AB'][a] = csd['c'+a]
                                except:
                                        try:
                                                self.csd['A'][a] = csd['c'+a[::-1]]
                                                self.csd['B'][a] = csd['c'+a[::-1]]
                                                self.csd['AB'][a] = csd['c'+a[::-1]]
                                        except:
                                                try:
                                                        self.csd['A'][a] = csd[a]
                                                        self.csd['B'][a] = csd[a]
                                                        self.csd['AB'][a] = csd[a]
                                                except:
                                                        self.csd['A'][a] = csd[a[::-1]]
                                                        self.csd['B'][a] = csd[a[::-1]]
                                                        self.csd['AB'][a] = csd[a[::-1]]
        copy = _logic_helpers.copy
        append = _logic_helpers.append
        __add__ = _logic_helpers.__add__
        _add_traces = _logic_helpers._add_traces
        __truediv__ = _logic_helpers.__truediv__
        _meta = _logic_helpers._meta
        
        ft = _signal_qa.ft
        psd = _signal_qa.psd
        Phase = _signal_qa.Phase
        Admittance = _signal_qa.Admittance
        Coherence = _signal_qa.Coherence
        CrossSpec = _signal_qa.CrossSpec
        Metrics = _signal_qa.Metrics
        ppsd = _signal_qa.calc_ppsd
        smooth = _signal_qa.smooth

        preparetraces = _plot_helpers.preparetraces
        get_arrivals = _plot_helpers.get_arrivals

        spectrogram = _plotters.spectrogram
        plottrace = _plotters.plottrace
        ppsd_plot = _plotters.obs_ppsd_plot

        _window = _signal_helpers._window
        _calc_csd = _signal_helpers._calc_csd
        _psd = _signal_helpers._psd
        _calc_phase = _signal_helpers._calc_phase
        _calc_admittance = _signal_helpers._calc_admittance
        _calc_coherence = _signal_helpers._calc_coherence
        _csd_helper = _signal_helpers._csd_helper
        _stft = _signal_helpers._stft
        librosa_stft = _signal_helpers.librosa_stft
        _updatespec = _signal_helpers._updatespec