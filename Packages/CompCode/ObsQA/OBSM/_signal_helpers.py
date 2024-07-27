# ---------------------------------------------------------------------------------------------------------
# Signal helpers ------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
import numpy as _np
from scipy.signal import stft, detrend

def _window(self,window=None,overlap=None):
        
        if overlap is None:
                overlap=self.overlap
        if window is None:
                window=self.window
        # Points in window
        ws = int(window/self.dt)
        # Number of points to overlap
        ss = int(window*overlap/self.dt)
        # hanning window
        hanning = _np.hanning(2*ss)
        wind = _np.ones(ws)
        wind[0:ss] = hanning[0:ss]
        wind[-ss:ws] = hanning[ss:ws]
        return wind,ws,ss
def _calc_csd(self,a_ft,b_ft):
        cab = a_ft*_np.conj(b_ft)
        return cab
def _psd(self,tr,window=None,overlap=None):
        f,_t,ft = self._stft(tr,scaling='psd',return_onesided=True,window=window,overlap=overlap)
        psd = _np.abs(ft)**2
        return f,psd
def _calc_phase(self,ab):
        ph = _np.angle(ab,deg=True)
        return ph
def _calc_admittance(self,ab,bb):
        ad = _np.abs(ab)/bb
        return ad
def _calc_coherence(self,ab,aa,bb):
        coh = ((abs(ab)**2)/(aa*bb))
        return coh
def _csd_helper(self,a,b,window=None,overlap=None):
        f,_t,a_ft = self._stft(a,window=window,overlap=overlap)
        f,_t,b_ft = self._stft(b,window=window,overlap=overlap)
        cab = _np.mean(self._calc_csd(a_ft,b_ft),axis=0)
        return f,cab
def _stft(self,tr,scaling='spectrum',window=None,overlap=None,return_onesided=False):
        wind,ws,ss = self._window(window=window,overlap=overlap)
        _f, _t, ft = stft(tr,self.fs, return_onesided=return_onesided, boundary=None,padded=False, window=wind, nperseg=ws, noverlap=ss,detrend='constant',scaling=scaling)
        _f = _f.reshape(-1)
        ft = _np.atleast_2d(ft)
        ft = ft*ws
        # ft = _np.mean(ft,axis=-1)
        return _f.T,_t.T,ft.T
def librosa_stft(self,x,window=120,hop_length_fraction=4):
        win_length_samples = _np.int(_np.round(2**_np.ceil(_np.log2(window*self.fs))))
        hop_length=win_length_samples // hop_length_fraction
        ft, phase = librosa.magphase(librosa.stft(x,n_fft=win_length_samples,hop_length=hop_length,win_length=win_length_samples))
        df = self.fs/win_length_samples
        _f = _np.arange(ft.shape[0]) * df
        _t = _np.arange(ft.shape[1]) * hop_length
        _t = _t/self.fs #Time axis in samples d
        return _f.T,_t.T,ft.T
def _updatespec(self,window=None,overlap=None):
        self.csd = dict()
        self.csd['A'] = dict()
        self.csd['B'] = dict()
        self.csd['AB'] = dict()
        for p in self.csd_pairs:
                self.csd['A'][p] = []
                self.csd['B'][p] = []
                self.csd['AB'][p] = []
        for i in range(len(self.traces.select(channel='*Z'))):
                for p in self.csd_pairs:
                        pp = p.replace('P','H')
                        A = self.A.select(channel='*'+pp[0])[i].copy()
                        B = self.B.select(channel='*'+pp[1])[i].copy()
                        f,spec_AB= self._csd_helper(A.data,B.data,window=window,overlap=overlap)
                        self.csd['AB'][p].append(spec_AB)
                        A = self.A.select(channel='*'+pp[0])[i].copy()
                        B = self.A.select(channel='*'+pp[1])[i].copy()
                        f,spec_AB= self._csd_helper(A.data,B.data,window=window,overlap=overlap)
                        self.csd['A'][p].append(spec_AB)
                        A = self.B.select(channel='*'+pp[0])[i].copy()
                        B = self.B.select(channel='*'+pp[1])[i].copy()
                        f,spec_AB= self._csd_helper(A.data,B.data,window=window,overlap=overlap)
                        self.csd['B'][p].append(spec_AB)
        self.f = f
