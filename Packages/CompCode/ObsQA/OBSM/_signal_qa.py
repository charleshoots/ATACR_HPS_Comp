# ---------------------------------------------------------------------------------------------------------
# Signal QA -----------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
import numpy as _np
from obspy import Stream
from obspy.core.util.attribdict import AttribDict
from obspy import read
from obspy.io.xseed import Parser
from obspy.signal import PPSD
def ft(self,r,window=None,overlap=None):
        tr = self.traces[r]
        ft = []
        for tr in self.traces[r]:
                f,_t,ft0 = self._stft(tr,window=window,overlap=overlap)
                ft.append(_np.mean(ft0,axis=-1))
        return f,_np.array(ft).squeeze()
def psd(self,r,window=None,overlap=None):
        tr = self.traces[r]
        ft = []
        f = self.f
        for tr in self.traces[r]:
                f,ft0 = self._psd(tr,window=window,overlap=overlap)
                ft.append(ft0)
        ft = _np.mean(_np.array(ft).squeeze(),axis=0)
        ft = ft[f>=0]
        f = f[f>=0]
        return f,ft
def Phase(self,r,s=False):
        '''
        Takes a single NDarray containing the cross-psd (complex) between two components and pulls the phase out of the imaginary component of the array. Output is in degrees for a quadrature circle (-180,180).
        '''
        # if r[0].upper=='Z':
        #         r = r[::-1]
        f = self.f[self.f>=0]
        ph = [self._calc_phase(ab) for ab in self.csd['AB'][r]]
        ph = _np.real(_np.array(ph).squeeze()[self.f>=0])
        if s:
                ph = self.smooth(ph)
        ph = ph[f>=0]
        f = f[f>=0]
        return f,ph
def Admittance(self,r,s=False):
        '''
        Takes two NDarrays of equal shape, the first (a) containing the cross-power-spectral density between the two components and the second (b) containing the auto-power-spectral density of the secondary component.
        ie, For the spectral admittance between Z and P, a is the cross psd between the two and b is the auto-psd of P.
        '''
        # r = ''.join(sorted(r))
        # if r[0].upper=='Z':
        #         r = r[::-1]
        f = self.f[self.f>=0]
        adm = [self._calc_admittance(self.csd['AB'][r][i],self.csd['B'][r[1] + r[1]][i]) for i in range(len(self.csd['AB'][r]))]
        adm = _np.real(_np.array(adm).squeeze()[self.f>=0])
        if s:
                adm = self.smooth(adm)
        adm = adm[f>=0]
        f = f[f>=0]
        return f,adm
def Coherence(self,r,s=False,plot=False,lw=1,ax=None,xlabel=None,ylabel=None,label=None,c='k',ls='-',alpha=0.3,ttl=None):
        '''
        Takes three NDarrays of equal shape, the first (ab) containing the cross-power-spectral density between the two components, the second and third (aa and bb) containing the auto-power-spectral density of these two components.
        ie, For the spectral coherence between Z and P, a is the cross psd between the two and b and c are the auto-psd of Z and P, respectively.
        '''
        # r = ''.join(sorted(r))
        # if r[0].upper=='Z':
        #         r = r[::-1]
        f = self.f[self.f>=0]
        coh = [self._calc_coherence(   self.csd['AB'][r][i],    self.csd['A'][r[0] + r[0]][i],   self.csd['B'][r[1] + r[1]][i])     for i in range(len(self.csd['AB'][r]))]
        coh = _np.abs(_np.array(coh).squeeze()[self.f>=0])
        if s:
                # coh = _np.array([self.smooth(c) for c in coh])
                coh = self.smooth(coh)
        coh = coh[f>=0]
        f = f[f>=0]
        if plot:
                if ax is None:
                        fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(10,10),height_ratios=[1],width_ratios=[1],layout='constrained',squeeze=False,sharey='row',sharex='col')
                        ax = axes[0,0]
                gax = ax
                gax.plot(f,coh,label=label,c=c,alpha=alpha,linewidth=lw,ls=ls)
                if ttl is not None:
                        gax.set_title(ttl,fontweight='bold')
                if ylabel is not None:
                        gax.set_ylabel(ylabel + ' Coherence',fontweight='bold')
                if xlabel:
                        gax.set_xlabel('Hz',fontweight='bold')
                gax.set_xscale('log')
                gax.set_ylim(0,1)
                gax.set_xlim(f[1],f[-1])
                return ax
        else:
                return f,coh
def CrossSpec(self,A,B,fs=None,return_onesided=False,window=None):
        if isinstance(A,obspy.core.trace.Trace):
                fs = 1/A.stats.delta
                A = A.data
        if isinstance(B,obspy.core.trace.Trace):
                fs = 1/B.stats.delta
                B = B.data
        f,spec_AB= self._csd_helper(A,B,window=window)
        f,spec_AA= self._csd_helper(A,A,window=window)
        f,spec_BB= self._csd_helper(B,B,window=window)
        self.f = f
        return f,spec_AB,_np.abs(spec_AA),_np.abs(spec_BB)
def Metrics(self,r):
        _,coh = self.Coherence(r)
        _,adm = self.Admittance(r)
        _,ph = self.Phase(r)
        return self.f,coh,adm,ph
def calc_ppsd(self,verbose=False,ppsd_length=None,db_bins=(-200, -50, 1.0),overlap=0.5,skip_on_gaps=True,rebuild=False):
        if (len(self.traces.keys())==len(self.traces_ppsd.keys())) and not rebuild:
                pass
        else:
                self.traces_ppsd = dict() 
                for key in list(self.traces.keys()):
                        st = Stream(self.traces[key])
                        if not ppsd_length:
                                fs = 1/st[0].stats.delta
                                # ppsd_length = int(round(2**np.ceil(np.log2(120*(fs))))) / st[0].stats.delta
                                ppsd_length = int(round(2**_np.ceil(_np.log2(26*(fs))))) / st[0].stats.delta
                        metadata = {'gain':1., 'sensitivity': 1., 'poles':[1-1j], 'zeros': [1-1j]}
                        stats = st[0].stats
                        ppsd_prefs = AttribDict()
                        ppsd_prefs.skip_on_gaps=skip_on_gaps
                        ppsd_prefs.ppsd_length=ppsd_length
                        ppsd_prefs.db_bins=db_bins
                        ppsd_prefs.overlap=overlap
                        ppsd = PPSD(stats = stats,metadata = metadata,**ppsd_prefs)
                        success = ppsd.add(st,verbose=verbose)
                        if success:
                                print(key + ' | Number of psd segments:', len(ppsd.times_processed))
                        else:
                                print('PPSD Failed')
                        self.traces_ppsd[key] = ppsd
        return self
def smooth(self,d,k=10):
        return _np.convolve(d, _np.ones(k) / k, mode='same')