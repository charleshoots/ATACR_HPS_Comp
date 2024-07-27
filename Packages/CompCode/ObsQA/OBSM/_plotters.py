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
# ---------------------------------------------------------------------------------------------------------
# Plotters ------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
def spectrogram(self,r,fig=None,ax=None,num=0,window=None,overlap=None,plot=False,yscale='log',cmap='magma',vmin=None,vmax=None,ylabel=True,xlabel=True,cbar_on=True):
        # f,t,s = self._stft(self.traces[r][num].data,scaling='psd',window=window,overlap=overlap,return_onesided=True)
        f,t,s = self.librosa_stft(self.traces[r][num].data)
        f,t,s = f.T,t.T,s.T
        # s = _np.abs(s)
        if plot:
                s = _np.abs(s)
                if ax is None:
                        fig, cax = plt.subplots(nrows=2,ncols=1,figsize=(17,10),sharex='all',layout='constrained')
                        x = self.traces['Z'][0].times()
                        y = self.traces['Z'][0].data
                        gax = cax[0]
                        gax.plot(x,y)
                        gax.set_xlim(x[0],x[-1])
                        gax.set_xlim(t[0],t[-1])
                        gax = cax[1]
                else:
                        gax = ax
                pc = gax.pcolormesh(t,f, 10*_np.log10(s), cmap = cmap, shading= 'auto')
                if (vmin is not None) and (vmax is not None):
                        pc.set_clim(vmin, vmax)
                if ylabel:
                        gax.set_ylabel('Frequency (Hz)',fontweight='bold')
                if xlabel:
                        gax.set_xlabel('Time (s)',fontweight='bold')
                gax.set_yscale(yscale)
                gax.set_ylim(f[1],f[-1])
                gax.set_xlim(t[0],t[-1])
                if cbar_on:
                        fig.colorbar(pc, ax=gax, pad=0.01, label='dB')
        if plot:
                return pc,gax
        else:
                return f,t,s
def plottrace(self,ax=None,fig=None,component='Z',num=0,band=None,color='k',linewidth=0.6,label=None,st_llaz=None,ev_llaz=None,ttl=None,ylabel=False,lg=False,phases=('P','S','SKS'),figsize=(15,5),trim=(None,None),max_percentage=0.01,max_length=500,type='hann',norm=True):
        if ax is None:
                figs, axes = plt.subplots(nrows=1, ncols=1,figsize=figsize,height_ratios=[1],width_ratios=[1],layout='constrained',squeeze=False,sharey='row',sharex='col',lg=False)
                ax = axes[0,0]
                ttl = self.traces[component][num].stats.network + '.' + self.traces[component][num].stats.station
        if band is not None:
                fmin = _np.min(band)
                fmax = _np.max(band)
        # label = label.replace('Raw','PreFilt')
        if label is None:
                label = component + '|' + self.traces[component][0].stats.location
        else:
                label = label
        if band is not None:
                ax.set_title(label,fontweight='bold')
                if ylabel and fig is None:
                        ax.set_ylabel('Bandpassed ' + str(int(1/fmax)) + '-' + str(int(1/fmin)) + 's',fontweight='bold')
        else:
                ax.set_title(label,fontweight='bold')
        if (ttl is not None) and (fig is not None):
                figs.suptitle(ttl,fontweight='bold',y=1.02)

        tr = self.traces[component][0].copy()
        if band is not None:
                tr = self.preparetraces(tr,trim=trim,band=[fmin,fmax],max_percentage=max_percentage,max_length=max_length,type=type,norm=norm)
                # tr.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=4, zerophase=True)

        x,y = tr.times(),tr.data
        ax.plot(x,y,c=color,linewidth=linewidth,label=label)

        ax.set_xlabel('Time (s)',fontweight='bold')
        ax.set_xlim(x[0],x[-1])
        yl = _np.abs(ax.get_ylim()).max()
        ax.set_ylim(-yl,yl)
        if lg:
                lg = ax.legend(loc='upper right',prop={'weight':'bold'})
                [legobj.set_linewidth(2.0) for legobj in lg.legend_handles]

        if (st_llaz is not None) and (ev_llaz is not None):
                arrivals = self.get_arrivals(st_llaz,ev_llaz,phases=phases)
                if len(arrivals)>0:
                        p_phases = [phase_t for phase_t in arrivals if (phase_t[0][0].upper()=='P')]
                        s_phases = [phase_t for phase_t in arrivals if (phase_t[0][0].upper()=='S')]
                        [ax.axvline(phase_t[1],c='r',alpha=0.5,linewidth=0.3,ls='-') for phase_t in arrivals if len(arrivals)>0]
                        [ax.text(phase_t[1],-(0.98- (len(s_phases)-ti)*0.1)*yl,phase_t[0],c='b',fontweight='bold',alpha=0.5) for ti,phase_t in enumerate(s_phases) if (len(s_phases)>0) and (phase_t[0][0].upper()=='S')]
                        [ax.text(phase_t[1],(0.98 - (len(p_phases)-ti)*0.1)*yl,phase_t[0],c='b',fontweight='bold',alpha=0.5) for ti,phase_t in enumerate(p_phases) if (len(p_phases)>0) and (phase_t[0][0].upper()=='P')]
        return ax
def obs_ppsd_plot(self,comp='Z',
        # --------Post plot tweaks
        figsize = [15,8],
        eq_lines = ['gray',0.9,'-'],
        suptitle = None,
        average_line = ['k',2,'-.'],
        earthquakes = False,
        period = [1,200],
        show_coverage = False,
        # cmap = obspy.imaging.cm.viridis
        cmap = obspy.imaging.cm.pqlx,
        xaxis_frequency = False,
        append_title = None,
        NM = ['k',0.3]
        # --------
        ):
        ppsd = self.traces_ppsd[comp]
        if earthquakes==True:
                earthquakes = (5,7,3000)
        plot_prefs = AttribDict()
        if average_line:
                plot_prefs.show_mean = True
        else:
                plot_prefs.show_mean = False
        plot_prefs.show_coverage=show_coverage
        plot_prefs.show_earthquakes=earthquakes
        plot_prefs.period_lim=period
        plot_prefs.cmap=cmap
        plot_prefs.show=False
        plot_prefs.xaxis_frequency=xaxis_frequency
        # average_line = None
        hevent = ppsd.plot(**plot_prefs)
        _ = hevent.set_figwidth(figsize[0])
        _ = hevent.set_figheight(figsize[1])
        hevent.suptitle(suptitle)
        # Assuming 'hevent' is your figure object
        fig = hevent
        # Access the first axes in the figure (adjust the index if needed)
        ax = fig.axes[0]
        if earthquakes:
                texts = [text for text in ax.get_children() if isinstance(text,plt.Text)][:-3]
                [t.set_position([t.get_position()[0]*(2.1+0.0*ti),t.get_position()[1]*(0.911-0.013*ti)]) for ti,t in enumerate(texts)]
                [t.set_text(t.get_text().replace('\n',' - ')) for t in texts]
                [t.set_rotation(20) for t in texts]
        # Identify line objects
        lines = [line for line in ax.get_children() if isinstance(line, plt.Line2D)]
        stack = 0
        if average_line:
                # averaging line
                lines[0].set_color(average_line[0])  # Change color of the first line to red
                lines[0].set_linewidth(average_line[1])
                lines[0].set_linestyle(average_line[2])
        if plot_prefs.show_mean:
                stack +=1
        # noise model lines
        lines[stack].set_color(NM[0])  # Change color of the first line to red
        lines[stack].set_alpha(NM[1])
        lines[stack].set_linestyle(':')
        lines[stack+1].set_color(NM[0])  # Change color of the second line to blue
        lines[stack+1].set_alpha(NM[1])
        lines[stack+1].set_linestyle(':')

        stack+=2
        if plot_prefs.show_earthquakes:
                _ = [l.set_color(eq_lines[0]) for l in lines[stack:]]
                _ = [l.set_linewidth(eq_lines[1]) for l in lines[stack:]]
                _ = [l.set_linestyle(eq_lines[2]) for l in lines[stack:]]
        N = int(len(ppsd.times_processed))
        NS = int(UTCDateTime(ppsd.times_processed[2]) - UTCDateTime(ppsd.times_processed[0]))
        ttl = fig.axes[0].get_title().replace('  ',' | ')
        ttl = ttl.replace('(','').replace(')','')
        ttl = ttl.replace('{N}/{N}'.format(N=N),str(N))
        ttl = ttl.replace('segments', '({NS}s) segments'.format(NS=NS))
        if append_title:
                ttl = ttl + '\n' + append_title
        fig.axes[0].set_title(ttl)
        fig.canvas.draw()
        plt.close('all')
        # Show the updated figure (if not already displayed)
        return fig
