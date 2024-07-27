# ===================================================================================================
# ============================================ IMPORTS ==============================================
# ===================================================================================================
import scipy
import numpy as np
import os
import matplotlib.pyplot as plt
import sys
import obspy
import pickle as pkl
import glob as g
from obspy.clients.fdsn import Client
import datetime
import re
import math
import matplotlib.colors as mcolors
import matplotlib.cm as cm2
from scipy.stats import norm
import scipy.stats as stats
import pandas as pd
from pathlib import Path
def save_tight(filename,fig=None,dpi=200):
        plt.margins(0.1,0.1)
        if fig is None:
                plt.savefig(filename, bbox_inches = 'tight',pad_inches = 0.05,dpi=dpi)
        else:
                fig.savefig(filename, bbox_inches = 'tight',pad_inches = 0.05,dpi=dpi)
        return 'Complete'
# ===================================================================================================
# ===================================================================================================
# ===================================================================================================
#
#
#
catalog = pd.read_pickle('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python/EVENTS/sta_catalog_proxima.pkl')
D = pd.read_excel('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/METHODS/ATaCR/ATaCR_Python/utilities/Janiszewski_etal_2023_StationList.xlsx')
Dstations = [str(e) for e in D.Station]
DExp = np.unique(['[' + str(n) + '] ' + ' ' + str(x) + ' ' for n,x in zip(D.Network,D.Experiment)])
DStart = [str(e) for e in D.Start]
DEnd = [str(e) for e in D.End]
n = 1
for i,s in enumerate(Dstations):
    s = ' '*(n*(15)) + s + ' '*(np.abs(n-1)*(15))
    n = np.abs(n-1)
    Dstations[i] = s
client = Client()
evts = client.get_events(starttime=np.min(D.Start), endtime=np.max(D.End),minmagnitude=6.0,maxmagnitude=7.0)
evtimes = [d.origins[0].time for d in evts]
days = np.array([np.min(D.Start) + datetime.timedelta(days=i) for i in range((np.max(D.End) - np.min(D.Start) + datetime.timedelta(seconds=1)).days)])
counts = np.zeros(days.shape)
for ev in evtimes:
    counts[np.where(days==datetime.datetime(year=ev.year,day=ev.day,month=ev.month))[0][0]] +=int(1)
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib as mpl
from matplotlib.ticker import AutoLocator
import matplotlib.dates as mdates
# from matplotlib import cm
def tmp_time_cdf(s_start,s_end,evtimes):
    e = [evts[i] for i in np.where(np.logical_and((s_end  > np.array(evtimes)),(s_start < np.array(evtimes))))[0]]
    e_times = [k.origins[0].time for k in e]
    e_times = [datetime.datetime(year=e_times[k].year,month=e_times[k].month,day=e_times[k].day,hour=e_times[k].hour,minute=e_times[k].minute) for k in range(len(e_times))]
    to_timestamp = np.vectorize(lambda x: x.timestamp())
    time_stamps = to_timestamp(e_times)
    hist,edges = np.histogram(time_stamps,density=False,bins=len(e_times))
    histc = hist
    # hist = hist / np.sum(hist)
    days = np.array([s_start + datetime.timedelta(days=i) for i in range((s_end - s_start).days+1)])
    prob = scipy.interpolate.interp1d(edges[0:-1], hist,kind='nearest',fill_value=0,bounds_error=False)
    prob = prob(to_timestamp(days))
    edges = [datetime.date.fromtimestamp(j) for j in edges][0:-1]
    days = days[prob>0]
    prob = prob[prob>0]
    return days,prob
def day_counts(s_start,s_end,days,counts,min=0,max=1000):
    s_start = datetime.datetime(year=s_start.year,day=s_start.day,month=s_start.month)
    s_end = datetime.datetime(year=s_end.year,day=s_end.day,month=s_end.month)
    d=days[np.where(np.logical_and(days>s_start,days<s_end))[0]]
    c=counts[np.where(np.logical_and(days>s_start,days<s_end))[0]]
    d = d[c>0]
    c = c[c>0]
    d = d[c>min]
    c = c[c>min]
    d = d[c<max]
    c = c[c<max]
    d = d[np.argsort(c)]
    c = c[np.argsort(c)]
    return d,c
fig, axes = plt.subplots(nrows=1, ncols=2,figsize=(45,45),width_ratios=[0.2,1],layout='constrained',squeeze=False,sharey='all')
ax = axes[0][1]
lk = []
yt = [(i*100) for i in range(len(Dstations))]
expyt = [yt[np.max(np.where(D.Network==j))]+1 for j in np.unique(D.Network)]
cyt = []
norm = mpl.colors.Normalize(vmin=1, vmax=int(np.max(counts[counts<np.max(counts)])))
for s in catalog.Station:
    k = (np.where(s ==D.Station)[0])
    if len(k)==0:
        continue
    elif len(k)==1:
        cyt.append(yt[k[0]])
[plt.axhline(y,color='r',linewidth=1.5,ls=':',zorder=-1,alpha=0.75) for y in cyt]
[plt.axhline(j,ls='-',c='w',zorder=-2) for j in expyt]
colors = [(0, .2, .6),(0, 1, 0)]  # R -> G -> B
n = len(np.unique(counts[counts<np.max(counts)]))
cmap = LinearSegmentedColormap.from_list('my_list', colors, N=100).resampled(len(np.unique(counts[counts<np.max(counts)])[1:]))
c = cmap(np.linspace(0,1,3)).tolist()
c.append((1, 0, 0))
cmap = LinearSegmentedColormap.from_list('my_list', c, N=100).resampled(len(np.unique(counts[counts<np.max(counts)])[1:]))
anomaly_bar = np.ceil(np.std(counts[counts>0])*2) + np.ceil(np.mean(counts[counts>0]))
for ii,(s_start,s_end,yy) in enumerate(zip(D.Start,D.End,yt)):
    for min,max in zip([0,6],[6,100]):
        dys,cts = day_counts(s_start,s_end,days,counts,min=min,max=max)
        if len(cts)>0:
            lk.append(np.max(cts))
        else:
            continue
        if np.max(cts)>anomaly_bar:
            plt.scatter(dys,cts*0+yy,c=cts,s=30, marker='s',cmap=cmap,norm=norm,edgecolor='k',linewidth=0.5)
        else:
            plt.scatter(dys,cts*0+yy,c=cts,s=20, marker='s',cmap=cmap,norm=norm,edgecolor='k',linewidth=0.5)
plt.yticks(ticks = yt, labels = Dstations,fontsize=7,fontweight='bold')
ax.xaxis.set_major_locator(AutoLocator())
ax.grid(visible=True,axis='x')
plt.xlim(np.min(D.Start),np.max(D.End))
xx = np.min(plt.gca().get_xlim()) + 50
v = [plt.text(xx,i,j,fontweight='bold',bbox=dict(facecolor='w', alpha=1),fontsize=9) for i,j in zip(expyt,DExp)]
yl = plt.gca().get_ylim()
plt.ylim(np.min(yt),np.max(yt)*1.01)
plt.gca().set_facecolor('dimgray')
plt.tick_params(labelbottom=True, labeltop=True, labelleft=False, labelright=False,bottom=True, top=True, left=False, right=False)
ax.set_xticklabels(labels=ax.get_xticklabels(),fontsize=10,fontweight='bold')
dtFmt = mdates.DateFormatter('%y/%m') # define the formatting
ax.xaxis.set_major_formatter(dtFmt)
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=4))
cb = plt.colorbar(location='top',pad=-0.02,fraction=0.15,aspect=80)
plt.title('OBS Experiments Deployment Coverage\n Colored by # of m6.0-7.0 events',fontweight='bold',y=1.02,fontsize=14)
plt.clim(vmin=1,vmax=int(np.max(counts[counts<np.max(counts)])))
ax = axes[0][0]
ax.xaxis.set_major_locator(AutoLocator())
x = D['Water Depth (m)']
ax.tick_params(labeltop=True, top=True, labelbottom=True, labelleft=False, labelright=True,bottom=True,left=True, right=True)
ax.barh(yt,x,height=60,color='k')
xlab = [str(int(e)) for e in ax.get_xticks()]
xt = [int(e) for e in ax.get_xticks()]
ax.set_yticks(ticks = yt, labels = Dstations,fontsize=7,fontweight='bold')
ax.set_xticklabels(labels=xlab,fontsize=10,fontweight='bold')
ax.set_xticks(xt)
ax.set_ylim(np.min(yt),np.max(yt)*1.01)
ax.set_xlim(np.min(x),np.max(x)*1.01)
ax.grid(visible=True)
ax.set_xlabel('OBS Depth (m)',fontweight='bold',fontsize=14)
ax.xaxis.set_label_coords(0.5,1.01)
save_tight('/Users/charlesh/Documents/Codes/OBS_Methods/NOISE/COMPS/FigureArchive/_GEN4/DeploymentCoverage.png',dpi=200)