"""Single-hit four-end plot in the user's supplied reference style.
R001 shaping; R007/R014/R015 nominal delay-line routing. Preset tau=5 ns.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
from strip_delay_readout import route
R=Path(__file__).resolve().parents[2];P=R/'deliverables/Single_Hit_Four_Ends_20261031';P.mkdir(parents=True,exist_ok=True)
trace=pd.read_csv(R/'rpc_response/examples/single_electron_trace.csv');top=json.loads((R/'rpc_response/config/pcb_delay_topology.json').read_text())
rows=[]
for net in ['11','Y21']:
 for r in trace.itertuples(index=False):rows.append(dict(net=net,t0_ns=r.t0_ns,t1_ns=r.t1_ns,current_nA=r.current_nA*.5))
raw=route(pd.DataFrame(rows),top,dt=.02,end=400);dt=.02;t0=raw['t_center_ns'][0];ht=np.arange(0,60+dt/2,dt);h=(ht/5)*np.exp(1-ht/5)
spec=[('X_net1','X1',40,'#0000ff'),('X_net81','X2',280,'#ff0000'),('Y_net1','Y1',80,'#008000'),('Y_net81','Y2',240,'#cc00cc')]
v={k:fftconvolve(raw[k],h,mode='full')*dt for k,*_ in spec};scale=max(np.max(a) for a in v.values());t=t0+np.arange(len(next(iter(v.values()))))*dt
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.linewidth':1.,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'legend.frameon':False})
fig,ax=plt.subplots(figsize=(9,4.9))
for k,label,delay,color in spec:ax.plot(t,v[k]/scale,c=color,lw=1.4,label=f'{label}, delay {delay:.0f} ns')
ax.set(xlabel='Time since electron birth (ns)',ylabel='Endpoint output / single-pulse peak',xlim=(20,350),ylim=(-.04,1.05));ax.minorticks_on();ax.tick_params(which='both',top=True,right=True);ax.legend(loc='upper center',bbox_to_anchor=(.53,1),fontsize=10)
fig.tight_layout();fig.savefig(P/'single_hit_four_ends.png',dpi=220);fig.savefig(P/'single_hit_four_ends.pdf')
pd.DataFrame({'time_ns':t,**{k:a/scale for k,a in v.items()}}).to_csv(P/'single_hit_four_ends.csv',index=False)
(P/'metadata.json').write_text(json.dumps({'hit':'A only: X net11 and Y net21; t0=0 ns','shape_tau_ns':5,'delay_ns':{k:delay for k,label,delay,col in spec},'normalization':'one common single-pulse peak','endpoint_names':'X1/X2,Y1/Y2 denote PCB network1/network81, not experimentally confirmed L/R/U/D','source':'existing single-electron 1mm linear avalanche template, seed20261021','coupling':'preset equal X/Y partition and equal endpoint split','references':['R001','R007','R014','R015'],'status':'simplified connection and preset shaping demonstration, uncalibrated voltage'},indent=2),encoding='utf-8')
print(P)
