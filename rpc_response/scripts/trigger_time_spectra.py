"""Standalone multistrip, shaped-waveform -> 20% timing -> position study.
No Geant4 rerun. Independent uniform positions; same existing avalanche template.
References R001/R007/R014/R015; 20% fraction is a user setting.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
from strip_charge_sharing import fractions
from electrode_signals import deposit
from trigger_readout import trigger,coordinate
args=argparse.ArgumentParser();args.add_argument('--events',type=int,default=5000);args.add_argument('--dt-ns',type=float,default=.1);args.add_argument('--output',default='deliverables/Trigger20_Spectra_20261031');args=args.parse_args()
if args.events<1 or args.dt_ns<=0:raise ValueError('Invalid scan size')
R=Path(__file__).resolve().parents[2];P=R/args.output;P.mkdir(parents=True,exist_ok=True)
g=json.loads((R/'rpc_response/config/pcb_geometry.json').read_text(encoding='utf-8'));trace=pd.read_csv(R/'rpc_response/examples/single_electron_trace.csv')
centres=[]
for axis in ['X','Y']:
 c=[]
 for n in range(1,82):
  if axis=='X':
   a=[q for q in g['tracks'] if q['net']==str(n) and abs(q['y2_mm']-q['y1_mm'])>250 and abs(q['x2_mm']-q['x1_mm'])<1e-6];c.append(np.mean([.5*(q['x1_mm']+q['x2_mm']) for q in a]))
  else:
   a=[q for q in g['fills'] if q['net']=='Y'+str(n)];c.append(np.average([.5*(q['y1_mm']+q['y2_mm']) for q in a],weights=[abs((q['x2_mm']-q['x1_mm'])*(q['y2_mm']-q['y1_mm'])) for q in a]))
 centres.append(np.array(c));assert np.all(np.diff(c)>0)
dt=args.dt_ns;t=(np.arange(int(round(440/dt)))+.5)*dt;base=np.zeros(int(round(100/dt)))
for r in trace.itertuples(index=False):deposit(base,r.t0_ns,r.t1_ns,r.current_nA*(r.t1_ns-r.t0_ns),dt)
ht=np.arange(0,60+dt/2,dt);h=(ht/5)*np.exp(1-ht/5);v=fftconvolve(base,h)*dt;vt=(np.arange(len(v))+.5)*dt
bank=np.array([np.interp(t-4*k,vt,v,left=0,right=0) for k in range(81)])
sigma=4/(2*np.sqrt(2*np.log(2)));rng=np.random.default_rng(20261031);N=args.events
truth=rng.uniform(-100,100,size=(N,2));records=[{'event':i,'t0_ns':0.,'x_true_mm':p[0],'y_true_mm':p[1]} for i,p in enumerate(truth)]
for ai,axis in enumerate(['X','Y']):
 c=centres[ai]
 for first in range(0,N,128):
  end=min(N,first+128);weights=np.array([fractions(c,p+c[40],sigma) for p in truth[first:end,ai]])
  waves=[weights@bank*.25,weights@bank[::-1]*.25]
  for j,i in enumerate(range(first,end)):
   valid=True
   for side,w in zip(['1','81'],waves):
    result=trigger(t,w[j],.2);records[i][f'T_{axis}{side}_ns']=result['time_ns'];records[i][f'{axis}{side}_valid']=result['valid'];records[i][f'{axis}{side}_status']=result['status'];valid &=result['valid']
   records[i][axis+'_valid']=valid
   records[i][axis.lower()+'_rec_mm']=coordinate(records[i][f'T_{axis}1_ns'],records[i][f'T_{axis}81_ns'],c) if valid else float('nan')
df=pd.DataFrame(records);df.to_csv(P/'trigger_times_and_positions.csv',index=False)
plt.rcParams.update({'font.size':13,'axes.linewidth':1.2,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'legend.frameon':False})
fig,axes=plt.subplots(1,2,figsize=(13,5.3))
summary={}
for ax,axis in zip(axes,['X','Y']):
 ok=df[axis+'_valid'];a=df.loc[ok,f'T_{axis}1_ns'];b=df.loc[ok,f'T_{axis}81_ns']
 hist=ax.hist2d(a,b,bins=[np.linspace(20,320,151)]*2,cmin=1,cmap='viridis');fig.colorbar(hist[3],ax=ax,label='Events / (2 ns x 2 ns)')
 ax.set(xlabel=f'T({axis} net1), 20% crossing (ns)',ylabel=f'T({axis} net81), 20% crossing (ns)',title=f'{axis}: two-end trigger-time spectrum');ax.set_aspect('equal');ax.minorticks_on()
 residual=df.loc[ok,axis.lower()+'_rec_mm']-df.loc[ok,axis.lower()+'_true_mm']
 summary[axis]={'valid_events':int(ok.sum()),'invalid_events':int((~ok).sum()),'mean_bias_mm':float(residual.mean()),'rms_error_mm':float(np.sqrt(np.mean(residual**2))),'max_abs_error_mm':float(residual.abs().max()),'time_sum_mean_ns':float((a+b).mean())}
fig.suptitle('Multi-strip readout, 5 ns shaping, 20%-of-own-peak leading edge; t0 = 0',fontsize=15);fig.tight_layout();fig.savefig(P/'trigger_time_2d.png',dpi=190);fig.savefig(P/'trigger_time_2d.pdf');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,4.6))
for ax,axis,col in zip(axes,['X','Y'],['blue','red']):
 ok=df[axis+'_valid'];ax.scatter(df.loc[ok,axis.lower()+'_true_mm'],df.loc[ok,axis.lower()+'_rec_mm']-df.loc[ok,axis.lower()+'_true_mm'],s=1,c=col);ax.set(xlabel=f'True {axis} relative to centre (mm)',ylabel=f'Reconstructed {axis} - true {axis} (mm)',title='20% timing coordinate bias (no noise)');ax.minorticks_on()
fig.tight_layout();fig.savefig(P/'position_bias.png',dpi=180)
meta={'events':N,'seed':20261031,'threshold_fraction':.2,'definition':'first rising 20%-of-own-positive-peak crossing, linearly interpolated; multi-peak ambiguity flagged','dt_ns':dt,'tau_ns':5,'charge_profile_fwhm_mm':4,'position_source':'uniform independent x,y in [-100,100] mm, not Geant4 or Gaisser-Tang','t0_ns':0,'noise':'none; measured sigma not applied','avalanche':'same existing single-electron template for all events, not independent full muon avalanches','coordinate':'tap=(T1-T81+320)/8; interpolate actual PCB centres; subtract centre channel coordinate','references':['R001','R007','R014','R015'],'summary':summary}
(P/'summary.json').write_text(json.dumps(meta,indent=2),encoding='utf-8');print(json.dumps(summary))
