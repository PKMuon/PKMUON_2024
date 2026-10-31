"""1GeV Heed events -> sampled linear avalanche histories -> strip timing.
Requires outputs from muon_waveform_reference.cc. Never re-labels guard failures
as valid detector response. No spatial sigma or artificial time jitter applied.
"""
from pathlib import Path
import json,hashlib,argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
from strip_charge_sharing import fractions
from trigger_readout import trigger,coordinate
args=argparse.ArgumentParser();args.add_argument('--waveforms',default='rpc_response/results/muon_trigger20_highmarker_20261031');args.add_argument('--output',default='deliverables/Muon1GeV_Trigger20_HighMarker_20261031');args=args.parse_args()
R=Path(__file__).resolve().parents[2];S=R/'rpc_response/results';W=R/args.waveforms;P=R/args.output;P.mkdir(parents=True,exist_ok=True)
ev=pd.read_csv(S/'muon_1GeV_normal_v1_events.csv');status=pd.read_csv(W/'currents_status.csv');N=len(ev);assert np.array_equal(ev.event,status.event)
status=status.rename(columns={'n_outside_linear_guard':'n_above_population_marker'})
generation=json.loads((W/'currents_generation.json').read_text()) if (W/'currents_generation.json').exists() else {'population_marker':10000000}
marker=generation['population_marker']
curr=np.fromfile(W/'currents.f64',dtype='<f8').reshape(N,600);dt=.1;shape_t=np.arange(601)*dt;h=(shape_t/5)*np.exp(1-shape_t/5)
g=json.loads((R/'rpc_response/config/pcb_geometry.json').read_text(encoding='utf-8'));centres={}
for axis in ['X','Y']:
 c=[]
 for n in range(1,82):
  if axis=='X':
   a=[q for q in g['tracks'] if q['net']==str(n) and abs(q['y2_mm']-q['y1_mm'])>250 and abs(q['x2_mm']-q['x1_mm'])<1e-6];c.append(np.mean([.5*(q['x1_mm']+q['x2_mm']) for q in a]))
  else:
   a=[q for q in g['fills'] if q['net']=='Y'+str(n)];c.append(np.average([.5*(q['y1_mm']+q['y2_mm']) for q in a],weights=[abs((q['x2_mm']-q['x1_mm'])*(q['y2_mm']-q['y1_mm'])) for q in a]))
 centres[axis]=np.asarray(c);assert np.all(np.diff(c)>0)
t=(np.arange(4500)+.5)*dt;rows=[];example={};sigma=4/2.354820045
for i,r in enumerate(ev.itertuples(index=False)):
 local=fftconvolve(curr[i],h)*dt;scale=np.max(local)
 if scale>0:local=local/scale # single common event scaling for numerical robustness, does not alter fraction times
 row={'max_family_population':int(status.iloc[i].max_family_population) if 'max_family_population' in status else None,'event':r.event,'energy_GeV':1.,'x_true_mm':r.x_true_mm,'y_true_mm':r.y_true_mm,'n_electrons':int(status.iloc[i].n_electrons),'n_above_population_marker':int(status.iloc[i].n_above_population_marker),'below_population_marker':bool(status.iloc[i].n_above_population_marker==0),'formal_electron_charge_fC':status.iloc[i].charge_fC}
 for axis,pos in [('X',r.x_true_mm),('Y',r.y_true_mm)]:
  c=centres[axis];weights=fractions(c,pos+c[40],sigma);times=[];ok=True
  for side in ['1','81']:
   wave=np.zeros(len(t))
   for k,w in enumerate(weights):
    if w==0:continue
    start=(k if side=='1' else 80-k)*40;wave[start:start+len(local)]+=local*w*.25
   result=trigger(t,wave,.2);row['T_'+axis+side+'_ns']=result['time_ns'];row[axis+side+'_status']=result['status'];ok &=result['valid'];times.append(result['time_ns'])
   if i==0:example[axis+'_net'+side]=wave
  row[axis+'_timing_valid']=ok;row[axis.lower()+'_rec_mm']=coordinate(*times,c) if ok else float('nan')
 rows.append(row)
df=pd.DataFrame(rows);df.to_csv(P/'muon_four_trigger_times.csv',index=False);pd.DataFrame({'time_ns':t,**example}).to_csv(P/'event0_shaped_waveforms.csv',index=False)
plt.rcParams.update({'font.size':12,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'axes.linewidth':1.2})
for guarded in [False,True]:
 fig,axes=plt.subplots(1,2,figsize=(13,5.3))
 for ax,axis in zip(axes,['X','Y']):
  ok=df[axis+'_timing_valid']
  if guarded:ok=ok&df.below_population_marker
  hist=ax.hist2d(df.loc[ok,'T_'+axis+'1_ns'],df.loc[ok,'T_'+axis+'81_ns'],bins=[np.linspace(20,340,161)]*2,cmin=1,cmap='viridis');fig.colorbar(hist[3],ax=ax,label='Events / (2 ns x 2 ns)')
  ax.set(xlabel=f'T({axis} net1), 20% crossing (ns)',ylabel=f'T({axis} net81), 20% crossing (ns)',title=f'{axis} direction: {int(ok.sum())} timed events');ax.set_aspect('equal');ax.minorticks_on()
 title=f'Below user population marker {marker:.0e}' if guarded else f'Unsaturated linear reference; user marker {marker:.0e}, no population-based rejection'
 fig.suptitle('1 GeV mu- Heed sample: '+title,fontsize=14);fig.tight_layout();name='below_marker_time_2d' if guarded else 'muon1GeV_time_2d';fig.savefig(P/(name+'.png'),dpi=190);fig.savefig(P/(name+'.pdf'));plt.close(fig)
summary={'source':'existing standalone TrackHeed mu- sample, kinetic energy 1GeV, normal +Z incidence, uniform x,y[-100,100]mm, t0=0','source_events':N,'ionization_seed':20261019,'avalanche_seed':20261031,'local_current':'sum of all starting electrons at actual Heed birth times and drift distances; independent linear birth/death histories, 0.01mm steps','new_run':'new avalanche and waveform digitization; reused existing Heed ionization, not a new Geant4 run','below_marker_events':int(df.below_population_marker.sum()),'above_marker_events':int((~df.below_population_marker).sum()),'population_marker':marker,'marker_meaning':'user-selected diagnostic marker, not a physical saturation/validity boundary; no clipping or resampling; raising marker does not remove space-charge limitations','timing_valid':{a:int(df[a+'_timing_valid'].sum()) for a in ['X','Y']},'assumptions':['4mm FWHM Gaussian sharing assigned around incident muon x/y; individual delta-electron transverse offsets not available in source CSV','5ns preset unipolar shaping','equal X/Y and endpoint splitting','no space charge/saturation, noise or empirical coordinate sigma'],'threshold_fraction':.2,'dt_ns':dt,'references':['R001','R002','R006','R007','R014','R015'],'source_sha256':hashlib.sha256((S/'muon_1GeV_normal_v1_electrons.csv').read_bytes()).hexdigest()}
(P/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8');print(json.dumps(summary))
