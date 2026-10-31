"""Revised multi-strip signal demo. R001/R007/R014/R015 in REFERENCES.md."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
from strip_charge_sharing import fractions
from strip_delay_readout import route
R=Path(__file__).resolve().parents[2];P=R/'deliverables/Shared_Strip_Four_Ends_20261031';P.mkdir(parents=True,exist_ok=True)
g=json.loads((R/'config/pcb_geometry.json').read_text(encoding='utf-8'));top=json.loads((R/'config/pcb_delay_topology.json').read_text(encoding='utf-8'))
trace=pd.read_csv(R/'rpc_response/examples/single_electron_trace.csv');dt=.02;tau=5.;fwhm=4.;sigma=fwhm/(2*np.sqrt(2*np.log(2)))
# Extract actual channel centres, independent of the old Geant4 idealized geometry.
xc=[];yc=[]
for n in range(1,82):
 ts=[a for a in g['tracks'] if a['net']==str(n) and abs(a['y2_mm']-a['y1_mm'])>250 and abs(a['x2_mm']-a['x1_mm'])<1e-6]
 assert ts;xc.append(np.mean([.5*(a['x1_mm']+a['x2_mm']) for a in ts]))
 fs=[a for a in g['fills'] if a['net']=='Y'+str(n)];assert fs
 areas=np.array([abs((a['x2_mm']-a['x1_mm'])*(a['y2_mm']-a['y1_mm'])) for a in fs]);yc.append(np.average([.5*(a['y1_mm']+a['y2_mm']) for a in fs],weights=areas))
xc=np.array(xc);yc=np.array(yc);assert np.all(np.diff(xc)>0) and np.all(np.diff(yc)>0)
px=xc[10]+.35*(xc[11]-xc[10]);py=yc[20]-.25*(yc[20]-yc[19]);wx=fractions(xc,px,sigma);wy=fractions(yc,py,sigma)
assert np.isclose(wx.sum(),1) and np.isclose(wy.sum(),1)
# Zero-width limit reproduces one channel; an edge case must preserve lost charge.
assert np.argmax(fractions(xc,px,0))==10
assert fractions(xc,xc[0]-.5*(xc[1]-xc[0]),sigma).sum()<.51
rows=[];old=[];weights=[]
for axis,c,pos,w,prefix,near in [('X',xc,px,wx,'',11),('Y',yc,py,wy,'Y',21)]:
 for n,weight in enumerate(w,1):
  weights.append(dict(axis=axis,net=prefix+str(n),centre_pcb_mm=c[n-1],hit_pcb_mm=pos,charge_fraction_within_axis=weight))
  if weight==0:continue
  for r in trace.itertuples(index=False):rows.append(dict(net=prefix+str(n),t0_ns=r.t0_ns,t1_ns=r.t1_ns,current_nA=r.current_nA*.5*weight))
 for r in trace.itertuples(index=False):old.append(dict(net=prefix+str(near),t0_ns=r.t0_ns,t1_ns=r.t1_ns,current_nA=r.current_nA*.5))
pd.DataFrame(weights).to_csv(P/'strip_charge_fractions.csv',index=False)
pd.DataFrame(rows).to_csv(P/'strip_input_currents.csv',index=False)
raw=route(pd.DataFrame(rows),top,dt=dt,end=400);single=route(pd.DataFrame(old),top,dt=dt,end=400)
ht=np.arange(0,12*tau+dt/2,dt);h=(ht/tau)*np.exp(1-ht/tau)
spec=[('X_net1','X1','#0000ff'),('X_net81','X2','#ff0000'),('Y_net1','Y1','#008000'),('Y_net81','Y2','#cc00cc')]
new={k:fftconvolve(raw[k],h)*dt for k,_,_ in spec};baseline={k:fftconvolve(single[k],h)*dt for k,_,_ in spec};scale=max(a.max() for a in baseline.values());t=dt/2+np.arange(len(next(iter(new.values()))))*dt
for axis,prefix,w in [('X','',wx),('Y','Y',wy)]:
 source=sum(raw['strip_'+prefix+str(n)].sum()*dt for n in range(1,82) if 'strip_'+prefix+str(n) in raw)
 dest=sum(raw[axis+'_net'+e].sum()*dt for e in ['1','81']);assert np.isclose(source,dest)
 for e in ['1','81']:assert np.isclose(new[axis+'_net'+e].sum(),baseline[axis+'_net'+e].sum(),rtol=1e-10)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.linewidth':1.,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'legend.frameon':False})
for compare in [False,True]:
 fig,ax=plt.subplots(figsize=(9,4.9))
 if compare:
  for i,(k,_,_) in enumerate(spec):ax.plot(t,baseline[k]/scale,c='.65',ls='--',lw=1.1,label='Previous single-strip' if i==0 else None)
 for k,lab,col in spec:ax.plot(t,new[k]/scale,c=col,lw=1.5,label=lab+' (multi-strip)')
 ax.set(xlabel='Time since electron birth (ns)',ylabel='Endpoint output / previous single-strip peak',xlim=(20,350),ylim=(-.035,1.05));ax.minorticks_on();ax.legend(loc='upper center',bbox_to_anchor=(.53,1),fontsize=9)
 fig.tight_layout();name='comparison' if compare else 'shared_strip_four_ends';fig.savefig(P/(name+'.png'),dpi=220);fig.savefig(P/(name+'.pdf'));plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10,4))
for ax,c,pos,w,lab in [(axes[0],xc,px,wx,'X'),(axes[1],yc,py,wy,'Y')]:
 ax.bar(c-pos,w,width=np.median(np.diff(c))*.8,facecolor='none',edgecolor='blue');ax.set(xlim=(-12,12),xlabel=f'{lab} channel centre - hit (mm)',ylabel='Charge fraction within axis',title=f'{lab} geometric sharing');ax.minorticks_on()
fig.tight_layout();fig.savefig(P/'strip_sharing.png',dpi=180)
pd.DataFrame({'time_ns':t,**{k:a/scale for k,a in new.items()},**{'single_'+k:a/scale for k,a in baseline.items()}}).to_csv(P/'four_end_waveforms.csv',index=False)
meta={'model':'Gaussian geometric charge sharing over actual PCB channel Voronoi cells + per-strip delays + endpoint sum + unipolar shaping','fwhm_mm':fwhm,'sigma_mm':sigma,'width_source':'R015 ~4mm measured avalanche charge profile; Gaussian shape is a preset approximation, not a fit to this PCB','hit_pcb_mm':[px,py],'position_presets':'X11 centre +0.35 pitch; Y21 centre -0.25 pitch','tau_ns':tau,'tap_delay_ns':4,'axis_coupling':.5,'endpoint_split':.5,'normalization':'same old single-strip peak for every curve, no independent peak renormalization','area_and_charge_conservation_pass':True,'references':['R001','R007','R014','R015'],'limitations':['not a weighting-field simulation','assumed common temporal shape across strips','no noise, reflection, saturation or calibrated voltage']}
(P/'metadata.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps({'output':str(P),'X fractions >1e-3':[(i+1,round(float(w),4)) for i,w in enumerate(wx) if w>1e-3],'Y fractions >1e-3':[(i+1,round(float(w),4)) for i,w in enumerate(wy) if w>1e-3]}))
