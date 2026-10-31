"""Causal unipolar shaping of PCB strip-delay endpoint currents.
R001 Garfield++ User Guide, Signals/Shaper example, h=(t/tau)exp(1-t/tau).
R007/R014/R015 for the upstream delay topology. Time constants are presets.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve
R=Path(__file__).resolve().parents[2]
SOURCE=R/'deliverables/Strip_Waveforms_20261030'
P=R/'deliverables/Strip_Waveforms_Shaped_20261030';P.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(SOURCE/'strip_and_endpoint_currents.csv');t=d.t_center_ns.to_numpy();dt=float(np.median(np.diff(t)))
assert np.allclose(np.diff(t),dt)
taus=[5.,10.,25.];ht=np.arange(int(round(12*max(taus)/dt))+1)*dt
endpoints=['X_net1','X_net81','Y_net1','Y_net81'];raw={};checks=[]
for tau in taus:
 h=(ht/tau)*np.exp(1-ht/tau)
 assert h[0]==0 and np.isclose(ht[h.argmax()],tau)
 for endpoint in endpoints:
  a=d[endpoint].to_numpy();v=fftconvolve(a,h,mode='full')*dt
  # FFT roundoff only: no clipping physical negative signals.
  bound=max(1.,np.max(np.abs(v)))*1e-12
  v[np.abs(v)<bound]=0
  assert np.isclose(v.sum()*dt,a.sum()*dt*h.sum()*dt,rtol=1e-9)
  first=np.flatnonzero(np.abs(a)>0)[0]
  assert np.max(np.abs(v[:first+1]),initial=0)<bound
  raw[f'{endpoint}_tau{tau:g}_ns']=v
  checks.append({'endpoint':endpoint,'tau_ns':tau,'causality_pass':True,'discrete_convolution_area_pass':True})
scale=max(np.max(np.abs(v)) for v in raw.values());out={k:v/scale for k,v in raw.items()}
tout=t[0]+np.arange(len(next(iter(out.values()))))*dt
pd.DataFrame({'t_center_ns':tout,**out}).to_csv(P/'shaped_endpoint_waveforms.csv',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.linewidth':1.2,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'legend.frameon':False})
fig,axes=plt.subplots(2,2,figsize=(14,9))
for ax,endpoint in zip(axes.flat,endpoints):
 for tau,col,ls in zip(taus,['#0000cc','#d60000','#008000'],['-','--','-.']):
  ax.plot(tout,out[f'{endpoint}_tau{tau:g}_ns'],c=col,ls=ls,lw=1.8,label=f'Preset shaping time {tau:g} ns')
 ax.set(title=endpoint.replace('_',' '),xlabel='Time from first electron birth (ns)',ylabel='Shaped output (common relative scale)',xlim=(0,480),ylim=(-.025,1.08));ax.legend(fontsize=10);ax.minorticks_on();ax.tick_params(which='both',top=True,right=True)
fig.suptitle('Strip delay-line outputs with causal unipolar front-end shaping',fontsize=18)
fig.text(.5,.015,'Same source pulses and PCB delays as before. One common normalization for all 12 curves; no calibrated voltage gain, noise or reflections.',ha='center',fontsize=10)
fig.tight_layout(rect=(0,.04,1,.96));fig.savefig(P/'shaped_strip_waveforms.png',dpi=190);fig.savefig(P/'shaped_strip_waveforms.pdf')
meta={'source':str(SOURCE/'strip_and_endpoint_currents.csv'),'source_metadata':json.loads((SOURCE/'metadata.json').read_text(encoding='utf-8')),'model':'V proportional to causal convolution I*h','h':'(t/tau)*exp(1-t/tau) for t>=0, otherwise 0','preset_tau_ns':taus,'dt_ns':dt,'impulse_response_tail_ns':float(ht[-1]),'normalization':'one common maximum across all endpoints and all tau settings; not independent curve normalization','relative_scale_divisor_nA_ns':scale,'output_end_ns':float(tout[-1]+dt/2),'references':['R001','R007','R014','R015'],'checks':checks,'limitations':['preset front-end response, not measured electronics','no calibrated mV scale','noise/reflections/dispersion absent','inherited simplified strip assignment, equal X/Y coupling and equal endpoint current split']}
(P/'metadata.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
print('PASS: 12 causal convolutions, discrete area identities, full tails retained. Output:',P)
