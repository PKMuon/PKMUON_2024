from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parents[2]/'.part1_build/assets'
d=pd.read_csv(P/'trace.csv'); meta=json.loads((P/'trace_metadata.json').read_text(encoding='utf-8-sig'))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':16,'axes.linewidth':1.3,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'legend.frameon':False,'figure.facecolor':'white'})
C=['#0000cc','#d60000','#008000','#990099']
def save(fig,name):
    for a in fig.axes: a.minorticks_on();a.tick_params(which='both',top=True,right=True)
    fig.tight_layout();fig.savefig(P/(name+'.png'),dpi=180);fig.savefig(P/(name+'.pdf'));plt.close(fig)
t=np.r_[0,d.t1_ns]; n=np.r_[d.n0.iloc[0],d.n1]
f,a=plt.subplots(figsize=(9,5.4));a.step(t,n,where='post',c=C[0],label='One birth/death realization');a.plot(t,np.exp((13.59775179838482-.177720119410425)*t*.10572738077531714),'--',c=C[1],label='Ensemble mean, not this realization');a.set(xlabel='Time since electron birth (ns)',ylabel='Number of electrons',yscale='log');a.legend(fontsize=12);save(f,'01_avalanche')
f,a=plt.subplots(figsize=(9,5.4));a.stairs(d.current_nA,np.r_[d.t0_ns.iloc[0],d.t1_ns],color=C[0],label='Electron-only planar Ramo current');a.set(xlabel='Time since electron birth (ns)',ylabel='Induced current magnitude (nA)');a.legend(fontsize=12);save(f,'02_current')
# Fine uniform grid: integrate each piecewise constant current bin exactly at sampling centres.
dt=.02;grid=np.arange(0,220+dt/2,dt);centre=grid+dt/2
idx=np.searchsorted(d.t1_ns.to_numpy(),centre,side='right');current=np.where(idx<len(d),d.current_nA.to_numpy()[np.minimum(idx,len(d)-1)],0.)
def shaped(tau):
    h=(grid/tau)*np.exp(1-grid/tau);return np.convolve(current,h)[:len(grid)]*dt
vs={tau:shaped(tau) for tau in [5,10,25]};base=vs[5];peak=base.max()
f,a=plt.subplots(figsize=(9,5.4))
for i,(tau,v) in enumerate(vs.items()): a.plot(grid,v/v.max(),color=C[i],label=f'Preset shaping time {tau} ns')
a.set(xlabel='Time since electron birth (ns)',ylabel='Output / its own peak',xlim=(0,150));a.legend(fontsize=12);save(f,'03_shaping')
delay={'XL':160+(50-1.062)/(2*.432125),'XR':160-(50-1.062)/(2*.432125),'YD':160-30/(2*.432875),'YU':160+30/(2*.432875)}
tg=np.arange(0,420+dt/2,dt)
f,a=plt.subplots(figsize=(9,5.4))
for i,(name,shift) in enumerate(delay.items()): a.plot(tg,np.interp(tg-shift,grid,base,left=0,right=0)/peak,c=C[i],label=f'{name}, delay {shift:.2f} ns')
a.set(xlabel='Time since electron birth (ns)',ylabel='Endpoint output / single-pulse peak',xlim=(80,300));a.legend(fontsize=12);save(f,'04_endpoints')
f,a=plt.subplots(figsize=(9,5.4))
for i,shift in enumerate([0,10,25]): a.plot(grid,(base+np.interp(grid-shift,grid,base,left=0,right=0))/peak,c=C[i],label=f'Two equal pulses, offset {shift} ns')
a.set(xlabel='Time since first electron birth (ns)',ylabel='Sum / one isolated pulse peak',xlim=(0,90));a.legend(fontsize=12);save(f,'05_overlap')
charge=np.sum(d.current_nA*(d.t1_ns-d.t0_ns))*.001
assert np.isclose(charge,meta['ramo_fC'],rtol=1e-12)
assert np.isclose(.432125*(delay['XL']-delay['XR'])+1.062,50)
assert np.isclose(.432875*(delay['YD']-delay['YU']),-30)
meta.update({'shaper':'h(t)=(t/tau)*exp(1-t/tau), t>=0','tau_presets_ns':[5,10,25],'endpoint_example_x_mm':50,'endpoint_example_y_mm':-30,'endpoint_delays_ns':delay,'waveforms':'model demonstration; no calibrated voltage scale, ion tail, diffusion, saturation, noise or reflection','integral_current_fC':charge,'endpoint_roundtrip_pass':True})
(P/'waveform_metadata.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
pd.DataFrame({'t_ns':grid,'electron_current_nA':current,'shaped_tau5_relative':base/peak}).to_csv(P/'waveforms.csv',index=False)
print(json.dumps(meta,indent=2))
