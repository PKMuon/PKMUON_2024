"""Waveform example for simplified PCB strip routing. References R007/R014/R015.
Uses actual existing single-electron linear avalanche trace; channel assignments
and equal X/Y coupling are explicit demonstration presets, not weighting fields.
"""
from pathlib import Path
import sys,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from strip_delay_readout import route
R=Path(__file__).resolve().parents[2]
P=R/'deliverables/Strip_Waveforms_20261030';P.mkdir(parents=True,exist_ok=True)
trace=pd.read_csv(R/'rpc_response/examples/single_electron_trace.csv')
meta=json.loads((R/'rpc_response/examples/single_electron_trace_metadata.json').read_text(encoding='utf-8-sig'))
top=json.loads((R/'rpc_response/config/pcb_delay_topology.json').read_text())
rows=[]
# Each template's total current is divided equally between the X/Y views.
# Hit B is a shifted/scaled copy, not an independent avalanche realization.
for label,nx,ny,shift,scale in [('A','11','Y21',0.,1.),('B','71','Y61',25.,.7)]:
 for net in [nx,ny]:
  for r in trace.itertuples(index=False):rows.append(dict(net=net,t0_ns=r.t0_ns+shift,t1_ns=r.t1_ns+shift,current_nA=r.current_nA*scale*.5))
inp=pd.DataFrame(rows);inp.to_csv(P/'strip_input_currents.csv',index=False)
dt=.02;out=route(inp,top,dt=dt,end=400);df=pd.DataFrame(out);df.to_csv(P/'strip_and_endpoint_currents.csv',index=False)
for axis in ['X','Y']:
 nets=['11','71'] if axis=='X' else ['Y21','Y61']
 src=sum(out['strip_'+n].sum()*dt for n in nets)
 dst=sum(out[axis+'_net'+e].sum()*dt for e in ['1','81'])
 assert np.isclose(src,dst,rtol=1e-12)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.linewidth':1.2,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,'legend.frameon':False})
fig,axes=plt.subplots(2,2,figsize=(14,9))
t=out['t_center_ns']
for a,nets,title in [(axes[0,0],['11','71'],'(a) X strip currents'),(axes[0,1],['Y21','Y61'],'(b) Y strip currents')]:
 for n,c,lab in zip(nets,['#0000cc','#d60000'],['A: t0=0 ns, scale=1','B: t0=25 ns, scale=0.7']):a.plot(t,out['strip_'+n],c=c,lw=1.8,label=f'Net {n}; {lab}')
 a.set(xlim=(0,45),xlabel='Time from first electron birth (ns)',ylabel='Strip current magnitude (nA)',title=title);a.legend(fontsize=10)
for a,axis,title in [(axes[1,0],'X','(c) X delay-chain outputs'),(axes[1,1],'Y','(d) Y delay-chain outputs')]:
 for e,c,ls in [('1','#0000cc','-'),('81','#d60000','--')]:a.plot(t,out[axis+'_net'+e],c=c,ls=ls,lw=1.8,label=f'{axis} chain: net {e} end')
 a.set(xlim=(0,350),xlabel='Time from first electron birth (ns)',ylabel='Endpoint current magnitude (nA)',title=title);a.legend(fontsize=10)
for a in axes.flat:a.minorticks_on();a.tick_params(which='both',top=True,right=True)
fig.suptitle('Simplified strip readout: separate strip inputs, PCB tap delays, endpoint sums',fontsize=17)
fig.text(.5,.015,'Linear avalanche template; preset strip assignment and 1/2 X/Y coupling; 1/2 split to each end; no amplifier shaping, noise or reflection.',ha='center',fontsize=10)
fig.tight_layout(rect=(0,.04,1,.96));fig.savefig(P/'strip_waveforms.png',dpi=190);fig.savefig(P/'strip_waveforms.pdf')
peaks=[]
for hit,nets,shift,scale in [('A',['11','Y21'],0,1),('B',['71','Y61'],25,.7)]:
 for net in nets:
  r=top['signal_routes'][net]
  for e in ['1','81']:peaks.append(dict(hit=hit,net=net,endpoint=r['axis']+'_net'+e,delay_ns=r['delay_to_net'+e+'_ns'],source_shift_ns=shift,source_scale=scale))
pd.DataFrame(peaks).to_csv(P/'delay_assignments.csv',index=False)
metadata={'source_trace':'rpc_response/examples/single_electron_trace.csv','source_avalanche':meta,'hits':[{'label':'A','X_net':'11','Y_net':'Y21','t0_ns':0,'scale':1},{'label':'B','X_net':'71','Y_net':'Y61','t0_ns':25,'scale':.7}],'coupling_fraction_each_axis':.5,'split_fraction_each_end':.5,'tap_delay_ns':4,'chain_delay_ns':320,'group_window_ns':50,'endpoint_labels':'PCB net1/net81; not experimentally identified L/R/U/D','sign':'current magnitude inherited from electron-only planar trace, not signed strip weighting-field prediction','model':'demonstration of channel partition and delay routing; not a full Heed muon event or calibrated detector voltage waveform','references':['R007','R014','R015'],'charge_conservation_pass':True}
(P/'metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
print(P)
