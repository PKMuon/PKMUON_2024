from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection,PatchCollection
from matplotlib.patches import Rectangle
P=Path(__file__).resolve().parents[1]/'results/strip_readout_20261030'
j=json.loads((P/'pcb_geometry.json').read_text(encoding='utf-8'))
edges=[]
for cid,c in enumerate(j['components']):
 if c['PATTERN']!='Stest2':continue
 pins={int(p['pin']):p['net'] for p in j['pads'] if p['component_id']==cid}
 assert set(pins)==set(range(1,15))
 assert all(pins[i]=='GND' for i in [1,2,14])
 for k in range(3,13):edges.append({'component_id':cid,'pin0':k,'pin1':k+1,'net0':pins[k],'net1':pins[k+1],'nominal_delay_ns':4})
assert len(edges)==160
routes={}
for prefix,axis in [('', 'X'),('Y','Y')]:
 for n in range(1,82):
  name=prefix+str(n)
  if n<81:assert any(e['net0']==name and e['net1']==prefix+str(n+1) for e in edges)
  routes[name]={'axis':axis,'tap_index':n-1,'delay_to_net1_ns':4*(n-1),'delay_to_net81_ns':4*(81-n),'split_fraction':None}
extra=sorted({f['net'] for f in j['fills']}|{t['net'] for t in j['tracks']},key=str)
extra=[n for n in extra if n not in routes and n is not None]
(P/'pcb_delay_topology.json').write_text(json.dumps({'source_sha256':j['sha256'],'status':'pin netlist decoded; 4 ns uses previously specified 1507-40A; copper routing propagation, impedances and transfer functions not calculated','edges':edges,'signal_routes':routes,'additional_copper_nets_not_on_delay_chain':extra},indent=2),encoding='utf-8')
plt.rcParams.update({'font.size':12,'xtick.direction':'in','ytick.direction':'in'})
fig,ax=plt.subplots(figsize=(10,10))
fills=[f for f in j['fills'] if f['layer'] in [1,32]]
pc=PatchCollection([Rectangle((min(f['x1_mm'],f['x2_mm']),min(f['y1_mm'],f['y2_mm'])),abs(f['x2_mm']-f['x1_mm']),abs(f['y2_mm']-f['y1_mm'])) for f in fills],color='#0c8b43',alpha=.6,linewidth=0);ax.add_collection(pc)
ts=[t for t in j['tracks'] if t['layer'] in [1,32]]
# Geometry overview only: true centre lines; width is encoded in the JSON, not to scale here.
lc=LineCollection([[(t['x1_mm'],t['y1_mm']),(t['x2_mm'],t['y2_mm'])] for t in ts],colors=['#c00000' if t['layer']==1 else '#0000cc' for t in ts],linewidths=.4);ax.add_collection(lc)
ax.autoscale();ax.set_aspect('equal');ax.set(xlabel='PCB X (mm)',ylabel='PCB Y (mm)',title='Decoded PCB: copper fills (green), track centre lines (red/blue)')
fig.tight_layout();fig.savefig(P/'pcb_copper_overview.png',dpi=180)
print('PASS: 16 delay footprints; 160 pin-to-pin intervals; both chains 81 taps / 320 ns. Extra copper nets:',extra)
