"""Read four shaped waveform columns and reconstruct positions by 20% timing.
CSV columns: time_ns,X_net1,X_net81,Y_net1,Y_net81; optional event column.
Output coordinates are relative to PCB centre channel, not global Geant4 axes.
"""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from trigger_readout import trigger,coordinate
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('input');p.add_argument('output');p.add_argument('--fraction',type=float,default=.2);p.add_argument('--geometry',default=str(Path(__file__).parents[1]/'config/pcb_geometry.json'));a=p.parse_args()
 dest=Path(a.output)
 if dest.exists():raise FileExistsError(dest)
 g=json.loads(Path(a.geometry).read_text(encoding='utf-8'));centres={}
 for axis in ['X','Y']:
  c=[]
  for n in range(1,82):
   if axis=='X':
    q=[x for x in g['tracks'] if x['net']==str(n) and abs(x['y2_mm']-x['y1_mm'])>250 and abs(x['x2_mm']-x['x1_mm'])<1e-6];c.append(np.mean([.5*(x['x1_mm']+x['x2_mm']) for x in q]))
   else:
    q=[x for x in g['fills'] if x['net']=='Y'+str(n)];c.append(np.average([.5*(x['y1_mm']+x['y2_mm']) for x in q],weights=[abs((x['x2_mm']-x['x1_mm'])*(x['y2_mm']-x['y1_mm'])) for x in q]))
  if not np.all(np.isfinite(c)) or not np.all(np.diff(c)>0):raise ValueError('Invalid PCB channel centres')
  centres[axis]=c
 data=pd.read_csv(a.input)
 if 'event' not in data:data['event']=0
 out=[]
 for event,df in data.groupby('event',sort=True):
  row={'event':event}
  for axis in ['X','Y']:
   valid=True;tt=[]
   for side in ['1','81']:
    r=trigger(df.time_ns.to_numpy(),df[axis+'_net'+side].to_numpy(),a.fraction);row['T_'+axis+side+'_ns']=r['time_ns'];row[axis+side+'_status']=r['status'];row[axis+side+'_peak']=r['peak'];row[axis+side+'_threshold']=r['threshold'];valid &=r['valid'];tt.append(r['time_ns'])
   value=coordinate(*tt,centres[axis]) if valid else float('nan');row[axis.lower()+'_rec_mm']=value;row[axis+'_valid']=bool(valid and np.isfinite(value))
  out.append(row)
 pd.DataFrame(out).to_csv(dest,index=False)
if __name__=='__main__':main()
