"""PCB-connected ideal strip readout, no weighting-field solver required.

Input CSV: net,t0_ns,t1_ns,current_nA. Each row is a constant current interval
on one strip. Multiple tracks/rows on a strip add linearly. Net names are the
actual PCB names: 1..81, Y1..Y81. R007/R014/R015 in REFERENCES.md support the
nominal 4 ns tap interval. Endpoint labels denote net1/net81, not measured L/R.
Equal, lossless 1/2 current splitting is an explicit circuit approximation.
"""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from electrode_signals import deposit
def route(rows,topology,dt=.1,end=500,fraction_to_1=.5):
    if not (np.isfinite([dt,end,fraction_to_1]).all() and dt>0 and end>0 and 0<=fraction_to_1<=1):raise ValueError('Invalid output grid or split')
    routes=topology['signal_routes'];nb=int(np.ceil(end/dt));out={k:np.zeros(nb) for k in ['X_net1','X_net81','Y_net1','Y_net81']};strips={}
    for row in rows.itertuples(index=False):
        net=str(row.net)
        if net not in routes:raise ValueError('Unknown or non-readout PCB net: '+net)
        r=routes[net];q=float(row.current_nA)*(row.t1_ns-row.t0_ns)
        local=strips.setdefault('strip_'+net,np.zeros(nb));deposit(local,row.t0_ns,row.t1_ns,q,dt)
        for side,fraction in [('1',fraction_to_1),('81',1-fraction_to_1)]:
            delay=r['delay_to_net'+side+'_ns'];deposit(out[r['axis']+'_net'+side],row.t0_ns+delay,row.t1_ns+delay,q*fraction,dt)
    return {'t_center_ns':(np.arange(nb)+.5)*dt,**strips,**out}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input');p.add_argument('output');p.add_argument('--topology',default=str(Path(__file__).parents[1]/'config/pcb_delay_topology.json'));p.add_argument('--dt-ns',type=float,default=.1);p.add_argument('--end-ns',type=float,default=500);p.add_argument('--fraction-to-net1',type=float,default=.5);a=p.parse_args()
    dest=Path(a.output)
    if dest.exists():raise FileExistsError(dest)
    topology=json.loads(Path(a.topology).read_text(encoding='utf-8'));rows=pd.read_csv(a.input,dtype={'net':str});data=route(rows,topology,a.dt_ns,a.end_ns,a.fraction_to_net1)
    pd.DataFrame(data).to_csv(dest,index=False)
    dest.with_suffix('.metadata.json').write_text(json.dumps({'model':'PCB net topology + ideal strip delay/sum','references':['R007','R014','R015'],'pcb_sha256':topology['source_sha256'],'tap_delay_ns':4,'chain_delay_ns':320,'fraction_to_net1':a.fraction_to_net1,'dt_ns':a.dt_ns,'current_unit':'nA','not_included':['interstrip charge-sharing calculation','impedance dispersion and reflections','electronics transfer function','noise and discriminator']},indent=2),encoding='utf-8')
if __name__=='__main__':main()
