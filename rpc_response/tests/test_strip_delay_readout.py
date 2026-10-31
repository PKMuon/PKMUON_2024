import sys,json
from pathlib import Path
import pandas as pd
import numpy as np
sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
from strip_delay_readout import route
t=json.loads((Path(__file__).parents[2]/'config/pcb_delay_topology.json').read_text())
# Two temporally coincident separated strips must retain two pulses at each end.
rows=pd.DataFrame([dict(net='11',t0_ns=0.,t1_ns=1.,current_nA=2.),dict(net='71',t0_ns=0.,t1_ns=1.,current_nA=4.)])
o=route(rows,t,dt=.1,end=400)
assert np.isclose(o['X_net1'][400],1) and np.isclose(o['X_net1'][2800],2)
assert np.isclose(o['X_net81'][400],2) and np.isclose(o['X_net81'][2800],1)
assert np.isclose((o['X_net1'].sum()+o['X_net81'].sum())*.1,6)
# Overlapping contributions on the same strip superpose, including fractional ns.
r=pd.DataFrame([dict(net='Y41',t0_ns=.13,t1_ns=1.17,current_nA=-3.)]*2)
o=route(r,t,dt=.1,end=400)
assert np.isclose((o['Y_net1'].sum()+o['Y_net81'].sum())*.1,-6*1.04)
try:route(pd.DataFrame([dict(net='GND',t0_ns=0.,t1_ns=1.,current_nA=1.)]),t)
except ValueError:pass
else:raise AssertionError('Ground net must not become a readout strip')
print('PASS: separated simultaneous strips, reversed endpoint ordering, same-strip overlap, signed charge conservation, grounded-net rejection')
