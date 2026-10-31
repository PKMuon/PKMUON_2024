import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parents[1]/'scripts'))
import numpy as np
from trigger_readout import trigger,coordinate
t=np.arange(0,30,.1);y=np.maximum(0,1-np.abs(t-10)/5)
a=trigger(t,y);assert a['valid'] and np.isclose(a['time_ns'],6)
b=trigger(t,y*7);assert np.isclose(a['time_ns'],b['time_ns'])
c=trigger(t+2.37,y);assert np.isclose(c['time_ns'],a['time_ns']+2.37)
assert trigger(t,np.zeros(len(t)))['status']=='no_positive_signal'
assert trigger(t,np.ones(len(t)))['status']=='left_truncated'
y2=y+np.maximum(0,1-np.abs(t-23)/3);assert not trigger(t,y2)['valid']
assert np.isclose(coordinate(50,290,np.arange(81)*3.457),-30*3.457)
print('PASS: analytic crossing, amplitude invariance, time translation, no signal, truncation, multi-pulse ambiguity, coordinate convention')
