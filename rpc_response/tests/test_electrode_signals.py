import importlib.util,tempfile,json
from pathlib import Path
import numpy as np
import pandas as pd
spec=importlib.util.spec_from_file_location('signals',Path(__file__).parents[1]/'scripts/electrode_signals.py');s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
with tempfile.TemporaryDirectory() as d:
 p=Path(d);g=np.linspace(0,1,5);x,y,z=np.meshgrid(g,g,g,indexing='ij')
 # Analytic Laplace solutions: x*z and (1-x)*z are both harmonic.
 np.savez(p/'map.npz',x_mm=g,y_mm=g,z_mm=g,names=np.array(['A','B']),phi=np.array([x*z,(1-x)*z]))
 pd.DataFrame([dict(x0_mm=.25,y0_mm=.5,z0_mm=0,x1_mm=.25,y1_mm=.5,z1_mm=1,t0_ns=.13,t1_ns=1.17,q_e=-1)]).to_csv(p/'paths.csv',index=False)
 w={'electrodes':{n:{'L':{'delay_ns':2.37,'fraction':.5},'R':{'delay_ns':4.12,'fraction':.5}} for n in ['A','B']}}
 out,m=s.calculate(p/'map.npz',p/'paths.csv',w,.1,10,.01)
 assert np.isclose(m['A']['signed_charge_fC'],-.25*.0001602176634)
 assert np.isclose(m['B']['signed_charge_fC'],-.75*.0001602176634)
 assert np.isclose(out['L'].sum(),out['R'].sum())
 assert np.isclose(out['L'].sum()+out['R'].sum(),out['electrode:A'].sum()+out['electrode:B'].sum())
 # Reverse trajectory changes the induced current polarity.
 r=pd.read_csv(p/'paths.csv');r['z0_mm']=1;r['z1_mm']=0;r.to_csv(p/'reverse.csv',index=False)
 _,rev=s.calculate(p/'map.npz',p/'reverse.csv',w,.1,10,.01)
 assert np.isclose(rev['A']['signed_charge_fC'],-m['A']['signed_charge_fC'])
 try:s.calculate(p/'map.npz',p/'paths.csv',w,.1,1,.01)
 except ValueError:pass
 else:raise AssertionError('Truncated output must fail')
 print('PASS: separate electrode responses, signed Ramo charge, delay/split conservation, reversed motion, time truncation rejection')
