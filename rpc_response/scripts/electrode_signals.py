"""Signed per-electrode Ramo currents from 3D weighting potentials and charge paths.
References: R001, R006, R008 in ../REFERENCES.md and ../references.json.
Map NPZ: x_mm,y_mm,z_mm,names,phi[electrode,x,y,z]. No planar fallback.
CSV paths: x0_mm,y0_mm,z0_mm,x1_mm,y1_mm,z1_mm,t0_ns,t1_ns,q_e.
q_e is SIGNED charge in elementary charges (electron=-1). An avalanche packet
may use the average population over a short segment; do not add an edep weight.
Wiring JSON: electrodes[name][endpoint]={delay_ns,fraction}; fractions must be
explicit, as current splitting is a circuit boundary condition, not Ramo theory.
Output current is nA. No voltage, amplification, noise or waveform normalization.
"""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
E_NS_TO_NA=.1602176634
def deposit(bins,t0,t1,charge,dt):
    if not (np.isfinite([t0,t1,charge,dt]).all() and t1>t0 and t0>=0):raise ValueError('Invalid charge interval')
    if t1>len(bins)*dt+1e-9:raise ValueError('Output time range truncates signal')
    for k in range(int(np.floor(t0/dt)),min(len(bins),int(np.ceil(t1/dt)))):
        overlap=max(0.,min(t1,(k+1)*dt)-max(t0,k*dt))
        bins[k]+=charge*overlap/((t1-t0)*dt)
def calculate(mapfile,paths,wiring,dt,end,substep):
    m=np.load(mapfile,allow_pickle=False);names=list(m['names'].astype(str));grid=tuple(m[k] for k in ['x_mm','y_mm','z_mm'])
    if m['phi'].shape!=(len(names),*[len(g) for g in grid]):raise ValueError('Invalid map shape')
    for g in grid:
        if len(g)<2 or not np.all(np.diff(g)>0):raise ValueError('Map axes must increase')
    if dt<=0 or end<=0 or substep<=0:raise ValueError('Invalid integration grid')
    rows=pd.read_csv(paths);nb=int(np.ceil(end/dt));out={};metadata={}
    if set(names)!=set(wiring['electrodes']):raise ValueError('Every mapped electrode needs explicit circuit routing, including unused electrodes')
    for i,name in enumerate(names):
        field=RegularGridInterpolator(grid,m['phi'][i],bounds_error=True)
        local=np.zeros(nb);routes=wiring['electrodes'][name]
        for r in rows.itertuples(index=False):
            p0=np.array([r.x0_mm,r.y0_mm,r.z0_mm]);p1=np.array([r.x1_mm,r.y1_mm,r.z1_mm])
            if r.t1_ns<=r.t0_ns or not np.isfinite(r.q_e):raise ValueError('Invalid transport segment')
            n=max(1,int(np.ceil(np.linalg.norm(p1-p0)/substep)))
            u=np.linspace(0,1,n+1);pos=p0[None,:]+u[:,None]*(p1-p0);phi=field(pos)
            times=r.t0_ns+u*(r.t1_ns-r.t0_ns)
            # i=-q*v.Ew = q*d(phi_w)/dt. Integrate before time binning.
            dq=r.q_e*E_NS_TO_NA*np.diff(phi) # nA ns
            for a,b,q in zip(times[:-1],times[1:],dq):
                deposit(local,a,b,q,dt)
                for endpoint,route in routes.items():
                    delay=float(route['delay_ns']);fraction=float(route['fraction'])
                    if delay<0 or not 0<=fraction<=1:raise ValueError('Invalid circuit delay/split')
                    dest=out.setdefault(endpoint,np.zeros(nb));deposit(dest,a+delay,b+delay,q*fraction,dt)
        out['electrode:'+name]=local
        metadata[name]={'signed_charge_fC':float(local.sum()*dt*.001)}
    return {'t_center_ns':(np.arange(nb)+.5)*dt,**out},metadata
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('map');p.add_argument('paths');p.add_argument('wiring');p.add_argument('output');p.add_argument('--dt-ns',type=float,default=.1);p.add_argument('--end-ns',type=float,default=500);p.add_argument('--substep-mm',type=float,default=.01);a=p.parse_args()
    out=Path(a.output)
    if out.exists():raise FileExistsError(out)
    w=json.loads(Path(a.wiring).read_text(encoding='utf-8'));data,meta=calculate(a.map,a.paths,w,a.dt_ns,a.end_ns,a.substep_mm)
    pd.DataFrame(data).to_csv(out,index=False)
    out.with_suffix('.metadata.json').write_text(json.dumps({'model':'static weighting potential line integral; ideal user-specified circuit delays/splits','current_unit':'nA','electrodes':meta},indent=2),encoding='utf-8')
if __name__=='__main__':main()
