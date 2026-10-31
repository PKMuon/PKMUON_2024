"""Layered-plane and sheet-mode diagnostics ONLY, not PCB weighting maps.
R008 provides layered dielectric/resistive theory. R010/R011/R012 values are
referenced in the manifest. All scan endpoints are research presets.
"""
from pathlib import Path
import itertools,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[2];P=R/'rpc_response/results/strip_readout_20261030/material_scan_pku';P.mkdir(parents=True,exist_ok=True)
c=json.loads((R/'config/strip_material_scan.json').read_text());v=c['parameters'];e0=8.8541878128e-12
def base(k):return v[k]['baseline']
def weighting(eg,ep,ec,tc):return 1/(2+2*2.7/eg+2*.1/ep+2*tc/ec)
ref=weighting(8,3.25,10,.2)
rows=[]
for eg,ep,ec,tc,rho,rs in itertools.product(*(v[k]['scan'] for k in ['glass_er','pet_er','graphite_er','graphite_thickness_mm','glass_bulk_resistivity_ohm_cm','graphite_sheet_ohm_sq'])):
 ew=weighting(eg,ep,ec,tc)
 # Maxwell relaxation of bulk glass. rho[ohm cm] -> ohm m.
 tau=rho*.01*e0*eg
 # Sheet diffusion diagnostic: sheet between two local capacitive loads,
 # each approximated by glass + film. NOT the actual RPC field solution.
 cap=2*e0/(.0027/eg+.0001/ep)
 wavelength=.003457;k=2*np.pi/wavelength;sheet_tau=rs*cap/k**2
 rows.append(dict(glass_er=eg,pet_er=ep,graphite_er=ec,graphite_thickness_mm=tc,glass_rho_ohm_cm=rho,graphite_rs_ohm_sq=rs,plane_Ew_per_mm=ew,relative_to_plane_baseline=ew/ref,glass_relaxation_s=tau,sheet_mode_diagnostic_tau_ns=sheet_tau*1e9))
d=pd.DataFrame(rows);d.to_csv(P/'material_diagnostics.csv',index=False)
assert np.all(d.plane_Ew_per_mm>0) and np.all(d.plane_Ew_per_mm<=.5)
plt.rcParams.update({'font.size':13,'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True})
fig,ax=plt.subplots(1,2,figsize=(12,4.7))
for ep,col in zip(v['pet_er']['scan'],['blue','red','green']):
 eg=np.linspace(6,10,101);ax[0].plot(eg,weighting(eg,ep,10,.2)/ref,c=col,label=f'PET er preset {ep}')
ax[0].set(xlabel='Glass relative permittivity (preset)',ylabel='Plane weighting field / baseline',title='Layered-plane cross-check, NOT PCB result');ax[0].legend(fontsize=10)
for eg,col in zip(v['glass_er']['scan'],['blue','red','green']):
 rho=np.logspace(10,12,80);ax[1].loglog(rho,rho*.01*e0*eg,c=col,label=f'Glass er preset {eg}')
ax[1].set(xlabel='Glass bulk resistivity (ohm cm)',ylabel='Bulk Maxwell relaxation (s)',title='Bulk glass timescale; not coating dynamics');ax[1].legend(fontsize=10)
fig.tight_layout();fig.savefig(P/'material_checks.png',dpi=180)
summary={'cases':len(d),'kind':'analytic diagnostics only, not real strip waveforms','references':sorted({r for p in v.values() for r in p.get('references',[])}),'graphite_sheet_baseline_ohm_sq':base('graphite_sheet_ohm_sq'),'baseline_plane_Ew_per_mm':ref,'plane_relative_range':[float(d.relative_to_plane_baseline.min()),float(d.relative_to_plane_baseline.max())],'glass_relaxation_s_range':[float(d.glass_relaxation_s.min()),float(d.glass_relaxation_s.max())],'sheet_mode_diagnostic_ns_range':[float(d.sheet_mode_diagnostic_tau_ns.min()),float(d.sheet_mode_diagnostic_tau_ns.max())],'warning':'Bulk relaxation cannot establish graphite sheet transparency. Sheet-mode estimate assumes local load and sinusoidal wavelength equal to pitch; it is not the actual transient weighting field.'}
(P/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
