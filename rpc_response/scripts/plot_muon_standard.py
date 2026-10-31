"""Plot exported ROOT/Heed samples with explicit model/selection labels.

Usage: python plot_muon_standard.py RESULTS/muon_1GeV_normal_v1
Writes PNG, vector PDF pages, multipage PDF, and JSON numerical summaries.
"""
import sys
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

prefix = Path(sys.argv[1]).resolve()
out = prefix.parent / (prefix.name + "_plots")
out.mkdir(exist_ok=True)
def read(suffix):
    return np.genfromtxt(str(prefix) + suffix + ".csv", delimiter=",", names=True)
events, electrons, clusters, kernel, transport = [read(s) for s in
    ("_events", "_electrons", "_clusters", "_kernel", "_transport")]
meta = json.loads(Path(str(prefix)+"_metadata.json").read_text(encoding="utf-8"))
plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":11,
    "axes.spines.top":False, "axes.spines.right":False, "axes.titlesize":12,
    "figure.facecolor":"#f8fafc", "axes.facecolor":"white", "axes.grid":True,
    "grid.alpha":.18, "savefig.facecolor":"#f8fafc"})
blue, teal, orange, red = "#2563eb", "#0d9488", "#d97706", "#dc2626"
summary = {}
def stats(a):
    return dict(mean=float(np.mean(a)), median=float(np.median(a)),
                std=float(np.std(a, ddof=1)), p05=float(np.quantile(a,.05)),
                p95=float(np.quantile(a,.95)))
for key in ("n_clusters", "n_electrons", "energy_transfer_ev", "max_drift_ns"):
    summary[key] = stats(events[key])
rx = events["x_rec_mm"]-events["x_true_mm"]
ry = events["y_rec_mm"]-events["y_true_mm"]
summary.update(events=len(events), electrons=len(electrons),
    invalid_electron_fraction=float(np.mean(electrons["valid"]==0)),
    invalid_muon_charge_fraction=float(np.mean(events["charge_valid"]==0)),
    residual_x_mm=stats(rx), residual_y_mm=stats(ry),
    muon_kinetic_ev=meta["muon_kinetic_ev"], angle_deg=meta["angle_deg"], sigma_input_mm=meta["sigma_input_mm"],
    model="Heed standalone mu- ionization; linear single-electron reference; empirical readout",
    caveat="No saturation/space charge/threshold. Invalid charge totals are NaN, never zero or resampled.")
assert len(events)==meta["events"] and np.all(np.isfinite(rx+ry))
assert np.all(events["n_invalid"]<=events["n_electrons"])
assert np.all(np.isnan(events["charge_fc"][events["charge_valid"]==0]))
assert abs(np.std(rx)-.749)<.03 and abs(np.std(ry)-.749)<.03
assert np.all((electrons["depth_mm"]>=0)&(electrons["depth_mm"]<=2))
assert len(kernel)==41000

pdf = PdfPages(out / "muon_standard_distributions.pdf")
def page(fig, name, title, footer):
    fig.suptitle(title, fontsize=19, fontweight="bold", x=.055, ha="left", y=.985)
    fig.text(.055,.015,footer,fontsize=9,color="#475569",va="bottom")
    fig.tight_layout(rect=(.025,.07,.99,.93),h_pad=2.3,w_pad=2.5)
    fig.savefig(out/(name+".png"),dpi=160)
    fig.savefig(out/(name+".pdf"))
    pdf.savefig(fig)
    plt.close(fig)

# Page 1: complete muon ionization sample; no charge-validity cut.
fig, ax = plt.subplots(2,3,figsize=(15,9))
a=ax[0,0]; n=events["n_clusters"]
a.hist(n,bins=np.arange(n.min()-.5,n.max()+1.5),color=blue,alpha=.85)
a.set(title=f"Primary clusters / muon | mean {n.mean():.2f}",xlabel="Number of clusters",ylabel="Muons")
a=ax[0,1]; n=events["n_electrons"]
cap=np.quantile(n,.995)
a.hist(n,bins=np.arange(0,int(cap)+3,2),color=teal,alpha=.85)
a.set(title=f"In-gap ionization electrons | median {np.median(n):.0f}",xlabel="Electrons / muon",ylabel="Muons")
a.text(.97,.93,f"Shown to 99.5th percentile\n{np.sum(n>cap)} events above range",ha="right",va="top",transform=a.transAxes,fontsize=9)
a=ax[0,2]; en=events["energy_transfer_ev"]/1000
positive=en[en>0]
a.hist(positive,bins=np.geomspace(positive.min()*.95,positive.max()*1.05,65),color=orange,alpha=.85)
a.set_xscale("log"); a.set_yscale("log")
a.set(title=f"Heed energy transfer | median {np.median(en):.3f} keV",xlabel="Transferred energy (keV)",ylabel="Muons / logarithmic bin")
a=ax[1,0]; a.hist(electrons["depth_mm"],bins=40,range=(0,2),color=blue,alpha=.85)
a.set(title="Electron starting depth",xlabel="Distance to electron-collection plane (mm)",ylabel="Electrons")
a=ax[1,1]; a.hist(events["max_drift_ns"],bins=60,color=teal,alpha=.85)
a.axvline(meta["gas_gap_mm"]/meta["velocity_mm_ns"],color=red,ls="--",label="2 mm / drift velocity")
a.set(title="Latest nominal electron arrival / muon",xlabel="Creation time + distance / velocity (ns)",ylabel="Muons"); a.legend(fontsize=9)
a=ax[1,2]; invalid=events["n_invalid"]/np.maximum(events["n_electrons"],1)
a.hist(invalid,bins=35,range=(0,1),color=red,alpha=.8)
a.set(title="Linear-kernel invalid fraction / muon",xlabel="Invalid electron draws / in-gap electrons",ylabel="Muons")
a.text(.97,.93,f"{100*summary['invalid_muon_charge_fraction']:.2f}% of muons have\nundefined TOTAL linear charge",ha="right",va="top",transform=a.transAxes,fontsize=10)
page(fig,"01_muon_ionization",f"Muon standard | {meta['muon_kinetic_ev']/1e9:g} GeV mu- at normal incidence | {len(events):,} events",
     "Pure R134a | 2 mm gas | 298.15 K | 760 Torr | 50.5 kV/cm. Heed includes delta cascades; energy transfer is not necessarily local deposition.\nNominal arrival neglects diffusion; invalid fraction is a model-guard statistic, NOT detector inefficiency.")

# Page 2: single-electron response; selection is deliberately visible.
fig,ax=plt.subplots(2,2,figsize=(13,9))
depths=np.unique(kernel["depth_mm"])
fractions=np.array([np.mean(kernel["valid"][kernel["depth_mm"]==d]==0) for d in depths])
a=ax[0,0]; a.plot(depths,100*fractions,"o-",color=red,ms=4)
a.set(title="Invalid fraction vs starting drift distance",xlabel="Drift distance (mm)",ylabel="Invalid draws (%)",ylim=(-2,102))
a.axhline(50,color="#94a3b8",ls=":")
a=ax[0,1]
for d,c in zip((.25,.5,.75,1.),(blue,teal,orange,"#7c3aed")):
    rows=kernel[np.isclose(kernel["depth_mm"],d)]
    q=rows["charge_fc"][(rows["valid"]==1)&(rows["charge_fc"]>0)]
    a.hist(q,bins=np.geomspace(1e-7,100,90),histtype="step",lw=1.7,color=c,
           label=f"{d:g} mm; valid {np.mean(rows['valid'])*100:.1f}%")
a.set_xscale("log"); a.set_yscale("log")
a.set(title="Single-electron induced charge | valid draws only",xlabel="Electron-only planar induced charge (fC)",ylabel="Samples / logarithmic bin"); a.legend(fontsize=9)
a=ax[1,0]
step=ax[1,0].hist2d(electrons["depth_mm"],electrons["nominal_arrival_ns"],bins=(60,60),cmap="Blues",cmin=1)
fig.colorbar(step[3],ax=a,label="Electrons")
a.set(title="Nominal drift-arrival mapping",xlabel="Drift distance (mm)",ylabel="Creation time + distance / velocity (ns)")
a=ax[1,1]; r=meta["alpha_per_mm"]-meta["eta_per_mm"]
xs=np.linspace(0,2,200)
a.plot(xs,np.exp(r*xs),color=blue,label="Unsaturated analytic mean exp[(alpha-eta)d]")
a.axhline(1e7,color=red,ls="--",label="Per-realization numerical guard: 10^7")
a.set_yscale("log"); a.set(title="Why the full-gap linear model fails",xlabel="Drift distance (mm)",ylabel="Mean electron population"); a.legend(fontsize=8)
page(fig,"02_linear_reference","Linear avalanche reference | validity limits are retained",
     "41 depth positions x 1,000 draws | step 0.01 mm | independent birth/death statistics | NO space charge or physical saturation.\nCharge spectra above are conditional on validity; they do not represent complete muon charge spectra or trigger efficiency.")

# Page 3: calibrated coordinate baseline, not a prediction of electronics performance.
fig,ax=plt.subplots(2,2,figsize=(13,9))
a=ax[0,0]
bins=np.linspace(-3.5,3.5,80)
a.hist(rx,bins=bins,histtype="step",lw=1.6,color=blue,label=f"X: sample sigma {np.std(rx,ddof=1):.4f} mm")
a.hist(ry,bins=bins,histtype="step",lw=1.6,color=teal,label=f"Y: sample sigma {np.std(ry,ddof=1):.4f} mm")
xs=np.linspace(-3.5,3.5,200)
a.plot(xs,len(events)*(bins[1]-bins[0])/(.749*np.sqrt(2*np.pi))*np.exp(-.5*(xs/.749)**2),"--",color=orange,label="Input Gaussian sigma = 0.749 mm")
a.set(title="Coordinate residuals | empirical calibration",xlabel="Reconstructed - true coordinate (mm)",ylabel="Muons / bin"); a.legend(fontsize=9)
a=ax[0,1]; im=a.hist2d(rx,ry,bins=60,range=((-3,3),(-3,3)),cmin=1,cmap="Blues")
fig.colorbar(im[3],ax=a,label="Muons")
a.set(title=f"X/Y residual correlation = {np.corrcoef(rx,ry)[0,1]:.3f}",xlabel="X residual (mm)",ylabel="Y residual (mm)")
a=ax[1,0]; dt=events["txl_ns"]-events["txr_ns"]
a.scatter(events["x_true_mm"][::10],dt[::10],s=5,alpha=.35,color=blue,rasterized=True)
xs=np.array([-100,100]); a.plot(xs,(xs-1.062)/(3.457/8),color=orange,lw=2,label="Nominal geometry / 4 ns tap spacing")
a.set(title="X delay difference vs known position",xlabel="True X (mm)",ylabel="T_XL - T_XR (ns)"); a.legend(fontsize=9)
a=ax[1,1]
ds=events["txl_ns"]+events["txr_ns"]-events["tyu_ns"]-events["tyd_ns"]
a.hist(ds,bins=70,color=teal,alpha=.85)
a.set(title=f"Time-sum difference | sample sigma {np.std(ds,ddof=1):.3f} ns",xlabel="(T_XL + T_XR) - (T_YU + T_YD) (ns)",ylabel="Muons")
summary["time_sum_difference_ns"]=stats(ds)
page(fig,"03_delay_coordinate","Delay-line coordinate baseline | measured spatial sigma imposed",
     "8 x 1507-40A per axis | nominal total delay 320 ns | equal independent endpoint jitter is assumed.\nThese residuals reproduce the supplied sigma by construction. Absolute timing, threshold losses and multi-hit pairing are NOT calibrated.")

# Page 4: gas transport parameters versus field, not event-by-event random distributions.
fig,ax=plt.subplots(2,2,figsize=(13,9))
f=transport["field_kv_cm"]
for a,key,label,color in [(ax[0,0],"velocity_mm_ns","Electron drift speed (mm/ns)",blue),
                         (ax[1,0],"eta_per_mm","Attachment coefficient (1/mm)",orange)]:
    a.plot(f,transport[key],"o-",color=color); a.set(xlabel="Electric field (kV/cm)",ylabel=label,title=label)
a=ax[0,1]; a.plot(f,transport["alpha_per_mm"],"o-",color=blue,label="Townsend alpha")
a.plot(f,transport["alpha_per_mm"]-transport["eta_per_mm"],"s--",color=teal,label="Effective alpha - eta")
a.set(xlabel="Electric field (kV/cm)",ylabel="Coefficient (1/mm)",title="Multiplication coefficients"); a.legend()
a=ax[1,1]; a.plot(f,transport["dl_sqrt_mm"],"o-",color=blue,label="Longitudinal")
a.plot(f,transport["dt_sqrt_mm"],"s-",color=teal,label="Transverse")
a.set(xlabel="Electric field (kV/cm)",ylabel="Diffusion coefficient (sqrt(mm))",title="Diffusion coefficients"); a.legend()
for a in ax.flat: a.axvline(50.5,color=red,ls=":",lw=1.5)
page(fig,"04_gas_transport","Gas transport scan | operating field marked at 50.5 kV/cm",
     "Magboltz: pure R134a, 298.15 K, 760 Torr | 10^7 collisions / field point (exploratory statistics).\n50.5 kV/cm is interpolated inside the table. Transport curves are coefficients vs field, not measured event distributions.")
pdf.close()
(out/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))
print("Plots:",out)
