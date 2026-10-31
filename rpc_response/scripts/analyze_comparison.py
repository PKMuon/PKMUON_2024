"""Matched Edep/current-baseline comparison and explicitly idealized two-hit study.
Usage: python analyze_comparison.py RESULTS/comparison_20261020
All stochastic response assumptions and selections are written into the report.
"""
import sys
import json
import itertools
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

root=Path(sys.argv[1]).resolve(); figures=root/"figures"; figures.mkdir(exist_ok=True)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.titlesize":12,
    "axes.spines.top":False,"axes.spines.right":False,"axes.grid":True,"grid.alpha":.17,
    "figure.facecolor":"#f8fafc","axes.facecolor":"white","savefig.facecolor":"#f8fafc"})
blue,teal,orange,red,purple="#2563eb","#0d9488","#d97706","#dc2626","#7c3aed"
sigma=.749; window=50.; axdelay=3.457/8; aydelay=3.463/8
rng=np.random.default_rng(20261020)
pdf=PdfPages(figures/"all_comparison_figures.pdf")
summary={"window_ns":window,"sigma_mm":sigma,"energy_scan":[],"datasets":{},"overlay":{}}
def save(fig,name,title,footer):
    fig.suptitle(title,x=.055,y=.985,ha="left",fontweight="bold",fontsize=17)
    fig.text(.055,.014,footer,fontsize=8.5,color="#475569",va="bottom")
    fig.tight_layout(rect=(.025,.075,.99,.935),h_pad=2,w_pad=2)
    fig.savefig(figures/(name+".png"),dpi=160); fig.savefig(figures/(name+".pdf")); pdf.savefig(fig); plt.close(fig)
def rms(v):
    a=np.asarray(v,dtype=float); a=a[np.isfinite(a)]
    return float(np.sqrt(np.mean(a*a))) if len(a) else None
def finite_stats(v):
    a=np.asarray(v,dtype=float); a=a[np.isfinite(a)]
    return {"n":len(a),"mean":float(a.mean()),"median":float(np.median(a)),"rms":rms(a),
            "std":float(a.std(ddof=1)),"p95_abs":float(np.quantile(np.abs(a),.95))} if len(a) else {"n":0}

# Full Heed standards, no selection on validity.
energy_frames=[]
for meta_path in (root/"energies").glob("*_metadata.json"):
    prefix=Path(str(meta_path).replace("_metadata.json","")); meta=json.loads(meta_path.read_text())
    e=pd.read_csv(str(prefix)+"_events.csv"); energy_frames.append((meta["muon_kinetic_ev"]/1e9,e))
    assert len(e)==meta["events"]==10000
    assert np.all(np.isnan(e.loc[e.charge_valid==0,"charge_fc"]))
    summary["energy_scan"].append({"energy_gev":meta["muon_kinetic_ev"]/1e9,"events":len(e),
      "clusters_mean":float(e.n_clusters.mean()),"electrons_mean":float(e.n_electrons.mean()),
      "transfer_ev_median":float(e.energy_transfer_ev.median()),"invalid_muon_fraction":float((e.charge_valid==0).mean()),
      "sigma_x_mm":float((e.x_rec_mm-e.x_true_mm).std()),"sigma_y_mm":float((e.y_rec_mm-e.y_true_mm).std())})
energy_frames.sort(key=lambda p:p[0]); summary["energy_scan"].sort(key=lambda d:d["energy_gev"])
scan=pd.DataFrame(summary["energy_scan"]); scan.to_csv(root/"energy_summary.csv",index=False)
fig,a=plt.subplots(2,3,figsize=(15,9))
for i,(energy,e) in enumerate(energy_frames):
    c=plt.cm.viridis(.08+.82*i/max(1,len(energy_frames)-1)); label=f"{energy:g} GeV"
    a[0,0].hist(e.n_clusters,bins=np.arange(-.5,65.5),density=True,histtype="step",color=c,label=label)
    a[0,1].hist(e.n_electrons,bins=np.arange(-.5,301.5,3),weights=np.ones(len(e))/len(e),histtype="step",color=c,label=label)
    a[0,2].hist(e.energy_transfer_ev/1000,bins=np.geomspace(.05,200,90),weights=np.ones(len(e))/len(e),histtype="step",color=c,label=label)
a[0,0].set(title="Heed primary-cluster multiplicity",xlabel="Clusters / muon",ylabel="Probability / unit bin"); a[0,0].legend(fontsize=8)
a[0,1].set(title="Heed in-gap electron multiplicity",xlabel="Electrons / muon (display 0-300)",ylabel="Probability / 3-electron bin")
a[0,2].set(xscale="log",yscale="log",title="Energy transfer; full tails retained in CSV",xlabel="Transferred energy (keV)",ylabel="Probability / log bin")
for key,label,c in [("clusters_mean","Mean clusters",blue),("electrons_mean","Mean electrons",teal)]:
    a[1,0].plot(scan.energy_gev,scan[key],"o-",label=label,color=c)
a[1,0].set(xscale="log",title="Energy dependence of ionization",xlabel="Muon kinetic energy (GeV)",ylabel="Mean multiplicity"); a[1,0].legend()
a[1,1].plot(scan.energy_gev,scan.invalid_muon_fraction*100,"o-",color=red)
a[1,1].set(xscale="log",title="Muon total linear charge undefined",xlabel="Muon kinetic energy (GeV)",ylabel="Events with >=1 invalid draw (%)",ylim=(98,100.05))
a[1,2].plot(scan.energy_gev,scan.sigma_x_mm,"o-",label="X",color=blue); a[1,2].plot(scan.energy_gev,scan.sigma_y_mm,"s-",label="Y",color=teal)
a[1,2].axhline(sigma,color=orange,ls="--",label="Imposed sigma 0.749 mm")
a[1,2].set(xscale="log",title="Empirical coordinate check",xlabel="Muon kinetic energy (GeV)",ylabel="Sample coordinate sigma (mm)"); a[1,2].legend(fontsize=8)
save(fig,"01_energy_standards","Six muon standard energies | 10,000 perpendicular mu- events per point",
     "Pure R134a, 2 mm, 50.5 kV/cm, 298.15 K, 760 Torr. Heed owns ionization; linear kernel remains unsaturated.\nCoordinate widths are imposed experimental calibration, not microscopic predictions. Display ranges do not remove tails from stored data.")

def prepare_dataset(name,directory):
    s=pd.read_csv(directory/"segments.csv"); p=pd.read_csv(directory/"primaries.csv")
    legacy=pd.read_csv(directory/"legacy_edep.csv"); baseline=pd.read_csv(directory/"baseline.csv")
    assert len(p)==100000 and p.event.nunique()==100000
    assert (s.root_track==1).all(),"Natural single-primary data must have one root ancestor"
    total_segments=len(s); total_gas_events=int(s.event.nunique()); baseline_all=baseline
    extra={}
    first_time=baseline_all.groupby("event").t_ns.transform("min")
    prompt=baseline_all[baseline_all.t_ns-first_time<=window].sort_values("t_ns").drop_duplicates(["event","layer"])
    per_event=prompt.groupby("event").agg(layers=("layer","nunique"),all_valid=("valid_coordinate","min"),any_ambiguous=("ambiguous","max"))
    extra["prompt_four_gas_layer_events"]=int((per_event.layers==4).sum())
    extra["prompt_four_layer_all_valid_events"]=int(((per_event.layers==4)&(per_event.all_valid==1)).sum())
    extra["events_with_ambiguous_prompt_window"]=int((per_event.any_ambiguous==1).sum())
    if (directory/"legacy_analysis_summary.csv").exists():
        extra["literal_analysis_events"]=len(pd.read_csv(directory/"legacy_analysis_summary.csv"))
    if "source_mode" in p and (p.source_mode==2).all():
        assert p.dz.min()>0 and p.z_mm.std()<1e-8
        assert p.x_mm.abs().max()<=200 and p.y_mm.abs().max()<=200
        planes=np.sort(s.gap_center_z_mm.unique())
        assert len(planes)==4
        for center in planes:
            for face in (center-1,center+1):
                length=(face-p.z_mm)/p.dz
                assert (length>0).all()
                assert (np.abs(p.x_mm+length*p.dx)<140+1e-6).all()
                assert (np.abs(p.y_mm+length*p.dy)<140+1e-6).all()
        extra.update(plane_z_mm=float(p.z_mm.iloc[0]),plane_size_mm=[400,400],
                   generated_proposals=int(p.source_proposals.sum()),
                   geometric_acceptance=float(len(p)/p.source_proposals.sum()),straight_four_gap_validation="PASS",
                   gas_plane_z_mm=[float(z) for z in planes])
    # Retain census counts; bound paired-analysis cost with a reproducible uniform sample.
    if len(baseline)>20000:
        baseline=baseline.sample(n=20000,random_state=20261020).sort_values(["event","layer","t_ns"])
        keys=pd.MultiIndex.from_frame(baseline[["event","layer"]]).unique()
        s=s[pd.MultiIndex.from_frame(s[["event","layer"]]).isin(keys)].copy()
    legacy["wx"]=legacy.energy_ev*legacy.x_mm; legacy["wy"]=legacy.energy_ev*legacy.y_mm
    old=legacy.groupby(["event","layer"])[["energy_ev","wx","wy","legacy_trigger"]].agg(
      {"energy_ev":"sum","wx":"sum","wy":"sum","legacy_trigger":"max"})
    old["old_x"]=old.wx/old.energy_ev; old["old_y"]=old.wy/old.energy_ev
    source=p.set_index("event")
    prim=s[(s.parent==0)&(s.pid==13)].copy()
    dz=prim.z1_mm-prim.z0_mm
    f=(prim.gap_center_z_mm-prim.z0_mm)/dz.replace(0,np.nan)
    crossed=prim[(f>=-1e-8)&(f<=1+1e-8)].copy(); f=f.loc[crossed.index].clip(0,1)
    crossed["truth_x"]=crossed.x0_mm+f*(crossed.x1_mm-crossed.x0_mm)
    crossed["truth_y"]=crossed.y0_mm+f*(crossed.y1_mm-crossed.y0_mm)
    crossed["truth_t"]=crossed.t0_ns+f*(crossed.t1_ns-crossed.t0_ns)
    truth=crossed.sort_values("truth_t").drop_duplicates(["event","layer"]).set_index(["event","layer"])
    positive=s.assign(ion=s.edep_ev-s.nonion_ev)
    positive=positive[positive.ion>0].copy(); positive["tmid"]=(positive.t0_ns+positive.t1_ns)/2
    groups={k:g.sort_values("tmid") for k,g in positive.groupby(["event","layer"],sort=True)}
    rows=[]; groups_by_row={}
    for b in baseline.itertuples(index=False):
        key=(int(b.event),int(b.layer)); g=groups[key]
        g=g[(g.tmid>=b.t_ns-1e-7)&(g.tmid<=b.t_ns+window+1e-7)]
        assert len(g)==b.n_segments and g.track.nunique()==b.n_tracks
        w=g.ion.to_numpy(); wsum=w.sum(); xm=(g.x0_mm+g.x1_mm).to_numpy()/2; ym=(g.y0_mm+g.y1_mm).to_numpy()/2
        cx=np.average(xm,weights=w); cy=np.average(ym,weights=w)
        rawx=np.average(g.x1_mm,weights=g.edep_ev); rawy=np.average(g.y1_mm,weights=g.edep_ev)
        wtracks=g.assign(wpx=w*xm,wpy=w*ym).groupby("track")[["ion","wpx","wpy"]].sum()
        txy=(wtracks[["wpx","wpy"]].to_numpy()/wtracks.ion.to_numpy()[:,None])
        spread=float(np.max(np.linalg.norm(txy[:,None,:]-txy[None,:,:],axis=2))) if len(txy)>1 else 0.
        tx=ty=tt=np.nan
        if key in truth.index:
            tr=truth.loc[key]
            if b.t_ns-1e-3<=tr.truth_t<=b.t_ns+window+1e-3:
                tx,ty,tt=tr.truth_x,tr.truth_y,tr.truth_t
        oldx=oldy=oldenergy=np.nan; triggered=False
        if key in old.index:
            o=old.loc[key]; oldx,oldy,oldenergy=o.old_x,o.old_y,o.energy_ev; triggered=bool(o.legacy_trigger)
        noise=np.array([b.x_mm-cx,b.y_mm-cy]) if b.valid_coordinate else rng.normal(0,sigma,2)
        row=dict(event=int(b.event),layer=int(b.layer),time_ns=b.t_ns,n_tracks=int(b.n_tracks),
          n_roots=int(g.root_track.nunique()),n_secondary_tracks=int(g.loc[g.parent>0,"track"].nunique()),
          n_segments=len(g),energy_ev=wsum,old_energy_ev=oldenergy,old_x=oldx,old_y=oldy,
          gas_post_x=rawx,gas_post_y=rawy,gas_mid_x=cx,gas_mid_y=cy,
          baseline_x=b.x_mm,baseline_y=b.y_mm,valid=bool(b.valid_coordinate),ambiguous=bool(b.ambiguous),
          old_noisy_x=oldx+noise[0],old_noisy_y=oldy+noise[1],legacy_trigger=triggered,
          truth_x=tx,truth_y=ty,truth_t=tt,track_spread_mm=spread,
          source_energy_gev=source.loc[key[0],"kinetic_ev"]/1e9)
        rows.append(row); groups_by_row[len(rows)-1]=g
    r=pd.DataFrame(rows)
    r["old_new_centroid_distance_mm"]=np.hypot(r.old_x-r.gas_mid_x,r.old_y-r.gas_mid_y)
    r["same_gas_post_mid_distance_mm"]=np.hypot(r.gas_post_x-r.gas_mid_x,r.gas_post_y-r.gas_mid_y)
    r["old_truth_error_mm"]=np.hypot(r.old_noisy_x-r.truth_x,r.old_noisy_y-r.truth_y)
    r["baseline_truth_error_mm"]=np.hypot(r.baseline_x-r.truth_x,r.baseline_y-r.truth_y)
    r.to_csv(root/(name+"_matched_windows.csv"),index=False)
    hastruth=np.isfinite(r.truth_x+r.truth_y); common=r.valid&hastruth&np.isfinite(r.old_x+r.old_y)
    stats=dict(generated=len(p),gas_events=total_gas_events,segments=total_segments,windows=len(baseline_all),
      analyzed_windows=len(r),sampling="all windows" if len(r)==len(baseline_all) else "uniform 20000-window sample; fixed seed",
      valid_windows=int(baseline_all.valid_coordinate.sum()),multi_track_windows=int(baseline_all.ambiguous.sum()),
      natural_multi_root_windows=int((r.n_roots>1).sum()),legacy_four_layer_events=int(legacy.loc[legacy.legacy_trigger==1,"event"].nunique()),
      truth_crossing_windows=int(hastruth.sum()),common_valid_truth_windows=int(common.sum()),
      old_new_centroid_distance=finite_stats(r.old_new_centroid_distance_mm),
      same_gas_post_mid_distance=finite_stats(r.same_gas_post_mid_distance_mm),
      paired_old_error=finite_stats(r.loc[common,"old_truth_error_mm"]),
      paired_baseline_error=finite_stats(r.loc[common,"baseline_truth_error_mm"]),
      source_mean_dy=float(p.dy.mean()),source_mean_dz=float(p.dz.mean()),
      secondary_spread_mm=finite_stats(r.loc[r.ambiguous,"track_spread_mm"]))
    stats.update(extra)
    summary["datasets"][name]=stats
    return s,p,r,groups_by_row

original_dir=root/"original_y_with_legacy" if (root/"original_y_with_legacy"/"baseline.csv").exists() else root/"csv"
datasets={"original_y":prepare_dataset("original_y",original_dir)}
if (root/"aligned_z_csv"/"baseline.csv").exists():
    datasets["aligned_z"]=prepare_dataset("aligned_z",root/"aligned_z_csv")
if (root/"gaisser_plane_100k_csv"/"baseline.csv").exists():
    datasets["plane_z"]=prepare_dataset("plane_z",root/"gaisser_plane_100k_csv")

# Source and selection: geometry acceptance must not be mistaken for electronics efficiency.
fig,a=plt.subplots(2,2,figsize=(13,9))
for name,(s,p,r,groups) in datasets.items():
    c=blue if name=="original_y" else teal
    label={"original_y":"Existing -Y sphere","aligned_z":"Aligned -Z sphere","plane_z":"Negative-Z plane, forward +Z"}[name]
    a[0,0].hist(p.kinetic_ev/1e9,bins=np.geomspace(.1,1000,80),weights=np.ones(len(p))/len(p),histtype="step",lw=1.5,color=c,label=label)
    a[0,1].hist(np.abs(p.dz),bins=np.linspace(0,1,50),weights=np.ones(len(p))/len(p),histtype="step",color=c,label=label)
    counts=summary["datasets"][name]
    vals=[counts["gas_events"],counts["legacy_four_layer_events"]]
    a[1,0].bar(np.arange(2)+(-.18 if name=="original_y" else .18),vals,width=.34,color=c,label=label)
    a[1,1].bar([label],[counts["multi_track_windows"]],color=c)
a[0,0].set(xscale="log",yscale="log",title="Generated Gaisser-Tang energy sample",xlabel="Kinetic energy (GeV)",ylabel="Probability / logarithmic bin"); a[0,0].legend(fontsize=9)
a[0,1].set(title="Direction relative to RPC normal",xlabel="Absolute Z direction cosine",ylabel="Probability / bin"); a[0,1].legend(fontsize=9)
a[1,0].set(yscale="log",title="Geometry / original selection counts",ylabel="Events out of 100,000",xticks=[0,1],xticklabels=["Any real-gas step","Original four-layer trigger"]); a[1,0].legend(fontsize=9)
a[1,1].set(title="Natural multi-track layer windows",ylabel="50 ns windows")
save(fig,"02_gaisser_source","Gaisser-Tang Geant4 run | source convention and selection",
     "Energy law and spherical source placement follow the existing generator; each event contains one incident muon.\nGas acceptance and original four-layer selection are not detector efficiencies. Root ancestry separates descendants from independent incident particles.")

if "plane_z" in datasets:
    s,p,rr,_=datasets["plane_z"]; info=summary["datasets"]["plane_z"]
    fig,a=plt.subplots(2,2,figsize=(13,9))
    zplanes=np.array(info["gas_plane_z_mm"]); zs=np.r_[info["plane_z_mm"],np.sort(zplanes),zplanes.max()+5]
    for aa,position,component in [(a[0,0],"x_mm","dx"),(a[0,1],"y_mm","dy")]:
        for row in p.sample(60,random_state=42).itertuples():
            coordinates=getattr(row,position)+(zs-row.z_mm)/row.dz*getattr(row,component)
            aa.plot(zs,coordinates,color=blue,alpha=.18,lw=.8)
        for z in zplanes: aa.vlines(z,-140,140,color="#334155",lw=2)
        aa.vlines(info["plane_z_mm"],-200,200,color=red,lw=3,label="400 mm source plane")
        aa.set(title="Accepted straight-ray "+position[0].upper()+"-Z projection",xlabel="Z (mm)",ylabel=position[0].upper()+" (mm)",ylim=(-210,210)); aa.legend(fontsize=8)
    im=a[1,0].hist2d(p.x_mm,p.y_mm,bins=60,range=((-200,200),(-200,200)),cmap="Blues",cmin=1)
    fig.colorbar(im[3],ax=a[1,0],label="Accepted source points")
    a[1,0].set(title=f"Conditional source positions | acceptance {100*info['geometric_acceptance']:.2f}%",xlabel="Source X (mm)",ylabel="Source Y (mm)",aspect="equal")
    prior=datasets["original_y"][1]
    a[1,1].hist(np.degrees(np.arccos(np.clip(-prior.dy,0,1))),bins=np.linspace(0,90,60),density=True,histtype="step",color=orange,label="Unselected existing Gaisser-Tang proposals")
    a[1,1].hist(np.degrees(np.arccos(np.clip(p.dz,0,1))),bins=np.linspace(0,90,60),density=True,histtype="step",color=blue,label="Plane sample after 4-gap ray selection")
    a[1,1].set(title="Geometric angular conditioning",xlabel="Angle from source forward axis (degrees)",ylabel="Probability density / degree"); a[1,1].legend(fontsize=8)
    save(fig,"00_plane_geometry","Requested planar source | negative-Z side, positive-Z flight, four-gap ray constraint",
         "Plane size matches the 400 x 400 mm RPC envelope; each accepted ray intersects both faces of all four 280 x 280 mm active gas regions.\nThis is geometry-conditioned source sampling, not absolute cosmic flux. Multiple scattering / stopping may still prevent actual four-layer passage.")

# Matched coordinate comparisons on identical simulated records.
for name,(s,p,r,groups) in datasets.items():
    fig,a=plt.subplots(2,3,figsize=(15,9)); csel=r.valid&np.isfinite(r.truth_x+r.old_x)
    for values,label,c in [(r.loc[csel,"old_noisy_x"]-r.loc[csel,"truth_x"],"Original scoring Edep",orange),
                          (r.loc[csel,"baseline_x"]-r.loc[csel,"truth_x"],"Current gas baseline",blue)]:
        a[0,0].hist(values,bins=np.linspace(-15,15,100),histtype="step",color=c,lw=1.4,label=label)
    a[0,0].set(title=f"Matched single-track X residual | n={csel.sum()}",xlabel="Reconstructed X - actual primary crossing (mm)",ylabel="Windows"); a[0,0].legend(fontsize=8)
    a[0,1].hist(r.old_new_centroid_distance_mm.dropna(),bins=np.geomspace(1e-4,300,90),histtype="step",color=orange,label="Original region vs real-gas midpoint")
    a[0,1].hist(r.same_gas_post_mid_distance_mm,bins=np.geomspace(1e-4,300,90),histtype="step",color=teal,label="Same gas: post-step vs midpoint")
    a[0,1].set(xscale="log",yscale="log",title="Centroid shifts before readout noise",xlabel="2D centroid difference (mm)",ylabel="Windows / log bin"); a[0,1].legend(fontsize=8)
    a[0,2].scatter(r.loc[csel,"source_energy_gev"],r.loc[csel,"old_new_centroid_distance_mm"],s=4,alpha=.25,color=blue,rasterized=True)
    a[0,2].set(xscale="log",yscale="log",title="Centroid shift vs generated muon energy",xlabel="Kinetic energy (GeV)",ylabel="Original-to-gas centroid difference (mm)")
    a[1,0].hist(r.n_tracks,bins=np.arange(.5,min(30,r.n_tracks.max())+1.5),color=blue,alpha=.8)
    a[1,0].set(yscale="log",title="Track multiplicity in a 50 ns window",xlabel="Tracks (includes secondaries)",ylabel="Windows")
    multi=r[r.ambiguous]
    a[1,1].hist(multi.track_spread_mm,bins=np.geomspace(1e-5,300,80),color=red,alpha=.8)
    a[1,1].set(xscale="log",title="Natural secondary-track spatial spread",xlabel="Maximum separation of track deposit centroids (mm)",ylabel="Multi-track windows / log bin")
    a[1,2].bar(["Original Edep\n(one coordinate)","Current valid\n(single-track)","Current ambiguous\n(no coordinate)"],
               [int(np.isfinite(r.old_x).sum()),int(r.valid.sum()),int(r.ambiguous.sum())],color=[orange,blue,red])
    a[1,2].set(title="Output coverage in analyzed window subset",ylabel="Layer windows")
    save(fig,"03_coordinates_"+name,"Matched coordinate comparison | "+name.replace("_"," "),
         "Common-noise comparison on identical valid windows. Original scoring region and real gas differ; same-gas post/midpoint contrast isolates step-position effects.\nLarge runs use a fixed-seed 20,000-window paired subset; selection counts use all records. No Garfield charge-weighted coordinate is claimed.")

# Natural multi-track examples: all descendants share one primary root.
chosen_name="plane_z" if "plane_z" in datasets else ("aligned_z" if "aligned_z" in datasets else "original_y")
s,p,r,groups=datasets[chosen_name]
multi=r[r.ambiguous&np.isfinite(r.truth_x)].sort_values("track_spread_mm",ascending=False)
examples=[]
if len(multi):
    candidates=[multi.index[0],multi.index[len(multi)//2]]
    fig,a=plt.subplots(1,2,figsize=(13,6.5))
    for aa,index in zip(a,candidates):
        row=r.loc[index]; g=groups[index]
        for track,tg in g.groupby("track"):
            pid=int(tg.pid.iloc[0]); aa.plot(np.r_[tg.x0_mm.to_numpy(),tg.x1_mm.iloc[-1]],np.r_[tg.y0_mm.to_numpy(),tg.y1_mm.iloc[-1]],".-",ms=3,lw=.8,label=f"track {track}; PDG {pid}")
        aa.scatter(row.truth_x,row.truth_y,s=100,marker="*",c="black",label="Primary midplane crossing",zorder=6)
        aa.scatter(row.old_x,row.old_y,s=65,marker="x",c=orange,label="Original Edep centroid",zorder=6)
        aa.scatter(row.gas_mid_x,row.gas_mid_y,s=60,facecolors="none",edgecolors=blue,label="Gas centroid (diagnostic)",zorder=6)
        aa.set(title=f"Event {int(row.event)}, layer {int(row.layer)} | {int(row.n_tracks)} tracks",xlabel="X (mm)",ylabel="Y (mm)",aspect="equal")
        aa.legend(fontsize=7,loc="best")
        examples.append(row.to_dict()); g.to_csv(root/f"natural_example_{int(row.event)}_layer{int(row.layer)}.csv",index=False)
    save(fig,"04_natural_multitrack","Natural multi-track examples | one incident muon, multiple descendants",
         "Both examples are inside 50 ns; the current program flags them ambiguous and emits NaN, not a reconstructed multi-hit result.\nTrack identities and the primary crossing shown here are simulation truth, unavailable to a real detector.")

# Controlled overlap: independent Geant4 events, same layer, each individually unambiguous.
pool=r[r.valid&np.isfinite(r.truth_x+r.truth_y)].copy()
assert len(pool)>100,"Insufficient pure windows for independent-event overlap"
by_layer={int(layer):g.reset_index(drop=True) for layer,g in pool.groupby("layer") if len(g)>1}
layers=np.array(list(by_layer)); probs=np.array([len(by_layer[l]) for l in layers],float); probs/=probs.sum()
n=5000; pair_rows=[]
for _ in range(n):
    layer=int(rng.choice(layers,p=probs)); g=by_layer[layer]
    i,j=rng.integers(0,len(g),2)
    while g.iloc[i].event==g.iloc[j].event: j=int(rng.integers(0,len(g)))
    aa,bb=g.iloc[i],g.iloc[j]
    pair_rows.append([aa.event,bb.event,layer,aa.gas_mid_x,aa.gas_mid_y,bb.gas_mid_x,bb.gas_mid_y,aa.energy_ev,bb.energy_ev])
pairs=np.array(pair_rows); xa,ya,xb,yb=pairs[:,3],pairs[:,4],pairs[:,5],pairs[:,6]
separation=np.hypot(xa-xb,ya-yb)
end_sx=sigma/(np.sqrt(2)*axdelay); end_sy=sigma/(np.sqrt(2)*aydelay)
noise=rng.normal(0,1,(n,2,4))*np.array([end_sx,end_sx,end_sy,end_sy])[None,None,:]
matching_radius=3*sigma
arrival_tolerance=3*np.sqrt((end_sx**2+end_sy**2)/2)
overlay=[]; case_frames=[]; example_pair=None
for dt in (0,10,25,49,50,51,75,100):
    t=np.array([0.,float(dt)])[None,:]
    xs=np.stack([xa,xb],axis=1); ys=np.stack([ya,yb],axis=1)
    xl=160+t+(xs-1.062)/(2*axdelay)+noise[:,:,0]; xr=160+t-(xs-1.062)/(2*axdelay)+noise[:,:,1]
    yd=160+t+ys/(2*aydelay)+noise[:,:,2]; yu=160+t-ys/(2*aydelay)+noise[:,:,3]
    minpulse=np.min(np.abs(np.stack([xl[:,0]-xl[:,1],xr[:,0]-xr[:,1],yd[:,0]-yd[:,1],yu[:,0]-yu[:,1]],axis=1)),axis=1)
    idealx=axdelay*(xl-xr)+1.062; idealy=aydelay*(yd-yu)
    idealdist=np.hypot(idealx-xs,idealy-ys)
    ncandidates=np.zeros(n,dtype=int); true_pass=np.zeros(n,dtype=bool); candidate_coordinates={}
    inferred_x=np.full((n,2),np.nan); inferred_y=np.full((n,2),np.nan)
    # All endpoint pulses assumed extracted. Eight assignments, no charge information.
    for px,py,xy in itertools.product((False,True),repeat=3):
        xrr=xr[:,::-1] if px else xr; yuu=yu[:,::-1] if py else yu
        cx=axdelay*(xl-xrr)+1.062; cy=aydelay*(yd-yuu)
        tx=(xl+xrr-320)/2; ty=(yd+yuu-320)/2
        if xy: cy=cy[:,::-1]; ty=ty[:,::-1]
        valid=np.all(np.abs(tx-ty)<=arrival_tolerance,axis=1)&np.all(np.abs(cx)<=140,axis=1)&np.all(np.abs(cy)<=140,axis=1)
        valid &= np.all((tx>=-arrival_tolerance)&(tx<=dt+arrival_tolerance),axis=1)
        valid &= np.all((ty>=-arrival_tolerance)&(ty<=dt+arrival_tolerance),axis=1)
        ncandidates+=valid
        inferred_x[valid]=cx[valid]; inferred_y[valid]=cy[valid]
        if not(px or py or xy): true_pass=valid.copy()
        candidate_coordinates[(px,py,xy)]=(valid,cx,cy)
    unique=ncandidates==1
    inferred_x[~unique]=np.nan; inferred_y[~unique]=np.nan
    direct=np.hypot(inferred_x-xs,inferred_y-ys)
    swapped=np.hypot(inferred_x[:,::-1]-xs,inferred_y[:,::-1]-ys)
    swap_assignment=np.nansum(swapped**2,axis=1)<np.nansum(direct**2,axis=1)
    inferred_errors=np.where(swap_assignment[:,None],swapped,direct)
    matched_both=unique&np.all(inferred_errors<=matching_radius,axis=1)
    same=dt<=window
    ea,eb=pairs[:,7],pairs[:,8]
    mergedx=(ea*xa+eb*xb)/(ea+eb); mergedy=(ea*ya+eb*yb)/(ea+eb)
    olderr=np.stack([np.hypot(mergedx-xa,mergedy-ya),np.hypot(mergedx-xb,mergedy-yb)],axis=1)
    row=dict(dt_ns=dt,pairs=n,same_50ns_window=same,
      original_coordinates_per_pair=1 if same else 2,current_coordinates_per_pair=0 if same else 2,
      ideal_truth_assisted_coordinates_per_pair=2,
      ideal_endpoint_minsep_gt2ns=float(np.mean(minpulse>2)),
      ideal_endpoint_minsep_gt5ns=float(np.mean(minpulse>5)),
      unique_correct_pairing=float(np.mean((ncandidates==1)&true_pass)),
      prototype_coordinates_per_pair=float(2*np.mean(unique)),
      prototype_both_truthmatched_fraction=float(np.mean(matched_both)),
      prototype_both_matched_and_endpoint_gt2ns=float(np.mean(matched_both&(minpulse>2))),
      prototype_unique_error_rms_mm=rms(inferred_errors[unique]),
      unique_correct_and_endpoint_gt2ns=float(np.mean(unique&true_pass&(minpulse>2))),
      ambiguous_pairing=float(np.mean(ncandidates>1)),
      no_pairing=float(np.mean(ncandidates==0)),
      false_unique_pairing=float(np.mean((ncandidates==1)&~true_pass)),
      merged_centroid_truth_rms_mm=rms(olderr) if same else 0,
      ideal_truth_assisted_error_rms_mm=rms(idealdist),
      matching_radius_mm=matching_radius,arrival_tolerance_ns=arrival_tolerance)
    overlay.append(row)
    case=pd.DataFrame({"event_a":pairs[:,0].astype(int),"event_b":pairs[:,1].astype(int),"layer":pairs[:,2].astype(int),
      "dt_ns":dt,"separation_mm":separation,"min_endpoint_separation_ns":minpulse,
      "pairing_candidates":ncandidates,"true_pairing_passes":true_pass,
      "old_centroid_error_a_mm":olderr[:,0],"old_centroid_error_b_mm":olderr[:,1],
      "ideal_error_a_mm":idealdist[:,0],"ideal_error_b_mm":idealdist[:,1]})
    for j in range(2):
        case[f"prototype_x{j+1}_mm"]=inferred_x[:,j]; case[f"prototype_y{j+1}_mm"]=inferred_y[:,j]
        case[f"prototype_matched_error{j+1}_mm"]=inferred_errors[:,j]
    case["prototype_both_truthmatched"]=matched_both
    case_frames.append(case)
    if dt==25:
        choice=np.where((ncandidates==1)&true_pass&(minpulse>5)&(separation>30))[0]
        if len(choice):
            k=int(choice[0]); example_pair=dict(index=k,dt=dt,xs=xs[k],ys=ys[k],xl=xl[k],xr=xr[k],yd=yd[k],yu=yu[k],
              idealx=idealx[k],idealy=idealy[k],mergedx=mergedx[k],mergedy=mergedy[k],details=case.iloc[k].to_dict())
over=pd.DataFrame(overlay); over.to_csv(root/"overlay_summary.csv",index=False)
pd.concat(case_frames,ignore_index=True).to_csv(root/"overlay_cases.csv",index=False)
summary["overlay"]={"source_dataset":chosen_name,"pairs_per_time_offset":n,"offsets_ns":over.dt_ns.tolist(),
  "pairing_time_tolerance_ns":arrival_tolerance,"pulse_thresholds_ns":[2,5],
  "caveat":"Artificial independent-event mixing; all endpoint pulses assumed extracted. Pulse-separation thresholds and pairing cut are sensitivity assumptions, not calibrated detector performance."}
fig,a=plt.subplots(2,3,figsize=(15,9))
same=over.same_50ns_window
a[0,0].plot(over.dt_ns,over.original_coordinates_per_pair,"o-",color=orange,label="Edep per 50 ns frame")
a[0,0].plot(over.dt_ns,over.current_coordinates_per_pair,"s-",color=blue,label="Current ambiguity policy")
a[0,0].plot(over.dt_ns,over.prototype_coordinates_per_pair,"^-",color=purple,label="Ideal-pulse pairing prototype")
a[0,0].plot(over.dt_ns,over.ideal_truth_assisted_coordinates_per_pair,":",color=teal,label="Truth-assisted ideal")
a[0,0].axvline(50,color=red,ls="--"); a[0,0].set(title="Number of coordinate outputs",xlabel="Injected arrival offset (ns)",ylabel="Coordinates / two-muon pair",ylim=(-.1,2.3)); a[0,0].legend(fontsize=8)
a[0,1].hist(separation,bins=50,color=blue,alpha=.8); a[0,1].set(title="Independent-pair position separation",xlabel="Distance between gas centroids (mm)",ylabel="Pairs")
a[0,2].scatter(separation,np.sqrt(np.mean(olderr**2,axis=1)),s=3,alpha=.25,color=orange,rasterized=True)
a[0,2].set(title="Single merged-centroid loss of information",xlabel="True hit separation (mm)",ylabel="RMS distance to both known hits (mm)")
a[1,0].plot(over.dt_ns,100*over.unique_correct_pairing,"o-",color=teal,label="Unique correct")
a[1,0].plot(over.dt_ns,100*over.ambiguous_pairing,"s-",color=orange,label="Multiple assignments")
a[1,0].plot(over.dt_ns,100*(over.no_pairing+over.false_unique_pairing),"^-",color=red,label="No / false unique")
a[1,0].set(title="Ideal extracted-pulse pairing toy model",xlabel="Injected arrival offset (ns)",ylabel="Pairs (%)"); a[1,0].legend(fontsize=8)
a[1,1].plot(over.dt_ns,100*over.ideal_endpoint_minsep_gt2ns,"o-",color=blue,label="Every endpoint gap > 2 ns")
a[1,1].plot(over.dt_ns,100*over.ideal_endpoint_minsep_gt5ns,"s-",color=purple,label="Every endpoint gap > 5 ns")
a[1,1].set(title="Necessary pulse-separation condition only",xlabel="Injected arrival offset (ns)",ylabel="Pairs meeting illustrative condition (%)"); a[1,1].legend(fontsize=8)
a[1,2].axis("off"); a[1,2].text(0,.95,"READING THIS STUDY\n\n50 ns is the supplied grouping window.\nIt does not specify electronics dead time.\n\nCurrent software rejects unresolved multi-track\nwindows; it does not reconstruct two hits.\n\nThe toy model assumes extracted endpoint pulses.\nTruth-assisted coordinates are an ideal reference,\nnot measured multi-hit efficiency.",va="top",fontsize=11)
save(fig,"05_independent_overlap","Independent two-muon overlap | 5,000 pairs per injected time offset",
     "Pairs are built from different Geant4 event IDs in the same layer; no incident rate is inferred.\nPairing enumerates endpoint assignments using nominal timing only. Values beyond 50 ns characterize separately grouped arrivals / a hypothetical wider buffer.")

if example_pair is not None:
    e=example_pair; fig,a=plt.subplots(1,2,figsize=(13,6.5))
    a[0].scatter(e["xs"],e["ys"],s=150,marker="*",color="black",label="Known separate gas centroids")
    a[0].scatter(e["idealx"],e["idealy"],s=80,facecolors="none",edgecolors=teal,label="Unique ideal-pulse assignment")
    a[0].scatter([e["mergedx"]],[e["mergedy"]],s=100,marker="x",color=orange,label="One merged Edep coordinate")
    for i in range(2): a[0].annotate(f"Muon {i+1}",(e["xs"][i],e["ys"][i]),xytext=(5,6),textcoords="offset points")
    a[0].set(title="Two positions inside one 50 ns window",xlabel="X (mm)",ylabel="Y (mm)",aspect="equal"); a[0].legend(fontsize=8)
    for j,(key,label) in enumerate([("xl","XL"),("xr","XR"),("yd","YD"),("yu","YU")]):
        for i,c in enumerate([blue,purple]): a[1].vlines(e[key][i],j-.25,j+.25,color=c,lw=2,label=f"Muon {i+1}" if j==0 else None)
    a[1].set(title="Ideal extracted endpoint pulses (not waveforms)",xlabel="Nominal endpoint time (ns)",yticks=range(4),yticklabels=["XL","XR","YD","YU"]); a[1].legend()
    save(fig,"06_overlap_example","A potentially separable two-muon example | injected offset 25 ns",
         "This example has one accepted assignment under the stated ideal pulse-extraction and timing-consistency assumptions.\nThe existing digitizer still flags this combined window ambiguous; this illustration is not an implemented electronics reconstruction.")
    summary["overlay_example"]=e["details"]

fig,a=plt.subplots(2,2,figsize=(13,9))
for offset,c in [(0,orange),(10,blue),(25,teal)]:
    case=case_frames[list(over.dt_ns).index(offset)]
    errors=case[["prototype_matched_error1_mm","prototype_matched_error2_mm"]].to_numpy().ravel()
    errors=errors[np.isfinite(errors)]
    a[0,0].hist(errors,bins=np.geomspace(.01,400,90),histtype="step",color=c,label=f"dt={offset} ns; {len(errors)//2} unique pairs")
a[0,0].set(xscale="log",yscale="log",title="Prototype coordinate errors | unique assignments",xlabel="Best two-to-two truth-matched error (mm)",ylabel="Coordinates / log bin"); a[0,0].legend(fontsize=8)
a[0,1].plot(over.dt_ns,100*over.prototype_both_truthmatched_fraction,"o-",color=teal,label="Both errors < 3 sigma")
a[0,1].plot(over.dt_ns,100*over.prototype_both_matched_and_endpoint_gt2ns,"s-",color=blue,label="Also every endpoint gap > 2 ns")
a[0,1].set(title="Toy-model truth-match fraction, not detector efficiency",xlabel="Injected time offset (ns)",ylabel="All injected pairs (%)"); a[0,1].legend(fontsize=8)
case=case_frames[list(over.dt_ns).index(25)]; unique_case=case[case.pairing_candidates==1]
a[1,0].scatter(unique_case.separation_mm,unique_case.prototype_matched_error1_mm,s=4,alpha=.3,color=blue,rasterized=True)
a[1,0].set(yscale="log",title="Unique-pair reconstruction at dt=25 ns",xlabel="Input hit separation (mm)",ylabel="Matched coordinate error (mm)")
a[1,1].axis("off"); a[1,1].text(0,.95,
    "PROTOTYPE ASSUMPTIONS\n\nAll two-particle pulses exist and are extracted.\nNo avalanche amplitude / threshold / dead time.\nAll endpoint permutations are enumerated.\nX/Y mean times agree within 3.674 ns.\nOnly unique hypotheses emit coordinates.\nTruth matching is evaluation only (2.247 mm cut).\n\nFalse unique assignments remain in the error tail.\nThis algorithm is separate from rpc_digitize.",va="top",fontsize=11)
save(fig,"07_pairing_prototype","Multi-hit coordinate pairing prototype | conditional ideal-input study",
     "The algorithm does not use truth to select a pairing; truth is used afterward to evaluate the emitted two coordinates.\nAssumed endpoint jitter comes from single-hit spatial sigma; real electronics pulse extraction and correlated timing are uncalibrated.")
pdf.close()
(root/"comparison_summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False),encoding="utf-8")
print(json.dumps(summary,indent=2,allow_nan=False))
print("Figures:",figures)
