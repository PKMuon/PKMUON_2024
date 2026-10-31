"""Spatial multi-position baseline and held-out distribution interpolation QA."""
import base64
import html
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from response_interpolation import MuonStandardInterpolator

base = Path(sys.argv[1]).resolve()
out = base/"spatial_interpolation"
figures = out/"figures"; figures.mkdir(exist_ok=True)
plt.rcParams.update({"font.family":"DejaVu Sans", "font.size":10, "axes.spines.top":False,
    "axes.spines.right":False, "axes.grid":True, "grid.alpha":.18,
    "figure.facecolor":"#f8fafc", "savefig.facecolor":"#f8fafc"})
blue, teal, orange, red = "#2563eb", "#0d9488", "#d97706", "#dc2626"
sigma = .749
pdf = PdfPages(figures/"spatial_interpolation_figures.pdf")

def save(fig, name, title, footer):
    fig.suptitle(title, x=.05, y=.985, ha="left", fontsize=17, fontweight="bold")
    fig.text(.05, .015, footer, fontsize=8.5, color="#475569", va="bottom")
    fig.tight_layout(rect=(.02,.14,.99,.94), h_pad=2, w_pad=2)
    fig.savefig(figures/(name+".png"), dpi=160)
    fig.savefig(figures/(name+".pdf")); pdf.savefig(fig); plt.close(fig)

old = pd.read_csv(base/"gaisser_plane_100k_csv"/"baseline.csv")
old = old.sort_values(["event", "layer", "t_ns"]).reset_index(drop=True)
old["window_index"] = old.groupby("event").cumcount()
keys = ["event", "layer", "window_index"]
spatial_rows = []; tables = {}
for factor in (3,5,10):
    table = pd.read_csv(out/f"gaisser_plane_spatial_{factor}sigma.csv")
    tables[factor] = table
    groups = table.groupby(keys, sort=True)
    windows = groups.agg(n_positions=("n_positions","first"), n_tracks=("window_n_tracks","first"),
        weight=("ionizing_edep_ev","sum"), n_segments=("n_segments","sum"),
        actual_positions=("position_index","size"), window_t_ns=("window_t_ns","first")).reset_index()
    check = windows.merge(old[keys+["ionizing_edep_ev","n_segments","n_tracks","t_ns"]],
                          on=keys, validate="one_to_one", suffixes=("","_old"))
    assert len(check) == len(old) == 397659
    assert np.array_equal(check.n_positions, check.actual_positions)
    assert np.array_equal(check.n_segments, check.n_segments_old)
    assert np.array_equal(check.n_tracks, check.n_tracks_old)
    assert np.allclose(check.weight, check.ionizing_edep_ev, rtol=1e-12, atol=1e-7)
    assert np.allclose(check.window_t_ns, check.t_ns, atol=1e-10, rtol=0)
    assert (table.valid_coordinate == 1).all() and (table.ideal_spatial_model == 1).all()
    assert (table.readout_pairing_validated == 0).all()
    # Check actual emitted component centers, not only the track-center graph.
    min_distances=[]
    for _, g in table[table.n_positions>1].groupby(keys, sort=False):
        xy=g[["centroid_x_mm","centroid_y_mm"]].to_numpy()
        distances=np.hypot((xy[:,None]-xy[None,:])[:,:,0],(xy[:,None]-xy[None,:])[:,:,1])
        distances[np.diag_indices(len(xy))]=np.inf
        min_distances.append(distances.min())
    assert min(min_distances)>factor*sigma
    multi = windows.n_tracks>1
    row={"factor_sigma":factor,"threshold_mm":factor*sigma,"windows":len(windows),
         "positions":len(table),"multi_track_windows":int(multi.sum()),
         "multi_position_windows":int((windows.n_positions>1).sum()),
         "two_position_windows":int((windows.n_positions==2).sum()),
         "three_or_more_position_windows":int((windows.n_positions>=3).sum()),
         "close_multi_track_single_position_windows":int((multi&(windows.n_positions==1)).sum()),
         "max_positions_per_window":int(windows.n_positions.max()),
         "sigma_x_check_mm":float((table.x_mm-table.centroid_x_mm).std()),
         "sigma_y_check_mm":float((table.y_mm-table.centroid_y_mm).std())}
    spatial_rows.append(row)
    windows.to_csv(out/f"spatial_{factor}sigma_windows.csv",index=False)
spatial = pd.DataFrame(spatial_rows); spatial.to_csv(out/"spatial_summary.csv",index=False)
five = tables[5]

fig, a = plt.subplots(2,2,figsize=(13,9))
a[0,0].plot(spatial.threshold_mm, spatial.multi_position_windows,"o-",color=blue)
a[0,0].set(xlabel="Assumed separation threshold (mm)",ylabel="Layer windows",title="Spatially separated multiple positions")
a[0,1].bar(["Old rejection\npolicy", "5-sigma spatial\napproximation"], [int((old.valid_coordinate==1).sum()),len(five)],color=[orange,teal])
a[0,1].set(ylabel="Output coordinates",title="Outputs from the same 397,659 windows")
window5 = five.drop_duplicates(keys)
counts=window5.n_positions.value_counts().sort_index()
a[1,0].bar(counts.index,counts.values,color=blue)
a[1,0].set(yscale="log",xlabel="Output positions / window",ylabel="Windows",title="5-sigma position multiplicity")
for col,label,c in [("x_mm","X",blue),("y_mm","Y",teal)]:
    center="centroid_"+col
    a[1,1].hist(five[col]-five[center],bins=100,range=(-3.5,3.5),histtype="step",density=True,color=c,label=label)
a[1,1].set(xlabel="Output - component centroid (mm)",ylabel="Density / mm",title="Imposed 0.749 mm coordinate noise")
a[1,1].legend()
save(fig,"01_spatial_scan","Same-window separated tracks now emit multiple spatial positions",
     "Ideal spatial-response approximation: simulated track centroids, close-component merging, independent calibrated Gaussian noise.\nCounts are not efficiency or validated delay-line pairing; no threshold/amplitude, dead time, or waveform extraction is included.")

# Plot an actual two-position same-window example and the old merged position.
examples=five[(five.n_positions==2)&(five.window_n_tracks==2)]
seps=examples.groupby(keys).agg(xmin=("centroid_x_mm","min"),xmax=("centroid_x_mm","max"),
                              ymin=("centroid_y_mm","min"),ymax=("centroid_y_mm","max"))
seps["separation_mm"]=np.hypot(seps.xmax-seps.xmin,seps.ymax-seps.ymin)
chosen=seps.iloc[(seps.separation_mm-30).abs().argsort()[:1]].reset_index().iloc[0]
select=np.logical_and.reduce([five[k]==chosen[k] for k in keys]); example=five[select]
example.to_csv(out/"natural_two_position_example.csv",index=False)
segments=pd.read_csv(base/"gaisser_plane_100k_csv"/"segments.csv")
segments=segments[(segments.event==chosen.event)&(segments.layer==chosen.layer)&
    ((segments.t0_ns+segments.t1_ns)/2>=example.window_t_ns.iloc[0])&
    ((segments.t0_ns+segments.t1_ns)/2<=example.window_t_ns.iloc[0]+50)]
raw=old[np.logical_and.reduce([old[k]==chosen[k] for k in keys])].iloc[0]
weights=example.ionizing_edep_ev.to_numpy()
merged=np.average(example[["centroid_x_mm","centroid_y_mm"]],axis=0,weights=weights)
fig,a=plt.subplots(1,2,figsize=(13,6))
for tid,g in segments.groupby("track"):
    a[0].plot(np.r_[g.x0_mm.to_numpy(),g.x1_mm.iloc[-1]],np.r_[g.y0_mm.to_numpy(),g.y1_mm.iloc[-1]],".-",label=f"Track {tid}; PDG {g.pid.iloc[0]}")
a[0].scatter(example.centroid_x_mm,example.centroid_y_mm,marker="*",s=140,c="black",label="Separate deposit centers")
a[0].scatter(example.x_mm,example.y_mm,s=80,facecolors="none",edgecolors=teal,label="Spatial-response outputs")
a[0].scatter(*merged,marker="x",s=100,c=orange,label="One gas Edep centroid")
a[0].set(xlabel="X (mm)",ylabel="Y (mm)",title=f"Event {int(chosen.event)}, layer {int(chosen.layer)}; 50 ns",aspect="equal")
a[0].legend(fontsize=8)
a[1].axis("off")
lines=["5-sigma spatial baseline",f"Threshold = {5*sigma:.3f} mm",f"Separation = {chosen.separation_mm:.3f} mm",
       f"Tracks = {int(example.window_n_tracks.iloc[0])}",f"Position outputs = {len(example)}","",
       "Original conservative baseline: NaN", "New spatial model: two measured positions", "",
       "This is a simulation-assisted approximation.", "Actual four-end pairing remains uncalibrated."]
a[1].text(0,.95,"\n".join(lines),va="top",fontsize=12)
save(fig,"02_natural_two_positions","Actual Geant4 secondary-track example | two positions in one window",
     "Track identities are used only by this ideal simulation response. Extended deposits on one track still have one centroid.\nThe original single-centroid and conservative-rejection models are retained for comparison.")
del segments

# Independent double-muon overlays: same input pairs as the previous report.
cases=pd.read_csv(base/"overlay_cases.csv")
overlay_rows=[]
for factor in (3,5,10):
    for dt,g in cases.groupby("dt_ns",sort=True):
        separated=(g.separation_mm>factor*sigma)|(g.dt_ns>50)
        overlay_rows.append({"factor_sigma":factor,"threshold_mm":factor*sigma,"dt_ns":dt,
          "pairs":len(g),"two_positions_fraction":float(separated.mean()),
          "spatial_coordinates_per_pair":float(1+separated.mean()),
          "old_edep_coordinates_per_pair":1 if dt<=50 else 2,
          "four_end_unique_fraction":float((g.pairing_candidates==1).mean())})
overlay=pd.DataFrame(overlay_rows);overlay.to_csv(out/"spatial_overlay_summary.csv",index=False)
fig,a=plt.subplots(1,2,figsize=(13,6))
for factor,c in [(3,blue),(5,teal),(10,orange)]:
    g=overlay[overlay.factor_sigma==factor]
    a[0].plot(g.dt_ns,g.spatial_coordinates_per_pair,"o-",color=c,label=f"Spatial {factor} sigma")
prototype=overlay[overlay.factor_sigma==5]
a[0].plot(prototype.dt_ns,2*prototype.four_end_unique_fraction,"s--",color=red,label="Previous four-end unique pairing")
a[0].set(xlabel="Injected offset (ns)",ylabel="Coordinates / pair",title="Spatial model vs endpoint pairing",ylim=(-.1,2.2));a[0].legend(fontsize=8)
g=cases[cases.dt_ns==0]
a[1].hist(g.separation_mm,bins=80,color=blue,alpha=.7)
for factor,c in [(3,blue),(5,teal),(10,orange)]:a[1].axvline(factor*sigma,color=c,label=f"{factor} sigma = {factor*sigma:.3f} mm")
a[1].set(xlabel="Independent-pair separation (mm)",ylabel="Pairs",title="Spatial separability at simultaneous arrival");a[1].legend(fontsize=8)
save(fig,"03_spatial_overlays","50 ns coincidence does not force separated spatial signals into one position",
     "User-authorized spatial approximation: far-separated components emit two coordinates even at dt=0.\nA four-end delay-line device can still have X/Y association degeneracies; the two curves answer different questions.")

# Distribution-level cross-validation. Never include a held-out energy in training.
train=MuonStandardInterpolator(base/"energies")
validation=MuonStandardInterpolator(out/"validation_energies")
features={"n_clusters":"Ionization clusters", "n_electrons":"In-gap electrons", "energy_transfer_ev":"Transferred energy (eV)"}
methods=("nearest","linear_mixture","log_mixture")
colors={"nearest":orange,"linear_mixture":blue,"log_mixture":teal}

def quantile(values,weights,p):
    order=np.argsort(values);v=np.asarray(values)[order];w=np.asarray(weights)[order]
    return np.interp(p,np.cumsum(w)/w.sum(),v)

def ks_distance(actual,predicted,weights):
    a=np.sort(np.asarray(actual)); order=np.argsort(predicted)
    b=np.asarray(predicted)[order];w=np.asarray(weights)[order];cw=np.r_[0.,np.cumsum(w)]
    grid=np.unique(np.r_[a,b]);fa=np.searchsorted(a,grid,side="right")/len(a)
    fb=cw[np.searchsorted(b,grid,side="right")]
    return float(np.max(np.abs(fa-fb)))

def compare(actual, predicted,weights,feature):
    a=actual[feature].to_numpy();p=predicted[feature].to_numpy()
    am=float(a.mean());pm=float(weights@p)
    aq=np.quantile(a,[.1,.5,.9,.99]);pq=quantile(p,weights,[.1,.5,.9,.99])
    iqr=float(np.quantile(a,.75)-np.quantile(a,.25))
    return {"actual_mean":am,"predicted_mean":pm,"mean_relative_error":(pm/am-1),
      "actual_median":float(aq[1]),"predicted_median":float(pq[1]),
      "median_relative_error":float(pq[1]/aq[1]-1),
      "actual_p90":float(aq[2]),"predicted_p90":float(pq[2]),
      "p90_relative_error":float(pq[2]/aq[2]-1),
      "actual_p99":float(aq[3]),"predicted_p99":float(pq[3]),
      "p99_relative_error":float(pq[3]/aq[3]-1),
      "ks_distance":ks_distance(a,p,weights),
      "wasserstein":float(wasserstein_distance(a,p,v_weights=weights)),
      "wasserstein_over_iqr":float(wasserstein_distance(a,p,v_weights=weights)/iqr)}

rows=[];joint_rows=[];examples_cv={}
tests=[("leave_one_out",e,train.frames[e],MuonStandardInterpolator(base/"energies",exclude_gev=e))
       for e in train.energies[1:-1]]
tests += [("independent_midpoint",e,f,train) for e,f in validation.frames.items()]
for kind,e,actual,library in tests:
    examples_cv[(kind,e)]=(actual,library)
    assert e not in library.frames
    for method in methods:
        predicted,weights=library.distribution(e,method)
        low,high,w=library.bracket(e,method)
        for feature in features:
            row={"validation":kind,"energy_gev":e,"feature":feature,"method":method,
                 "low_gev":low,"high_gev":high,"high_weight":w,"actual_events":len(actual)}
            row.update(compare(actual,predicted,weights,feature)); rows.append(row)
        # Evaluate correlations at the distribution level as joint CDFs on a
        # held-out quantile grid. Atomic event sampling alone does not guarantee
        # that interpolated joint statistics match the target energy.
        for fx,fy in (("n_clusters","n_electrons"),("n_electrons","energy_transfer_ev")):
            bins=[np.r_[-np.inf,np.unique(np.quantile(actual[f],np.linspace(.05,.95,19))),np.inf]
                  for f in (fx,fy)]
            observed=np.histogram2d(actual[fx],actual[fy],bins=bins)[0]/len(actual)
            expected=np.histogram2d(predicted[fx],predicted[fy],bins=bins,weights=weights)[0]
            joint_rows.append({"validation":kind,"energy_gev":e,"method":method,
              "feature_x":fx,"feature_y":fy,
              "joint_cdf_max_distance":float(np.max(np.abs(observed.cumsum(0).cumsum(1)-expected.cumsum(0).cumsum(1))))})
metrics=pd.DataFrame(rows);metrics.to_csv(out/"interpolation_cross_validation.csv",index=False)
joint=pd.DataFrame(joint_rows);joint.to_csv(out/"interpolation_joint_cross_validation.csv",index=False)

# Within-standard split differences show finite-sample scale, not a p-value.
noise=[]
for e,actual in {**train.frames,**validation.frames}.items():
    for feature in features:
        a=actual[feature].iloc[::2].to_numpy();b=actual[feature].iloc[1::2].to_numpy()
        noise.append({"energy_gev":e,"feature":feature,"half_sample_events":len(a),
          "split_ks_distance":ks_distance(a,b,np.full(len(b),1/len(b)))})
pd.DataFrame(noise).to_csv(out/"finite_sample_split_reference.csv",index=False)

fig,a=plt.subplots(2,3,figsize=(15,9))
for col,(feature,label) in enumerate(features.items()):
    for row,kind in enumerate(("leave_one_out","independent_midpoint")):
        for method in methods:
            g=metrics[(metrics.feature==feature)&(metrics.validation==kind)&(metrics.method==method)].sort_values("energy_gev")
            a[row,col].plot(g.energy_gev,g.ks_distance,"o-",color=colors[method],label=method.replace("_"," "))
        a[row,col].axhline(.05,color=red,ls=":",lw=1,label="Provisional KS=0.05 guide")
        a[row,col].set(xscale="log",xlabel="Held-out energy (GeV)",ylabel="CDF maximum distance (KS)",title=label+" | "+("leave one out" if row==0 else "fresh intermediate samples"))
        if col==0:a[row,col].legend(fontsize=7)
save(fig,"04_interpolation_ks","Interpolation cross-validation | complete distributions, no invalid-charge filtering",
     "Leave-one-out removes each interior standard; fresh samples at 0.2, 0.6, 2, 6, 30 GeV use independent seeds.\nKS is a descriptive distance, not a significance test. The 0.05 line is a diagnostic guide, not an experimentally validated tolerance.")

fig,a=plt.subplots(2,3,figsize=(15,9))
for col,(feature,label) in enumerate(features.items()):
    for row,kind in enumerate(("leave_one_out","independent_midpoint")):
        for method in methods:
            g=metrics[(metrics.feature==feature)&(metrics.validation==kind)&(metrics.method==method)].sort_values("energy_gev")
            a[row,col].plot(g.energy_gev,100*g.median_relative_error,"o-",color=colors[method],label=method.replace("_"," "))
        a[row,col].axhline(0,color="gray",lw=.8)
        a[row,col].set(xscale="log",xlabel="Held-out energy (GeV)",ylabel="Median relative error (%)",title=label+" | "+("leave one out" if row==0 else "fresh intermediate samples"))
        if col==0:a[row,col].legend(fontsize=7)
save(fig,"05_interpolation_bias","Median prediction errors | separate from heavy-tail mean and charge",
     "All ionization events are retained. Mean, p90, p99, Wasserstein/IQR and medians are saved in interpolation_cross_validation.csv.\nMeasured spatial sigma is imposed in every energy library; its apparent interpolation agreement is not an independent physics validation.")

# Worst blind point for log mixture: plot all three distributions.
blind=metrics[(metrics.validation=="independent_midpoint")&(metrics.method=="log_mixture")]
worst=blind.sort_values("ks_distance",ascending=False).iloc[0]
e=float(worst.energy_gev);actual,library=examples_cv[("independent_midpoint",e)]
fig,a=plt.subplots(1,3,figsize=(15,5.5))
for axis,(feature,label) in zip(a,features.items()):
    if feature=="energy_transfer_ev":
        bins=np.geomspace(50,200000,75);axis.set_xscale("log")
    else:bins=np.arange(-.5,65.5) if feature=="n_clusters" else np.arange(-.5,301.5,3)
    axis.hist(actual[feature],bins=bins,weights=np.full(len(actual),1/len(actual)),histtype="step",color="black",lw=2,label=f"Direct Heed {e:g} GeV")
    for method in methods:
        predicted,weights=library.distribution(e,method)
        axis.hist(predicted[feature],bins=bins,weights=weights,histtype="step",color=colors[method],label=method.replace("_"," "))
    axis.set(xlabel=label,ylabel="Probability / bin",title=label)
    if feature=="n_clusters":axis.legend(fontsize=7)
save(fig,"06_interpolation_distributions",f"Direct simulation vs interpolated libraries | worst blind log-mixture point: {e:g} GeV",
     "Plot display limits do not remove events from validation metrics. Mixture interpolation preserves integer counts and joint event records.\nUncalibrated total avalanche charge is not validated; no energy extrapolation or incidence-angle extrapolation is provided.")

# Reproducible deployable example and checks of atomic event interpolation.
sample_path=out/"interpolated_2GeV_example.csv"
train.sample(2.,10000,20261021).to_csv(sample_path,index=False)
primaries=pd.read_csv(base/"gaisser_plane_100k_csv"/"primaries.csv")
outside=(primaries.kinetic_ev<train.energies.min()*1e9)|(primaries.kinetic_ev>train.energies.max()*1e9)
summary={"window_ns":50,"sigma_mm":sigma,"selected_spatial_factor":5,
    "spatial":spatial_rows,"interpolation_methods":list(methods),
    "training_energies_gev":train.energies.tolist(),"independent_validation_energies_gev":sorted(validation.frames),
    "interpolation_domain":"normal-incidence mu-; same gas, gap, field and imposed sigma; 0.1 to 100 GeV",
    "gaisser_primary_outside_energy_range":int(outside.sum()),
    "worst_independent_log_mixture_ks":float(worst.ks_distance),
    "worst_independent_log_mixture_energy_gev":e,"worst_independent_log_mixture_feature":worst.feature,
    "linear_charge_validated":False,"delay_pairing_validated":False}
by_method=metrics.groupby(["validation","method"]).agg(mean_ks=("ks_distance","mean"),
    max_ks=("ks_distance","max"),max_abs_median_error=("median_relative_error",lambda s:s.abs().max())).reset_index()
by_method.to_csv(out/"interpolation_method_summary.csv",index=False)
summary["cross_validation_summary"]=by_method.to_dict("records")
(out/"spatial_interpolation_summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False),encoding="utf-8")
pdf.close()

def md_table(frame):
    head="| "+" | ".join(frame.columns)+" |\n|"+"|".join(["---"]*len(frame.columns))+"|\n"
    return head+"\n".join("| "+" | ".join(str(v) for v in row)+" |" for row in frame.itertuples(index=False,name=None))

spatial_table=spatial[["factor_sigma","threshold_mm","positions","multi_position_windows","two_position_windows","three_or_more_position_windows"]].copy()
spatial_table.columns=["分离门限/σ","分离门限/mm","输出位置数","多位置窗口","双位置窗口","三位置及以上窗口"]
five_summary=spatial[spatial.factor_sigma==5].iloc[0]
method_table=by_method.copy()
for k in ("mean_ks","max_ks","max_abs_median_error"):method_table[k]=method_table[k].map(lambda x:f"{x:.5f}")
method_table["validation"]=method_table.validation.map({"independent_midpoint":"独立中间点","leave_one_out":"标准点留一"})
method_table["method"]=method_table.method.map({"nearest":"最近邻","linear_mixture":"E线性混合","log_mixture":"log(E)混合"})
method_table.columns=["验证方式","插值方法","平均KS距离","最大KS距离","最大中位数相对偏差"]
report=f"""# 同窗多位置读出与标准库插值交叉验证

本报告更新于2026-10-21，替代旧报告中“所有多轨迹窗口均无坐标”的唯一处理策略；旧结果保留为对照。
用户确认50 ns为事例分组/符合窗；空间分离扫描3σ、5σ、10σ，5σ作为本轮理想空间读出基线。

## 1. 空间多位置模型

σ=0.749 mm；三种假设分离门限为2.247、3.745、7.490 mm。门限是用户授权的模型选择，不是由单命中σ自动推导的实测双命中能力。
每个事件、每层先按首段时间建立固定50 ns组；每条模拟轨迹形成一个电离沉积加权中心，中心间距不超过门限的连通分量合成一个位置，分离分量各输出一个位置，分别添加实测σ约束的名义延迟端点噪声。
归并后若输出中心再次过近则继续合并，保证噪声前输出中心两两间距大于门限。靠近的次级电子与主粒子可形成一个响应，远处次级轨迹可以形成另一响应。

该模型采用模拟轨迹身份，并假设相隔足够远的位置能够被同时读出。它没有从真实波形寻找峰，也没有验证延迟线X/Y配对；几乎同时到达时的四端配对退化仍可能存在。空间可分与四端唯一配对分别列图。
一个轨迹内的延展沉积仍只取一个中心，不能把本模型理解为任意连续电离云的自动峰识别。没有脉冲幅度阈值时，极低沉积轨迹也可能产生理想响应，计数不能当作探测效率。

## 2. 同一批100000次Geant4平面源事件的重数字化

不重新抽取Geant4事件，使用原有1051625段真实气隙记录。三种空间门限均重处理全量397659个层时间窗；验证每窗段数、轨迹数、电离沉积总量及时间边界保持一致。

{md_table(spatial_table)}

5σ方案得到{int(five_summary.positions)}个坐标，{int(five_summary.multi_position_windows)}个窗口含多个位置；其中两个位置窗口{int(five_summary.two_position_windows)}个，三个及以上{int(five_summary.three_or_more_position_windows)}个。
原保守拒绝策略只有350620个有效坐标。新增坐标既包含近距离多轨迹合成的单位置，也包含远距离多轨迹的多个位置，不能直接把坐标数增量解释为真实效率提升。
三种门限的全部ROOT、逐位置CSV、逐窗口CSV均保留。输出以(event,layer,window_index,position_index)区分；n_positions为同组位置数，window_n_tracks为同组轨迹数，n_tracks为当前位置包含的轨迹数。
unresolved_component表示当前位置内合并了多条轨迹；valid_coordinate=true仅表示在本理想空间模型中坐标有效，readout_pairing_validated始终为false。

独立双缪子复用先前5000对输入，对0、10、25、49、50、51、75、100 ns各偏移比较。空间分开且同窗时现在允许两个坐标；间距不足则仍合成一个位置。51 ns以上作为分窗对照。没有由该受控叠加推算宇宙线偶然符合率。

## 3. 实际实现的插值算法

默认在相邻能量的log(E)上做完整事例分布的凸组合：w=ln(E/E低)/ln(E高/E低)，P(E)=(1-w)P低+wP高。
查询时先按w选高/低能量库，再抽取其中完整的一条事例，保留电离簇、整数电子数、转移能量之间的事例相关性和原始电荷无效标志。不是把两个随机事例的整数电子数直接取加权平均。
同时比较最近邻、E线性混合和log(E)混合三种方法。它是经验分布插值，并非由Ramo定理推导出的能量缩放；Ramo定理本身不能决定初级电离随缪子能量的变化。
代码response_interpolation.py可独立使用；输出interpolated_2GeV_example.csv含10000条插值查询结果及来源标准点/事例编号。
该示例的绝对坐标仍属于标准库内的均匀取点；用于指定Geant4位置时，应使用相对库内真值的响应偏移平移到实际入射位置，不能直接把抽到的库内绝对坐标当作Geant4测量坐标。本轮空间多位置重数字化使用实际Geant4沉积中心，不使用这份示例的绝对位置。

## 4. 交叉验证设计与结果

训练点为0.1、0.3、1、3、10、100 GeV，各10000个垂直mu-。
第一组留一验证：逐一完全移除0.3、1、3、10 GeV，再用剩余库相邻点预测被移除点。端点不参加插值留一验证，因为移除端点需要外推，而当前程序拒绝外推。
第二组独立验证：新生成0.2、0.6、2、6、30 GeV各10000个Heed事例，使用各不相同且区别于训练库的随机种子。这五点不加入训练库。
分别检验电离簇数、气隙电子数、转移能量的完整分布；保存KS距离、Wasserstein距离/IQR、均值、中位数、p90、p99及相对偏差。
额外验证“簇数—电子数”和“电子数—转移能量”二维累计分布，在留出样本分位网格上计算最大差，保存为interpolation_joint_cross_validation.csv；它是描述性距离，不是二维显著性检验。
KS为累计分布最大绝对差，不报告独立性不满足时的p值。训练库具有共同随机数，留一结果可能带相关性；新种子中间点用于补充检验。另附每库奇偶半样本的KS差异作为有限样本量参照，半样本参照不等于全样本检验阈值。

{md_table(method_table)}

独立中间点中log(E)混合最差KS={float(worst.ks_distance):.5f}，对应{e:g} GeV的{worst.feature}。图中的KS=0.05仅作临时诊断参考；本轮没有事先规定实验接受门限，不能据此宣称所有能量已验证通过。
留一验证最差点在0.3 GeV附近的低电离区：移除该点后，相邻0.1与1 GeV库的混合无法恢复电离谷，应保留0.3 GeV标准点，并优先在0.1—0.3 GeV加密。不能将log(E)本轮平均较好解释为每个点都优于其他方法。
CSV同时保留长尾均值与p99，不用核心中位数掩盖大能量转移尾部。应按实际精度目标选择更密的标准点，尤其检查留一验证暴露的非线性区间。

## 5. 适用边界

插值仅验证相同气体、298.15 K、760 Torr、2 mm、10.1 kV低计数率稳态、垂直mu-的能量维度，不能把垂直参考直接用于任意倾角、次级电子或所有Geant4步的局部动能。
Gaisser平面源100000次初级入射中有{int(outside.sum())}次超出0.1–100 GeV标准范围，查询必须明确返回越界，不能静默用端点替代；入射角仍需另建验证。
所有能量点的0.749 mm展宽是外部输入，所以坐标宽度一致不是插值物理正确性的独立证据。
线性雪崩完整总电荷仍大量无效，本轮没有用少量有效电荷来验证完整电荷插值，也没有据此预测触发效率、饱和电荷或多脉冲幅度。

Heed簇和δ电子的定义参考Garfield++官方手册：https://garfieldpp.web.cern.ch/doxygen/UserGuide.pdf 。本轮凸组合算法与理想空间规则为本项目明确声明的近似，未声称它们是手册规定的RPC多命中模型。

## 6. 文件与复现

run_spatial_interpolation.sh：在Ubuntu构建、运行内核验证、三门限重数字化、五个独立能量验证库。
analyze_spatial_interpolation.py：生成本报告、全量一致性检查、空间与插值图表。
gaisser_plane_spatial_5sigma.root及同名CSV：主要多位置基线。
interpolation_cross_validation.csv：全部81组“验证点×方法×物理量”指标。
interpolation_method_summary.csv：方法汇总；finite_sample_split_reference.csv：样本统计参照。
validation_energies/：五套独立Heed验证库及逐电子/簇/事例导出。
figures/spatial_interpolation_figures.pdf：全部六组科学图表；每组另存PNG、PDF。
"""
(out/"SPATIAL_INTERPOLATION_REPORT.md").write_text(report,encoding="utf-8")
parts=[]
for paragraph in report.split("\n\n"):
    if paragraph.startswith("#"):
        level=len(paragraph.split(" ",1)[0]);parts.append(f"<h{level}>"+html.escape(paragraph[level:].strip())+f"</h{level}>")
    elif paragraph.startswith("|"):
        rows=[s.strip().strip("|").split("|") for s in paragraph.splitlines()]
        parts.append("<table>"+"".join("<tr>"+"".join("<td>"+html.escape(c.strip())+"</td>" for c in r)+"</tr>" for r in [rows[0]]+rows[2:])+"</table>")
    else:parts.append("<p>"+html.escape(paragraph).replace("\n","<br>")+"</p>")
for path in sorted(figures.glob("*.png")):
    parts.append("<h2>"+html.escape(path.stem)+"</h2><img src='data:image/png;base64,"+base64.b64encode(path.read_bytes()).decode()+"'>")
style="body{font:16px/1.75 'Microsoft YaHei',sans-serif;max-width:1200px;margin:40px auto;padding:0 25px;color:#0f172a;background:#f8fafc}h1,h2{color:#1e3a8a}img{width:100%}table{border-collapse:collapse;font-size:14px}td{border:1px solid #cbd5e1;padding:8px}tr:first-child{background:#dbeafe}"
(out/"spatial_interpolation_report.html").write_text("<!doctype html><html lang='zh-CN'><meta charset='utf-8'><title>RPC多位置与插值验证</title><style>"+style+"</style><body>"+"\n".join(parts)+"</body></html>",encoding="utf-8")
print(spatial.to_string(index=False))
print(by_method.to_string(index=False))
print("Report:",out/"spatial_interpolation_report.html")
