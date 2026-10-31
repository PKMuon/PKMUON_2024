"""Redraw existing result figures with the user's root_macro/TDR conventions.

All numerical work is delegated unchanged to the original analysis scripts in an
isolated staging directory. No detector rerun and no edits to root_macro inputs.
Matplotlib is the renderer; this is a visual adaptation, not a claim of ROOT use.
"""
import gc
import hashlib
import json
import os
import runpy
import sys
import textwrap
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.colors import to_hex, to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Polygon
from matplotlib.ticker import AutoMinorLocator

PROJECT=Path(__file__).resolve().parents[2]
SOURCE=PROJECT/"rpc_response/results"
DEST=PROJECT/"rpc_response/results/presentation_tdr_20261021"
DEST.mkdir(exist_ok=True)
scripts=PROJECT/"rpc_response/scripts"
reference=Path("D:/muon/root_macro")
COLORS=["#000000","#d00000","#008000","#0000cc","#a000a0","#e06000"]
known={"#2563eb":COLORS[3],"#0d9488":COLORS[1],"#d97706":COLORS[2],
       "#dc2626":COLORS[4],"#7c3aed":COLORS[5]}
panels=[];grouped=[]
original_save=Figure.savefig

def link_inputs(src,dst):
    dst.mkdir(parents=True,exist_ok=True)
    for f in src.iterdir():
        if f.is_file() and f.suffix in (".csv",".json"):
            target=dst/f.name
            if not target.exists():os.link(f,target)

initial=SOURCE/"muon_1GeV_normal_v1"
for f in SOURCE.glob("muon_1GeV_normal_v1_*"):
    if f.is_file() and f.suffix in (".csv",".json"):
        target=DEST/f.name
        if not target.exists():os.link(f,target)
base=DEST/"comparison_20261020";base.mkdir(exist_ok=True)
source_base=SOURCE/"comparison_20261020"
for folder in ("energies","original_y_with_legacy","gaisser_plane_100k_csv"):
    link_inputs(source_base/folder,base/folder)
sp=base/"spatial_interpolation";sp.mkdir(exist_ok=True)
for f in (source_base/"spatial_interpolation").glob("gaisser_plane_spatial_*sigma.csv"):
    target=sp/f.name
    if not target.exists():os.link(f,target)
link_inputs(source_base/"spatial_interpolation/validation_energies",sp/"validation_energies")

def root_style(fig):
    if getattr(fig,"_rpc_tdr_styled",False):return
    fig._rpc_tdr_styled=True
    fig.set_facecolor("white")
    mapping=dict(known);mapping.update({c:c for c in COLORS});extra=0
    def color(value):
        nonlocal extra
        try:
            rgba=to_rgba(value)
            if rgba[3]==0:return value
            key=to_hex(rgba)
            if max(rgba[:3])-min(rgba[:3])<.025:return key
            if key not in mapping:
                mapping[key]=COLORS[extra%len(COLORS)];extra+=1
            return mapping[key]
        except (ValueError,TypeError):return value
    if fig._suptitle is not None:fig._suptitle.set_visible(False)
    for t in fig.texts:
        t.set_fontfamily("Arial");t.set_color("black")
    for ax in fig.axes:
        ax.set_facecolor("white");ax.grid(False,which="both")
        if not ax.axison:continue
        for spine in ax.spines.values():spine.set_visible(True);spine.set_color("black");spine.set_linewidth(1)
        ax.tick_params(axis="both",which="major",direction="in",top=True,right=True,length=6,width=.9,labelsize=11,pad=5)
        ax.tick_params(axis="both",which="minor",direction="in",top=True,right=True,length=3,width=.7)
        for axis in (ax.xaxis,ax.yaxis):
            axis.label.set_fontsize(13);axis.label.set_fontfamily("Arial");axis.label.set_color("black")
            for label in axis.get_ticklabels():label.set_fontfamily("Arial");label.set_color("black")
        if ax.get_xscale()=="linear" and "Category" not in type(ax.xaxis.get_major_locator()).__name__:
            ax.xaxis.set_minor_locator(AutoMinorLocator(5))
        if ax.get_yscale()=="linear":ax.yaxis.set_minor_locator(AutoMinorLocator(5))
        title=ax.get_title()
        ax.set_title("\n".join(textwrap.wrap(title,46,break_long_words=False)),fontsize=11.5,fontfamily="Arial",fontweight="normal",color="black",pad=11)
        ax.xaxis.label.set_horizontalalignment("right");ax.xaxis.label.set_x(1.)
        for line in ax.lines:
            line.set_color(color(line.get_color()));line.set_linewidth(1.6)
            if line.get_marker() not in (None,"None","none",""):line.set_markersize(4.8)
        # TH1::Draw("HIST") has a continuous step outline, not a collection
        # of outlined bars with a vertical line from zero at every bin edge.
        for container in list(ax.containers):
            patches=getattr(container,"patches",[])
            if len(patches)<12 or not all(isinstance(p,Rectangle) for p in patches):continue
            if not all(abs(p.get_y())<1e-12 for p in patches):continue
            ordered=sorted(patches,key=lambda p:p.get_x())
            if not all(abs(a.get_x()+a.get_width()-b.get_x())<1e-7 for a,b in zip(ordered,ordered[1:])):continue
            values=[p.get_height() for p in ordered]
            edges=[p.get_x() for p in ordered]+[ordered[-1].get_x()+ordered[-1].get_width()]
            c=color(ordered[0].get_facecolor())
            for p in patches:p.set_visible(False)
            ax.stairs(values,edges,baseline=0,fill=False,color=c,linewidth=1.6)
        for p in ax.patches:
            edge=p.get_edgecolor();face=p.get_facecolor()
            if isinstance(p,(Rectangle,Polygon)):
                use=face if face[3]>0 and sum(face[:3])<2.98 else edge
                p.set_edgecolor(color(use));p.set_facecolor("none");p.set_linewidth(1.25);p.set_alpha(1)
        for collection in ax.collections:
            if collection.get_array() is not None:continue  # keep scalar density maps and their colorbar
            fc=collection.get_facecolors();ec=collection.get_edgecolors()
            if len(fc)==1:collection.set_facecolor(color(fc[0]))
            if len(ec)==1:collection.set_edgecolor(color(ec[0]))
        leg=ax.get_legend()
        if leg:
            leg.set_frame_on(False)
            for text in leg.get_texts():text.set_fontfamily("Arial");text.set_fontsize(9.5);text.set_color("black")
            for line in leg.get_lines():
                line.set_color(color(line.get_color()));line.set_linewidth(1.6)
                line.set_markerfacecolor(color(line.get_markerfacecolor()))
                line.set_markeredgecolor(color(line.get_markeredgecolor()))
            for handle in leg.legend_handles:
                if hasattr(handle,"get_facecolors"):
                    fc=handle.get_facecolors();ec=handle.get_edgecolors()
                    if len(fc)==1:handle.set_facecolor(color(fc[0]))
                    if len(ec)==1:handle.set_edgecolor(color(ec[0]))
            for p in leg.get_patches():
                face=p.get_facecolor();edge=p.get_edgecolor();use=face if face[3]>0 else edge
                p.set_edgecolor(color(use));p.set_facecolor("none");p.set_linewidth(1.3)
        for text in ax.texts:text.set_fontfamily("Arial");text.set_color("black")
    fig.tight_layout(rect=(.025,.085,.99,.97),h_pad=3.0,w_pad=3.0)

def savefig(fig,fname,*args,**kwargs):
    root_style(fig)
    kwargs["facecolor"]="white"
    result=original_save(fig,fname,*args,**kwargs)
    if isinstance(fname,(str,Path)) and str(fname).endswith(".png"):
        filename=Path(fname)
        grouped.append(str(filename.relative_to(DEST)))
        fig.canvas.draw();renderer=fig.canvas.get_renderer()
        axes=[a for a in fig.axes if a.axison and a.get_position().width>.12]
        for i,ax in enumerate(axes,1):
            # Retain the axis, labels, title and legend. A colorbar belongs to its
            # density panel; include it when Matplotlib records one on a mappable.
            boxes=[ax.get_tightbbox(renderer)]
            for collection in ax.collections:
                cb=getattr(collection,"colorbar",None)
                if cb is not None:boxes.append(cb.ax.get_tightbbox(renderer))
            from matplotlib.transforms import Bbox
            bounds=Bbox.union(boxes).expanded(1.045,1.065).transformed(fig.dpi_scale_trans.inverted())
            target=filename.with_name(filename.stem+f"__panel_{i:02d}.png")
            original_save(fig,target,dpi=260,bbox_inches=bounds,facecolor="white",pad_inches=.04)
            panels.append({"group":str(filename.relative_to(DEST)),"panel":i,"title":ax.get_title(),
                "xlabel":ax.get_xlabel(),"ylabel":ax.get_ylabel(),"path":str(target.relative_to(DEST))})
    return result

Figure.savefig=savefig
plt.rcParams["font.family"]="Arial"
sys.path.insert(0,str(scripts))
for filename,arg in [("plot_muon_standard.py",DEST/"muon_1GeV_normal_v1"),
                     ("analyze_comparison.py",base),
                     ("analyze_spatial_interpolation.py",base)]:
    sys.argv=[str(scripts/filename),str(arg)]
    print("REDRAW",filename,flush=True)
    runpy.run_path(str(scripts/filename),run_name="__main__")
    gc.collect();plt.close("all")

# The appearance must not change numerical conclusions.
comparisons=[(SOURCE/"muon_1GeV_normal_v1_plots/summary.json",DEST/"muon_1GeV_normal_v1_plots/summary.json"),
 (source_base/"comparison_summary.json",base/"comparison_summary.json"),
 (source_base/"spatial_interpolation/spatial_interpolation_summary.json",sp/"spatial_interpolation_summary.json")]
for a,b in comparisons:
    assert json.loads(a.read_text())==json.loads(b.read_text()),f"Numerical summary changed: {a.name}"
style_files=["tdrstyle.C","plotOffsetStandardComparison.cc","plotMaterialsThicknessRatio.cc","plotEnergyDistributions.cc"]
manifest={"renderer":"Matplotlib; visual conventions adapted from user ROOT/TDR macros",
  "references":[{"path":str(reference/f),"sha256":hashlib.sha256((reference/f).read_bytes()).hexdigest()} for f in style_files],
  "style":{"canvas":"white","font":"Arial / ROOT font42 analogue","grid":False,"ticks":"inward on all four sides",
           "histograms":"unfilled outlines","colors":COLORS,"legend":"no frame"},
  "numeric_summary_checks":"all three summaries exactly unchanged","grouped_figures":sorted(set(grouped)),"panels":panels}
(DEST/"root_style_manifest.json").write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding="utf-8")
print("REDRAW COMPLETE groups=",len(set(grouped)),"panels=",len(panels),flush=True)
