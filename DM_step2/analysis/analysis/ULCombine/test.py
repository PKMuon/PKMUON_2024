import pyhf
import numpy as np
import matplotlib.pyplot as plt
from pyhf.contrib.viz import brazil

import ROOT
from ROOT import TGraph,TCanvas,TLegend

def plotUpperLimits(values,limits,XaxisTitle,YaxisTitle,theoryprediction = None,logy=False,logx=False):
    # see CMS plot guidelines: https://ghm.web.cern.ch/ghm/plots/
 
    N = len(values)
    yellow = TGraph(2*N)    # yellow band
    green = TGraph(2*N)     # green band
    median = TGraph(N)      # median line
    if theoryprediction:
        theory = TGraph(N)      # theory line
 
    up2s = [ ]
    dn2s = [ ]
    for i in range(N):
        limit = limits[i]
        up2s.append(limit[4])
        dn2s.append(limit[0])
        yellow.SetPoint(    i,    values[i], limit[4] ) # + 2 sigma
        green.SetPoint(     i,    values[i], limit[3] ) # + 1 sigma
        median.SetPoint(    i,    values[i], limit[2] ) # median
        if theoryprediction:
            theory.SetPoint(    i,    values[i], theoryprediction[i] ) # median
        green.SetPoint(  2*N-1-i, values[i], limit[1] ) # - 1 sigma
        yellow.SetPoint( 2*N-1-i, values[i], limit[0] ) # - 2 sigma
 
    W = 800
    H  = 800
    T = 0.08*H
    B = 0.12*H
    L = 0.12*W
    R = 0.04*W
    c = TCanvas("c","c",100,100,W,H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin( L/W )
    c.SetRightMargin( R/W )
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )
    c.SetTickx(0)
    c.SetTicky(0)
    c.SetGrid()
    if logy:
        c.SetLogy()
    if logx:
        c.SetLogx()
    c.cd()
    
    frame = c.DrawFrame(min(values),min(dn2s), max(values), max(up2s))
    frame.GetYaxis().CenterTitle()
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.03)
    frame.GetYaxis().SetLabelSize(0.03)
    frame.GetYaxis().SetTitleOffset(1.1)
    frame.GetXaxis().SetNdivisions(508)
    frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetTitle(YaxisTitle)
#    frame.GetYaxis().SetTitle("95% upper limit on #sigma #times BR / (#sigma #times BR)_{SM}")
    frame.GetXaxis().SetTitle(XaxisTitle)
    if logy:
        frame.SetMinimum(min(dn2s))
        frame.SetMaximum(max(up2s))
    else:
        frame.SetMinimum(0)
        frame.SetMaximum(max(up2s))
    frame.GetXaxis().SetLimits(min(values),max(values))
 
    yellow.SetFillColor(ROOT.kOrange)
    yellow.SetLineColor(ROOT.kOrange)
    yellow.SetFillStyle(1001)
    yellow.Draw('F')
 
    green.SetFillColor(ROOT.kGreen+1)
    green.SetLineColor(ROOT.kGreen+1)
    green.SetFillStyle(1001)
    green.Draw('Fsame')
 
    median.SetLineColor(1)
    median.SetLineWidth(2)
    median.SetLineStyle(2)
    median.Draw('Lsame')

    if theoryprediction:
        theory.SetLineColor(2)
        theory.SetLineWidth(2)
        theory.SetLineStyle(1)
        theory.Draw('Lsame')
 
#     CMS_lumi.CMS_lumi(c,14,11)
    ROOT.gPad.SetTicks(1,1)
    frame.Draw('sameaxis')
 
    x1 = 0.1
    x2 = 0.3
    y2 = 0.9
    y1 = 0.7
    legend = TLegend(0.1,0.7,0.3,0.9,"brNDC")
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)
    legend.SetTextFont(42)
    legend.AddEntry(median, "Asymptotic CL_{s} expected",'L')
    legend.AddEntry(green, "#pm 1 std. deviation",'f')
#    legend.AddEntry(green, "Asymptotic CL_{s} #pm 1 std. deviation",'f')
    legend.AddEntry(yellow,"#pm 2 std. deviation",'f')
    if theoryprediction:
        legend.AddEntry(theory, "theory prediciotn",'L')
#    legend.AddEntry(green, "Asymptotic CL_{s} #pm 2 std. deviation",'f')
    legend.Draw()

    c.Update()
    #c.Print("plot11.eps")
     
    Graphs = [yellow,median,green,]
    if theoryprediction:
        Graphs = [yellow,median,green,theory,legend]

    return c,Graphs

sigmalimits_const_air = [[6.274541219075522e-09,
  8.337606986363731e-09,
  1.1555989583333336e-08,
  1.602406899134318e-08,
  2.1325774987538656e-08],
 [2.7926127115885414e-09,
  3.710822264353434e-09,
  5.143229166666667e-09,
  7.13183879852295e-09,
  9.491471449534098e-09],
 [3.2708740234374996e-09,
  4.3633874257405605e-09,
  6.067708333333334e-09,
  8.413764635721842e-09,
  1.119753360748291e-08],
 [3.9477539062500015e-09,
  5.570621490478516e-09,
  8.020833333333336e-09,
  1.150571425755819e-08,
  1.553211053212484e-08],
 [6.375757853190105e-09,
  9.20810600121816e-09,
  1.3216145833333335e-08,
  1.9379721085230512e-08,
  2.645141084988912e-08],
 [9.4024658203125e-09,
  1.359138488769531e-08,
  2.0572916666666666e-08,
  3.049546480178833e-08,
  4.2143630981445315e-08],
 [5.526224772135417e-08,
  8.428225914637248e-08,
  1.3346354166666666e-07,
  2.0368639628092448e-07,
  2.89760947227478e-07],
 [3.9947509765625003e-07,
  7.445573806762696e-07,
  1.2174479166666665e-06,
  1.8337533871332806e-06,
  2.738724748293559e-06]]

MDM=[0.005, 0.05, 0.1, 0.2, 0.5, 1, 10, 100]

c,_ = plotUpperLimits(MDM,sigmalimits_const_air,"DM mass (GeV)","95% upper limit on #sigma_{#mu,DM}",logy=True,logx=True)
c.Draw()

'''
values=MDM
limits=sigmalimits_const_air
logy=True
logx=True
XaxisTitle="DM mass (GeV)"
YaxisTitle="95% upper limit on #sigma_{#mu,DM}"

N = len(values)
yellow = TGraph(2*N)    # yellow band
green = TGraph(2*N)     # green band
median = TGraph(N)      # median line

up2s = [ ]
dn2s = [ ]
for i in range(N):
    limit = limits[i]
    up2s.append(limit[4])
    dn2s.append(limit[0])
    yellow.SetPoint(    i,    values[i], limit[4] ) # + 2 sigma
    green.SetPoint(     i,    values[i], limit[3] ) # + 1 sigma
    median.SetPoint(    i,    values[i], limit[2] ) # median
    green.SetPoint(  2*N-1-i, values[i], limit[1] ) # - 1 sigma
    yellow.SetPoint( 2*N-1-i, values[i], limit[0] ) # - 2 sigma

W = 800
H  = 800
T = 0.08*H
B = 0.12*H
L = 0.12*W
R = 0.04*W
c = TCanvas("c","c",100,100,W,H)
c.SetFillColor(0)
c.SetBorderMode(0)
c.SetFrameFillStyle(0)
c.SetFrameBorderMode(0)
c.SetLeftMargin( L/W )
c.SetRightMargin( R/W )
c.SetTopMargin( T/H )
c.SetBottomMargin( B/H )
c.SetTickx(0)
c.SetTicky(0)
c.SetGrid()
if logy:
    c.SetLogy()
if logx:
    c.SetLogx()
c.cd()

frame = c.DrawFrame(min(values)*0.8,min(dn2s)*0.1, max(values)*2, max(up2s)*10)
frame.GetYaxis().CenterTitle()
frame.GetYaxis().SetTitleSize(0.05)
frame.GetXaxis().SetTitleSize(0.05)
frame.GetXaxis().SetLabelSize(0.03)
frame.GetYaxis().SetLabelSize(0.03)
frame.GetYaxis().SetTitleOffset(1.1)
frame.GetXaxis().SetNdivisions(508)
frame.GetYaxis().CenterTitle(True)
frame.GetYaxis().SetTitle(YaxisTitle)
frame.GetXaxis().SetTitle(XaxisTitle)
if logy:
    frame.SetMinimum(min(dn2s)*0.1)
    frame.SetMaximum(max(up2s)*10)
else:
    frame.SetMinimum(0)
    frame.SetMaximum(max(up2s)*1.2)
frame.GetXaxis().SetLimits(min(values),max(values))

yellow.SetFillColor(ROOT.kOrange)
yellow.SetLineColor(ROOT.kOrange)
yellow.SetFillStyle(1001)
yellow.Draw('F')

green.SetFillColor(ROOT.kGreen+1)
green.SetLineColor(ROOT.kGreen+1)
green.SetFillStyle(1001)
green.Draw('Fsame')

median.SetLineColor(1)
median.SetLineWidth(2)
median.SetLineStyle(2)
median.Draw('Lsame')

ROOT.gPad.SetTicks(1,1)
frame.Draw('sameaxis')

x1 = 0.1
x2 = 0.3
y2 = 0.9
y1 = 0.7
legend = TLegend(0.2,0.7,0.8,0.9)
legend.SetFillStyle(0)
legend.SetBorderSize(0)
legend.SetTextSize(0.03)
legend.SetTextFont(42)
legend.AddEntry(median, "Asymptotic CL_{s} expected",'L')
legend.AddEntry(green, "#pm 1 std. deviation",'f')
legend.AddEntry(yellow,"#pm 2 std. deviation",'f')
legend.Draw()
 
# 将图形保存到文件
c.Print("histogram.png")
'''

# 等待用户输入
input("Press Enter to continue...")

print("Done")
