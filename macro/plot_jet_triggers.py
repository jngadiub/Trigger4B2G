import ROOT
from ROOT import *
import time
ROOT.gROOT.SetBatch(True)
import CMS_lumi, tdrstyle
tdrstyle.setTDRStyle()
#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_13TeV = "SingleMuon 2018A, 10.3 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPeriod = 4
iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
 
def get_text():

 pt = ROOT.TPaveText(0.445,0.73,0.80,0.89,"NDC")
 pt.SetTextFont(42)
 pt.SetTextSize(0.04)
 pt.SetTextAlign(12)
 pt.SetFillColor(0)
 pt.SetBorderSize(0)
 pt.SetFillStyle(0)
 
 return pt

def get_canvas(cname):

 H_ref = 600; 
 W_ref = 800; 
 W = W_ref
 H  = H_ref

 # references for T, B, L, R
 T = 0.08*H_ref
 B = 0.12*H_ref 
 L = 0.12*W_ref
 R = 0.04*W_ref

 canvas = ROOT.TCanvas(cname,cname,50,50,W,H)
 canvas.SetFillColor(0)
 canvas.SetBorderMode(0)
 canvas.SetFrameFillStyle(0)
 canvas.SetFrameBorderMode(0)
 canvas.SetLeftMargin( L/W )
 canvas.SetRightMargin( R/W )
 canvas.SetTopMargin( T/H )
 canvas.SetBottomMargin( B/H )
 canvas.SetTickx(0)
 canvas.SetTicky(0)

 return canvas

colors = {}
triggers = {}

colors["HLT_PFHT1050"] = ROOT.kAzure+7
colors["HLT_AK8PFJet500"] = ROOT.kOrange+1
colors["HLT_AK8PFJet400_TrimMass30"] = 210
colors["HLT_AK8PFHT800_TrimMass50"] = ROOT.kPink-1
colors["All triggers"] = ROOT.kBlack

f = ROOT.TFile.Open('jet_triggers.root','R')
mjj_num = {}
mjj_den = ROOT.TH1F()
mjet1_num = {}
mjet1_den = ROOT.TH1F()
mjet2_num = {}
mjet2_den = ROOT.TH1F()

for k in f.GetListOfKeys():
 if 'den' in str(k.GetName()) and 'mjj' in str(k.GetName()): mjj_den = f.Get(k.GetName())
 if 'den' in str(k.GetName()) and 'mjet1' in str(k.GetName()): mjet1_den = f.Get(k.GetName()) 
 if 'den' in str(k.GetName()) and 'mjet2' in str(k.GetName()): mjet2_den = f.Get(k.GetName())
 if 'num' in str(k.GetName()) and 'mjj' in str(k.GetName()): mjj_num[k.GetName().replace('mjj_num_','').replace('all_triggers','All triggers')] = f.Get(k.GetName())
 if 'num' in str(k.GetName()) and 'mjet1' in str(k.GetName()): mjet1_num[k.GetName().replace('mjet1_num_','').replace('all_triggers','All triggers')] = f.Get(k.GetName())
 if 'num' in str(k.GetName()) and 'mjet2' in str(k.GetName()): mjet2_num[k.GetName().replace('mjet2_num_','').replace('all_triggers','All triggers')] = f.Get(k.GetName())

outf = ROOT.TFile.Open('effs_triggers.root','RECREATE')
effs_mjj = {}
effs_mjet1 = {}
effs_mjet2 = {}
for k in colors.keys():

 c = get_canvas("c_mjj_%s"%k) 
 c.cd()  
 
 pEff = ROOT.TEfficiency(mjj_num[k],mjj_den)
 pEff.SetLineColor(colors[k])
 pEff.SetMarkerColor(colors[k])
 pEff.SetTitle("mjj_%s;m_{jj} [GeV]; Efficiency"%k)
 pEff.SetName("mjj_eff_%s"%k)
 pEff.Draw()
 ROOT.gPad.Update()
 pEff.GetPaintedGraph().SetMinimum(0.)
 pEff.GetPaintedGraph().SetMaximum(1.4)
 effs_mjj[k] = pEff
 
 outf.cd()
 pEff.Write()
 
 pt = get_text()
 text = pt.AddText(k)
 text.SetTextFont(62)
 pt.AddText("AK8 jets, p_{T} > 200 GeV, |#Delta#eta_{jj}| < 1.3")
 pt.AddText("m_{jet}^{SD} > 55 GeV")
 pt.Draw()
   
 CMS_lumi.CMS_lumi(c, iPeriod, iPos) 
 c.cd()
 c.Update()
 c.RedrawAxis()
 frame = c.GetFrame()
 frame.Draw()
 c.SaveAs(c.GetName()+".png")
 
 for j in range(1,3):
  c = get_canvas("c_mjet%i_%s"%(j,k)) 
  c.cd()  
 
  if j==1: pEff = ROOT.TEfficiency(mjet1_num[k],mjet1_den)
  else: pEff = ROOT.TEfficiency(mjet2_num[k],mjet2_den)
  pEff.SetLineColor(colors[k])
  pEff.SetMarkerColor(colors[k])
  pEff.SetTitle("mjet%i_%s;m_{jet%i} [GeV]; Efficiency"%(j,k,j))
  pEff.SetName("mjet%i_eff_%s"%(j,k))
  pEff.Draw()
  ROOT.gPad.Update()
  pEff.GetPaintedGraph().SetMinimum(0.)
  pEff.GetPaintedGraph().SetMaximum(1.4)
  if j == 1: effs_mjet1[k] = pEff
  else: effs_mjet2[k] = pEff

  outf.cd()
  pEff.Write()
  
  pt = get_text()
  text = pt.AddText(k)
  text.SetTextFont(62)
  pt.AddText("AK8 jets, p_{T} > 200 GeV, |#Delta#eta_{jj}| < 1.3")
  pt.AddText("m_{jj} > 1.2 TeV")
  pt.Draw()
   
  CMS_lumi.CMS_lumi(c, iPeriod, iPos) 
  c.cd()
  c.Update()
  c.RedrawAxis()
  frame = c.GetFrame()
  frame.Draw()
  c.SaveAs(c.GetName()+".png") 

f.Close()
outf.Write()
outf.Close()

print effs_mjj.keys()
print effs_mjet1.keys()
print effs_mjet1.keys()
print colors.keys()

c = get_canvas("c_mjj_all") 
c.cd()

leg = ROOT.TLegend(0.55,0.20,0.75,0.40)
leg.SetBorderSize(0)
leg.SetTextSize(0.031)
leg.SetLineColor(1)
leg.SetLineStyle(1)
leg.SetShadowColor(0)
leg.SetLineWidth(1)
leg.SetFillColor(0)
leg.SetTextFont(42)

effs_mjj[effs_mjj.keys()[0]].Draw()
for k in effs_mjj.keys():
 effs_mjj[k].Draw("SAME")
 leg.AddEntry(effs_mjj[k],k,"LP")
leg.Draw()

pt = get_text()
#text = pt.AddText("All triggers")
#text.SetTextFont(62)
pt.AddText("AK8 jets, p_{T} > 200 GeV, |#Delta#eta_{jj}| < 1.3")
pt.AddText("m_{jet}^{SD} > 55 GeV")
pt.Draw()

CMS_lumi.CMS_lumi(c, iPeriod, iPos) 
c.cd()
c.Update()
c.RedrawAxis()
frame = c.GetFrame()
frame.Draw()
c.SaveAs(c.GetName()+".png")  

for j in range(1,3):
 c = get_canvas("c_mjet%i_all"%j) 
 c.cd()

 leg = ROOT.TLegend(0.55,0.20,0.75,0.40)
 leg.SetBorderSize(0)
 leg.SetTextSize(0.031)
 leg.SetLineColor(1)
 leg.SetLineStyle(1)
 leg.SetShadowColor(0)
 leg.SetLineWidth(1)
 leg.SetFillColor(0)
 leg.SetTextFont(42)

 if j==1:
  effs_mjet1[effs_mjet1.keys()[0]].Draw()
  for k in effs_mjet1.keys():
   effs_mjet1[k].Draw("SAME")
   leg.AddEntry(effs_mjet1[k],k,"LP")
 else:
  effs_mjet2[effs_mjet2.keys()[0]].Draw()
  for k in effs_mjet2.keys():
   effs_mjet2[k].Draw("SAME")
   leg.AddEntry(effs_mjet2[k],k,"LP")

 leg.Draw()

 pt = get_text()
 #text = pt.AddText("All triggers")
 #text.SetTextFont(62)
 pt.AddText("AK8 jets, p_{T} > 200 GeV, |#Delta#eta_{jj}| < 1.3")
 pt.AddText("m_{jj} > 1.2 TeV")
 pt.Draw()

 CMS_lumi.CMS_lumi(c, iPeriod, iPos) 
 c.cd()
 c.Update()
 c.RedrawAxis()
 frame = c.GetFrame()
 frame.Draw()
 c.SaveAs(c.GetName()+".png")