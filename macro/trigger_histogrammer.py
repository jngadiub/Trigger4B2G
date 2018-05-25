#! /usr/bin/env python

import os, multiprocessing
import copy
import math
import numpy as np
from array import array
from ROOT import ROOT, gROOT, gStyle, gRandom, TSystemDirectory
from ROOT import TFile, TChain, TTree, TCut, TH1, TH1F, TH2F, THStack, TGraph, TGraphAsymmErrors, TF1
from ROOT import TStyle, TCanvas, TPad
from ROOT import TLegend, TLatex, TText, TLine, TBox


#from Analysis.ALPHA.drawUtils import *
from variables_trigger_longlabel import *
#from Analysis.ALPHA.selections import *
#from Analysis.ALPHA.samples import sample, samples

########## SETTINGS ##########

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-v", "--variable", action="store", type="string", dest="variable", default="met_pt_nomu_L")
parser.add_option("-c", "--cut", action="store", type="string", dest="cut", default="")
parser.add_option("-l", "--lepton", action="store", type="string", dest="lepton", default="mu3")
parser.add_option("-r", "--run", action="store", type="string", dest="run", default="")
parser.add_option("-g", "--goodonlinejec", action="store_true", default=False, dest="goodonlinejec")
parser.add_option("-a", "--all", action="store_true", default=False, dest="all")
parser.add_option("-b", "--bash", action="store_true", default=False, dest="bash")
parser.add_option("-B", "--blind", action="store_true", default=False, dest="blind")
parser.add_option("-f", "--final", action="store_true", default=False, dest="final")
parser.add_option("-R", "--rebin", action="store_true", default=False, dest="rebin")
parser.add_option("-p", "--public", action="store_true", default=False, dest="public")
(options, args) = parser.parse_args()
if options.bash: gROOT.SetBatch(True)

########## SETTINGS ##########

gStyle.SetOptStat(0)

NTUPLEDIR   = "$CMSSW_BASE/src/TrigAnalyzer/Trigger4B2G/v1/"
#LUMI        = 10090#4889.443#3785.218#Golden JSON # in pb-1
#LUMI        = 3819.072 + 676.870 + 1264.178 + 3992.295 + 2246.983#at 5 sept, prod v8
#LUMI        = 3819.072 + 676.870 + 1264.178 + 3992.295 + 3978.673#at 13 sept, prod v9
#LUMI        = 15337#at 15 sept, prod v10
#LUMI        = 9869.850 #ONLY Goodonline jec run b-c-d
LUMI        = 7.93
SIGNAL      = 1.
RATIO       = 4 # 0: No ratio plot; !=0: ratio between the top and bottom pads
BLIND       = False
POISSON     = False
jobs        = []
REBIN       = True#options.rebin#False#True
NICE        = False

########## SAMPLES ##########

#sign = ["XZZ_M600", "XZZ_M650", "XZZ_M700", "XZZ_M750", "XZZ_M800", "XZZ_M1000", "XZZ_M1200", "XZZ_M1400", "XZZ_M1800", "XZZ_M2000", "XZZ_M2500", "XZZ_M3000", "XZZ_M3500", "XZZ_M4000", "XZZ_M4500"]
### WprimeToWZToWhadZlep
#sign = ["XWZ_M600", "XWZ_M800", "XWZ_M1000", "XWZ_M1200", "XWZ_M1400", "XWZ_M1600", "XWZ_M1800", "XWZ_M2000", "XWZ_M2500", "XWZ_M3000", "XWZ_M3500", "XWZ_M4000", "XWZ_M4500"]
#ZhadZinv
sign = ["XZZInv_M600", "XZZInv_M800", "XZZInv_M1000", "XZZInv_M1200", "XZZInv_M1400", "XZZInv_M1600", "XZZInv_M1800", "XZZInv_M2000", "XZZInv_M2500", "XZZInv_M3000", "XZZInv_M3500", "XZZInv_M4500"]#"XZZInv_M4000", 
colors = [634, 410, 856, 2, 401, 418, 881, 798, 602, 921]
########## ######## ##########

#gROOT.SetBatch(True)

sign_sampl = {}
sign_sampl_zb = {
    'ZeroBiasRunC' : {
        'order' : 5,
        'files' : ['ZeroBiasRun2016C-23Sep2016-v1'],
        'howmany' : "RunC",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "ZeroBias",
        'weight': 1.,
        'plot': True,
    },
}
sign_sampl_mu = {
    'test' : {
        'order' : 5,
        'files' : ['output'],
        'howmany' : "RunA",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleMu",
        'nice_label' : "Single Muon RunA",
        'weight': 1.,
        'plot': True,
    },
    'SingleMuRunA' : {
        'order' : 5,
        'files' : ['SingleMuon_Run2017A-PromptReco-v2', 'SingleMuon_Run2017A-PromptReco-v3'],
        'howmany' : "RunA",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleMu",
        'nice_label' : "Single Muon RunA",
        'weight': 1.,
        'plot': True,
    },
    'SingleMuRunB' : {
        'order' : 5,
        'files' : ['SingleMuon_Run2017B-PromptReco-v1','SingleMuon_Run2017B-PromptReco-v2'],
        'howmany' : "RunB",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleMu",
        'nice_label' : "Single Muon RunB",
        'weight': 1.,
        'plot': True,
    },
    'SingleMuRunB-v1' : {
        'order' : 5,
        'files' : ['SingleMuon_Run2017B-PromptReco-v1'],
        'howmany' : "RunB-v1",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleMu",
        'weight': 1.,
        'plot': True,
    },
    'SingleMuRunB-v2' : {
        'order' : 5,
        'files' : ['SingleMuon_Run2017B-PromptReco-v2'],
        'howmany' : "RunB-v2",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleMu",
        'weight': 1.,
        'plot': True,
    },
    'SingleMuRunC-v1' : {
        'order' : 5,
        'files' : ['SingleMuon_Run2017C-PromptReco-v1'],
        'howmany' : "RunC-v1",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleMu",
        'weight': 1.,
        'plot': True,
    },
    'SingleMuRunC-v2' : {
        'order' : 5,
        'files' : ['SingleMuon_Run2017C-PromptReco-v2'],
        'howmany' : "RunC-v2",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleMu",
        'weight': 1.,
        'plot': True,
    },
    'SingleMuRunC-v3' : {
        'order' : 5,
        'files' : ['SingleMuon_Run2017C-PromptReco-v3'],
        'howmany' : "RunC-v3",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleMu",
        'weight': 1.,
        'plot': True,
    },
    'SingleMuRunC' : {
        'order' : 5,
        'files' : ['SingleMuon_Run2017C-PromptReco-v1','SingleMuon_Run2017C-PromptReco-v2','SingleMuon_Run2017C-PromptReco-v3'],
        'howmany' : "RunC",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleMu",
        'nice_label' : "Single Muon RunC",
        'weight': 1.,
        'plot': True,
    },
    'SingleMuRunD-v1' : {
        'order' : 5,
        'files' : ['SingleMuon_Run2017D-PromptReco-v1'],
        'howmany' : "RunD-v1",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleMu",
        'weight': 1.,
        'plot': True,
    },
    'SingleMuRunD' : {
        'order' : 5,
        'files' : ['SingleMuon_Run2017D-PromptReco-v1'],
        'howmany' : "RunD",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleMu",
        'nice_label' : "Single Muon RunD",
        'weight': 1.,
        'plot': True,
    },
    'SingleMuAll' : {
        'order' : 5,
        'files' : ['SingleMuon_Run2017B-PromptReco-v1', 'SingleMuon_Run2017B-PromptReco-v2', 'SingleMuon_Run2017C-PromptReco-v1','SingleMuon_Run2017C-PromptReco-v2','SingleMuon_Run2017C-PromptReco-v3','SingleMuon_Run2017D-PromptReco-v1'],
        'howmany' : "",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleMu",
        'nice_label' : "Single Muon",
        'weight': 1.,
        'plot': True,
    },
}

sign_sampl_ele = {
    'SingleEleRunB' : {
        'order' : 5,
        'files' : ['SingleElectron_Run2016B-03Feb2017_ver2-v2'],
        'howmany' : "RunB",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleEle",
        'weight': 1.,
        'plot': True,
    },
    'SingleEleRunC' : {
        'order' : 5,
        'files' : ['SingleElectron_Run2016C-03Feb2017-v1'],
        'howmany' : "RunC",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleEle",
        'weight': 1.,
        'plot': True,
    },
    'SingleEleRunD' : {
        'order' : 5,
        'files' : ['SingleElectron_Run2016D-03Feb2017-v1'],
        'howmany' : "RunD",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleEle",
        'weight': 1.,
        'plot': True,
    },
    'SingleEleRunE' : {
        'order' : 5,
        'files' : ['SingleElectron_Run2016E-03Feb2017-v1'],
        'howmany' : "RunE",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleEle",
        'weight': 1.,
        'plot': True,
    },
    'SingleEleRunF' : {
        'order' : 5,
        'files' : ['SingleElectron_Run2016F-03Feb2017-v1'],
        'howmany' : "RunF",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleEle",
        'weight': 1.,
        'plot': True,
    },
    'SingleEleRunG' : {
        'order' : 5,
        'files' : ['SingleElectron_Run2016G-03Feb2017-v1'],
        'howmany' : "RunG",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleEle",
        'weight': 1.,
        'plot': True,
    },
    'SingleEleRunH' : {
        'order' : 5,
        'files' : ['SingleElectron_Run2016H-03Feb2017_ver2-v1','SingleElectron_Run2016H-03Feb2017_ver3-v1'],
        'howmany' : "RunH",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleEle",
        'weight': 1.,
        'plot': True,
    },
    'SingleEleAll' : {
        'order' : 5,
        'files' : ['SingleElectron_Run2016B-03Feb2017_ver2-v2','SingleElectron_Run2016C-03Feb2017-v1','SingleElectron_Run2016D-03Feb2017-v1','SingleElectron_Run2016E-03Feb2017-v1','SingleElectron_Run2016F-03Feb2017-v1','SingleElectron_Run2016G-03Feb2017-v1','SingleElectron_Run2016H-03Feb2017_ver2-v1','SingleElectron_Run2016H-03Feb2017_ver3-v1'],
        'howmany' : "AllRuns",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "SingleEle",
        'weight': 1.,
        'plot': True,
    },

}

#JetHT
sign_sampl_jetHT = {
    'JetHTRunB' : {
        'order' : 5,
        'files' : ['JetHT_Run2016B-03Feb2017_ver2-v2'],
        'howmany' : "RunB",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "JetHT",
        'weight': 1.,
        'plot': True,
    },
    'JetHTRunC' : {
        'order' : 5,
        'files' : ['JetHT_Run2016C-03Feb2017-v1'],
        'howmany' : "RunC",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "JetHT",
        'weight': 1.,
        'plot': True,
    },
    'JetHTRunD' : {
        'order' : 5,
        'files' : ['JetHT_Run2016D-03Feb2017-v1'],
        'howmany' : "RunD",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "JetHT",
        'weight': 1.,
        'plot': True,
    },
    'JetHTRunE' : {
        'order' : 5,
        'files' : ['JetHT_Run2016E-03Feb2017-v1'],
        'howmany' : "RunE",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "JetHT",
        'weight': 1.,
        'plot': True,
    },
    'JetHTRunF' : {
        'order' : 5,
        'files' : ['JetHT_Run2016F-03Feb2017-v1'],
        'howmany' : "RunF",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "JetHT",
        'weight': 1.,
        'plot': True,
    },
    'JetHTRunG' : {
        'order' : 5,
        'files' : ['JetHT_Run2016G-03Feb2017-v1'],
        'howmany' : "RunG",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "JetHT",
        'weight': 1.,
        'plot': True,
    },
    'JetHTRunH' : {
        'order' : 5,
        'files' : ['JetHT_Run2016H-03Feb2017_ver2-v1','JetHT_Run2016H-03Feb2017_ver3-v1'],
        'howmany' : "RunH",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "JetHT",
        'weight': 1.,
        'plot': True,
    },
    'JetHTAll' : {
        'order' : 5,
        'files' : ['JetHT_Run2016B-03Feb2017_ver2-v2','JetHT_Run2016C-03Feb2017-v1','JetHT_Run2016D-03Feb2017-v1','JetHT_Run2016E-03Feb2017-v1','JetHT_Run2016F-03Feb2017-v1','JetHT_Run2016G-03Feb2017-v1','JetHT_Run2016H-03Feb2017_ver2-v1','JetHT_Run2016H-03Feb2017_ver3-v1'],
        'howmany' : "AllRuns",
        'fillcolor' : 60,
        'fillstyle' : 1001,
        'linecolor' : 60,
        'linewidth' : 2,
        'linestyle' : 1,
        'marker' : 25,
        'label' : "JetHT",
        'weight': 1.,
        'plot': True,
    },

}

chain = {}
hist = {}
graph = {}
graph_200 = {}
graph_250 = {}
graph_300 = {}
hist_met = {}
hist_num = {}
hist_num110 = {}
hist_num120 = {}
hist_num120_200 = {}
hist_num120_250 = {}
hist_num120_300 = {}
hist_den_200 = {}
hist_den_250 = {}
hist_den_300 = {}
hist_num130 = {}
hist_num140 = {}
hist_numOR = {}
hist_num170 = {}
hist_den = {}
chain_num = {}
chain_num110 = {}
chain_num120 = {}
chain_num130 = {}
chain_num140 = {}
chain_numOR = {}
chain_num170 = {}
chain_den = {}
goodstring = ''
goodlabel = ''
var = "met_pt_nomu"#"nFatJets"#"Lepton1.pt"#"FatJet1.pt"#"MEt.phi"#
#cut_num = "isZtoNN && FatJet1.pt>170 && HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v"
#cut_den = "isZtoNN && FatJet1.pt>170"
can = TCanvas("can","can", 1000, 800)
can.SetGrid()
can.cd()
if options.variable == "nPV":
    if NICE:
        leg = TLegend(0.12, 0.12, 0.48, 0.50)
    else:
        leg = TLegend(0.12, 0.12, 0.52, 0.52)#0.12, 0.12, 0.56, 0.50
else:
    if NICE:
        leg = TLegend(0.52, 0.12, 0.88, 0.3)#0.45)
    else:
        leg = TLegend(0.38, 0.12, 0.88, 0.3)#0.45)
leg.SetTextSize(0.03)
#set of cuts:
# 0 Alejandro
# 1 Alejandro
# 2 VZ
# 3 VH-Monojet
# 4 Veto
if options.lepton == "mu0":#Alejandro 0
    sign_sampl = sign_sampl_mu
    cut_den = "nFatJets>0 && FatJet1.isTight && FatJet1.pt>200"
    cut_num90 = "nFatJets>0 && FatJet1.isTight && FatJet1.pt>200 && (HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v)"
    cut_num110 = "nFatJets>0 && FatJet1.isTight && FatJet1.pt>200 && (HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v)"
    cut_num120 = "nFatJets>0 && FatJet1.isTight && FatJet1.pt>200 && (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v)"
    cut_num170 = "nFatJets>0 && FatJet1.isTight && FatJet1.pt>200 && (HLT_PFMET170_NoiseCleaned_v || HLT_PFMET170_JetIdCleaned_v || HLT_PFMET170_HBHECleaned_v)"
    cut_numOR = "nFatJets>0 && FatJet1.isTight && FatJet1.pt>200 && (HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v || HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v)"

elif options.lepton == "mu3":#VH
    sign_sampl = sign_sampl_mu
    cut_den = "nTightFatJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_globalSuperTightHalo2016Filter && Flag_BadChCand && Flag_BadPFMuon)"
    cut_num110 = cut_den + " && HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v"#"nTightJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_globalSuperTightHalo2016Filter && Flag_BadChCand && Flag_BadPFMuon) && (HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v)"
    cut_num120 = cut_den + " && (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v)"#"nTightJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_globalSuperTightHalo2016Filter && Flag_BadChCand && Flag_BadPFMuon) && (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v)"
    cut_num130 = cut_den + " && HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v"#"nTightJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_globalSuperTightHalo2016Filter && Flag_BadChCand && Flag_BadPFMuon) && (HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v)"
    cut_num140 = cut_den + " && HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v"#"nTightJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_globalSuperTightHalo2016Filter && Flag_BadChCand && Flag_BadPFMuon) && (HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v)"
    cut_numOR = cut_den + " && ((HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v) || HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v || HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v)"#"nTightJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_globalSuperTightHalo2016Filter && Flag_BadChCand && Flag_BadPFMuon) && (HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v || HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v || HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v)"

elif options.lepton == "mu3nPV":
    sign_sampl = sign_sampl_mu
    cut_den = "nTightFatJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_globalSuperTightHalo2016Filter && Flag_BadChCand && Flag_BadPFMuon) && met_pt_nomu_L>200"
    cut_den_200 = "nTightFatJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_globalSuperTightHalo2016Filter && Flag_BadChCand && Flag_BadPFMuon) && met_pt_nomu_L>200"
    cut_den_250 = "nTightFatJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_globalSuperTightHalo2016Filter && Flag_BadChCand && Flag_BadPFMuon) && met_pt_nomu_L>250"
    cut_den_300 = "nTightFatJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_globalSuperTightHalo2016Filter && Flag_BadChCand && Flag_BadPFMuon) && met_pt_nomu_L>300"
    cut_num110 = cut_den + " && HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v"
    cut_num120 = cut_den + " && (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v)"
    cut_num120_200 = cut_den_200 + " && (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v)"
    cut_num120_250 = cut_den_250 + " && (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v)"
    cut_num120_300 = cut_den_300 + " && (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v)"
    cut_num130 = cut_den + " && HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v"
    cut_num140 = cut_den + " && HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v"
    cut_numOR = cut_den + " && ((HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v) || HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v || HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v)"#"nTightJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_globalSuperTightHalo2016Filter && Flag_BadChCand && Flag_BadPFMuon) && (HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v || HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v || HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v)"


elif options.lepton == "mu3nofilt":#VH
    sign_sampl = sign_sampl_mu
    cut_den = "nTightFatJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v)"
    cut_num110 = "nTightJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v) && (HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v)"
    cut_num120 = "nTightJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v) && ((HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v))"
    cut_num130 = "nTightJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v) && (HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v)"
    cut_num140 = "nTightJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v) && (HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v)"
    cut_num170 = "nTightJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v)"
    cut_numOR = "nTightJets>0 && nTightMuons==1 && Muon1_pt>35 && Muon1_pfIso04<0.15 && Muon1_isTight && FatJet1_isTight && FatJet1_pt>170 && (HLT_IsoMu24_v) && (HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v || (HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v) || HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v || HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v)"



s=''
if (options.run=="A" and ("mu" in str(options.lepton))):
    s+="SingleMuRunA"
elif (options.run=="test" and ("mu" in str(options.lepton))):
    s+="test"
elif (options.run=="B-v1" and ("mu" in str(options.lepton))):
    s+="SingleMuRunB-v1"
elif (options.run=="B-v2" and ("mu" in str(options.lepton))):
    s+="SingleMuRunB-v2"
elif (options.run=="B" and ("mu" in str(options.lepton))):
    s+="SingleMuRunB"
elif (options.run=="C-v1" and ("mu" in str(options.lepton))):
    s+="SingleMuRunC-v1"
elif (options.run=="C-v2" and ("mu" in str(options.lepton))):
    s+="SingleMuRunC-v2"
elif (options.run=="C-v3" and ("mu" in str(options.lepton))):
    s+="SingleMuRunC-v3"
elif (options.run=="C" and ("mu" in str(options.lepton))):
    s+="SingleMuRunC"
elif (options.run=="D" and ("mu" in str(options.lepton))):
    s+="SingleMuRunD"
elif (options.run=="E" and ("mu" in str(options.lepton))):
    s+="SingleMuRunE"
elif (options.run=="F" and ("mu" in str(options.lepton))):
    s+="SingleMuRunF"
elif (options.run=="G" and ("mu" in str(options.lepton))):
    s+="SingleMuRunG"
elif (options.run=="H" and ("mu" in str(options.lepton))):
    s+="SingleMuRunH"
elif (options.run=="All" and ("mu" in str(options.lepton))):
    s+="SingleMuAll"
elif (options.run=="B" and ("ele" in str(options.lepton))):
    s+="SingleEleRunB"
elif (options.run=="C" and ("ele" in str(options.lepton))):
    s+="SingleEleRunC"
elif (options.run=="D" and ("ele" in str(options.lepton))):
    s+="SingleEleRunD"
elif (options.run=="E" and ("ele" in str(options.lepton))):
    s+="SingleEleRunE"
elif (options.run=="F" and ("ele" in str(options.lepton))):
    s+="SingleEleRunF"
elif (options.run=="G" and ("ele" in str(options.lepton))):
    s+="SingleEleRunG"
elif (options.run=="H" and ("ele" in str(options.lepton))):
    s+="SingleEleRunH"
elif (options.run=="All" and ("ele" in str(options.lepton))):
    s+="SingleEleAll"
elif (options.run=="C" and ("zb" in str(options.lepton))):
    s+="ZeroBiasRunC"
elif (options.run=="B" and ("jetHT" in str(options.lepton))):
    s+="JetHTRunB"
elif (options.run=="C" and ("jetHT" in str(options.lepton))):
    s+="JetHTRunC"
elif (options.run=="D" and ("jetHT" in str(options.lepton))):
    s+="JetHTRunD"
elif (options.run=="E" and ("jetHT" in str(options.lepton))):
    s+="JetHTRunE"
elif (options.run=="F" and ("jetHT" in str(options.lepton))):
    s+="JetHTRunF"
elif (options.run=="G" and ("jetHT" in str(options.lepton))):
    s+="JetHTRunG"
elif (options.run=="H" and ("jetHT" in str(options.lepton))):
    s+="JetHTRunH"
elif (options.run=="All" and ("jetHT" in str(options.lepton))):
    s+="JetHTAll"
print "faccio l'istogramma!", s

chain[s] = TChain("trigger/tree")
print "s: ", s
for j,ss in enumerate(sign_sampl[s]['files']):
    print "ss: ", ss
    chain[s].Add(NTUPLEDIR + ss + ".root")
if variable[options.variable]['nbins']>0:
    hist[s] = TH1F(s, ";"+variable[options.variable]['title'], variable[options.variable]['nbins'], variable[options.variable]['min'], variable[options.variable]['max'])
    if options.variable == "nPV":
        hist_den_200[s] = TH1F(s+"_den_200", ";"+variable[options.variable]['title'], variable[options.variable]['nbins'], variable[options.variable]['min'], variable[options.variable]['max'])
        hist_den_250[s] = TH1F(s+"_den_250", ";"+variable[options.variable]['title'], variable[options.variable]['nbins'], variable[options.variable]['min'], variable[options.variable]['max'])
        hist_den_300[s] = TH1F(s+"_den_300", ";"+variable[options.variable]['title'], variable[options.variable]['nbins'], variable[options.variable]['min'], variable[options.variable]['max'])
        hist_num120_200[s] = TH1F(s+"_num120_200", ";"+variable[options.variable]['title'], variable[options.variable]['nbins'], variable[options.variable]['min'], variable[options.variable]['max'])
        hist_num120_250[s] = TH1F(s+"_num120_250", ";"+variable[options.variable]['title'], variable[options.variable]['nbins'], variable[options.variable]['min'], variable[options.variable]['max'])
        hist_num120_300[s] = TH1F(s+"_num120_300", ";"+variable[options.variable]['title'], variable[options.variable]['nbins'], variable[options.variable]['min'], variable[options.variable]['max'])
    else:
        hist_num120[s] = TH1F(s+"_num120", ";"+variable[options.variable]['title'], variable[options.variable]['nbins'], variable[options.variable]['min'], variable[options.variable]['max'])
        hist_num120[s].Sumw2()

hist[s].Sumw2()

##hist_met[s].Sumw2()
#hist_num110[s].Sumw2()
if options.variable == "nPV":
    hist_den_200[s].Sumw2()
    hist_den_250[s].Sumw2()
    hist_den_300[s].Sumw2()
    hist_num120_200[s].Sumw2()
    hist_num120_250[s].Sumw2()
    hist_num120_300[s].Sumw2()

print "cut denominator: ", cut_den
chain[s].Project(s, options.variable, cut_den)           
hist[s].SetOption("%s" % chain[s].GetTree().GetEntriesFast())

if options.variable == "nPV":
    chain[s].Project(s+"_den_200", options.variable, cut_den_200)
    hist_den_200[s].SetOption("%s" % chain[s].GetTree().GetEntriesFast())
    chain[s].Project(s+"_den_250", options.variable, cut_den_250)
    hist_den_250[s].SetOption("%s" % chain[s].GetTree().GetEntriesFast())
    chain[s].Project(s+"_den_300", options.variable, cut_den_300)
    hist_den_300[s].SetOption("%s" % chain[s].GetTree().GetEntriesFast())
    
    chain[s].Project(s+"_num120_200", options.variable, cut_num120_200)
    hist_num120_200[s].SetOption("%s" % chain[s].GetTree().GetEntriesFast())
    chain[s].Project(s+"_num120_250", options.variable, cut_num120_250)
    hist_num120_250[s].SetOption("%s" % chain[s].GetTree().GetEntriesFast())
    chain[s].Project(s+"_num120_300", options.variable, cut_num120_300)
    hist_num120_300[s].SetOption("%s" % chain[s].GetTree().GetEntriesFast())

else:
    chain[s].Project(s+"_num120", options.variable, cut_num120)
    hist_num120[s].SetOption("%s" % chain[s].GetTree().GetEntriesFast())
    hist_num120[s].SetMarkerSize(1.)
    hist_num120[s].SetMarkerStyle(21)#(sign_sampl[s]['marker'])
    hist_num120[s].SetMarkerColor(418)#(2)
    hist_num120[s].SetFillColor(418)#(2) 
    hist_num120[s].SetLineColor(418)#(2)
    hist_num120[s].SetLineWidth(2)
    hist_num120[s].GetYaxis().SetTitleOffset(1.2)
    hist_num120[s].GetYaxis().SetTitle("Efficiency")
#hist[s].GetYaxis().SetTitle("Efficiency")
#hist[s].GetXaxis().SetTitle(variable[options.variable]['title'])

if options.variable == "nPV":
    hist_num120_200[s].SetMarkerSize(1.)
    hist_num120_200[s].SetMarkerStyle(21)#(sign_sampl[s]['marker'])
    hist_num120_200[s].SetMarkerColor(418)#(2)
    hist_num120_200[s].SetFillColor(418)#(2) 
    hist_num120_200[s].SetLineColor(418)#(2)
    hist_num120_200[s].SetLineWidth(2)
    
    hist_num120_250[s].SetMarkerSize(1.)
    hist_num120_250[s].SetMarkerStyle(21)#(sign_sampl[s]['marker'])
    hist_num120_250[s].SetMarkerColor(801)#(2)
    hist_num120_250[s].SetFillColor(801)#(2) 
    hist_num120_250[s].SetLineColor(801)#(2)
    hist_num120_250[s].SetLineWidth(2)
    
    hist_num120_300[s].SetMarkerSize(1.)
    hist_num120_300[s].SetMarkerStyle(21)#(sign_sampl[s]['marker'])
    hist_num120_300[s].SetMarkerColor(866)#(2)
    hist_num120_300[s].SetFillColor(866)#(2) 
    hist_num120_300[s].SetLineColor(866)#(2)
    hist_num120_300[s].SetLineWidth(2)

    hist_num120_200[s].GetYaxis().SetTitle("Efficiency")
    hist_num120_250[s].GetYaxis().SetTitle("Efficiency")
    hist_num120_300[s].GetYaxis().SetTitle("Efficiency")


if options.variable == "nPV":
    bins = np.array([0.,10.,15.,20.,25.,30.,35.,40.,50.,70.])
else:
    bins = np.array([0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,225.,250.,275.,300.,325.,350.,375.,400.,425.,450.,475.,500.,525.,550.,575.,600.,625.,650.,675.,700.,725.,750.,775.,800.,850.,900.,1000.])

if REBIN:
    if options.variable == "nPV":
        hist_den_200[s] = hist_den_200[s].Rebin(len(bins)-1,s+"_den_200",bins)
        hist_num120_200[s] = hist_num120_200[s].Rebin(len(bins)-1,s+"_num120_200",bins)
        hist_den_250[s] = hist_den_250[s].Rebin(len(bins)-1,s+"_den_250",bins)
        hist_num120_250[s] = hist_num120_250[s].Rebin(len(bins)-1,s+"_num120_250",bins)
        hist_den_300[s] = hist_den_300[s].Rebin(len(bins)-1,s+"_den_300",bins)
        hist_num120_300[s] = hist_num120_300[s].Rebin(len(bins)-1,s+"_num120_300",bins)

        graph_300[s] = TGraphAsymmErrors()
        graph_300[s].BayesDivide(hist_num120_300[s],hist_den_300[s])
        graph_300[s].SetMarkerSize(1.)
        graph_300[s].SetMarkerStyle(21)#(sign_sampl[s]['marker'])
        graph_300[s].SetMarkerColor(866)#(2)
        graph_300[s].SetFillColor(866)#(2) 
        graph_300[s].SetLineColor(866)#(2)
        graph_300[s].SetLineWidth(2)
        graph_300[s].GetYaxis().SetTitleOffset(1.2)
        graph_300[s].GetYaxis().SetTitle("Efficiency")
        graph_300[s].GetYaxis().SetRangeUser(0.6,1.01)
        graph_300[s].SetMinimum(0.8)
        graph_300[s].SetMaximum(1.05)
        graph_300[s].GetXaxis().SetTitle(variable[options.variable]['title'])
        graph_300[s].GetXaxis().SetRangeUser(variable[options.variable]['min'], variable[options.variable]['max'])
        graph_300[s].Draw("AP,sames")

        graph_250[s] = TGraphAsymmErrors()
        graph_250[s].BayesDivide(hist_num120_250[s],hist_den_250[s])
        graph_250[s].SetMarkerSize(1.)
        graph_250[s].SetMarkerStyle(21)#(sign_sampl[s]['marker'])
        graph_250[s].SetMarkerColor(801)#(2)
        graph_250[s].SetFillColor(801)#(2) 
        graph_250[s].SetLineColor(801)#(2)
        graph_250[s].SetLineWidth(2)
        graph_250[s].GetYaxis().SetTitleOffset(1.2)
        graph_250[s].GetYaxis().SetTitle("Efficiency")
        graph_250[s].GetYaxis().SetRangeUser(0.6,1.01)
        graph_250[s].SetMinimum(0.8)
        graph_250[s].SetMaximum(1.05)
        graph_250[s].GetXaxis().SetTitle(variable[options.variable]['title'])
        graph_250[s].GetXaxis().SetRangeUser(variable[options.variable]['min'], variable[options.variable]['max'])
        graph_250[s].Draw("P,sames")

        graph_200[s] = TGraphAsymmErrors()
        graph_200[s].BayesDivide(hist_num120_200[s],hist_den_200[s])
        graph_200[s].SetMarkerSize(1.)
        graph_200[s].SetMarkerStyle(21)#(sign_sampl[s]['marker'])
        graph_200[s].SetMarkerColor(418)#(2)
        graph_200[s].SetFillColor(418)#(2) 
        graph_200[s].SetLineColor(418)#(2)
        graph_200[s].SetLineWidth(2)
        graph_200[s].GetYaxis().SetTitleOffset(1.2)
        graph_200[s].GetYaxis().SetTitle("Efficiency")
        graph_200[s].GetYaxis().SetRangeUser(0.6,1.01)
        graph_200[s].SetMinimum(0.8)
        graph_200[s].SetMaximum(1.05)
        graph_200[s].GetXaxis().SetTitle(variable[options.variable]['title'])
        graph_200[s].GetXaxis().SetRangeUser(variable[options.variable]['min'], variable[options.variable]['max'])
        graph_200[s].Draw("P,sames")

        if NICE:
            leg.AddEntry(graph_200[s], "HLT #slash{E}_{T}^{no #mu}, #slash{H}_{T}^{no #mu} > 120 GeV",'PL')
            leg.AddEntry("", "offline #slash{E}_{T}^{no #mu} > 200 GeV",'')
            leg.AddEntry(graph_250[s], "HLT #slash{E}_{T}^{no #mu}, #slash{H}_{T}^{no #mu} > 120 GeV",'PL')
            leg.AddEntry("", "offline #slash{E}_{T}^{no #mu} > 250 GeV",'')
            leg.AddEntry(graph_300[s], "HLT #slash{E}_{T}^{no #mu}, #slash{H}_{T}^{no #mu} > 120 GeV",'PL')
            leg.AddEntry("", "offline #slash{E}_{T}^{no #mu} > 300 GeV",'')
        else:
            leg.SetHeader("HLT E_{T}^{miss}(no #mu), H_{T}^{miss}(no #mu) > 120 GeV","C")#("","HLT E_{T}^{miss}(no #mu), H_{T}^{miss}(no #mu) > 120 GeV","")
            leg.AddEntry(graph_200[s], "Offline E_{T}^{miss}(no #mu) > 200 GeV",'PL')
                #leg.AddEntry("", "offline E_{T}^{miss}(no #mu) > 200 GeV",'')
            leg.AddEntry(graph_250[s], "Offline E_{T}^{miss}(no #mu) > 250 GeV",'PL')
                #leg.AddEntry("", "offline E_{T}^{miss}(no #mu) > 250 GeV",'')
            leg.AddEntry(graph_300[s], "Offline E_{T}^{miss}(no #mu) > 300 GeV",'PL')
                #leg.AddEntry("", "offline E_{T}^{miss}(no #mu) > 300 GeV",'')
    else:
        hist[s] = hist[s].Rebin(len(bins)-1,s+"_num90",bins)
        hist_num120[s] = hist_num120[s].Rebin(len(bins)-1,s+"_num120",bins)

        graph[s] = TGraphAsymmErrors()
        graph[s].BayesDivide(hist_num120[s],hist[s])
        graph[s].SetMarkerSize(1.)
        graph[s].SetMarkerStyle(21)#(sign_sampl[s]['marker'])
        graph[s].SetMarkerColor(418)#(2)
        graph[s].SetFillColor(418)#(2) 
        graph[s].SetLineColor(418)#(2)
        graph[s].SetLineWidth(2)
        graph[s].GetYaxis().SetTitle("Efficiency")
        graph[s].GetYaxis().SetRangeUser(0.,1.01)
        graph[s].GetXaxis().SetRangeUser(variable[options.variable]['min'], variable[options.variable]['max'])
        graph[s].GetXaxis().SetTitle(variable[options.variable]['title'])
        graph[s].Draw("AP,sames")
        #print "HERE!"
        #hist[s].Draw("PL")
        #hist_num120[s].Draw("APL,SAMES")
        if NICE:
            leg.AddEntry(graph[s], "HLT #slash{E}_{T}^{no #mu}, #slash{H}_{T}^{no #mu} > 120 GeV",'PL')
        else:
            leg.AddEntry(graph[s], "HLT E_{T}^{miss}(no #mu), H_{T}^{miss}(no #mu) > 120 GeV",'PL')
                

leg.SetBorderSize(0)
leg.Draw()

etichetta = TLatex()
etichetta.SetNDC()
etichetta.SetTextSize(0.04)
etichetta.SetTextColor(1)
if not options.public:
    etichetta.DrawLatex(0.3, 0.4, sign_sampl[s]['nice_label'])
if options.goodonlinejec and not options.public:
    good_etichetta = TLatex()
    good_etichetta.SetNDC()
    good_etichetta.SetTextSize(0.03)
    good_etichetta.SetTextColor(1)
    good_etichetta.DrawLatex(0.3, 0.3, goodlabel)
latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.SetTextColor(1)
latex.SetTextFont(42)
latex.SetTextAlign(33)
latex.DrawLatex(0.9, 0.96, "%.1f fb^{-1}  (13 TeV, 2018)" % (float(LUMI)))
latex.SetTextSize(0.04)
latex.SetTextFont(62)
latex.DrawLatex(0.20, 0.96, "CMS")
latex.SetTextFont(52)
latex.DrawLatex(0.36, 0.96, "Preliminary")

can.Update()
outpath = ""
if options.public:
    outpath += "$CMSSW_BASE/src/TrigAnalyzer/Trigger4B2G/macro/plots_public/"
else:
    outpath += "$CMSSW_BASE/src/TrigAnalyzer/Trigger4B2G/macro/plots/"

if options.goodonlinejec:
    goodstring = "_runnumber_299504"


can.Print(outpath + "TriggerTurnOn_" + sign_sampl[s]['label'] + sign_sampl[s]['howmany'] + "_" + str(options.variable)  + "_" + options.lepton + goodstring+"_v11_longlabel.png")
can.Print(outpath + "TriggerTurnOn_" + sign_sampl[s]['label'] + sign_sampl[s]['howmany'] + "_" + str(options.variable)  +  "_" + options.lepton + goodstring+"_v11_longlabel.pdf")

if not gROOT.IsBatch(): raw_input("Press Enter to continue...") 
