import ROOT
from ROOT import *
import os,sys,time
import math
from optparse import OptionParser
from array import array



def getBinning(binsMVV,minx,maxx,bins):
    l=[]
    if binsMVV=="":
        for i in range(0,bins+1):
            l.append(minx + i* (maxx - minx)/bins)
    else:
        s = binsMVV.split(",")
        for w in s:
            l.append(int(w))
    return l

def passMETfilters(tree):

 if not tree.Flag_HBHENoiseFilter: return False
 if not tree.Flag_HBHENoiseIsoFilter: return False
 if not tree.Flag_EcalDeadCellTriggerPrimitiveFilter: return False
 if not tree.Flag_goodVertices: return False
 if not tree.Flag_eeBadScFilter: return False
 if not tree.Flag_globalSuperTightHalo2016Filter: return False
 if not tree.Flag_BadChCand: return False
 if not tree.Flag_BadPFMuon: return False
 return True


parser = OptionParser()
parser.add_option('-f','--filelist',action='store',type='string',dest='filelist',default='list.txt', help='file list')
parser.add_option('-o','--outfname',action='store',type='string',dest='outfname',default='jet_triggers.root', help='outfile name')
(options,args) = parser.parse_args()
    
files = []
f = open(options.filelist,'r')
for l in f.readlines():
  files.append(l.replace('\n',''))

print files

bins = "0,3,6,10,16,23,31,40,50,61,74,88,103,119,137,156,176,197,220,244,270,296,325,354,386,419,453,489,526,565,606,649,693,740,788,838,890,944,1000,1058,1118,1181,1246,1313,1383,1455,1530,1607,1687,1770,1856,1945,2037,2132,2231,2332,2438,2546,2659,2775,2895,3019,3147,3279,3416,3558,3704,3854,4010,4171,4337,4509,4686,4869,5058,5253,5455,5663,5877,6099,6328,6564,7000"
binning = getBinning(bins,0,7000,82)

outf = ROOT.TFile(options.outfname,'RECREATE')

mjj_den = ROOT.TH1F('mjj_den','mjj_den',82,array('f',binning))

mjj_num_pfht1050 = ROOT.TH1F('mjj_num_HLT_PFHT1050','mjj_num_HLT_PFHT1050',82,array('f',binning))
mjj_num_pfjet500 = ROOT.TH1F('mjj_num_HLT_AK8PFJet500','mjj_num_HLT_AK8PFJet500',82,array('f',binning))
mjj_num_pfjet400_trimmass30 = ROOT.TH1F('mjj_num_HLT_AK8PFJet400_TrimMass30','mjj_num_HLT_AK8PFJet400_TrimMass300',82,array('f',binning))    
mjj_num_pfht800_trimmass50 = ROOT.TH1F('mjj_num_HLT_AK8PFHT800_TrimMass50','mjj_num_HLT_AK8PFHT800_TrimMass50',82,array('f',binning))
mjj_num_all_triggers = ROOT.TH1F('mjj_num_all_triggers','mjj_num_all_triggers',82,array('f',binning))    

mjet1_den = ROOT.TH1F('mjet1_den','mjet1_den',50,0,250)
mjet1_num_pfht1050 = ROOT.TH1F('mjet1_num_HLT_PFHT1050','mjet1_num_HLT_PFHT1050',50,0,250)
mjet1_num_pfjet500 = ROOT.TH1F('mjet1_num_HLT_AK8PFJet500','mjet1_num_HLT_AK8PFJet500',50,0,250)
mjet1_num_pfjet400_trimmass30 = ROOT.TH1F('mjet1_num_HLT_AK8PFJet400_TrimMass30','mjet1_num_HLT_AK8PFJet400_TrimMass30',50,0,250)
mjet1_num_pfht800_trimmass50 = ROOT.TH1F('mjet1_num_HLT_AK8PFHT800_TrimMass50','mjet1_num_HLT_AK8PFHT800_TrimMass50',50,0,250)
mjet1_num_all_triggers = ROOT.TH1F('mjet1_num_all_triggers','mjet1_num_all_triggers',50,0,250)    

mjet2_den = ROOT.TH1F('mjet2_den','mjet2_den',50,0,250)
mjet2_num_pfht1050 = ROOT.TH1F('mjet2_num_HLT_PFHT1050','mjet2_num_HLT_PFHT1050',50,0,250)
mjet2_num_pfjet500 = ROOT.TH1F('mjet2_num_HLT_AK8PFJet500','mjet2_num_HLT_AK8PFJet500',50,0,250)
mjet2_num_pfjet400_trimmass30 = ROOT.TH1F('mjet2_num_HLT_AK8PFJet400_TrimMass30','mjet2_num_HLT_AK8PFJet400_TrimMass30',50,0,250)
mjet2_num_pfht800_trimmass50 = ROOT.TH1F('mjet2_num_HLT_AK8PFHT800_TrimMass50','mjet2_num_HLT_AK8PFHT800_TrimMass50',50,0,250)
mjet2_num_all_triggers = ROOT.TH1F('mjet2_num_all_triggers','mjet2_num_all_triggers',50,0,250)    

nfiles = len(files)
countf = 0 
for f in files:
 
 if countf%10 == 0: print "Analyzing file",countf+1,"of",nfiles
 countf+=1
 
 tf = ROOT.TFile.Open(f,'READ')
 t = tf.Get('trigger/tree')
 
 for e in range(t.GetEntries()):
  t.GetEntry(e)

  if not passMETfilters(t): continue
  
  passReference = False
  if t.HLT_IsoMu24_v or t.HLT_IsoMu27_v or t.HLT_Mu50_v: passReference = True
  if not passReference: continue

  if t.AK8jets_N < 2: continue
  
  dijet = []
  j_index = []
  for j in range(t.AK8jets_N):
  
   if len(dijet) >= 2: break
   if not bool(t.AK8jets_isTight[j]): continue
   if t.AK8jets_pt[j] < 200: continue
   if math.fabs(t.AK8jets_eta[j]) > 2.4: continue
   
   tlv_j = ROOT.TLorentzVector()
   tlv_j.SetPtEtaPhiE(t.AK8jets_pt[j],t.AK8jets_eta[j],t.AK8jets_phi[j],t.AK8jets_e[j])
   
   tlv_m = 0   
   for m in range(t.muons_N):
   
    if t.muons_pt[m] < 30: continue
    if math.fabs(t.muons_eta[m]) > 2.4: continue
    if not bool(t.muons_isHighPt[m]): continue
    if t.muons_trkIso[m] > 0.1: continue
    tlv_m = ROOT.TLorentzVector()
    tlv_m.SetPtEtaPhiE(t.muons_pt[m],t.muons_eta[m],t.muons_phi[m],t.muons_e[m])
    break
  
   if tlv_m and tlv_m.DeltaR(tlv_j) < 0.8: continue
   
   dijet.append(tlv_j)
   j_index.append(j)
  
  if len(dijet) < 2: continue  
  detajj = math.fabs(dijet[0].Eta()-dijet[1].Eta())
  if detajj > 1.3: continue  
  
  mjj = (dijet[0]+dijet[1]).M()
  
  #fill mjj histos
  if t.AK8jets_softdrop_mass[j_index[0]] > 55 and t.AK8jets_softdrop_mass[j_index[1]] > 55:
    
   mjj_den.Fill(mjj)
    
   if t.HLT_PFHT1050_v: mjj_num_pfht1050.Fill(mjj) 
   if t.HLT_AK8PFJet500_v: mjj_num_pfjet500.Fill(mjj)
   if t.HLT_AK8PFJet400_TrimMass30_v: mjj_num_pfjet400_trimmass30.Fill(mjj)
   if t.HLT_AK8PFHT800_TrimMass50_v: mjj_num_pfht800_trimmass50.Fill(mjj)
   if t.HLT_PFHT1050_v or t.HLT_AK8PFJet500_v or t.HLT_AK8PFJet400_TrimMass30_v or t.HLT_AK8PFHT800_TrimMass50_v:
    mjj_num_all_triggers.Fill(mjj)
            
  #fill mjet histos     
  if mjj>1200:
  
   mjet1_den.Fill(t.AK8jets_softdrop_mass[j_index[0]])
   mjet2_den.Fill(t.AK8jets_softdrop_mass[j_index[1]])
    
   if t.HLT_PFHT1050_v:
    mjet1_num_pfht1050.Fill(t.AK8jets_softdrop_mass[j_index[0]]) 
    mjet2_num_pfht1050.Fill(t.AK8jets_softdrop_mass[j_index[1]])   
 
   if t.HLT_AK8PFJet500_v:
    mjet1_num_pfjet500.Fill(t.AK8jets_softdrop_mass[j_index[0]])
    mjet2_num_pfjet500.Fill(t.AK8jets_softdrop_mass[j_index[1]])     

   if t.HLT_AK8PFJet400_TrimMass30_v:
    mjet1_num_pfjet400_trimmass30.Fill(t.AK8jets_softdrop_mass[j_index[0]])
    mjet2_num_pfjet400_trimmass30.Fill(t.AK8jets_softdrop_mass[j_index[1]])

   if t.HLT_AK8PFHT800_TrimMass50_v:
    mjet1_num_pfht800_trimmass50.Fill(t.AK8jets_softdrop_mass[j_index[0]])
    mjet2_num_pfht800_trimmass50.Fill(t.AK8jets_softdrop_mass[j_index[1]])

   if t.HLT_PFHT1050_v or t.HLT_AK8PFJet500_v or t.HLT_AK8PFJet400_TrimMass30_v or t.HLT_AK8PFHT800_TrimMass50_v:
    mjet1_num_all_triggers.Fill(t.AK8jets_softdrop_mass[j_index[0]])
    mjet2_num_all_triggers.Fill(t.AK8jets_softdrop_mass[j_index[1]])
   

 tf.Close()
 
print "Tot. passed:",mjj_den.Integral()
print "Passed HLT_PFHT1050:",mjj_num_pfht1050.Integral()
print "Passed HLT_AK8PFJet500:",mjj_num_pfjet500.Integral()
print "Passed HLT_AK8PFJet400_TrimMass30:",mjj_num_pfjet400_trimmass30.Integral()
print "Passed HLT_AK8PFHT800_TrimMass50:",mjj_num_pfht800_trimmass50.Integral()
print "Passed all triggers:",mjj_num_all_triggers.Integral() 

#outf.cd()
#pEff = ROOT.TEfficiency(mjj_num,mjj_den)
#pEff.Draw()
#pEff.Write('mjj_eff_all_triggers')
outf.Write()
outf.Close()
#time.sleep(1000)
      