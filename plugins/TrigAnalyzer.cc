// -*- C++ -*-
//
// Package:    Trigger2017/TrigAnalyzer
// Class:      TrigAnalyzer
// 
/**\class TrigAnalyzer TrigAnalyzer.cc Trigger2017/TrigAnalyzer/plugins/TrigAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Lisa Benato
//         Created:  Sat, 08 Jul 2017 16:26:05 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include <string>
//#include "TrigAnalyzer.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TrigAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TrigAnalyzer(const edm::ParameterSet&);
      ~TrigAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      virtual bool isLooseJet(pat::Jet&);
      virtual bool isTightJet(pat::Jet&);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      bool passIDWP(std::string, bool, float, float, float, float, float, float, float, bool, int);
      void reset(void);

      // ----------member data ---------------------------
    edm::EDGetTokenT<edm::TriggerResults> trigResultsToken;
    edm::EDGetTokenT<edm::TriggerResults> filterResultsToken;
    edm::EDGetTokenT<trigger::TriggerEvent> trigobjectsRAWToken_;//for getting L1 bit
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;//for getting HLT in miniaod??
    edm::EDGetTokenT<bool> badChCandFilterToken;
    edm::EDGetTokenT<bool> badPFMuonFilterToken;
    edm::EDGetTokenT<pat::METCollection> metToken;
    edm::EDGetTokenT<pat::JetCollection> jetToken;
    edm::EDGetTokenT<pat::JetCollection> fatjetToken;
    edm::EDGetTokenT<pat::MuonCollection> muonToken;
    edm::EDGetTokenT<reco::VertexCollection> vertexToken;
    edm::EDGetTokenT<edm::View<pat::Electron> > electronToken;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleVetoIdMapToken;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleLooseIdMapToken;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleHEEPIdMapToken;

    std::vector<std::string> filterNames;

    TTree* tree;
    
    //global variables
    bool isVerbose;
    bool isMC;
    int EventNumber, LumiNumber, RunNumber, nPV;
    
    //muons
    int muons_N;
    std::vector<float> muons_pt, muons_eta, muons_phi, muons_e;
    std::vector<float> muons_pfIso04, muons_trkIso;
    std::vector<int> muons_isLoose, muons_isTight, muons_isHighPt;
    
    //electrons
    int eles_N;
    std::vector<float> eles_pt, eles_eta, eles_phi, eles_e;
    std::vector<int> eles_isVeto, eles_isLoose, eles_isMedium, eles_isTight, eles_isHeep;
    
    //ak4 jets
    float jet1_pt;
    int nLooseJets, nTightJets;
    
    //ak8 jets
    int AK8jets_N;
    std::vector<float> AK8jets_pt, AK8jets_eta, AK8jets_phi, AK8jets_e, AK8jets_m;
    std::vector<float> AK8jets_softdrop_mass, AK8jets_tau1, AK8jets_tau2, AK8jets_tau3;
    std::vector<int> AK8jets_isLoose, AK8jets_isTight;
    
    //HT (AK4 CHS tight jets with pT > 30 GeV and |eta| < 3.0)
    float HT;
    
    //met/mht
    float met_pt, met_pt_nomu_L, met_pt_nomu_T, m_ht, m_ht_nomu_L, m_ht_nomu_T, min_met_mht, min_met_mht_nomu_L, min_met_mht_nomu_T, met_phi, met_phi_nomu_L, met_phi_nomu_T;

    //trigger bits
    int   trig_bit_pfmet110_pfmht110;
    int   trig_bit_pfmet120_pfmht120;
    int   trig_bit_pfmet120_pfmht120_PFHT60;
    int   trig_bit_pfmet130_pfmht130;
    int   trig_bit_pfmet140_pfmht140;
    int   trig_bit_pfmetTypeOne110_pfmht110;
    int   trig_bit_pfmetTypeOne120_pfmht120;
    int   trig_bit_pfmetTypeOne120_pfmht120_PFHT60;
    int   trig_bit_pfmetTypeOne130_pfmht130;
    int   trig_bit_pfmetTypeOne140_pfmht140;
    int   trig_bit_pfmetnomu110_pfmhtnomu110;
    int   trig_bit_pfmetnomu120_pfmhtnomu120;
    int   trig_bit_pfmetnomu120_pfmhtnomu120_PFHT60;
    int   trig_bit_pfmetnomu130_pfmhtnomu130;
    int   trig_bit_pfmetnomu140_pfmhtnomu140;
    int   trig_bit_ele27_wptight_gsf;
    int   trig_bit_isomu24;
    int   trig_bit_isomu27;
    int   trig_bit_mu50;    
    int   trig_bit_pfht1050;
    int   trig_bit_ak8pfjet400;
    int   trig_bit_ak8pfjet450;
    int   trig_bit_ak8pfjet500;
    int   trig_bit_ak8pfjet550;
    int   trig_bit_ak8pfht750_trimmass50;
    int   trig_bit_ak8pfht800_trimmass50;
    int   trig_bit_ak8pfht850_trimmass50;
    int   trig_bit_ak8pfht900_trimmass50;
    int   trig_bit_ak8pfjet360_trimmass30;
    int   trig_bit_ak8pfjet380_trimmass30;
    int   trig_bit_ak8pfjet400_trimmass30;
    int   trig_bit_ak8pfjet420_trimmass30;        
    //L1 bits
    int trig_bit_hltMHT90;
    int trig_bit_hltL1sAllETMHFSeeds;
    int trig_bit_hltL1sAllETMHadSeeds;
    int trig_bit_hltMETClean80;
    int trig_bit_hltMET90;
    int trig_bit_hltPFMHTNoMuTightID120;
    int trig_bit_hltPFMETNoMu120;
    int trig_bit_hltL1sAllETMHFHTT60Seeds;
    int trig_bit_hltPFHTJet30;
    int trig_bit_hltPFHT60Jet30;
    //MET filters
    int trig_bit_flag_HBHENoiseFilter;
    int trig_bit_flag_HBHENoiseIsoFilter;
    int trig_bit_flag_EcalDeadCellTriggerPrimitiveFilter;
    int trig_bit_flag_goodVertices;
    int trig_bit_flag_eeBadScFilter;
    int trig_bit_flag_globalSuperTightHalo2016Filter;
    int flag_BadChCand;
    int flag_BadPFMuon;

    //muon trigger objects
    float mu_trigObj_pt_IsoMu24, mu_trigObj_eta_IsoMu24, mu_trigObj_phi_IsoMu24, mu_trigObj_e_IsoMu24;
    float mu_trigObj_pt_IsoMu27, mu_trigObj_eta_IsoMu27, mu_trigObj_phi_IsoMu27, mu_trigObj_e_IsoMu27;
    float mu_trigObj_pt_Mu50, mu_trigObj_eta_Mu50, mu_trigObj_phi_Mu50, mu_trigObj_e_Mu50;
    
};


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrigAnalyzer::TrigAnalyzer(const edm::ParameterSet& iConfig)

{

    filterNames.push_back("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07"); //IsoMu24
    filterNames.push_back("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07"); //IsoMu27
    filterNames.push_back("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"); //Mu50
    
    //Input tags
    edm::InputTag IT_trigResults = edm::InputTag("TriggerResults::HLT");
    trigResultsToken= consumes<edm::TriggerResults>(IT_trigResults);
    edm::InputTag IT_filterResults = edm::InputTag("TriggerResults::RECO");
    filterResultsToken= consumes<edm::TriggerResults>(IT_filterResults);

    edm::InputTag IT_badChCandFilter = edm::InputTag("BadChargedCandidateFilter");
    badChCandFilterToken= consumes<bool>(IT_badChCandFilter);
    edm::InputTag IT_badPFMuonFilter = edm::InputTag("BadPFMuonFilter");
    badPFMuonFilterToken= consumes<bool>(IT_badPFMuonFilter);

    edm::InputTag IT_met = edm::InputTag("slimmedMETs");
    metToken = consumes<pat::METCollection>(IT_met);
    edm::InputTag IT_jets = edm::InputTag("slimmedJets");
    jetToken = consumes<pat::JetCollection>(IT_jets);
    edm::InputTag IT_fatjets = edm::InputTag("slimmedJetsAK8");
    fatjetToken = consumes<pat::JetCollection>(IT_fatjets);
    edm::InputTag IT_vertices = edm::InputTag("offlineSlimmedPrimaryVertices");
    vertexToken = consumes<reco::VertexCollection>(IT_vertices);
    edm::InputTag IT_muons = edm::InputTag("slimmedMuons");
    muonToken = consumes<pat::MuonCollection>(IT_muons);
    edm::InputTag IT_electrons = edm::InputTag("slimmedElectrons");
    electronToken = consumes<edm::View<pat::Electron> >(IT_electrons);   
    
    edm::InputTag IT_trigObj = edm::InputTag("slimmedPatTrigger");
    triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(IT_trigObj);
;
    
    /*edm::InputTag IT_eleVetoId = edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-veto");
    eleVetoIdMapToken = consumes<edm::ValueMap<bool> >(IT_eleVetoId);
    edm::InputTag IT_eleLooseId = edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-loose");
    eleLooseIdMapToken = consumes<edm::ValueMap<bool> >(IT_eleLooseId);
    edm::InputTag IT_eleMediumId = edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-medium");
    eleMediumIdMapToken = consumes<edm::ValueMap<bool> >(IT_eleMediumId);
    edm::InputTag IT_eleTightId = edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight");
    eleTightIdMapToken = consumes<edm::ValueMap<bool> >(IT_eleTightId);
    edm::InputTag IT_eleHeepId = edm::InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70");
    eleHEEPIdMapToken = consumes<edm::ValueMap<bool> >(IT_eleHeepId);*/
    
    isVerbose = iConfig.getParameter<bool> ("verbose");

    //now do what ever initialization is needed
    usesResource("TFileService");

    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree", "tree");
    
    //global variables
    tree -> Branch("isMC" , &isMC, "isMC/O");
    tree -> Branch("EventNumber" , &EventNumber , "EventNumber/I");
    tree -> Branch("LumiNumber" , &LumiNumber , "LumiNumber/I");
    tree -> Branch("RunNumber" , &RunNumber , "RunNumber/I");
    tree -> Branch("nPV" , &nPV , "nPV/I");
    
    //muons
    tree -> Branch("muons_N", &muons_N, "muons_N/I");
    tree -> Branch("muons_pt", &muons_pt);
    tree -> Branch("muons_eta", &muons_eta);
    tree -> Branch("muons_phi", &muons_phi);
    tree -> Branch("muons_e", &muons_e);
    tree -> Branch("muons_pfIso04", &muons_pfIso04);
    tree -> Branch("muons_trkIso", &muons_trkIso);
    tree -> Branch("muons_isLoose", &muons_isLoose);
    tree -> Branch("muons_isTight", &muons_isTight);
    tree -> Branch("muons_isHighPt", &muons_isHighPt);
    
    //electrons
    tree -> Branch("eles_N", &eles_N, "eles_N/I");
    tree -> Branch("eles_pt", &eles_pt);
    tree -> Branch("eles_eta", &eles_eta);
    tree -> Branch("eles_phi", &eles_phi);
    tree -> Branch("eles_e", &eles_e);
    tree -> Branch("eles_isVeto", &eles_isVeto);
    tree -> Branch("eles_isLoose", &eles_isLoose);
    tree -> Branch("eles_isMedium", &eles_isMedium);
    tree -> Branch("eles_isTight", &eles_isTight);
    tree -> Branch("eles_isHeep", &eles_isHeep);
        
    //ak8jets
    tree -> Branch("AK8jets_N", &AK8jets_N, "AK8jets_N/I");
    tree -> Branch("AK8jets_isLoose", &AK8jets_isLoose);
    tree -> Branch("AK8jets_isTight", &AK8jets_isTight);
    tree -> Branch("AK8jets_pt", &AK8jets_pt);
    tree -> Branch("AK8jets_eta", &AK8jets_eta);
    tree -> Branch("AK8jets_phi", &AK8jets_phi);
    tree -> Branch("AK8jets_m", &AK8jets_m);
    tree -> Branch("AK8jets_e", &AK8jets_e);
    tree -> Branch("AK8jets_softdrop_mass", &AK8jets_softdrop_mass);
    tree -> Branch("AK8jets_tau1", &AK8jets_tau1);
    tree -> Branch("AK8jets_tau2", &AK8jets_tau2);
    tree -> Branch("AK8jets_tau3", &AK8jets_tau3);
    
    //HT (AK4 CHS tight jets with pT > 30 GeV and |eta| < 3.0)
    tree -> Branch("HT", &HT, "HT/F");
    
    //met/mht
    tree -> Branch("MEt_pt", &met_pt, "MEt_pt/F");
    tree -> Branch("MEt_phi", &met_phi, "MEt_phi/F");    
    tree -> Branch("m_ht", &m_ht, "m_ht/F");
    tree -> Branch("m_ht_nomu_L", &m_ht_nomu_L, "m_ht_nomu_L/F");
    tree -> Branch("m_ht_nomu_T", &m_ht_nomu_T, "m_ht_nomu_T/F");
    tree -> Branch("min_met_mht", &min_met_mht, "min_met_mht/F");
    tree -> Branch("met_pt_nomu_L", &met_pt_nomu_L, "met_pt_nomu_L/F");
    tree -> Branch("met_pt_nomu_T", &met_pt_nomu_T, "met_pt_nomu_T/F");
    tree -> Branch("min_met_mht_nomu_L", &min_met_mht_nomu_L, "min_met_mht_nomu_L/F");
    tree -> Branch("min_met_mht_nomu_T", &min_met_mht_nomu_T, "min_met_mht_nomu_T/F");

    //ak4 jets
    tree -> Branch("Jet1_pt", &jet1_pt, "Jet1_pt/F");
    tree -> Branch("nTightJets" , &nTightJets , "nTightJets/I");
    //tree -> Branch("nLooseJets" , &nLooseJets , "nLooseJets/L");
    
    //trigger bits    
    tree -> Branch("HLT_PFMET110_PFMHT110_IDTight_v", &trig_bit_pfmet110_pfmht110, "HLT_PFMET110_PFMHT110_IDTight_v/I");
    tree -> Branch("HLT_PFMET120_PFMHT120_IDTight_v", &trig_bit_pfmet120_pfmht120, "HLT_PFMET120_PFMHT120_IDTight_v/I");
    tree -> Branch("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v", &trig_bit_pfmet120_pfmht120_PFHT60, "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v/I");
    tree -> Branch("HLT_PFMET130_PFMHT130_IDTight_v", &trig_bit_pfmet130_pfmht130, "HLT_PFMET130_PFMHT130_IDTight_v/I");
    tree -> Branch("HLT_PFMET140_PFMHT140_IDTight_v", &trig_bit_pfmet140_pfmht140, "HLT_PFMET140_PFMHT140_IDTight_v/I");
    tree -> Branch("HLT_PFMETTypeOne110_PFMHT110_IDTight_v", &trig_bit_pfmetTypeOne110_pfmht110, "HLT_PFMETTypeOne110_PFMHT110_IDTight_v/I");
    tree -> Branch("HLT_PFMETTypeOne120_PFMHT120_IDTight_v", &trig_bit_pfmetTypeOne120_pfmht120, "HLT_PFMETTypeOne120_PFMHT120_IDTight_v/I");
    tree -> Branch("HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v", &trig_bit_pfmetTypeOne120_pfmht120_PFHT60, "HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v/I");
    tree -> Branch("HLT_PFMETTypeOne130_PFMHT130_IDTight_v", &trig_bit_pfmetTypeOne130_pfmht130, "HLT_PFMETTypeOne130_PFMHT130_IDTight_v/I");
    tree -> Branch("HLT_PFMETTypeOne140_PFMHT140_IDTight_v", &trig_bit_pfmetTypeOne140_pfmht140, "HLT_PFMETTypeOne140_PFMHT140_IDTight_v/I");
    tree -> Branch("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v", &trig_bit_pfmetnomu110_pfmhtnomu110, "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v/I");
    tree -> Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v", &trig_bit_pfmetnomu120_pfmhtnomu120, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v/I");
    tree -> Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v", &trig_bit_pfmetnomu120_pfmhtnomu120_PFHT60, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v/I");
    tree -> Branch("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v", &trig_bit_pfmetnomu130_pfmhtnomu130, "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v/I");
    tree -> Branch("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v", &trig_bit_pfmetnomu140_pfmhtnomu140, "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v/I");
    tree -> Branch("HLT_Ele27_WPTight_Gsf_v", &trig_bit_ele27_wptight_gsf, "HLT_Ele27_WPTight_Gsf_v/I");
    tree -> Branch("HLT_IsoMu24_v", &trig_bit_isomu24, "HLT_IsoMu24_v/I");
    tree -> Branch("HLT_IsoMu27_v", &trig_bit_isomu27, "HLT_IsoMu27_v/I");
    tree -> Branch("HLT_Mu50_v", &trig_bit_mu50, "HLT_Mu50_v/I");
    tree -> Branch("HLT_PFHT1050_v", &trig_bit_pfht1050, "HLT_PFHT1050_v/I");
    tree -> Branch("HLT_AK8PFJet400_v", &trig_bit_ak8pfjet400, "HLT_AK8PFJet400_v/I");
    tree -> Branch("HLT_AK8PFJet450_v", &trig_bit_ak8pfjet450, "HLT_AK8PFJet450_v/I");
    tree -> Branch("HLT_AK8PFJet500_v", &trig_bit_ak8pfjet500, "HLT_AK8PFJet500_v/I");
    tree -> Branch("HLT_AK8PFJet550_v", &trig_bit_ak8pfjet400, "HLT_AK8PFJet550_v/I");
    tree -> Branch("HLT_AK8PFHT750_TrimMass50_v", &trig_bit_ak8pfht750_trimmass50, "HLT_AK8PFHT750_TrimMass50_v/I");
    tree -> Branch("HLT_AK8PFHT800_TrimMass50_v", &trig_bit_ak8pfht800_trimmass50, "HLT_AK8PFHT800_TrimMass50_v/I");
    tree -> Branch("HLT_AK8PFHT850_TrimMass50_v", &trig_bit_ak8pfht850_trimmass50, "HLT_AK8PFHT850_TrimMass50_v/I");
    tree -> Branch("HLT_AK8PFHT900_TrimMass50_v", &trig_bit_ak8pfht900_trimmass50, "HLT_AK8PFHT9000_TrimMass50_v/I");
    tree -> Branch("HLT_AK8PFJet360_TrimMass30_v", &trig_bit_ak8pfjet360_trimmass30, "HLT_AK8PFJet360_TrimMass30_v/I");
    tree -> Branch("HLT_AK8PFJet380_TrimMass30_v", &trig_bit_ak8pfjet380_trimmass30, "HLT_AK8PFJet380_TrimMass30_v/I");
    tree -> Branch("HLT_AK8PFJet400_TrimMass30_v", &trig_bit_ak8pfjet400_trimmass30, "HLT_AK8PFJet400_TrimMass30_v/I");
    tree -> Branch("HLT_AK8PFJet420_TrimMass30_v", &trig_bit_ak8pfjet420_trimmass30, "HLT_AK8PFJet420_TrimMass30_v/I");
    
    //met filters
    tree -> Branch("Flag_HBHENoiseFilter", &trig_bit_flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/I");
    tree -> Branch("Flag_HBHENoiseIsoFilter", &trig_bit_flag_HBHENoiseIsoFilter, "Flag_HBHENoiseIsoFilter/I");
    tree -> Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &trig_bit_flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/I");
    tree -> Branch("Flag_goodVertices", &trig_bit_flag_goodVertices, "Flag_goodVertices/I");
    tree -> Branch("Flag_eeBadScFilter", &trig_bit_flag_eeBadScFilter, "Flag_eeBadScFilter/I");
    tree -> Branch("Flag_globalSuperTightHalo2016Filter", &trig_bit_flag_globalSuperTightHalo2016Filter, "Flag_globalSuperTightHalo2016Filter/I");
    tree -> Branch("Flag_BadChCand", &flag_BadChCand, "Flag_BadChCand/I");
    tree -> Branch("Flag_BadPFMuon", &flag_BadPFMuon, "Flag_BadPFMuon/I");
    
    //muon trigger objects
    tree -> Branch("mu_trigObj_pt_IsoMu24", &mu_trigObj_pt_IsoMu24, "mu_trigObj_pt_IsoMu24/F"); 
    tree -> Branch("mu_trigObj_eta_IsoMu24", &mu_trigObj_eta_IsoMu24, "mu_trigObj_eta_IsoMu24/F"); 
    tree -> Branch("mu_trigObj_phi_IsoMu24", &mu_trigObj_phi_IsoMu24, "mu_trigObj_phi_IsoMu24/F"); 
    tree -> Branch("mu_trigObj_e_IsoMu24", &mu_trigObj_e_IsoMu24, "mu_trigObj_e_IsoMu24/F"); 
    tree -> Branch("mu_trigObj_pt_IsoMu27", &mu_trigObj_pt_IsoMu27, "mu_trigObj_pt_IsoMu27/F"); 
    tree -> Branch("mu_trigObj_eta_IsoMu27", &mu_trigObj_eta_IsoMu27, "mu_trigObj_eta_IsoMu27/F"); 
    tree -> Branch("mu_trigObj_phi_IsoMu27", &mu_trigObj_phi_IsoMu27, "mu_trigObj_phi_IsoMu27/F"); 
    tree -> Branch("mu_trigObj_e_IsoMu27", &mu_trigObj_e_IsoMu27, "mu_trigObj_e_IsoMu27/F"); 
    tree -> Branch("mu_trigObj_pt_Mu50", &mu_trigObj_pt_Mu50, "mu_trigObj_pt_Mu50/F"); 
    tree -> Branch("mu_trigObj_eta_Mu50", &mu_trigObj_eta_Mu50, "mu_trigObj_eta_Mu50/F"); 
    tree -> Branch("mu_trigObj_phi_Mu50", &mu_trigObj_phi_Mu50, "mu_trigObj_phi_Mu50/F"); 
    tree -> Branch("mu_trigObj_e_Mu50", &mu_trigObj_e_Mu50, "mu_trigObj_e_Mu50/F");         

}


TrigAnalyzer::~TrigAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrigAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace reco;
    using namespace std;

    reset();

    //Accessing trigger bits (same as AOD); thanks to Owen Long (SUSY)
    edm::Handle<edm::TriggerResults> trigResults;
    iEvent.getByToken(trigResultsToken, trigResults);

    if( !trigResults.failedToGet() ) {
        int N_Triggers = trigResults->size();
        const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);

        for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
            if (trigResults.product()->accept(i_Trig)) {
              TString TrigPath =trigName.triggerName(i_Trig);

              if ( TrigPath.Contains("HLT_PFMET110_PFMHT110_IDTight_v") ) trig_bit_pfmet110_pfmht110 = true;
              if ( TrigPath.Contains("HLT_PFMET120_PFMHT120_IDTight_v") ) trig_bit_pfmet120_pfmht120 = true;
              if ( TrigPath.Contains("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v") ) trig_bit_pfmet120_pfmht120_PFHT60 = true;
              if ( TrigPath.Contains("HLT_PFMET130_PFMHT130_IDTight_v") ) trig_bit_pfmet130_pfmht130 = true;
              if ( TrigPath.Contains("HLT_PFMET140_PFMHT140_IDTight_v") ) trig_bit_pfmet140_pfmht140 = true;

              if ( TrigPath.Contains("HLT_PFMETTypeOne110_PFMHT110_IDTight_v") ) trig_bit_pfmetTypeOne110_pfmht110 = true;
              if ( TrigPath.Contains("HLT_PFMETTypeOne120_PFMHT120_IDTight_v") ) trig_bit_pfmetTypeOne120_pfmht120 = true;
              if ( TrigPath.Contains("HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v") ) trig_bit_pfmetTypeOne120_pfmht120_PFHT60 = true;
              if ( TrigPath.Contains("HLT_PFMETTypeOne130_PFMHT130_IDTight_v") ) trig_bit_pfmetTypeOne130_pfmht130 = true;
              if ( TrigPath.Contains("HLT_PFMETTypeOne140_PFMHT140_IDTight_v") ) trig_bit_pfmetTypeOne140_pfmht140 = true;

              if ( TrigPath.Contains("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v") ) {
		  trig_bit_pfmetnomu110_pfmhtnomu110 = true;
		  std::cout << TrigPath << std::endl;
	      }
              if ( TrigPath.Contains("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") ) trig_bit_pfmetnomu120_pfmhtnomu120 = true;
              if ( TrigPath.Contains("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v") ) trig_bit_pfmetnomu120_pfmhtnomu120_PFHT60 = true;
              if ( TrigPath.Contains("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v") ) trig_bit_pfmetnomu130_pfmhtnomu130 = true;
              if ( TrigPath.Contains("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v") ) trig_bit_pfmetnomu140_pfmhtnomu140 = true;

              if ( TrigPath.Contains("HLT_Ele27_WPTight_Gsf_v") ) trig_bit_ele27_wptight_gsf = true;
              if ( TrigPath.Contains("HLT_IsoMu24_v") ) trig_bit_isomu24 = true;
              if ( TrigPath.Contains("HLT_IsoMu27_v") ) trig_bit_isomu27 = true;
              if ( TrigPath.Contains("HLT_Mu50_v") ) trig_bit_mu50 = true;
 
              if ( TrigPath.Contains("HLT_PFHT1050_v") ) trig_bit_pfht1050 = true;             
              if ( TrigPath.Contains("HLT_AK8PFJet400_v") ) trig_bit_ak8pfjet400 = true;
              if ( TrigPath.Contains("HLT_AK8PFJet450_v") ) trig_bit_ak8pfjet450 = true;
              if ( TrigPath.Contains("HLT_AK8PFJet500_v") ) trig_bit_ak8pfjet500 = true;
              if ( TrigPath.Contains("HLT_AK8PFJet550_v") ) trig_bit_ak8pfjet550 = true;
              if ( TrigPath.Contains("HLT_AK8PFHT750_TrimMass50_v") ) trig_bit_ak8pfht750_trimmass50 = true;
              if ( TrigPath.Contains("HLT_AK8PFHT800_TrimMass50_v") ) trig_bit_ak8pfht800_trimmass50 = true;
              if ( TrigPath.Contains("HLT_AK8PFHT850_TrimMass50_v") ) trig_bit_ak8pfht850_trimmass50 = true;
              if ( TrigPath.Contains("HLT_AK8PFHT900_TrimMass50_v") ) trig_bit_ak8pfht900_trimmass50 = true;
              if ( TrigPath.Contains("HLT_AK8PFJet360_TrimMass30_v") ) trig_bit_ak8pfjet360_trimmass30 = true;
              if ( TrigPath.Contains("HLT_AK8PFJet380_TrimMass30_v") ) trig_bit_ak8pfjet380_trimmass30 = true;
              if ( TrigPath.Contains("HLT_AK8PFJet400_TrimMass30_v") ) trig_bit_ak8pfjet400_trimmass30 = true;
              if ( TrigPath.Contains("HLT_AK8PFJet420_TrimMass30_v") ) trig_bit_ak8pfjet420_trimmass30 = true;
    
           }
        }
    }

    //trigger objects
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    iEvent.getByToken(triggerObjects_,triggerObjects);
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);
    
    if( trig_bit_isomu24 || trig_bit_isomu27 || trig_bit_mu50 ){
     //std::cout << "************* PASSED ISOMU24 " << trig_bit_isomu24 << " OR ISOMU27 " << trig_bit_isomu27 << " OR MU50 " << trig_bit_mu50 << std::endl;
     
     for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    
      obj.unpackPathNames(trigNames);
      obj.unpackFilterLabels(iEvent, *trigResults);
      std::vector<std::string> pathNamesAll  = obj.pathNames(false);
     
      if( trig_bit_isomu24 ){
       bool pathExist = false;
       for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h){
        if( std::string(pathNamesAll[h]).find("HLT_IsoMu24_v") != std::string::npos ){ /*std::cout << " FOUND ISOMU24 " << pathNamesAll[h] << std::endl;*/ pathExist = true; break; }
       }//close loop on pathNames;
       if( pathExist ){
        //bool filterExist = false;
        for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){
	 if( filterNames[0] == obj.filterLabels()[hh] ){
	  //std::cout << " ALSO FOUND FILTER FOR ISOMU24 " << obj.filterLabels()[hh] << std::endl;
	  mu_trigObj_pt_IsoMu24 = obj.pt();
	  mu_trigObj_eta_IsoMu24 = obj.eta();
	  mu_trigObj_phi_IsoMu24 = obj.phi();
	  mu_trigObj_e_IsoMu24 = obj.energy();
	  //filterExist = true;
	  break;
	 }
	}//close loop on filter labels of the object
	//if( filterExist ) break;
       }//close found path for IsoMu24
      }//close if Isomu24  
      


      if( trig_bit_isomu27 ){
       bool pathExist = false;
       for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h){
        if( std::string(pathNamesAll[h]).find("HLT_IsoMu27_v") != std::string::npos ){ /*std::cout << " FOUND ISOMU27 " << pathNamesAll[h] << std::endl;*/ pathExist = true; break; }
       }//close loop on pathNames;
       if( pathExist ){
        //bool filterExist = false;
        for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){
	 if( filterNames[1] == obj.filterLabels()[hh] ){
	  //std::cout << " ALSO FOUND FILTER FOR ISOMU27 " << obj.filterLabels()[hh] << std::endl;
	  mu_trigObj_pt_IsoMu27 = obj.pt();
	  mu_trigObj_eta_IsoMu27 = obj.eta();
	  mu_trigObj_phi_IsoMu27 = obj.phi();
	  mu_trigObj_e_IsoMu27 = obj.energy();
	  //filterExist = true;
	  break;
	 }
	}//close loop on filter labels of the object
	//if( filterExist ) break;
       }//close found path for IsoMu27
      }//close if Isomu27  
            


      if( trig_bit_mu50 ){
       bool pathExist = false;
       for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h){
        if( std::string(pathNamesAll[h]).find("HLT_Mu50_v") != std::string::npos ){ /*std::cout << " FOUND Mu50 " << pathNamesAll[h] << std::endl;*/ pathExist = true; break; }
       }//close loop on pathNames;
       if( pathExist ){
        //bool filterExist = false;
        for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){
	 if( filterNames[2] == obj.filterLabels()[hh] ){
	  //std::cout << " ALSO FOUND FILTER FOR Mu50 " << obj.filterLabels()[hh] << std::endl;
	  mu_trigObj_pt_Mu50 = obj.pt();
	  mu_trigObj_eta_Mu50 = obj.eta();
	  mu_trigObj_phi_Mu50 = obj.phi();
	  mu_trigObj_e_Mu50 = obj.energy();
	  //filterExist = true;
	  break;
	 }
	}//close loop on filter labels of the object
	//if( filterExist ) break;
       }//close found path for Mu50
      }//close if Mu50  
               
     }//end loop on objects     
    }// if one of the muon trigger paths fired --> to save time
           
    //MET filters
    edm::Handle<edm::TriggerResults> filterResults; 
    iEvent.getByToken(filterResultsToken, filterResults);

    if( !filterResults.failedToGet() ) { 
        int N_Filters = filterResults->size();
        const edm::TriggerNames & filterName = iEvent.triggerNames(*filterResults);

        for( int i_Trig = 0; i_Trig < N_Filters; ++i_Trig ) { 
	    if (filterResults.product()->accept(i_Trig)) {
	        TString TrigPath =filterName.triggerName(i_Trig);

	        if ( TrigPath.Contains("Flag_HBHENoiseFilter") ) trig_bit_flag_HBHENoiseFilter = true;
	        if ( TrigPath.Contains("Flag_HBHENoiseIsoFilter") ) trig_bit_flag_HBHENoiseIsoFilter = true;
	        if ( TrigPath.Contains("Flag_EcalDeadCellTriggerPrimitiveFilter") ) trig_bit_flag_EcalDeadCellTriggerPrimitiveFilter = true;
	        if ( TrigPath.Contains("Flag_goodVertices") ) trig_bit_flag_goodVertices = true;
	        if ( TrigPath.Contains("Flag_eeBadScFilter") ) trig_bit_flag_eeBadScFilter = true;
	        if ( TrigPath.Contains("Flag_globalSuperTightHalo2016Filter") ) trig_bit_flag_globalSuperTightHalo2016Filter = true;
	    }
        }
    }

    //BadChCand and BadPFMuon filters
    edm::Handle<bool> filterBadChCand; 
    iEvent.getByToken(badChCandFilterToken, filterBadChCand);
    flag_BadChCand = *filterBadChCand;

    edm::Handle<bool> filterBadPFMuon; 
    iEvent.getByToken(badPFMuonFilterToken, filterBadPFMuon);
    flag_BadPFMuon = *filterBadPFMuon;

    //Event info
    isMC = !iEvent.isRealData();
    EventNumber = iEvent.id().event();
    LumiNumber = iEvent.luminosityBlock();
    RunNumber = iEvent.id().run();
    
    /*std::cout << "============================= DONE WITH EVENT " << EventNumber << " LUMI " << LumiNumber << " RUN " << RunNumber << std::endl;
    std::cout << "CHECK ISOMU24 " << mu_trigObj_pt_IsoMu24 << " " <<  mu_trigObj_eta_IsoMu24 << " " <<  mu_trigObj_phi_IsoMu24 << " " <<  mu_trigObj_e_IsoMu24 << std::endl;
    std::cout << "CHECK ISOMU27 " << mu_trigObj_pt_IsoMu27 << " " <<  mu_trigObj_eta_IsoMu27 << " " <<  mu_trigObj_phi_IsoMu27 << " " <<  mu_trigObj_e_IsoMu27 << std::endl;
    std::cout << "CHECK MU50 " << mu_trigObj_pt_Mu50 << " " <<  mu_trigObj_eta_Mu50 << " " <<  mu_trigObj_phi_Mu50 << " " <<  mu_trigObj_e_Mu50 << std::endl;*/
   
    //std::cout << " " << mu_trigObj_pt << " " << mu_trigObj_eta << " " << mu_trigObj_phi << " " << mu_trigObj_e << std::endl;

    //Initialize met no mu
    float met_pt_nomu_x_L(0.), met_pt_nomu_y_L(0.), met_pt_nomu_x_T(0.), met_pt_nomu_y_T(0.);

    //Vertices
    edm::Handle<reco::VertexCollection> VertexColl;
    iEvent.getByToken( vertexToken, VertexColl);
    nPV = VertexColl->size();
    const reco::Vertex* vertex=&VertexColl->front();
    reco::TrackBase::Point vtxPoint(0,0,0);
    if(  VertexColl->size() >= 1 ) {
        vtxPoint = VertexColl->at(0).position();
    }

    /////////////////////////////////////////////////////////////////
    //Fill muons collection
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken( muonToken, muons );
    
    for(std::vector<pat::Muon>::const_iterator it=muons->begin(); it!=muons->end(); it++) {
        pat::Muon m=*it;
        //if ( m.pt() < 30 ) continue; //this causes a jump at ~30 GeV, investigate
        if ( fabs( m.eta() ) > 2.4 ) continue; //this selection is necessary
        muons_N++;
        muons_pt.push_back(m.pt());
        muons_eta.push_back(m.eta());
        muons_phi.push_back(m.phi());
        muons_e.push_back(m.energy());
        muons_isLoose.push_back(m.isLooseMuon());
        muons_isTight.push_back(m.isTightMuon(*vertex));
        muons_isHighPt.push_back(m.isTightMuon(*vertex));
	    float pfIso04 = (m.pfIsolationR04().sumChargedHadronPt + std::max(m.pfIsolationR04().sumNeutralHadronEt + m.pfIsolationR04().sumPhotonEt - 0.5*m.pfIsolationR04().sumPUPt, 0.) ) / m.pt();
        muons_pfIso04.push_back(pfIso04);
        muons_trkIso.push_back(m.trackIso());
    }

    /////////////////////////////////////////////////////////////////
    //Fill electrons collection 
    /*edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
    edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
    edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
    edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
    edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
    
    iEvent.getByToken(eleVetoIdMapToken  , veto_id_decisions  );
    iEvent.getByToken(eleLooseIdMapToken , loose_id_decisions );
    iEvent.getByToken(eleMediumIdMapToken, medium_id_decisions);
    iEvent.getByToken(eleTightIdMapToken , tight_id_decisions );
    iEvent.getByToken(eleHEEPIdMapToken  , heep_id_decisions  );*/

    edm::Handle<edm::View<pat::Electron> > electrons;
    iEvent.getByToken( electronToken, electrons );
    
    for (const pat::Electron &e : *electrons) {

        if ( e.pt() < 30 ) continue;
        if ( fabs( e.eta() ) > 2.5 ) continue;
        
        //const auto el = electrons->ptrAt(eles_N);
        
        eles_N++;
        eles_pt.push_back(e.pt());
        eles_eta.push_back(e.eta());
        eles_phi.push_back(e.phi());
        eles_e.push_back(e.energy());
        /*eles_isVeto.push_back((*veto_id_decisions)[el]);
        eles_isLoose.push_back((*loose_id_decisions)[el]);
        eles_isMedium.push_back((*medium_id_decisions)[el]);
        eles_isTight.push_back((*tight_id_decisions)[el]);
        eles_isHeep.push_back((*heep_id_decisions)[el]);*/
    }
               
    /////////////////////////////////////////////////////////////////        
    //Fill AK8 jets collection
    edm::Handle<pat::JetCollection> fatjets;
    iEvent.getByToken( fatjetToken, fatjets );
    
    for(std::vector<pat::Jet>::const_iterator it=fatjets->begin(); it!=fatjets->end(); it++) {
        pat::Jet f=*it;
        if ( f.pt() < 170 ) continue;
        if ( fabs( f.eta() ) > 2.5 ) continue;
        AK8jets_N++;
        AK8jets_isLoose.push_back(isLooseJet(f));
        AK8jets_isTight.push_back(isTightJet(f));
	    AK8jets_pt.push_back(f.pt());
	    AK8jets_eta.push_back(f.eta());
	    AK8jets_phi.push_back(f.phi());
	    AK8jets_m.push_back(f.mass());
	    AK8jets_e.push_back(f.energy());
	    AK8jets_softdrop_mass.push_back(f.userFloat("ak8PFJetsPuppiSoftDropMass"));
	    AK8jets_tau1.push_back(f.userFloat("NjettinessAK8Puppi:tau1"));
	    AK8jets_tau2.push_back(f.userFloat("NjettinessAK8Puppi:tau2"));
	    AK8jets_tau3.push_back(f.userFloat("NjettinessAK8Puppi:tau3"));
    }

    /////////////////////////////////////////////////////////////////
    //Fill met/mht variables
    edm::Handle<pat::METCollection> MetColl;
    iEvent.getByToken( metToken, MetColl);
    pat::MET met = MetColl->front();
    met_pt = met.pt();
    met_phi = met.phi();
    met_pt_nomu_x_L = met_pt_nomu_x_T = met.px();//before summing up muons
    met_pt_nomu_y_L = met_pt_nomu_y_T = met.py();//before summing up muons

    //Loop on AK4 jets
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken( jetToken, jets );
    std::vector<pat::Jet> JetVect;

    for(std::vector<pat::Jet>::const_iterator it=jets->begin(); it!=jets->end(); it++) {
        pat::Jet j=*it;
        //std::cout << "jet info" << std::endl; // AK4 are always PF jets!!
	//std::cout << j.isPFJet() << std::endl;//if this works...
	//std::cout << j.hasPFSpecific() << std::endl;//if this works...
	//if ( !isLooseJet(j) ) continue;
        if ( !isTightJet(j) ) continue;
        if ( j.pt() < 30 ) continue;//this causes a jump at ~30? investigate!
        if ( fabs( j.eta() ) < 3.0 ) HT+=j.pt();
        if ( fabs( j.eta() ) > 2.5 ) continue;
        //nLooseJets++;
        nTightJets++;
	jet1_pt = j.pt();
	JetVect.push_back(j);
    }

    float m_ht_x(0.), m_ht_y(0.), m_ht_nomu_x_L(0.), m_ht_nomu_y_L(0.), m_ht_nomu_x_T(0.), m_ht_nomu_y_T(0.);
    
    for(unsigned int a=0; a<JetVect.size(); a++){
        m_ht_x -= JetVect.at(a).px();
        m_ht_y -= JetVect.at(a).py();
        m_ht_nomu_x_L -= JetVect.at(a).px();
        m_ht_nomu_y_L -= JetVect.at(a).py();
        m_ht_nomu_x_T -= JetVect.at(a).px();
        m_ht_nomu_y_T -= JetVect.at(a).py();
    }
    
    m_ht = sqrt( pow(m_ht_x,2) + pow(m_ht_y,2)  );
    min_met_mht = std::min(met_pt,m_ht);

    //for ( const pat::Muon &m : *muons) {
    for(std::vector<pat::Muon>::const_iterator it=muons->begin(); it!=muons->end(); it++) {
        pat::Muon m=*it;
        //if ( m.pt() < 30 ) continue; //this causes a jump at ~30 GeV, investigate
        if ( fabs( m.eta() ) > 2.4 ) continue; //this selection is necessary
        if(m.isLooseMuon()){
         met_pt_nomu_x_L += m.px();
         met_pt_nomu_y_L += m.py();
	     m_ht_nomu_x_L += m.px();
	     m_ht_nomu_y_L += m.py();
	    }
	    if (m.isTightMuon(*vertex)){
         met_pt_nomu_x_T += m.px();
         met_pt_nomu_y_T += m.py();
	     m_ht_nomu_x_T += m.px();
	     m_ht_nomu_y_T += m.py();
	    }

	//MuonVect.push_back(m);
    } // loop over muons, saving only tight muons

    met_pt_nomu_L = sqrt( pow(met_pt_nomu_x_L,2) + pow(met_pt_nomu_y_L,2) );
    met_pt_nomu_T = sqrt( pow(met_pt_nomu_x_T,2) + pow(met_pt_nomu_y_T,2) );
    m_ht_nomu_L = sqrt( pow(m_ht_nomu_x_L,2) + pow(m_ht_nomu_y_L,2)  );
    m_ht_nomu_T = sqrt( pow(m_ht_nomu_x_T,2) + pow(m_ht_nomu_y_T,2)  );
    min_met_mht_nomu_L = std::min(met_pt_nomu_L,m_ht_nomu_L);
    min_met_mht_nomu_T = std::min(met_pt_nomu_T,m_ht_nomu_T);


    //Loop on electrons ---> fix, missing Electron IDs
    /*
    edm::InputTag convLabel = edm::InputTag("reducedEgamma:reducedConversions");
    edm::Handle<reco::ConversionCollection> conversions;
    iEvent.getByLabel( convLabel, conversions );

    //edm::InputTag bsLabel = edm::InputTag("offlineBeamSpot");
    //edm::Handle<reco::BeamSpot> bsHandle;
    //iEvent.getByLabel( bsLabel, bsHandle );
    //const reco::BeamSpot &beamspot = *bsHandle.product();

    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken( electronToken, electrons );
    //std::vector<pat::Electron> ElectronVect;
    for(std::vector<pat::Electron>::const_iterator it=electrons->begin(); it!=electrons->end(); it++) {
        pat::Electron e=*it;
	** //
	GsfElectron::PflowIsolationVariables pfIso = e.pfIsolationVariables();
	bool isEB = e.isEB() ? true : false;
	float dEtaIn = e.deltaEtaSuperClusterTrackAtVtx();
	float dPhiIn = e.deltaPhiSuperClusterTrackAtVtx();
	float full5x5 = e.full5x5_sigmaIetaIeta();
	float hoe = e.hadronicOverEm();
	float absiso = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
	float relIsoWithDBeta_ = absiso/e.pt();
	float ooEmooP_; 
	if( e.ecalEnergy() == 0 ){
	    printf("Electron energy is zero!\n");
	    ooEmooP_ = 999;
	}
        else if( !std::isfinite(e.ecalEnergy())){
	    printf("Electron energy is not finite!\n");
	    ooEmooP_ = 998;
	}
        else{
	    ooEmooP_ = fabs(1.0/e.ecalEnergy() - e.eSuperClusterOverP()/e.ecalEnergy() );
	}
	std::cout << ooEmooP_ << std::endl;
	float d0 = (-1) * e.gsfTrack()->dxy(vtxPoint);
	float dz = e.gsfTrack()->dz(vtxPoint);
	float missHits = e.gsfTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS);
	bool hasMatchConv = ConversionTools::hasMatchedConversion(e, conversions, beamspot.position());
	bool isVeto = passIDWP("VETO",isEB, dEtaIn, dPhiIn, full5x5, hoe, d0, dz, ooEmooP_, hasMatchConv, missHits);
	bool isLoose = passIDWP("LOOSE",isEB, dEtaIn, dPhiIn, full5x5, hoe, d0, dz, ooEmooP_, hasMatchConv, missHits);
	bool isMedium = passIDWP("MEDIUM",isEB, dEtaIn, dPhiIn, full5x5, hoe, d0, dz, ooEmooP_, hasMatchConv, missHits);
	bool isTight = passIDWP("TIGHT",isEB, dEtaIn, dPhiIn, full5x5, hoe, d0, dz, ooEmooP_, hasMatchConv, missHits);
	// Look up the ID decision for this electron in 	
	if ( isLoose ) {
  	    std::cout << "LOOSE ele" << std::endl;
	}
        /**
        if ( e.pt() < 30 ) continue;
        if ( fabs( e.eta() ) > 2.5 ) continue;
	//std::cout << "Electron pt: " << e.pt() << std::endl;
	electron1_pt = e.pt();
	electron1_eta = e.eta();
	//ElectronVect.push_back(e);
    } // loop over electrons
    */
    
    tree -> Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
TrigAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrigAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrigAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//Method to define loose jet ID (2017 data)
bool TrigAnalyzer::isLooseJet(pat::Jet& jet){
    if(fabs(jet.eta())<=2.7){/// |eta| < 2.7
        if(jet.neutralHadronEnergyFraction()>=0.99) return false;
        if(jet.neutralEmEnergyFraction()>=0.99) return false;
        if((jet.chargedMultiplicity()+jet.neutralMultiplicity())<=1) return false;
            if(fabs(jet.eta())<=2.4) { /// |eta| < 2.4
                if(jet.chargedHadronEnergyFraction()<=0.) return false;
                if(jet.chargedMultiplicity()<=0) return false;
                if(jet.chargedEmEnergyFraction()>=0.99) return false;
            }
    }
    else{ /// |eta| > 2.7
        if(jet.neutralEmEnergyFraction()>=0.90) return false;
        if (fabs(jet.eta())<=3.0) { /// 2.7 < |eta| < 3.0
            if(jet.neutralMultiplicity()<=2) return false;
        }
        else{ /// |eta| > 3.0
           if(jet.neutralMultiplicity()<=10) return false;
        }
    }

    return true;
}

//Method to define tight jet ID (2017 data)
bool TrigAnalyzer::isTightJet(pat::Jet& jet){
    if(fabs(jet.eta())<=2.7){/// |eta| < 2.7
        if(jet.neutralHadronEnergyFraction()>=0.90) return false;
        if(jet.neutralEmEnergyFraction()>=0.90) return false;
        if((jet.chargedMultiplicity()+jet.neutralMultiplicity())<=1) return false;
            if(fabs(jet.eta())<=2.4) { /// |eta| < 2.4
                if(jet.chargedHadronEnergyFraction()<=0.) return false;
                if(jet.chargedMultiplicity()<=0) return false;
                if(jet.chargedEmEnergyFraction()>=0.99) return false;
            }
    }
    else{ /// |eta| > 2.7
        if (fabs(jet.eta())<=3.0) { /// 2.7 < |eta| < 3.0
            if(jet.neutralMultiplicity()<=2) return false;
            if(jet.neutralEmEnergyFraction()>=0.99) return false;
            if(jet.neutralEmEnergyFraction()<=0.02) return false;
        }
        else{ /// |eta| > 3.0
           if(jet.neutralMultiplicity()<=10) return false;
           if(jet.neutralHadronEnergyFraction()<=0.02) return false;
           if(jet.neutralEmEnergyFraction()>=0.90) return false;
        }
    }

    return true;
}


bool TrigAnalyzer::passIDWP(std::string WP, bool isEB, float dEtaIn, float dPhiIn, float full5x5, float hoe, float d0, float dz, float ooemoop, bool conv, int missHits){
  bool pass = false;

  if(WP == "VETO"){
    if(isEB){
      pass = (fabs(dEtaIn) <  0.0126 ) && (fabs(dPhiIn) <  0.107 ) && (full5x5 < 0.012 ) && (hoe <  0.186 ) && (fabs(d0) < 0.0621 ) && (fabs(dz) <  0.613 ) && (fabs(ooemoop) <  0.239 ) && !conv && (missHits <= 2);
    }
    else{
      pass = (fabs(dEtaIn) <  0.0109 ) && (fabs(dPhiIn) <  0.219 ) && (full5x5 < 0.0339 ) && (hoe <  0.0962 ) && (fabs(d0) < 0.279 ) && (fabs(dz) < 0.947 ) && (fabs(ooemoop) < 0.141 ) && !conv && (missHits <= 3);
    }
  }
  if(WP == "LOOSE"){
    if(isEB){
      pass = (fabs(dEtaIn) < 0.00976 ) && (fabs(dPhiIn) < 0.0929 ) && (full5x5 <  0.0105 ) && (hoe < 0.0765 ) && (fabs(d0) < 0.0227 ) && (fabs(dz) < 0.379 ) && (fabs(ooemoop) <  0.184 ) && !conv && (missHits <= 2);
    }
    else{
      pass = (fabs(dEtaIn) < 0.00952 ) && (fabs(dPhiIn) < 0.181 ) && (full5x5 < 0.0318 ) && (hoe < 0.0824 ) && (fabs(d0) < 0.242 ) && (fabs(dz) < 0.921 ) && (fabs(ooemoop) < 0.125 ) && !conv && (missHits <= 1);
    }
  }

  if(WP == "MEDIUM"){
    if(isEB){
      pass = (fabs(dEtaIn) <  0.0094 ) && (fabs(dPhiIn) <  0.0296 ) && (full5x5 <  0.0101 ) && (hoe <  0.0372 ) && (fabs(d0) <  0.0151 ) && (fabs(dz) <  0.238 ) && (fabs(ooemoop) <  0.118 ) && !conv && (missHits <= 2);
    }
    else{
      pass = (fabs(dEtaIn) <  0.00773 ) && (fabs(dPhiIn) <  0.148 ) && (full5x5 <  0.0287 ) && (hoe <  0.0546 ) && (fabs(d0) <  0.0535 ) && (fabs(dz) <  0.572 ) && (fabs(ooemoop) <  0.104 ) && !conv && (missHits <= 1);
    }
  }

  if(WP == "TIGHT"){
    if(isEB){
      pass = (fabs(dEtaIn) <  0.0095 ) && (fabs(dPhiIn) <  0.0291 ) && (full5x5 <  0.0101 ) && (hoe <  0.0372 ) && (fabs(d0) <  0.0144 ) && (fabs(dz) <  0.323 ) && (fabs(ooemoop) <  0.0174 ) && !conv && (missHits <= 2);
    }
    else{
      pass = (fabs(dEtaIn) <  0.00762 ) && (fabs(dPhiIn) <  0.0439 ) && (full5x5 <  0.0287 ) && (hoe <  0.0544 ) && (fabs(d0) <  0.0377 ) && (fabs(dz) <  0.571 ) && (fabs(ooemoop) <  0.01 ) && !conv && (missHits <= 1);
    }
  }
  return pass;
}

////////////////////////////reset all variables before filling tree/////////////////////////////////////////
void TrigAnalyzer::reset(void){

    //global variables
    isMC = false;
    EventNumber = LumiNumber = RunNumber = nPV = 0;
    
    //muons
    muons_N = 0;
    muons_pt.clear();
    muons_eta.clear();
    muons_phi.clear();
    muons_e.clear();
    muons_pfIso04.clear();
    muons_trkIso.clear();
    muons_isLoose.clear();
    muons_isTight.clear();
    muons_isHighPt.clear();
    
    //electrons
    eles_N = 0;
    eles_pt.clear();
    eles_eta.clear();
    eles_phi.clear();
    eles_e.clear();
    eles_isVeto.clear();
    eles_isLoose.clear();
    eles_isMedium.clear();
    eles_isTight.clear();
    eles_isHeep.clear();
    
    //ak4 jets
    jet1_pt = -9999;
    nLooseJets = nTightJets = 0;
    
    //ak8 jets
    AK8jets_N = 0;
    AK8jets_isLoose.clear();
    AK8jets_isTight.clear();
    AK8jets_pt.clear();
    AK8jets_eta.clear();
    AK8jets_phi.clear();
    AK8jets_e.clear();
    AK8jets_m.clear();
    AK8jets_softdrop_mass.clear();
    AK8jets_tau1.clear();
    AK8jets_tau2.clear();
    AK8jets_tau3.clear();
    
    //HT (AK4 CHS tight jets with pT > 30 GeV and |eta| < 3.0)
    HT = 0;
    
    //met/mht
    met_pt = met_pt_nomu_L = met_pt_nomu_T = m_ht = m_ht_nomu_L = m_ht_nomu_T = min_met_mht = min_met_mht_nomu_L = -9999;
    min_met_mht_nomu_T = met_phi = met_phi_nomu_L = met_phi_nomu_T = -9999;

    //trigger bits
    trig_bit_pfmet110_pfmht110 = false;
    trig_bit_pfmet120_pfmht120 = false;
    trig_bit_pfmet120_pfmht120_PFHT60 = false;
    trig_bit_pfmet130_pfmht130 = false;
    trig_bit_pfmet140_pfmht140 = false;
    trig_bit_pfmetTypeOne110_pfmht110 = false;
    trig_bit_pfmetTypeOne120_pfmht120 = false;
    trig_bit_pfmetTypeOne120_pfmht120_PFHT60 = false;
    trig_bit_pfmetTypeOne130_pfmht130 = false;
    trig_bit_pfmetTypeOne140_pfmht140 = false;
    trig_bit_pfmetnomu110_pfmhtnomu110 = false;
    trig_bit_pfmetnomu120_pfmhtnomu120 = false;
    trig_bit_pfmetnomu120_pfmhtnomu120_PFHT60 = false;
    trig_bit_pfmetnomu130_pfmhtnomu130 = false;
    trig_bit_pfmetnomu140_pfmhtnomu140 = false;
    trig_bit_ele27_wptight_gsf = false;
    trig_bit_isomu24 = false;
    trig_bit_isomu27 = false;
    trig_bit_mu50 = false;    
    trig_bit_pfht1050 = false;
    trig_bit_ak8pfjet400 = false;
    trig_bit_ak8pfjet450 = false;
    trig_bit_ak8pfjet500 = false;
    trig_bit_ak8pfjet550 = false;
    trig_bit_ak8pfht750_trimmass50 = false;
    trig_bit_ak8pfht800_trimmass50 = false;
    trig_bit_ak8pfht850_trimmass50 = false;
    trig_bit_ak8pfht900_trimmass50 = false;
    trig_bit_ak8pfjet360_trimmass30 = false;
    trig_bit_ak8pfjet380_trimmass30 = false;
    trig_bit_ak8pfjet400_trimmass30 = false;
    trig_bit_ak8pfjet420_trimmass30 = false;        
    //L1 bits
    trig_bit_hltMHT90 = false;
    trig_bit_hltL1sAllETMHFSeeds = false;
    trig_bit_hltL1sAllETMHadSeeds = false;
    trig_bit_hltMETClean80 = false;
    trig_bit_hltMET90 = false;
    trig_bit_hltPFMHTNoMuTightID120 = false;
    trig_bit_hltPFMETNoMu120 = false;
    trig_bit_hltL1sAllETMHFHTT60Seeds = false;
    trig_bit_hltPFHTJet30 = false;
    trig_bit_hltPFHT60Jet30 = false;
    //MET filters
    trig_bit_flag_HBHENoiseFilter = false;
    trig_bit_flag_HBHENoiseIsoFilter = false;
    trig_bit_flag_EcalDeadCellTriggerPrimitiveFilter = false;
    trig_bit_flag_goodVertices = false;
    trig_bit_flag_eeBadScFilter = false;
    trig_bit_flag_globalSuperTightHalo2016Filter = false;
    flag_BadChCand = false;
    flag_BadPFMuon = false;

    //muon trigger objects
    mu_trigObj_pt_IsoMu24 = mu_trigObj_eta_IsoMu24 = mu_trigObj_phi_IsoMu24 = mu_trigObj_e_IsoMu24 = -9999;
    mu_trigObj_pt_IsoMu27 = mu_trigObj_eta_IsoMu27 = mu_trigObj_phi_IsoMu27 = mu_trigObj_e_IsoMu27 = -9999;
    mu_trigObj_pt_Mu50 = mu_trigObj_eta_Mu50 = mu_trigObj_phi_Mu50 = mu_trigObj_e_Mu50 = -9999;
            
}


//define this as a plug-in
DEFINE_FWK_MODULE(TrigAnalyzer);
