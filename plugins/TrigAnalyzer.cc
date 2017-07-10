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

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"

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

      // ----------member data ---------------------------
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigobjectToken_;
    edm::EDGetTokenT<edm::TriggerResults> trgresultsToken_;
    edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<pat::JetCollection> jetToken_;
    edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;
    edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
    edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
    edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken_ ;


    TTree* tree ;
    bool verbose_ ;
    bool isMC;
    long int EventNumber, LumiNumber, RunNumber;
    float muon1_pt, muon1_pfIso04, electron1_pt, fatjet1_pt, jet1_pt;
    float met_pt, met_pt_nomu_L, met_pt_nomu_T, m_ht, min_met_mht, min_met_mht_nomu_L, min_met_mht_nomu_T, met_phi, met_phi_nomu_L, met_phi_nomu_T;
    bool fatjet1_isloose, fatjet1_istight;
  long int nTightMuons, nTightElectrons, nTightFatJets, nLooseMuons, nLooseElectrons, nLooseFatJets, nLooseJets, nTightJets;
    bool   trig_bit_pfmet110_pfmht110 ;
    bool   trig_bit_pfmet120_pfmht120 ;
    bool   trig_bit_pfmet130_pfmht130 ;
    bool   trig_bit_pfmet140_pfmht140 ;
    bool   trig_bit_pfmetnomu110_pfmhtnomu110 ;
    bool   trig_bit_pfmetnomu120_pfmhtnomu120 ;
    bool   trig_bit_pfmetnomu130_pfmhtnomu130 ;
    bool   trig_bit_pfmetnomu140_pfmhtnomu140 ;
    bool   trig_bit_ele27_wptight_gsf ;
    bool   trig_bit_isomu24 ;

    //Owen

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

    //-- In dump of event content, this appears as
    //  vector<pat::TriggerObjectStandAlone>    "slimmedPatTrigger"        ""                "PAT"
    edm::InputTag IT_triggerObjects= edm::InputTag("slimmedPatTrigger");
    trigobjectToken_ = consumes<pat::TriggerObjectStandAloneCollection>(IT_triggerObjects);

    //-- In dump of event content, this appears as
    //  edm::TriggerResults                   "TriggerResults"            ""                "HLT
    edm::InputTag IT_hltresults = edm::InputTag("TriggerResults::HLT");
    trgresultsToken_= consumes<edm::TriggerResults>(IT_hltresults);

    //-- In dump of event content, this appears as
    //  vector<pat::MET>                      "slimmedMETs"               ""                "PAT"
    edm::InputTag IT_met = edm::InputTag("slimmedMETs") ;
    metToken_ = consumes<pat::METCollection>(IT_met) ;

    //-- In dump of event content, this appears as
    //  vector<pat::Jet>                      "slimmedJets"               ""                "PAT"
    edm::InputTag IT_jets = edm::InputTag("slimmedJets") ;
    jetToken_ = consumes<pat::JetCollection>(IT_jets) ;

    //-- In dump of event content, this appears as
    //  vector<pat::Jet>                      "slimmedJets"               ""                "PAT"
    edm::InputTag IT_fatjets = edm::InputTag("slimmedJetsAK8") ;
    fatjetToken_ = consumes<pat::JetCollection>(IT_fatjets) ;

    //-- In dump of event content, this appears as
    //  vector<reco::Vertex>                      "slimmedMuons"               ""                "RECO"
    edm::InputTag IT_vertices = edm::InputTag("offlineSlimmedPrimaryVertices") ;
    vertexToken_ = consumes<reco::VertexCollection>(IT_vertices) ;

    //-- In dump of event content, this appears as
    //  vector<pat::Muon>                      "slimmedMuons"               ""                "PAT"
    edm::InputTag IT_muons = edm::InputTag("slimmedMuons") ;
    muonToken_ = consumes<pat::MuonCollection>(IT_muons) ;

    //-- In dump of event content, this appears as
    //  vector<pat::Electron>                      "slimmedElectrons"               ""                "PAT"
    edm::InputTag IT_electrons = edm::InputTag("slimmedElectrons") ;
    electronToken_ = consumes<pat::ElectronCollection>(IT_electrons) ;
    
    //-- In dump of event content, this appears as
    //  vector<pat::PackedCandidate>          "packedPFCandidates"        ""                "RECO"
    edm::InputTag IT_pfcands = edm::InputTag("packedPFCandidates") ;
    pfcandsToken_ = consumes<pat::PackedCandidateCollection>(IT_pfcands) ;

    verbose_ = iConfig.getParameter<bool> ("verbose") ;

   //now do what ever initialization is needed
    usesResource("TFileService");

    edm::Service<TFileService> fs ;
    tree = fs->make<TTree>( "tree", "tree" ) ;
    tree -> Branch( "isMC" , &isMC, "isMC/O");
    tree -> Branch( "EventNumber" , &EventNumber , "EventNumber/L");
    tree -> Branch( "LumiNumber" , &LumiNumber , "LumiNumber/L");
    tree -> Branch( "RunNumber" , &RunNumber , "RunNumber/L");
    tree -> Branch( "nLooseMuons" , &nLooseMuons , "nLooseMuons/L");
    tree -> Branch( "nLooseElectrons" , &nLooseElectrons , "nLooseElectrons/L");
    tree -> Branch( "nLooseFatJets" , &nLooseFatJets , "nLooseFatJets/L");
    tree -> Branch( "nLooseJets" , &nLooseJets , "nLooseJets/L");
    tree -> Branch( "nTightMuons" , &nTightMuons , "nTightMuons/L");
    tree -> Branch( "nTightElectrons" , &nTightElectrons , "nTightElectrons/L");
    tree -> Branch( "nTightFatJets" , &nTightFatJets , "nTightFatJets/L");
    tree -> Branch( "nTightJets" , &nTightJets , "nTightJets/L");
    tree -> Branch("Muon1_pt", &muon1_pt, "Muon1_pt/F");
    tree -> Branch("Muon1_pfIso04", &muon1_pfIso04, "Muon1_pfIso04/F");
    tree -> Branch("Electron1_pt", &electron1_pt, "Electron1_pt/F");
    tree -> Branch("FatJet1_pt", &fatjet1_pt, "FatJet1_pt/F");
    tree -> Branch("Jet1_pt", &jet1_pt, "Jet1_pt/F");
    tree -> Branch("MEt_pt", &met_pt, "MEt_pt/F");
    tree -> Branch("MEt_phi", &met_phi, "MEt_phi/F");
    tree -> Branch("m_ht", &m_ht, "m_ht/F");
    tree -> Branch("min_met_mht", &min_met_mht, "min_met_mht/F");
    tree -> Branch("met_pt_nomu_L", &met_pt_nomu_L, "met_pt_nomu_L/F");
    tree -> Branch("met_pt_nomu_T", &met_pt_nomu_T, "met_pt_nomu_T/F");
    tree -> Branch("min_met_mht_nomu_L", &min_met_mht_nomu_L, "min_met_mht_nomu_L/F");
    tree -> Branch("min_met_mht_nomu_T", &min_met_mht_nomu_T, "min_met_mht_nomu_T/F");
    tree -> Branch( "HLT_PFMET110_PFMHT110_IDTight_v", &trig_bit_pfmet110_pfmht110, "HLT_PFMET110_PFMHT110_IDTight_v/B" ) ;
    tree -> Branch( "HLT_PFMET120_PFMHT120_IDTight_v", &trig_bit_pfmet120_pfmht120, "HLT_PFMET120_PFMHT120_IDTight_v/B" ) ;
    tree -> Branch( "HLT_PFMET130_PFMHT130_IDTight_v", &trig_bit_pfmet130_pfmht130, "HLT_PFMET130_PFMHT130_IDTight_v/B" ) ;
    tree -> Branch( "HLT_PFMET140_PFMHT140_IDTight_v", &trig_bit_pfmet140_pfmht140, "HLT_PFMET140_PFMHT140_IDTight_v/B" ) ;
    tree -> Branch( "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v", &trig_bit_pfmetnomu110_pfmhtnomu110, "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v/B" ) ;
    tree -> Branch( "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v", &trig_bit_pfmetnomu120_pfmhtnomu120, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v/B" ) ;
    tree -> Branch( "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v", &trig_bit_pfmetnomu130_pfmhtnomu130, "HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v/B" ) ;
    tree -> Branch( "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v", &trig_bit_pfmetnomu140_pfmhtnomu140, "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v/B" ) ;
    tree -> Branch( "HLT_Ele27_WPTight_Gsf_v", &trig_bit_ele27_wptight_gsf, "HLT_Ele27_WPTight_Gsf_v/B" ) ;
    tree -> Branch( "HLT_IsoMu24_v", &trig_bit_isomu24, "HLT_IsoMu24_v/B" ) ;
    //Owen

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

    //Owen

    isMC = false;
    EventNumber = LumiNumber = RunNumber = 0;
    nTightMuons = nTightElectrons = nTightFatJets = nLooseMuons = nLooseElectrons = nLooseFatJets = 0;

    trig_bit_pfmet110_pfmht110 = false ;
    trig_bit_pfmet120_pfmht120 = false ;
    trig_bit_pfmet130_pfmht130 = false ;
    trig_bit_pfmet140_pfmht140 = false ;
    trig_bit_pfmetnomu110_pfmhtnomu110 = false ;
    trig_bit_pfmetnomu120_pfmhtnomu120 = false ;
    trig_bit_pfmetnomu130_pfmhtnomu130 = false ;
    trig_bit_pfmetnomu140_pfmhtnomu140 = false ;
    trig_bit_ele27_wptight_gsf = false ;
    trig_bit_isomu24 = false ;

    muon1_pt = 0.;
    muon1_pfIso04 = -1.;
    electron1_pt = 0.;
    met_pt = met_pt_nomu_L = met_pt_nomu_T = m_ht = min_met_mht = min_met_mht_nomu_L = min_met_mht_nomu_T = 0.;
    met_phi = met_phi_nomu_L = met_phi_nomu_T = -10.;

    //bool dump_event(false) ;
    //bool dump_event(true) ;


    //Accessing trigger bits (same as AOD):
    edm::Handle<edm::TriggerResults> trigResults;
    iEvent.getByToken(trgresultsToken_, trigResults);

    if( !trigResults.failedToGet() ) {
        int N_Triggers = trigResults->size();
        const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
        // put unpackPathNames( trigName ) here ?


        for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
            if (trigResults.product()->accept(i_Trig)) {
              TString TrigPath =trigName.triggerName(i_Trig);

              //if ( TrigPath.Contains("HLT_PFHT1050_v") ) trig_bit_hlt_pfht1050 = true ;
              //if ( TrigPath.Contains("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v") ) trig_bit_hlt_pfht500_pfmet100_pfmht100 = true ;
              //if ( TrigPath.Contains("HLT_PFHT500_PFMET110_PFMHT110_IDTight_v") ) trig_bit_hlt_pfht500_pfmet110_pfmht110 = true ;
              //if ( TrigPath.Contains("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v") ) trig_bit_hlt_pfht700_pfmet85_pfmht85 = true ;
              //if ( TrigPath.Contains("HLT_PFHT700_PFMET95_PFMHT95_IDTight_v") ) trig_bit_hlt_pfht700_pfmet95_pfmht95 = true ;
              //if ( TrigPath.Contains("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v") ) trig_bit_hlt_pfht800_pfmet75_pfmht75 = true ;
              //if ( TrigPath.Contains("HLT_PFHT800_PFMET85_PFMHT85_IDTight_v") ) trig_bit_hlt_pfht800_pfmet85_pfmht85 = true ;

              if ( TrigPath.Contains("HLT_PFMET110_PFMHT110_IDTight_v") ) trig_bit_pfmet110_pfmht110 = true ;
              if ( TrigPath.Contains("HLT_PFMET120_PFMHT120_IDTight_v") ) trig_bit_pfmet120_pfmht120 = true ;
              if ( TrigPath.Contains("HLT_PFMET130_PFMHT130_IDTight_v") ) trig_bit_pfmet130_pfmht130 = true ;
              if ( TrigPath.Contains("HLT_PFMET140_PFMHT140_IDTight_v") ) trig_bit_pfmet140_pfmht140 = true ;

              if ( TrigPath.Contains("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v") ) trig_bit_pfmetnomu110_pfmhtnomu110 = true ;
              if ( TrigPath.Contains("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") ) trig_bit_pfmetnomu120_pfmhtnomu120 = true ;
              if ( TrigPath.Contains("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v") ) trig_bit_pfmetnomu130_pfmhtnomu130 = true ;
              if ( TrigPath.Contains("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v") ) trig_bit_pfmetnomu140_pfmhtnomu140 = true ;

              //if ( TrigPath.Contains("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight_v") ) trig_bit_hlt_monocentralpfjet80_pfmetnomu120_pfmhtnomu120 = true ;

              if ( TrigPath.Contains("HLT_Ele27_WPTight_Gsf_v") ) trig_bit_ele27_wptight_gsf = true ;
              if ( TrigPath.Contains("HLT_IsoMu24_v") ) trig_bit_isomu24 = true ;
           }
        }
    }


    //Event info
    isMC = !iEvent.isRealData();
    EventNumber = iEvent.id().event();
    LumiNumber = iEvent.luminosityBlock();
    RunNumber = iEvent.id().run();

    //Vertices
    edm::Handle<reco::VertexCollection> VertexColl;
    iEvent.getByToken( vertexToken_, VertexColl);
    const reco::Vertex* vertex=&VertexColl->front();

    //Loop on MET
    edm::Handle<pat::METCollection> MetColl;
    iEvent.getByToken( metToken_, MetColl) ;
    pat::MET met = MetColl->front() ;
    met_pt = met.pt();
    met_phi = met.phi();

    //Loop on AK8 jets
    edm::Handle<pat::JetCollection> fatjets;
    iEvent.getByToken( fatjetToken_, fatjets );
    std::vector<pat::Jet> FatJetVect;

    for(std::vector<pat::Jet>::const_iterator it=fatjets->begin(); it!=fatjets->end(); it++) {
        pat::Jet f=*it;
	if ( !isLooseJet(f) ) continue;
        if ( f.pt() < 170 ) continue;
        if ( fabs( f.eta() ) > 2.5 ) continue;
	nLooseFatJets++;
	if ( !isTightJet(f) ) continue;
	nTightFatJets++;
	fatjet1_pt = f.pt();
	FatJetVect.push_back(f);
    }

    //Loop on AK4 jets
    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken( jetToken_, jets );
    std::vector<pat::Jet> JetVect;

    for(std::vector<pat::Jet>::const_iterator it=jets->begin(); it!=jets->end(); it++) {
        pat::Jet j=*it;
	if ( !isLooseJet(j) ) continue;
        if ( j.pt() < 30 ) continue;
        if ( fabs( j.eta() ) > 2.5 ) continue;
        nLooseJets++;
	jet1_pt = j.pt();
	JetVect.push_back(j);
    }

    float m_ht_x(0.), m_ht_y(0.);
    
    for(unsigned int a=0; a<JetVect.size(); a++){
        m_ht_x += JetVect.at(a).px();
        m_ht_y += JetVect.at(a).py();
    }
    
    m_ht = sqrt( pow(m_ht_x,2) + pow(m_ht_y,2)  );
    min_met_mht = std::min(met_pt,m_ht);

    //Loop on muons ---> fix
    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken( muonToken_, muons );
    std::vector<pat::Muon> MuonVect;

    float met_pt_nomu_x_L(0.), met_pt_nomu_y_L(0.), met_pt_nomu_x_T(0.), met_pt_nomu_y_T(0.);

    for ( const pat::Muon &m : *muons) {
        if ( m.pt() < 30 ) continue;
        if ( fabs( m.eta() ) > 2.4 ) continue;
	if (!m.isLooseMuon()) continue;
	float pfIso04 = (m.pfIsolationR04().sumChargedHadronPt + std::max(m.pfIsolationR04().sumNeutralHadronEt + m.pfIsolationR04().sumPhotonEt - 0.5*m.pfIsolationR04().sumPUPt, 0.) ) / m.pt();
        if (pfIso04>0.25) continue; //at least loose isolation
        met_pt_nomu_x_L += m.px();
        met_pt_nomu_y_L += m.py();
	nLooseMuons++;
        //muon1_isLoose = true;
	if (!m.isTightMuon(*vertex)) continue;
        met_pt_nomu_x_T += m.px();
        met_pt_nomu_y_T += m.py();
	nTightMuons++;
        //muon1_isTight = true;
	muon1_pt = m.pt();
        muon1_pfIso04 = pfIso04;
	MuonVect.push_back(m);
    } // loop over muons, saving only tight muons

    met_pt_nomu_L = sqrt( pow(met_pt_nomu_x_L,2) + pow(met_pt_nomu_y_L,2) );
    met_pt_nomu_T = sqrt( pow(met_pt_nomu_x_T,2) + pow(met_pt_nomu_y_T,2) );
    min_met_mht_nomu_L = std::min(met_pt_nomu_L,m_ht);
    min_met_mht_nomu_T = std::min(met_pt_nomu_T,m_ht);


    //Loop on electrons ---> fix, missing Electron IDs
    /*
    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken( electronToken_, electrons );
    std::vector<pat::Electron> ElectronVect;
    for ( const pat::Electron &e : *electrons) {
        if ( e.pt() < 30 ) continue;
        if ( fabs( e.eta() ) > 2.5 ) continue;
	std::cout << "Electron pt: " << e.pt() << std::endl;
	electron1_pt = e.pt();
	ElectronVect.push_back(e);
    } // loop over muons
    */
    tree -> Fill() ;




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

//define this as a plug-in
DEFINE_FWK_MODULE(TrigAnalyzer);
