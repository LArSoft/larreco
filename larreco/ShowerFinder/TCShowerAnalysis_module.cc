// -------------------------------------------------
// tcshower analysis tree
//
// Author: Rory Fitzpatrick (roryfitz@umich.edu)
// Created: 7/16/18
// -------------------------------------------------

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/MCBase/MCDataHolder.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <TROOT.h>
#include <TStyle.h>

const int kMaxShowers = 1000;  //maximum number of showers   

namespace shower {

  class TCShowerAnalysis : public art::EDAnalyzer {

  public:

    explicit TCShowerAnalysis(fhicl::ParameterSet const& pset);
    virtual ~TCShowerAnalysis();

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void analyze(const art::Event& evt);
   
  protected:

    void reset();

    void truthMatcher(std::vector<art::Ptr<recob::Hit>>all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet);

  private:
    
    TTree* fTree;
    int run;
    int subrun;
    int event;
    int nuPDG_truth;
    int ccnc_truth;
    int mode_truth;
    int nshws;                         //number of showers
    int shwid[kMaxShowers];             //recob::Shower::ID()
    float shwdcosx[kMaxShowers];      //shower direction cosin                                                        
    float shwdcosy[kMaxShowers];
    float shwdcosz[kMaxShowers];
    float shwstartx[kMaxShowers];     //shower start position (cm)
    float shwstarty[kMaxShowers];
    float shwstartz[kMaxShowers];
    double shwdedx[kMaxShowers][2];    //shower dE/dx of the initial track measured on the 3 plane (MeV/cm)
    int shwbestplane[kMaxShowers];      //recommended plane for energy and dE/dx information

    int highestHitsPDG;
    double highestHitsFrac;

    std::string fHitModuleLabel;
    std::string fShowerModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fDigitModuleLabel;

    calo::CalorimetryAlg fCalorimetryAlg;

  }; // class TCShowerAnalysis

} // shower

// -------------------------------------------------

shower::TCShowerAnalysis::TCShowerAnalysis(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset),
  fHitModuleLabel           (pset.get< std::string >("HitModuleLabel", "trajcluster" ) ),
  fShowerModuleLabel        (pset.get< std::string >("ShowerModuleLabel", "tcshower" ) ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel", "generator")  ),
  fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel", "largeant")  ),
  fCalorimetryAlg           (pset.get< fhicl::ParameterSet >("CalorimetryAlg") ) {
  this->reconfigure(pset);
} // TCShowerTemplateMaker

// -------------------------------------------------

shower::TCShowerAnalysis::~TCShowerAnalysis() {
} // ~TCShowerTemplateMaker

// -------------------------------------------------

void shower::TCShowerAnalysis::reconfigure(fhicl::ParameterSet const& pset) {
} // reconfigure

// -------------------------------------------------

void shower::TCShowerAnalysis::beginJob() {

  art::ServiceHandle<art::TFileService> tfs;

  fTree = tfs->make<TTree>("tcshowerana", "tcshowerana");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");

  fTree->Branch("nuPDG_truth",&nuPDG_truth,"nuPDG_truth/I");
  fTree->Branch("ccnc_truth",&ccnc_truth,"ccnc_truth/I");
  fTree->Branch("mode_truth",&mode_truth,"mode_truth/I");

  fTree->Branch("nshws",&nshws,"nshws/I");
  fTree->Branch("shwid",shwid,"shwid[nshws]/I");
  fTree->Branch("shwdcosx",shwdcosx,"shwdcosx[nshws]/F");
  fTree->Branch("shwdcosy",shwdcosy,"shwdcosy[nshws]/F");
  fTree->Branch("shwdcosz",shwdcosz,"shwdcosz[nshws]/F");
  fTree->Branch("shwstartx",shwstartx,"shwstartx[nshws]/F");
  fTree->Branch("shwstarty",shwstarty,"shwstarty[nshws]/F");
  fTree->Branch("shwstartz",shwstartz,"shwstartz[nshws]/F");
  fTree->Branch("shwdedx",shwdedx,"shwdedx[nshws][2]/D");
  fTree->Branch("shwbestplane",shwbestplane,"shwbestplane[nshws]/I");
  
  fTree->Branch("highestHitsPDG",&highestHitsPDG,"highestHitsPDG/I");
  fTree->Branch("highestHitsFrac",&highestHitsFrac,"highestHitsFrac/D");

} // beginJob

// -------------------------------------------------

void shower::TCShowerAnalysis::analyze(const art::Event& evt) {

  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  art::Handle< std::vector<sim::SimChannel> > scListHandle;
  std::vector<art::Ptr<sim::SimChannel> > simchanlist;
  if (evt.getByLabel(fDigitModuleLabel,scListHandle))
    art::fill_ptr_vector(simchanlist, scListHandle);

  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower> > showerlist;
  if (evt.getByLabel(fShowerModuleLabel,showerListHandle))
    art::fill_ptr_vector(showerlist, showerListHandle);
  
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);

  art::FindManyP<recob::Hit> shwfm(showerListHandle, evt, fShowerModuleLabel);

  //shower information                                                                                                     
  if (showerListHandle.isValid()) {
    nshws = showerlist.size();

    for (int i = 0; i < std::min(int(showerlist.size()),kMaxShowers); ++i) {
      shwid[i] = showerlist[i]->ID();
      shwdcosx[i] = showerlist[i]->Direction().X();
      shwdcosy[i] = showerlist[i]->Direction().Y();
      shwdcosz[i] = showerlist[i]->Direction().Z();
      shwstartx[i] = showerlist[i]->ShowerStart().X();
      shwstarty[i] = showerlist[i]->ShowerStart().Y();
      shwstartz[i] = showerlist[i]->ShowerStart().Z();
      for (size_t j = 0; j<(showerlist[i]->dEdx()).size(); ++j){
        shwdedx[i][j] = showerlist[i]->dEdx()[j];
      }
      shwbestplane[i] = showerlist[i]->best_plane();
    }

  } // shower info                                                                                                         


  if (mclist.size()) {
    art::Ptr<simb::MCTruth> mctruth = mclist[0];
    if (mctruth->NeutrinoSet()) {

      nuPDG_truth = mctruth->GetNeutrino().Nu().PdgCode();
      ccnc_truth = mctruth->GetNeutrino().CCNC();
      mode_truth = mctruth->GetNeutrino().Mode();

      if (showerlist.size()) { // only looks at the first shower since this is for tcshower
	std::vector< art::Ptr<recob::Hit> > showerhits = shwfm.at(0);
	// get shower truth info
	const simb::MCParticle *particle;
	double tmpEfrac = 0.0;
	double tmpEcomplet = 0;
	truthMatcher(hitlist, showerhits, particle, tmpEfrac, tmpEcomplet);
	if (particle) {
	  std::cout << "shower pdg: "<< particle->PdgCode() << " efrac " << tmpEfrac << std::endl;

	  highestHitsPDG = particle->PdgCode();
	  highestHitsFrac = tmpEfrac;
	}
      }
    }
  }

  fTree->Fill();

} // analyze

// -------------------------------------------------

void shower::TCShowerAnalysis::reset() {

  run = -99999;
  subrun = -99999;
  event = -99999;

  nuPDG_truth = -99999;
  ccnc_truth = -99999;
  mode_truth = -99999;

  nshws = 0;
  for (int i = 0; i < kMaxShowers; ++i){
    shwid[i] = -99999;
    shwdcosx[i] = -99999;
    shwdcosy[i] = -99999;
    shwdcosz[i] = -99999;
    shwstartx[i] = -99999;
    shwstarty[i] = -99999;
    shwstartz[i] = -99999;
    for (int j = 0; j < 2; ++j){
      shwdedx[i][j] = -99999;
    }
    shwbestplane[i] = -99999;
  }

  highestHitsPDG = -99999;
  highestHitsFrac = -99999;

  return;

} // reset

// -------------------------------------------------

void shower::TCShowerAnalysis::truthMatcher(std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> shower_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet){

  MCparticle=0;
  Efrac=1.0;
  Ecomplet=0;

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  std::map<int,double> trkID_E;
  for(size_t j = 0; j < shower_hits.size(); ++j){
    art::Ptr<recob::Hit> hit = shower_hits[j];
    //For know let's use collection plane to look at the shower reconstruction                                            
    if( hit->View() != 1) continue;                                                                                     
    std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToEveTrackIDEs(hit);
    for(size_t k = 0; k < TrackIDs.size(); k++){
      if (trkID_E.find(std::abs(TrackIDs[k].trackID))==trkID_E.end()) trkID_E[std::abs(TrackIDs[k].trackID)] = 0;
      trkID_E[std::abs(TrackIDs[k].trackID)] += TrackIDs[k].energy;
    }
  }
  double max_E = -999.0;
  double total_E = 0.0;
  int TrackID = -999;
  double partial_E=0.0;
  //double noEM_E = 0.0;  //non electromagnetic energy is defined as energy from charged pion and protons                 
  if( !trkID_E.size() ) return; //Ghost shower???                                                                         
  for(std::map<int,double>::iterator ii = trkID_E.begin(); ii!=trkID_E.end(); ++ii){
    total_E += ii->second;
    if((ii->second)>max_E){
      partial_E = ii->second;
      max_E = ii->second;
      TrackID = ii->first;
    }
    //int ID = ii->first;                                                                                                 
    // const simb::MCParticle *particle = pi_serv->TrackIDToParticle(ID);                                                 
    //if( abs(particle->PdgCode()) == 211 || particle->PdgCode() == 2212 ){                                               
    //if( particle->PdgCode() != 22 && abs(particle->PdgCode()) != 11){                                                   
    //noEM_E += ii->second;                                                                                               
    //}                                                                                                                   

  }
  MCparticle = pi_serv->TrackIdToParticle_P(TrackID);

  //  Efrac = 1-(partial_E/total_E);
  Efrac = partial_E/total_E;

  //completeness                                                                                                          
  double totenergy =0;
  for(size_t k = 0; k < all_hits.size(); ++k){
    art::Ptr<recob::Hit> hit = all_hits[k];
    std::vector<sim::TrackIDE> TrackIDs = bt_serv->HitToEveTrackIDEs(hit);
    for(size_t l = 0; l < TrackIDs.size(); ++l){
      if(std::abs(TrackIDs[l].trackID)==TrackID) {
	totenergy += TrackIDs[l].energy;
      }
    }
  }
  Ecomplet = partial_E/totenergy;

} // truthMatcher

// -------------------------------------------------

DEFINE_ART_MODULE(shower::TCShowerAnalysis)
