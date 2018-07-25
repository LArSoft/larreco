// -------------------------------------------------
// shower analysis module
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
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"

namespace shower {

  class TCShowerAnalysis : public art::EDAnalyzer {

  public:

    explicit TCShowerAnalysis(fhicl::ParameterSet const& pset);
    virtual ~TCShowerAnalysis();

    void reconfigure(fhicl::ParameterSet const& pset);
    void beginJob();
    void analyze(const art::Event& evt);
   
    void showerProfile(std::vector< art::Ptr<recob::Hit> > showerhits, TVector3 shwvtx, TVector3 shwdir);
    void showerProfileTrue(std::vector< art::Ptr<recob::Hit> > allhits);
 
  protected:

  private: 

    void resetVars();

    TTree* fTree;

    TProfile* fShowerProfile;

    std::string fClusterModuleLabel;
    std::string fTrackModuleLabel;
    std::string fHitModuleLabel;
    std::string fShowerModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fGenieGenModuleLabel;

    calo::CalorimetryAlg fCalorimetryAlg;

  }; // class TCShowerAnalysis

} // shower

// -------------------------------------------------

shower::TCShowerAnalysis::TCShowerAnalysis(fhicl::ParameterSet const& pset) :
  EDAnalyzer(pset),
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel", "trajcluster" ) ),
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel", "pmtrack" ) ),
  fHitModuleLabel           (pset.get< std::string >("HitModuleLabel", "trajcluster" ) ),
  fShowerModuleLabel        (pset.get< std::string >("ShowerModuleLabel", "tcshower" ) ),
  fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel", "calo")  ), 
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel", "generator")  ),
  fCalorimetryAlg           (pset.get< fhicl::ParameterSet >("CalorimetryAlg") ) {
  this->reconfigure(pset);
} // TCShowerAnalysis

// -------------------------------------------------

shower::TCShowerAnalysis::~TCShowerAnalysis() {
} // ~TCShowerAnalysis

// -------------------------------------------------

void shower::TCShowerAnalysis::reconfigure(fhicl::ParameterSet const& pset) {
} // reconfigure

// -------------------------------------------------

void shower::TCShowerAnalysis::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tcshowerana", "tcshowerana");
  fShowerProfile = tfs->make<TProfile>("fShowerProfile", "fShowerProfile", 15, 0, 5);

} // beginJob

// -------------------------------------------------

void shower::TCShowerAnalysis::analyze(const art::Event& evt) {

  resetVars();

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);
  
  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
    art::fill_ptr_vector(clusterlist, clusterListHandle);

  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower> > showerlist;
  if (evt.getByLabel(fShowerModuleLabel,showerListHandle))
    art::fill_ptr_vector(showerlist, showerListHandle);
  
  // mc info
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);

  showerProfileTrue(hitlist);

  /*
  art::FindManyP<recob::Hit> shwfm(showerListHandle, evt, fShowerModuleLabel);

  if (showerlist.size()) {
    std::vector< art::Ptr<recob::Hit> > showerhits = shwfm.at(0);

    // TODO: shower profile for E = 4-5 GeV
    if (mclist.size()) {
      art::Ptr<simb::MCTruth> mctruth = mclist[0];
      if (mctruth->NeutrinoSet()) {
	if (std::abs(mctruth->GetNeutrino().Nu().PdgCode()) == 12 && mctruth->GetNeutrino().CCNC() == 0) {
	  double elep =  mctruth->GetNeutrino().Lepton().E();
	  std::cout << "ELECTRON ENERGY: " << elep << std::endl;
	  if (elep > 4 && elep < 5) {
	    showerProfile(showerhits, showerlist[0]->ShowerStart(), showerlist[0]->Direction());
	  }
	}
      }
    }
  }
  */
  fTree->Fill();

} // analyze

// -------------------------------------------------

void shower::TCShowerAnalysis::resetVars() {

} // resetVars

// -------------------------------------------------

void shower::TCShowerAnalysis::showerProfile(std::vector< art::Ptr<recob::Hit> > showerhits, TVector3 shwvtx, TVector3 shwdir) {

  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<geo::Geometry> geom;

  auto collectionPlane = geo::PlaneID(0, 0, 1);

  double shwVtxTime = detprop->ConvertXToTicks(shwvtx[0], collectionPlane);
  double shwVtxWire = geom->WireCoordinate(shwvtx[1], shwvtx[2], collectionPlane);

  TVector3 shwort = shwdir.Orthogonal().Unit();
  double shwTwoTime = detprop->ConvertXToTicks(shwvtx[0]+shwort[0], collectionPlane);
  double shwTwoWire = geom->WireCoordinate(shwvtx[1]+shwort[1], shwvtx[2]+shwort[2], collectionPlane);

  TH1F* temp = new TH1F("temp", "temp", 15, 0, 5);

  for (size_t i = 0; i < showerhits.size(); ++i) {
    if (showerhits[i]->WireID().Plane != collectionPlane.Plane) continue;

    double wirePitch = geom->WirePitch(showerhits[i]->WireID());
    double tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
    tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns

    double xvtx = shwVtxTime * tickToDist;
    double yvtx = shwVtxWire * wirePitch;

    double xtwo = shwTwoTime * tickToDist;
    double ytwo = shwTwoWire * wirePitch;

    double xhit = showerhits[i]->PeakTime() * tickToDist;
    double yhit = showerhits[i]->WireID().Wire * wirePitch;

    double dist = std::abs((ytwo-yvtx)*xhit - (xtwo-xvtx)*yhit + xtwo*yvtx - ytwo*xvtx)/std::sqrt( pow((ytwo-yvtx), 2) + pow((xtwo-xvtx), 2) );

    double to3D = 1. / sqrt( pow(xvtx-xtwo,2) + pow(yvtx-ytwo,2) ) ; // distance between two points in 3D space is one 
    dist *= to3D;
    double Q = showerhits[i]->Integral() * fCalorimetryAlg.LifetimeCorrection(showerhits[i]->PeakTime());
    double t = dist / 14; // convert to radiation lengths    
    int bin = floor(t*3);

    //    fShowerProfile->Fill(t, Q);

    temp->SetBinContent(bin, temp->GetBinContent(bin) + Q);

  } // loop through showerhits

  for (int i = 0; i < 15; ++i) {
    fShowerProfile->Fill(temp->GetBinCenter(i), temp->GetBinContent(i));
  }

} // showerProfile

// -------------------------------------------------

void shower::TCShowerAnalysis::showerProfileTrue(std::vector< art::Ptr<recob::Hit> > allhits) {

  art::ServiceHandle<cheat::BackTrackerService> btserv;
  art::ServiceHandle<cheat::ParticleInventoryService> piserv;
  std::map<int,double> trkID_E;  

  for (size_t i = 0; i < allhits.size(); ++i) {

    art::Ptr<recob::Hit> hit = allhits[i];
    std::vector<sim::TrackIDE> trackIDs = btserv->HitToEveTrackIDEs(hit);

    for (size_t j = 0; j < trackIDs.size(); ++j) {
      // only want energy associated with the electron and electron must have neutrino mother
      if ( std::abs((piserv->TrackIdToParticle_P(trackIDs[j].trackID))->PdgCode()) != 11) continue;
      if ( std::abs((piserv->TrackIdToParticle_P(trackIDs[j].trackID))->Mother()) != 0) continue;

      trkID_E[std::abs(trackIDs[j].trackID)] += trackIDs[j].energy;
    } // loop through track IDE

  } // loop through all hits

  if (!trkID_E.size()) return; 

  for (std::map<int,double>::iterator ii = trkID_E.begin(); ii != trkID_E.end(); ++ii) {
    const simb::MCParticle* mcpart = piserv->TrackIdToParticle_P(ii->first);

    std::cout << "PARTICLE ID " << mcpart->PdgCode() << " " << mcpart->Vx() << " " << mcpart->Vy() << " " << mcpart->Vz() << " " << mcpart->Mother() << std::endl;

  } // loop through trkID_E

  return;
} // showerProfileTrue

// -------------------------------------------------

DEFINE_ART_MODULE(shower::TCShowerAnalysis)

