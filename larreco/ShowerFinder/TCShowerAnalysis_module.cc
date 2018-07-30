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
#include "TF1.h"
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
    void showerProfileTrue(std::vector< art::Ptr<sim::SimChannel> > allchan, simb::MCParticle electron);
 
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
    std::string fDigitModuleLabel;

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
  fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel", "largeant")  ),
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
  fShowerProfile = tfs->make<TProfile>("fShowerProfile", "fShowerProfile", 20, 0, 5);
  //  fShowerProf = tfs->make<TH2F>("fShowerProf", "fShowerProf", 50, 0, 5, 50, 0, 100000);

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

  art::Handle< std::vector<sim::SimChannel> > scListHandle;
  std::vector<art::Ptr<sim::SimChannel> > simchanlist;
  if (evt.getByLabel(fDigitModuleLabel,scListHandle))
    art::fill_ptr_vector(simchanlist, scListHandle);

  art::Handle< std::vector<recob::Shower> > showerListHandle;
  std::vector<art::Ptr<recob::Shower> > showerlist;
  if (evt.getByLabel(fShowerModuleLabel,showerListHandle))
    art::fill_ptr_vector(showerlist, showerListHandle);
  
  // mc info
  art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
  std::vector<art::Ptr<simb::MCTruth> > mclist;
  if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);

  /*
  art::FindManyP<recob::Hit> shwfm(showerListHandle, evt, fShowerModuleLabel);

  if (showerlist.size()) {
    std::vector< art::Ptr<recob::Hit> > showerhits = shwfm.at(0);
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

  /*
  if (mclist.size()) {
    art::Ptr<simb::MCTruth> mctruth = mclist[0];
    if (mctruth->NeutrinoSet()) {
      if (std::abs(mctruth->GetNeutrino().Nu().PdgCode()) == 12 && mctruth->GetNeutrino().CCNC() == 0) {
	double elep =  mctruth->GetNeutrino().Lepton().E();
	std::cout << "ELECTRON ENERGY: " << elep << std::endl;
	if (elep > 4 && elep < 5) {
	  showerProfileTrue(hitlist);
	}
      }
    }
  }
  */

  if (mclist.size()) {
    art::Ptr<simb::MCTruth> mctruth = mclist[0];
    if (mctruth->NeutrinoSet()) { 
      if (std::abs(mctruth->GetNeutrino().Nu().PdgCode()) == 12 && mctruth->GetNeutrino().CCNC() == 0) { 
	double elep =  mctruth->GetNeutrino().Lepton().E();
	std::cout << "ELECTRON ENERGY: " << elep << std::endl;
        if (elep > 4 && elep < 5) {
	  showerProfileTrue(simchanlist, mclist[0]->GetNeutrino().Lepton());
	}
      }
    }
  }

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

  double shwTwoTime = detprop->ConvertXToTicks(shwvtx[0]+shwdir[0], collectionPlane);
  double shwTwoWire = geom->WireCoordinate(shwvtx[1]+shwdir[1], shwvtx[2]+shwdir[2], collectionPlane);

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

    double xtwoorth = (ytwo - yvtx) + xvtx;
    double ytwoorth = -(xtwo - xvtx) + yvtx;
    
    double xhit = showerhits[i]->PeakTime() * tickToDist;
    double yhit = showerhits[i]->WireID().Wire * wirePitch;

    double dist = std::abs((ytwoorth-yvtx)*xhit - (xtwoorth-xvtx)*yhit + xtwoorth*yvtx - ytwoorth*xvtx)/std::sqrt( pow((ytwoorth-yvtx), 2) + pow((xtwoorth-xvtx), 2) );

    double to3D = 1. / sqrt( pow(xvtx-xtwo,2) + pow(yvtx-ytwo,2) ) ; // distance between two points in 3D space is one 
    dist *= to3D;
    double Q = showerhits[i]->Integral() * fCalorimetryAlg.LifetimeCorrection(showerhits[i]->PeakTime());
    double t = dist / 14; // convert to radiation lengths    
    std::cout << t << std::endl;
    int bin = floor(t*3);

    temp->SetBinContent(bin+1, temp->GetBinContent(bin+1) + Q);

  } // loop through showerhits

  for (int i = 0; i < 15; ++i) {
    fShowerProfile->Fill(temp->GetBinCenter(i+1), temp->GetBinContent(i+1));
  }

} // showerProfile

// -------------------------------------------------

void shower::TCShowerAnalysis::showerProfileTrue(std::vector< art::Ptr<recob::Hit> > allhits) {

  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<geo::Geometry> geom;
  auto collectionPlane = geo::PlaneID(0, 0, 1);
  art::ServiceHandle<cheat::BackTrackerService> btserv;
  art::ServiceHandle<cheat::ParticleInventoryService> piserv;
  std::map<int,double> trkID_E;  

  TH1F* temp = new TH1F("temp", "temp", 15, 0, 5);

  //  std::vector< art::Ptr<recob::Hit> > showerhits;

  double xvtx = -999;
  double yvtx = -999;
  double zvtx = -999;
  double xtwo = -999;
  double ytwo = -999;
  double ztwo = -999;
  double shwvtxT = -999;
  double shwvtxW = -999;
  double shwtwoT = -999;
  double shwtwoW = -999;

  double shwvtxx = -999;
  double shwvtxy = -999;
  double shwtwox = -999;
  double shwtwoy = -999;
  double xtwoorth = -999;
  double ytwoorth = -999;

  double wirePitch = -999;
  double tickToDist = -999;

  bool foundParent = false;

  for (size_t i = 0; i < allhits.size(); ++i) {
    if (allhits[i]->WireID().Plane != collectionPlane.Plane) continue;

    art::Ptr<recob::Hit> hit = allhits[i];
    std::vector<sim::TrackIDE> trackIDs = btserv->HitToEveTrackIDEs(hit);

    for (size_t j = 0; j < trackIDs.size(); ++j) {
      // only want energy associated with the electron and electron must have neutrino mother
      if ( std::abs((piserv->TrackIdToParticle_P(trackIDs[j].trackID))->PdgCode()) != 11) continue;
      if ( std::abs((piserv->TrackIdToParticle_P(trackIDs[j].trackID))->Mother()) != 0) continue;

      if (!foundParent) {
	xvtx = (piserv->TrackIdToParticle_P(trackIDs[j].trackID))->Vx();
	yvtx = (piserv->TrackIdToParticle_P(trackIDs[j].trackID))->Vy();
	zvtx = (piserv->TrackIdToParticle_P(trackIDs[j].trackID))->Vz();

	xtwo = (piserv->TrackIdToParticle_P(trackIDs[j].trackID))->EndX();
	ytwo = (piserv->TrackIdToParticle_P(trackIDs[j].trackID))->EndY();
	ztwo = (piserv->TrackIdToParticle_P(trackIDs[j].trackID))->EndZ();

	shwvtxT = detprop->ConvertXToTicks(xvtx, collectionPlane);
	shwvtxW = geom->WireCoordinate(yvtx, zvtx, collectionPlane);

	shwtwoT = detprop->ConvertXToTicks(xtwo, collectionPlane);
	shwtwoW = geom->WireCoordinate(ytwo, ztwo, collectionPlane);

	wirePitch = geom->WirePitch(hit->WireID());
	tickToDist = detprop->DriftVelocity(detprop->Efield(),detprop->Temperature());
	tickToDist *= 1.e-3 * detprop->SamplingRate(); // 1e-3 is conversion of 1/us to 1/ns

	shwvtxx = shwvtxT * tickToDist;
	shwvtxy = shwvtxW * wirePitch;

	shwtwox = shwtwoT * tickToDist;
	shwtwoy = shwtwoW * wirePitch;

	xtwoorth = (shwtwoy - shwvtxy) + shwvtxx;
	ytwoorth = -(shwtwox - shwvtxx) + shwvtxy;

	foundParent = true;
      }
      double xhit = hit->PeakTime() * tickToDist;
      double yhit = hit->WireID().Wire * wirePitch;

      double dist = abs((ytwoorth-shwvtxy)*xhit - (xtwoorth-shwvtxx)*yhit + xtwoorth*shwvtxy - ytwoorth*shwvtxx)/sqrt( pow((ytwoorth-shwvtxy), 2) + pow((xtwoorth-shwvtxx), 2) );

      double to3D = sqrt( pow(xvtx-xtwo , 2) + pow(yvtx-ytwo , 2) + pow(zvtx-ztwo , 2) ) / sqrt( pow(shwvtxx-shwtwox,2) + pow(shwvtxy-shwtwoy,2) ) ; // distance between two points in 3D space is one 
      dist *= to3D;
      double E = trackIDs[j].energy;
      double t = dist / 14; // convert to radiation lengths
      //      std::cout << t << " " << E << std::endl;
      int bin = floor(t*3);

      temp->SetBinContent(bin+1, temp->GetBinContent(bin+1) + E);

      break;
      //      trkID_E[std::abs(trackIDs[j].trackID)] += trackIDs[j].energy;
    } // loop through track IDE

  } // loop through all hits

  for (int i = 0; i < 15; ++i) {
    fShowerProfile->Fill(temp->GetBinCenter(i+1), temp->GetBinContent(i+1));
  }

  //  std::cout << xvtx << " " << yvtx << " " << zvtx << " " << xtwo << " " << ytwo << " " << ztwo << std::endl;

  /*
  if (!trkID_E.size()) return; 

  for (std::map<int,double>::iterator ii = trkID_E.begin(); ii != trkID_E.end(); ++ii) {
    const simb::MCParticle* mcpart = piserv->TrackIdToParticle_P(ii->first);

    std::cout << "PARTICLE ID " << mcpart->PdgCode() << " " << mcpart->Vx() << " " << mcpart->Vy() << " " << mcpart->Vz() << " " << mcpart->Mother() << std::endl;

  } // loop through trkID_E
  */
  return;
} // showerProfileTrue

// -------------------------------------------------

void shower::TCShowerAnalysis::showerProfileTrue(std::vector< art::Ptr<sim::SimChannel> > allchan, simb::MCParticle electron) {
  art::ServiceHandle<cheat::ParticleInventoryService> piserv;
  art::ServiceHandle<geo::Geometry> geom;

  std::vector<sim::MCEnDep> alledep;

  TH1D* temp = new TH1D("temp", "temp", 20, 0, 5);

  for (size_t i = 0; i < allchan.size(); ++i) {
    art::Ptr<sim::SimChannel> simchan = allchan[i];
    if ( geom->View(simchan->Channel()) != geo::kV ) continue;
    auto tdc_ide_map = simchan->TDCIDEMap();

    for(auto const& tdc_ide_pair : tdc_ide_map) {
      for (auto const& ide : tdc_ide_pair.second) {
	if (piserv->TrackIdToMotherParticle_P(ide.trackID) == NULL) continue;
	if ( std::abs(piserv->TrackIdToMotherParticle_P(ide.trackID)->PdgCode()) != 11 ) continue;

	sim::MCEnDep edep;
	edep.SetVertex(ide.x,ide.y,ide.z);
	edep.SetEnergy(ide.energy);
	edep.SetTrackId(ide.trackID);
	
	alledep.push_back(edep);
	
      } // loop through ide
    } // loop through tdc_ide

  } // loop through channels

  double x0 = electron.Vx();
  double y0 = electron.Vy();
  double z0 = electron.Vz();

  double charge4cm = 0;
 
  for (size_t i = 0; i < alledep.size(); ++i) {
    double x = (double)alledep[i].Vertex()[0];
    double y = (double)alledep[i].Vertex()[1];
    double z = (double)alledep[i].Vertex()[2];

    double dist = std::sqrt( std::pow(x0-x, 2) + std::pow(y0-y, 2) + std::pow(z0-z, 2) );

    if (dist < 3.5) charge4cm += alledep[i].Energy();

    double E = alledep[i].Energy();
    double t = dist / 14; // convert to radiation lengths                                                               

    int bin = floor(t*4);
    
    temp->SetBinContent(bin+1, temp->GetBinContent(bin+1) + E);

  }
  
  for (int i = 0; i < 20; ++i) {
    if (temp->GetBinContent(i+1) == 0) continue;
    if (temp->GetBinContent(i+2) == 0) continue;
    fShowerProfile->Fill(temp->GetBinCenter(i+1), temp->GetBinContent(i+1));
  }

  return;

} // showerProfileTrue 

// -------------------------------------------------

DEFINE_ART_MODULE(shower::TCShowerAnalysis)

