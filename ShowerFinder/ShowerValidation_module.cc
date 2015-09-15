//////////////////////////////////////////////////////////////////////////////////////
// Class:       ShowerValidation
// Module type: analyser
// File:        ShowerValidation_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffied.ac.uk), July 2015
//
// ...
//////////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/FindManyP.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Shower.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RawData/ExternalTrigger.h"
#include "MCCheater/BackTracker.h"
#include "AnalysisBase/ParticleID.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"

// ROOT & STL includes
#include "TTree.h"
#include "TTimeStamp.h"
#include "TLorentzVector.h"
#include "TH2F.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TLegend.h"
#include "TFolder.h"
#include "TStyle.h"
#include "TAttMarker.h"

#include <map>
#include <iostream>
#include <algorithm>
#include <iostream>
#include <fstream>

namespace ShowerValidation {
  class ShowerValidation;
  class ShowerCounter;
}

class ShowerValidation::ShowerCounter {
public:

  explicit ShowerCounter(int plane);

  void AddRawHit(int trackID);
  void AddSignalShowerHit(int showerID);
  void AddNoiseShowerHit(int showerID);
  void AssociateShowerAndTrack(int showerID, int trackID);
  double GetCompleteness(int showerID);
  double GetCleanliness(int showerID);
  int GetNumShowerHits(int showerID);
  std::vector<int> GetShowerIDs();

private:

  int fPlane;

  // Maps to hold hit information
  std::map<int,int> signalHits;
  std::map<int,int> noiseHits;
  std::map<int,int> rawHits;
  std::map<int,int> showerToTrackID;
  std::map<int,std::vector<int> > trackToShowerIDs;

};

class ShowerValidation::ShowerValidation : public art::EDAnalyzer {
public:

  explicit ShowerValidation(fhicl::ParameterSet const &p);

  void analyze(art::Event const &evt);
  void beginJob();
  void endJob();
  void reconfigure(fhicl::ParameterSet const &p);
  int FindTrackID(art::Ptr<recob::Hit> &hit);
  int FindTrueTrack(std::vector<art::Ptr<recob::Hit> > &showerHits);

private:

  // Showerings to compare and the hits which the clustering was run over
  std::string fShowerModuleLabel;
  std::string fHitsModuleLabel;

  // Histograms to fill
  TH1D *hCompleteness, *hCleanliness, *hMiddleCompletenessCleanliness;
  TH2D *hEnergyVCharge;

  // Canvas on which to save histograms
  TCanvas *fCanvas;

  // Map to hold all shower information
  std::map<int,std::unique_ptr<ShowerCounter> > showerMap;

  std::map<int,const simb::MCParticle*> trueParticles;
  art::ServiceHandle<cheat::BackTracker> backtracker;

};

ShowerValidation::ShowerCounter::ShowerCounter(int plane) {
  fPlane = plane;
}

void ShowerValidation::ShowerCounter::AddRawHit(int trackID) { ++rawHits[trackID]; }

void ShowerValidation::ShowerCounter::AddSignalShowerHit(int showerID) { ++signalHits[showerID]; }

void ShowerValidation::ShowerCounter::AddNoiseShowerHit(int showerID) { ++noiseHits[showerID]; }

void ShowerValidation::ShowerCounter::AssociateShowerAndTrack(int showerID, int trackID) { showerToTrackID[showerID] = trackID; trackToShowerIDs[trackID].push_back(showerID); }

double ShowerValidation::ShowerCounter::GetCompleteness(int showerID) { return (double)signalHits[showerID]/(double)rawHits[showerToTrackID[showerID]]; }

double ShowerValidation::ShowerCounter::GetCleanliness(int showerID) { return (double)signalHits[showerID]/GetNumShowerHits(showerID); }

int ShowerValidation::ShowerCounter::GetNumShowerHits(int showerID) { return signalHits[showerID] + noiseHits[showerID]; }

std::vector<int> ShowerValidation::ShowerCounter::GetShowerIDs() { std::vector<int> v; for (std::map<int,int>::iterator it = showerToTrackID.begin(); it != showerToTrackID.end(); ++it) v.push_back(it->first); return v; }

ShowerValidation::ShowerValidation::ShowerValidation(fhicl::ParameterSet const &pset) :  EDAnalyzer(pset) {
  this->reconfigure(pset);
  fCanvas = new TCanvas("fCanvas","",800,600);
  gStyle->SetOptStat(0);
}

void ShowerValidation::ShowerValidation::reconfigure(fhicl::ParameterSet const& p) {

}

void ShowerValidation::ShowerValidation::beginJob() {
  hCompleteness = new TH1D("Completeness",";Completeness;",101,0,1.01);
  hCleanliness = new TH1D("Cleanliness",";Cleanliness;",101,0,1.01);
  hMiddleCompletenessCleanliness = new TH1D("MiddleCompletenessCleanliness",";Completeness * Cleanliness for the second best/worst plane;",101,0,1.01);
  hEnergyVCharge = new TH2D("EnergyVCharge",";True Energy Deposited;Charge;",1000,0,2000,1000,0,100000);
}

void ShowerValidation::ShowerValidation::analyze(art::Event const &evt) {

  // Get the showers from the event
  art::Handle<std::vector<recob::Shower> > showerHandle;
  std::vector<art::Ptr<recob::Shower> > showers;
  if (evt.getByLabel("testshower",showerHandle))
    art::fill_ptr_vector(showers, showerHandle);

  // Get the hits out of the event
  art::Handle<std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if (evt.getByLabel("dcheat",hitHandle))
    art::fill_ptr_vector(hits, hitHandle);

  // All the associations
  art::FindManyP<recob::Hit> fmh(showerHandle, evt, "testshower");
  art::FindManyP<recob::SpacePoint> fmsp(showerHandle, evt, "testshower");
  art::FindManyP<recob::Cluster> fmc(showerHandle, evt, "testshower");
  art::FindManyP<recob::Track> fmt(showerHandle, evt, "testshower");

  // Make a map of ShowerCounter s, one for each plane
  for (int plane = 0; plane < 3; ++plane)
    showerMap[plane] = (std::unique_ptr<ShowerCounter>) new ShowerCounter(plane);

  // Save raw hits
  for (size_t hitIt = 0; hitIt < hits.size(); ++hitIt) {
    art::Ptr<recob::Hit> hit = hits.at(hitIt);
    int trackID = FindTrackID(hit);
    showerMap[hit->View()]->AddRawHit(trackID);
  }

  // Save true tracks
  trueParticles.clear();
  const sim::ParticleList& particles = backtracker->ParticleList();
  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt) {
    const simb::MCParticle *particle = particleIt->second;
    trueParticles[particle->TrackId()] = particle;
  }

  // Look at the showers
  for (auto &shower : showers) {

    int id = shower->ID();

    // All the hits associated with the shower
    std::vector<art::Ptr<recob::Hit> > showerHits = fmh.at(id);

    // Find the true track this shower is associated with
    int trueTrackID = FindTrueTrack(showerHits);

    // Associate shower and track
    for (int plane = 0; plane < 3; ++plane)
      showerMap[plane]->AssociateShowerAndTrack(id, trueTrackID);

    // Look at charge and energy of this shower/particle track
    std::vector<double> charges = shower->Energy();
    std::sort(charges.begin(), charges.end(), std::greater<double>());
    std::vector<sim::IDE> ides = backtracker->TrackIDToSimIDE(trueTrackID);
    double energyDeposited = 0;
    for (auto &ide : ides)
      energyDeposited += ide.energy;
    std::cout << "Energy... " << energyDeposited << " and charge " << charges.at(0) << std::endl;
    hEnergyVCharge->Fill(energyDeposited, charges.at(0));

    // Look at all the hits associated with the shower
    for (auto &hit : showerHits) {

      int trackID = FindTrackID(hit);
      if (trackID == trueTrackID) showerMap[hit->View()]->AddSignalShowerHit(id);
      else                        showerMap[hit->View()]->AddNoiseShowerHit(id);

    }

  }

  // Fill histograms
  for (unsigned int plane = 0; plane < 3; ++plane) {

    std::vector<int> showerIDs = showerMap[plane]->GetShowerIDs();

    for (auto &showerID : showerIDs) {

      hCompleteness->Fill(showerMap[plane]->GetCompleteness(showerID));
      hCleanliness->Fill(showerMap[plane]->GetCleanliness(showerID));

    }
  }

  // Histogram for showing second best/worst plane!
  for (auto &shower : showers) {
    std::vector<double> vec;
    vec.push_back(showerMap[0]->GetCompleteness(shower->ID()) * showerMap[0]->GetCleanliness(shower->ID()));
    vec.push_back(showerMap[1]->GetCompleteness(shower->ID()) * showerMap[1]->GetCleanliness(shower->ID()));
    vec.push_back(showerMap[2]->GetCompleteness(shower->ID()) * showerMap[2]->GetCleanliness(shower->ID()));
    std::sort(vec.begin(), vec.end());
    hMiddleCompletenessCleanliness->Fill(vec.at(1));
  }

}

int ShowerValidation::ShowerValidation::FindTrackID(art::Ptr<recob::Hit> &hit) {

  /// Finds the track which this hit is most associated with

  double particleEnergy = 0;
  int likelyTrackID = 0;
  std::vector<sim::TrackIDE> trackIDs = backtracker->HitToTrackID(hit);
  for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
    if (trackIDs.at(idIt).energy > particleEnergy) {
      particleEnergy = trackIDs.at(idIt).energy;
      likelyTrackID = std::abs(trackIDs.at(idIt).trackID);
    }
  }

  return likelyTrackID;

}

int ShowerValidation::ShowerValidation::FindTrueTrack(std::vector<art::Ptr<recob::Hit> > &showerHits) {

  /// Takes all the hits associated with a shower and determines which particle corresponds most

  std::map<int,double> trackMap;
  for (std::vector<art::Ptr<recob::Hit> >::iterator clusHitIt = showerHits.begin(); clusHitIt != showerHits.end(); ++clusHitIt) {
    art::Ptr<recob::Hit> hit = *clusHitIt;
    int trackID = FindTrackID(hit);
    trackMap[trackID] += hit->Integral();
  }

  double highestCharge = 0;
  int clusterTrack = 0;
  for (std::map<int,double>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt)
    if (trackIt->second > highestCharge) {
      highestCharge = trackIt->second;
      clusterTrack  = trackIt->first;
    }

  return clusterTrack;

}

void ShowerValidation::ShowerValidation::endJob() {

  fCanvas->cd();
  hCompleteness->Draw();
  fCanvas->SaveAs("Completeness.png");

  fCanvas->Clear();
  hCleanliness->Draw();
  fCanvas->SaveAs("Cleanliness.png");

  fCanvas->Clear();
  hEnergyVCharge->SetMarkerStyle(8);
  hEnergyVCharge->SetMarkerSize(1);
  hEnergyVCharge->Draw();
  fCanvas->SaveAs("EnergyVCharge.png");

  fCanvas->Clear();
  hMiddleCompletenessCleanliness->Draw();
  fCanvas->SaveAs("MiddleCompletenessCleanliness.png");

  delete fCanvas;

}

DEFINE_ART_MODULE(ShowerValidation::ShowerValidation)
