//////////////////////////////////////////////////////////////////////////////
// Class:       ClusteringValidation
// Module type: analyser
// File:        ClusteringValidation_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffied.ac.uk), May 2015
//
// A module to validate clustering algorithms.
// Compares the output of the clustering to the output of the hit finding
// (specified in fcl files) and makes all validation plots.
// Written to work over a pi0 sample.
//////////////////////////////////////////////////////////////////////////////

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

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
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

#include <map>
#include <iostream>
#include <algorithm>
#include <iostream>
#include <fstream>

namespace ClusteringValidation {
  class ClusteringValidation;
  class ClusterCounter;
}

enum class ClusterID : int { };
enum class TrackID   : int { };
typedef std::vector<ClusterID> ClusterIDs;
typedef std::vector<TrackID> TrackIDs;

class ClusteringValidation::ClusterCounter {
public:

  explicit ClusterCounter(unsigned int tpc, unsigned int plane);
  ~ClusterCounter();

  void				              AddHitPreClustering       (TrackID id);
  void				              AddSignalHitPostClustering(ClusterID id);
  void				              AddNoiseHitPostClustering (ClusterID id);
  void				              AssociateClusterAndTrack  (ClusterID clusID, TrackID trackID);
  double			              GetCompleteness           (ClusterID id);
  double			              GetCleanliness            (ClusterID id);
  double			              GetEfficiency             (TrackID id);
  ClusterIDs	                              GetListOfClusterIDs       ();
  TrackIDs		                      GetListOfTrackIDs         ();
  int				              GetNumberHitsFromTrack    (TrackID id);
  int	                                      GetNumberHitsInCluster    (ClusterID id);
  std::vector<std::pair<TrackID,ClusterIDs> > GetPhotons                ();
  TrackID                                     GetTrack                  (ClusterID id);
  bool		                              IsNoise                   (ClusterID id);
  bool                       	              IsNoise                   (TrackID id);
  bool                                        PassesCut                 ();

private:

  unsigned int tpc, plane;

  std::map<TrackID,int>                           numHitsPreClustering;
  std::map<ClusterID,int>                         numSignalHitsPostClustering;
  std::map<ClusterID,int>                         numNoiseHitsPostClustering;
  std::map<ClusterID,TrackID>                     clusterToTrackID;
  std::map<TrackID,ClusterIDs>                    trackToClusterIDs;
  std::map<TrackID,std::map<std::string,double> > particleProperties;
  std::map<TrackID,simb::MCParticle>              trueParticles;

  art::ServiceHandle<geo::Geometry> geometry;
  art::ServiceHandle<cheat::BackTracker> backtracker;

};

class ClusteringValidation::ClusteringValidation : public art::EDAnalyzer {
public:

  explicit ClusteringValidation(fhicl::ParameterSet const &p);
  virtual ~ClusteringValidation();

  void                    analyze(art::Event const &evt);
  TrackID                 FindTrackID(art::Ptr<recob::Hit> &hit);
  TrackID                 FindTrueTrack(std::vector<art::Ptr<recob::Hit> > &clusterHits);
  double                  FindPhotonAngle();
  const simb::MCParticle* GetPi0();
  double                  GetEndTrackDistance(TrackID id1, TrackID id2);
  void                    beginJob();
  void                    endJob();
  void                    reconfigure(fhicl::ParameterSet const &p);

private:

  // Clustering to validate and the hits which the clustering was run over
  std::string fClusterModuleLabel, fHitsModuleLabel;

  // hists
  TH1 *hCompleteness, *hCleanliness;
  TH1 *hPi0Angle, *hPi0Energy, *hPi0ConversionDistance, *hPi0ConversionSeparation, *hPi0AngleCut, *hPi0EnergyCut, *hPi0ConversionDistanceCut, *hPi0ConversionSeparationCut;
  TH2 *hNumHitsCompleteness, *hNumHitsEnergy;
  TProfile *hCompletenessEnergy, *hCompletenessAngle, *hCompletenessConversionDistance, *hCompletenessConversionSeparation;
  TProfile *hCleanlinessEnergy, *hCleanlinessAngle, *hCleanlinessConversionDistance, *hCleanlinessConversionSeparation;
  TEfficiency *hEfficiencyAngle, *hEfficiencyEnergy, *hEfficiencyConversionDistance, *hEfficiencyConversionSeparation;
  TCanvas *fCanvas;
  TObjArray fHistArray;

  std::map<unsigned int,std::map<unsigned int,ClusterCounter*> > clusterMap;
  std::map<TrackID,const simb::MCParticle*>                      trueParticles;

  art::ServiceHandle<geo::Geometry>      geometry;
  art::ServiceHandle<cheat::BackTracker> backtracker;

  // Average completenesses and cleanlinesses for all events
  double totCompleteness = 0, totCleanliness = 0;
  int nClusters = 0;

};

ClusteringValidation::ClusterCounter::ClusterCounter(unsigned int t, unsigned int p)  {
  tpc   = t;
  plane = p;
}
ClusteringValidation::ClusterCounter::~ClusterCounter() {
}
void ClusteringValidation::ClusterCounter::AddHitPreClustering(TrackID trackID) { ++numHitsPreClustering[trackID]; }
void ClusteringValidation::ClusterCounter::AddSignalHitPostClustering(ClusterID clusID) { ++numSignalHitsPostClustering[clusID]; }
void ClusteringValidation::ClusterCounter::AddNoiseHitPostClustering(ClusterID clusID) { ++numNoiseHitsPostClustering[clusID]; }
void ClusteringValidation::ClusterCounter::AssociateClusterAndTrack(ClusterID clusID, TrackID trackID) { clusterToTrackID[clusID] = trackID; trackToClusterIDs[trackID].push_back(clusID); }
double ClusteringValidation::ClusterCounter::GetCompleteness(ClusterID clusID) { return (double)numSignalHitsPostClustering[clusID]/(double)numHitsPreClustering[clusterToTrackID[clusID]]; }
double ClusteringValidation::ClusterCounter::GetCleanliness(ClusterID clusID) { return (double)numSignalHitsPostClustering[clusID]/(double)(GetNumberHitsInCluster(clusID)); }
double ClusteringValidation::ClusterCounter::GetEfficiency(TrackID trackID) { return 1/(double)trackToClusterIDs.at(trackID).size(); }
int ClusteringValidation::ClusterCounter::GetNumberHitsFromTrack(TrackID trackID) { return numHitsPreClustering[trackID]; }
int ClusteringValidation::ClusterCounter::GetNumberHitsInCluster(ClusterID clusID) { return numSignalHitsPostClustering[clusID] + numNoiseHitsPostClustering[clusID]; }
ClusterIDs ClusteringValidation::ClusterCounter::GetListOfClusterIDs() { ClusterIDs v; for (std::map<ClusterID,TrackID>::iterator i = clusterToTrackID.begin(); i != clusterToTrackID.end(); i++) v.push_back(i->first); return v; }
TrackIDs ClusteringValidation::ClusterCounter::GetListOfTrackIDs() { TrackIDs v; for (std::map<TrackID,ClusterIDs>::iterator i = trackToClusterIDs.begin(); i != trackToClusterIDs.end(); i++) v.push_back(i->first); return v; }
std::vector<std::pair<TrackID,ClusterIDs> > ClusteringValidation::ClusterCounter::GetPhotons() {
  std::vector<std::pair<TrackID,ClusterIDs> > photonVector;
  for (unsigned int track = 0; track < GetListOfTrackIDs().size(); ++track)
    if (!IsNoise(GetListOfTrackIDs().at(track)) && backtracker->TrackIDToParticle((int)GetListOfTrackIDs().at(track))->PdgCode() == 22)
      photonVector.push_back(std::pair<TrackID,ClusterIDs>(GetListOfTrackIDs().at(track),trackToClusterIDs.at(GetListOfTrackIDs().at(track))));
  return photonVector;
}
TrackID ClusteringValidation::ClusterCounter::GetTrack(ClusterID id) { return clusterToTrackID.at(id); }
bool ClusteringValidation::ClusterCounter::IsNoise(ClusterID clusID) { return IsNoise(clusterToTrackID.at(clusID)); }
bool ClusteringValidation::ClusterCounter::IsNoise(TrackID trackID) { return (int)trackID == 0 ? true : false; }
bool ClusteringValidation::ClusterCounter::PassesCut() {
  if (GetPhotons().size() > 2 || GetPhotons().size() == 0) return false;
  TrackIDs goodPhotons;
  for (unsigned int photon = 0; photon < GetPhotons().size(); ++photon)
    for (unsigned int cluster = 0; cluster < GetPhotons().at(photon).second.size(); ++cluster)
      if (GetCompleteness(GetPhotons().at(photon).second.at(cluster)) > 0.5) goodPhotons.push_back(GetPhotons().at(photon).first);
  bool pass = ( (GetPhotons().size() == 2 && goodPhotons.size() == 2) || (GetPhotons().size() == 1 && goodPhotons.size() == 1) );
  return pass;
}

ClusteringValidation::ClusteringValidation::ClusteringValidation(fhicl::ParameterSet const &pset) :  EDAnalyzer(pset) {
  this->reconfigure(pset);
  fCanvas = new TCanvas("fCanvas","",800,600);
}

ClusteringValidation::ClusteringValidation::~ClusteringValidation() { }

void ClusteringValidation::ClusteringValidation::reconfigure(fhicl::ParameterSet const& p) {
  fClusterModuleLabel = p.get<std::string>("ClusterModuleLabel");
  fHitsModuleLabel    = p.get<std::string>("HitsModuleLabel");
}

void ClusteringValidation::ClusteringValidation::analyze(art::Event const &evt)
{
  // Get the hits from the event
  art::Handle<std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if (evt.getByLabel(fHitsModuleLabel,hitHandle))
    art::fill_ptr_vector(hits, hitHandle);

  // Get the clusters from the event
  art::Handle<std::vector<recob::Cluster> > clusterHandle;
  std::vector<art::Ptr<recob::Cluster> > clusters;
  if (evt.getByLabel(fClusterModuleLabel,clusterHandle))
    art::fill_ptr_vector(clusters, clusterHandle);

  // Find the associations (the hits) for the clusters
  art::FindManyP<recob::Hit> fmh(clusterHandle,evt,fClusterModuleLabel);

  // Make a map of cluster counters in TPC/plane space
  for (unsigned int tpc = 0; tpc < geometry->NTPC(0); ++tpc) {
    for (unsigned int plane = 0; plane < geometry->Nplanes(tpc,0); ++plane) {
      clusterMap[tpc][plane] = new ClusterCounter(tpc, plane);
    }
  }

  // Loop over the preclustered hits first to save the info
  for (size_t hitIt = 0; hitIt < hits.size(); ++hitIt) {
    art::Ptr<recob::Hit> hit = hits.at(hitIt);
    TrackID trackID = FindTrackID(hit);
    clusterMap[hit->WireID().TPC][hit->WireID().Plane]->AddHitPreClustering(trackID);
  }

  // Look at true particles in the event
  trueParticles.clear();
  const sim::ParticleList& particles = backtracker->ParticleList();
  for (sim::ParticleList::const_iterator particleIt = particles.begin(); particleIt != particles.end(); ++particleIt) {
    const simb::MCParticle *particle = particleIt->second;
    trueParticles[(TrackID)particle->TrackId()] = particle;
    //std::cout << "True particle " << particle->TrackId() << " with pdg " << particle->PdgCode() << " which came from a " << particle->Mother() << " and process " << particle->Process() << std::endl;
  }

  // Loop over clusters
  for (size_t clusIt = 0; clusIt < clusters.size(); ++clusIt) {

    // Get cluster information
    unsigned int tpc   = clusters.at(clusIt)->Plane().TPC;
    unsigned int plane = clusters.at(clusIt)->Plane().Plane;
    ClusterID    id    = (ClusterID)clusters.at(clusIt)->ID();

    // Get the hits from the cluster
    std::vector<art::Ptr<recob::Hit> > clusterHits = fmh.at(clusIt);

    // Find which track this cluster belongs to
    TrackID trueTrackID = FindTrueTrack(clusterHits);

    // Save the info for this cluster
    clusterMap[tpc][plane]->AssociateClusterAndTrack(id, trueTrackID);
    for (std::vector<art::Ptr<recob::Hit> >::iterator clusHitIt = clusterHits.begin(); clusHitIt != clusterHits.end(); ++clusHitIt) {
      art::Ptr<recob::Hit> hit = *clusHitIt;
      TrackID trackID = FindTrackID(hit);
      if (trackID == trueTrackID) clusterMap[tpc][plane]->AddSignalHitPostClustering(id);
      else                        clusterMap[tpc][plane]->AddNoiseHitPostClustering (id);
    }

  } // cluster loop

  // Make plots
  for (unsigned int tpc = 0; tpc < geometry->NTPC(0); ++tpc) {
    for (unsigned int plane = 0; plane < geometry->Nplanes(tpc,0); ++plane) {
      ClusterIDs clusterIDs = clusterMap[tpc][plane]->GetListOfClusterIDs();
      if (clusterMap[tpc][plane]->GetPhotons().size() == 2) {
	hPi0Angle               ->Fill(FindPhotonAngle());
	hPi0Energy              ->Fill(GetPi0()->Momentum().E());
	hPi0ConversionDistance  ->Fill(std::min(GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(0).first, (TrackID)GetPi0()->TrackId()),
						GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(1).first, (TrackID)GetPi0()->TrackId())));
	hPi0ConversionSeparation->Fill(GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(0).first, clusterMap[tpc][plane]->GetPhotons().at(1).first));
	if (clusterMap[tpc][plane]->PassesCut()) {
	  hPi0AngleCut               ->Fill(FindPhotonAngle());
	  hPi0EnergyCut              ->Fill(GetPi0()->Momentum().E());
	  hPi0ConversionDistanceCut  ->Fill(std::min(GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(0).first, (TrackID)GetPi0()->TrackId()),
						     GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(1).first, (TrackID)GetPi0()->TrackId())));
	  hPi0ConversionSeparationCut->Fill(GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(0).first, clusterMap[tpc][plane]->GetPhotons().at(1).first));
	}
      }
      for (unsigned int cluster = 0; cluster < clusterIDs.size(); ++cluster) {
	ClusterID clusID = clusterIDs.at(cluster); ++nClusters;
	hCompleteness                  ->Fill(clusterMap[tpc][plane]->GetCompleteness(clusID)); totCompleteness += clusterMap[tpc][plane]->GetCompleteness(clusID);
	hCleanliness                   ->Fill(clusterMap[tpc][plane]->GetCleanliness (clusID)); totCleanliness  += clusterMap[tpc][plane]->GetCleanliness (clusID);
	hNumHitsCompleteness           ->Fill(clusterMap[tpc][plane]->GetCompleteness(clusID), clusterMap[tpc][plane]->GetNumberHitsInCluster(clusID));
	if (clusterMap[tpc][plane]->IsNoise(clusID)) continue;
	hCompletenessEnergy            ->Fill(GetPi0()->Momentum().E(),                                                                    clusterMap[tpc][plane]->GetCompleteness       (clusID));
	hCompletenessAngle             ->Fill(FindPhotonAngle(),                                                                           clusterMap[tpc][plane]->GetCompleteness       (clusID));
	hCompletenessConversionDistance->Fill(GetEndTrackDistance(clusterMap[tpc][plane]->GetTrack(clusID), (TrackID)GetPi0()->TrackId()), clusterMap[tpc][plane]->GetCompleteness       (clusID));
	hCleanlinessEnergy             ->Fill(GetPi0()->Momentum().E(),                                                                    clusterMap[tpc][plane]->GetCleanliness        (clusID));
	hCleanlinessAngle              ->Fill(FindPhotonAngle(),                                                                           clusterMap[tpc][plane]->GetCleanliness        (clusID));
	hCleanlinessConversionDistance ->Fill(GetEndTrackDistance(clusterMap[tpc][plane]->GetTrack(clusID), (TrackID)GetPi0()->TrackId()), clusterMap[tpc][plane]->GetCleanliness        (clusID));
	hNumHitsEnergy                 ->Fill(GetPi0()->Momentum().E(),                                                                    clusterMap[tpc][plane]->GetNumberHitsInCluster(clusID));
	if (clusterMap[tpc][plane]->GetPhotons().size() != 2) continue;
	hCompletenessConversionSeparation->Fill(GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(0).first, clusterMap[tpc][plane]->GetPhotons().at(1).first), clusterMap[tpc][plane]->GetCompleteness(clusID));
	hCleanlinessConversionSeparation ->Fill(GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(0).first, clusterMap[tpc][plane]->GetPhotons().at(1).first), clusterMap[tpc][plane]->GetCleanliness (clusID));
      }
    }
  }
}

TrackID ClusteringValidation::ClusteringValidation::FindTrackID(art::Ptr<recob::Hit> &hit) {
  double particleEnergy = 0;
  TrackID likelyTrackID = (TrackID)0;
  std::vector<sim::TrackIDE> trackIDs = backtracker->HitToTrackID(hit);
  for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
    if (trackIDs.at(idIt).energy > particleEnergy) {
      particleEnergy = trackIDs.at(idIt).energy;
      likelyTrackID = (TrackID)TMath::Abs(trackIDs.at(idIt).trackID);
    }
  }
  return likelyTrackID;
}

TrackID ClusteringValidation::ClusteringValidation::FindTrueTrack(std::vector<art::Ptr<recob::Hit> > &clusterHits) {
  std::map<TrackID,double> trackMap;
  for (std::vector<art::Ptr<recob::Hit> >::iterator clusHitIt = clusterHits.begin(); clusHitIt != clusterHits.end(); ++clusHitIt) {
    art::Ptr<recob::Hit> hit = *clusHitIt;
    TrackID trackID = FindTrackID(hit);
    trackMap[trackID] += hit->Integral();
  }
  //return std::max_element(trackMap.begin(), trackMap.end(), [](const std::pair<int,double>& p1, const std::pair<int,double>& p2) {return p1.second < p2.second;} )->first;
  double highestCharge = 0;
  TrackID clusterTrack = (TrackID)0;
  for (std::map<TrackID,double>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt)
    if (trackIt->second > highestCharge) {
      highestCharge = trackIt->second;
      clusterTrack  = trackIt->first;
    }
  return clusterTrack;
}

double ClusteringValidation::ClusteringValidation::FindPhotonAngle() {
  const simb::MCParticle* pi0 = GetPi0();
  if (pi0->NumberDaughters() != 2) return -999;
  double angle = (trueParticles.at((TrackID)pi0->Daughter(0))->Momentum().Angle(trueParticles.at((TrackID)pi0->Daughter(1))->Momentum().Vect()) * (180/TMath::Pi()));
  return angle;
}

const simb::MCParticle* ClusteringValidation::ClusteringValidation::GetPi0() {
  const simb::MCParticle* pi0 = nullptr;
  for (std::map<TrackID,const simb::MCParticle*>::iterator particleIt = trueParticles.begin(); particleIt != trueParticles.end(); ++particleIt)
    if (particleIt->second->PdgCode() == 111)
      pi0 = particleIt->second;
  return pi0;
}

double ClusteringValidation::ClusteringValidation::GetEndTrackDistance(TrackID id1, TrackID id2) {
  return TMath::Sqrt(TMath::Power(trueParticles.at(id1)->EndPosition().X() - trueParticles.at(id2)->EndPosition().X(),2)+
                     TMath::Power(trueParticles.at(id1)->EndPosition().Y() - trueParticles.at(id2)->EndPosition().Y(),2)+
		     TMath::Power(trueParticles.at(id1)->EndPosition().Z() - trueParticles.at(id2)->EndPosition().Z(),2));
}

void ClusteringValidation::ClusteringValidation::beginJob() {
  hCompleteness                     = new TH1D("Completeness",";Completeness;",101,0,1.01);
  hCompletenessEnergy               = new TProfile("CompletenessEnergy",";True Energy (GeV);Completeness",100,0,2);
  hCompletenessAngle                = new TProfile("CompletenessAngle",";True Angle (deg);Completeness;",100,0,180);
  hCompletenessConversionDistance   = new TProfile("CompletenessConversionDistance",";True Distance from Vertex (cm);Completeness",100,0,200);
  hCompletenessConversionSeparation = new TProfile("CompletenessConversionSeparation",";True Conversion Separation (cm);Completeness",100,0,200);
  hCleanliness                      = new TH1D("Cleanliness",";Cleanliness;",101,0,1.01);
  hCleanlinessEnergy                = new TProfile("CleanlinessEnergy",";True Energy (GeV);Cleanliness",100,0,2);
  hCleanlinessAngle                 = new TProfile("CleanlinessAngle",";True Angle (deg);Cleanliness;",100,0,180);
  hCleanlinessConversionDistance    = new TProfile("CleanlinessConversionDistance",";True Distance from Vertex (cm);Cleanliness",100,0,200);
  hCleanlinessConversionSeparation  = new TProfile("CleanlinessConversionSeparation",";True Conversion Separation (cm);Cleanliness",100,0,200);
  hPi0Energy                        = new TH1D("Pi0EnergyCut",";True Energy (GeV);",25,0,2);                                          hPi0Energy                 ->Sumw2();
  hPi0Angle                         = new TH1D("Pi0AngleCut",";True Angle (deg);",25,0,180);                                          hPi0Angle                  ->Sumw2();
  hPi0ConversionDistance            = new TH1D("Pi0ConversionDistanceCut",";True Distance from Vertex (cm);",25,0,200);               hPi0ConversionDistance     ->Sumw2();
  hPi0ConversionSeparation          = new TH1D("Pi0ConversionSeparationCut",";True Separation from Vertex (cm);",25,0,200);           hPi0ConversionSeparation   ->Sumw2();
  hPi0EnergyCut                     = new TH1D("Pi0EnergyCut",";True Energy (GeV);Efficiency",25,0,2);                                hPi0EnergyCut              ->Sumw2();
  hPi0AngleCut                      = new TH1D("Pi0AngleCut",";True Angle (deg);Efficiency",25,0,180);                                hPi0AngleCut               ->Sumw2();
  hPi0ConversionDistanceCut         = new TH1D("Pi0ConversionDistanceCut",";True Distance from Vertex (cm);Efficiency",25,0,200);     hPi0ConversionDistanceCut  ->Sumw2();
  hPi0ConversionSeparationCut       = new TH1D("Pi0ConversionSeparationCut",";True Separation from Vertex (cm);Efficiency",25,0,200); hPi0ConversionSeparationCut->Sumw2();
  hNumHitsCompleteness              = new TH2D("NumHitsCompleteness",";Completeness;Size of Cluster",101,0,1.01,100,0,100);
  hNumHitsEnergy                    = new TH2D("NumHitsEnergy",";True Energy (GeV);Size of Cluster",100,0,2,100,0,100);
}

void ClusteringValidation::ClusteringValidation::endJob()
{
  hEfficiencyEnergy               = new TEfficiency(*hPi0EnergyCut,*hPi0Energy);
  hEfficiencyAngle                = new TEfficiency(*hPi0AngleCut,*hPi0Angle);
  hEfficiencyConversionDistance   = new TEfficiency(*hPi0ConversionDistanceCut,*hPi0ConversionDistance);
  hEfficiencyConversionSeparation = new TEfficiency(*hPi0ConversionSeparationCut,*hPi0ConversionSeparation);

  hEfficiencyEnergy              ->SetName("EfficiencyEnergy");
  hEfficiencyAngle               ->SetName("EnergyAngle");
  hEfficiencyConversionDistance  ->SetName("EfficiencyConversionDistance");
  hEfficiencyConversionSeparation->SetName("EfficiencyConversionSeparation");

  fHistArray.Add(hCompleteness);     fHistArray.Add(hCompletenessEnergy); fHistArray.Add(hCompletenessAngle); fHistArray.Add(hCompletenessConversionDistance); fHistArray.Add(hCompletenessConversionSeparation);
  fHistArray.Add(hCleanliness);      fHistArray.Add(hCleanlinessEnergy);  fHistArray.Add(hCleanlinessAngle);  fHistArray.Add(hCleanlinessConversionDistance);  fHistArray.Add(hCleanlinessConversionSeparation);
  fHistArray.Add(hEfficiencyEnergy); fHistArray.Add(hEfficiencyAngle);    fHistArray.Add(hEfficiencyConversionDistance); fHistArray.Add(hEfficiencyConversionSeparation);
  fHistArray.Add(hNumHitsCompleteness); fHistArray.Add(hNumHitsEnergy);

  for (int histIt = 0; histIt < fHistArray.GetEntriesFast(); ++histIt) {
    fCanvas->cd();
    TH1* h = (TH1*)fHistArray.At(histIt);
    h->Draw();
    fCanvas->SaveAs(h->GetName()+TString(".png"));
  }

  // Average completeness/cleanliness
  double avCompleteness = totCompleteness / nClusters;
  double avCleanliness  = totCleanliness  / nClusters;

  // Write file
  ofstream outFile("effpur");
  outFile << avCompleteness << " " << avCleanliness;
  outFile.close();

}

DEFINE_ART_MODULE(ClusteringValidation::ClusteringValidation)
