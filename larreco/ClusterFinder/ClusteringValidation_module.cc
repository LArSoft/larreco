//////////////////////////////////////////////////////////////////////////////////////
// Class:       ClusteringValidation
// Module type: analyser
// File:        ClusteringValidation_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffied.ac.uk), May 2015
//
// A module to validate clustering algorithms.
// Compares the output of different clustering algorithms run over a pi0 sample.
//
// Usage: Specify the hit finder (HitsModuleLabel) and the clustering outputs
// to validate (ClusterModuleLabels) in the fhicl file.
// Module will make validation plots for all clusterings specified and also
// produce comparison plots. Number of clusterings to analyse can be one or more.
// Saves everything in the file validationHistograms.root.
//////////////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "fhiclcpp/ParameterSet.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ROOT & STL includes
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TStyle.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>

namespace ClusteringValidation {
  class ClusteringValidation;
  class ClusterAnalyser;
  class ClusterCounter;
}

enum class ClusterID : int {};
enum class TrackID : int {};
typedef std::vector<ClusterID> ClusterIDs;
typedef std::vector<TrackID> TrackIDs;

class ClusteringValidation::ClusterCounter {
public:
  explicit ClusterCounter(unsigned int& tpc, unsigned int& plane);

  void AddHitPreClustering(TrackID id);
  void AddSignalHitPostClustering(ClusterID id);
  void AddNoiseHitPostClustering(ClusterID id);
  void AssociateClusterAndTrack(ClusterID clusID, TrackID trackID);
  double GetCompleteness(ClusterID id);
  double GetCleanliness(ClusterID id);
  double GetEfficiency(TrackID id);
  ClusterIDs GetListOfClusterIDs();
  TrackIDs GetListOfTrackIDs();
  int GetNumberHitsFromTrack(TrackID id);
  int GetNumberHitsInCluster(ClusterID id);
  int GetNumberHitsInPlane();
  std::vector<std::pair<TrackID, ClusterIDs>> GetPhotons();
  TrackID GetTrack(ClusterID id);
  bool IsNoise(ClusterID id);
  bool IsNoise(TrackID id);
  bool PassesCut();

private:
  unsigned int tpc, plane;

  std::map<TrackID, int> numHitsPreClustering;
  std::map<ClusterID, int> numSignalHitsPostClustering;
  std::map<ClusterID, int> numNoiseHitsPostClustering;
  std::map<ClusterID, TrackID> clusterToTrackID;
  std::map<TrackID, ClusterIDs> trackToClusterIDs;
  std::map<TrackID, std::map<std::string, double>> particleProperties;
  std::map<TrackID, simb::MCParticle> trueParticles;

  art::ServiceHandle<geo::Geometry const> geometry;
  art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
};

ClusteringValidation::ClusterCounter::ClusterCounter(unsigned int& t, unsigned int& p)
{
  tpc = t;
  plane = p;
}

void ClusteringValidation::ClusterCounter::AddHitPreClustering(TrackID trackID)
{
  ++numHitsPreClustering[trackID];
}

void ClusteringValidation::ClusterCounter::AddSignalHitPostClustering(ClusterID clusID)
{
  ++numSignalHitsPostClustering[clusID];
}

void ClusteringValidation::ClusterCounter::AddNoiseHitPostClustering(ClusterID clusID)
{
  ++numNoiseHitsPostClustering[clusID];
}

void ClusteringValidation::ClusterCounter::AssociateClusterAndTrack(ClusterID clusID,
                                                                    TrackID trackID)
{
  clusterToTrackID[clusID] = trackID;
  trackToClusterIDs[trackID].push_back(clusID);
}

double ClusteringValidation::ClusterCounter::GetCompleteness(ClusterID clusID)
{
  return (double)numSignalHitsPostClustering[clusID] /
         (double)numHitsPreClustering[clusterToTrackID[clusID]];
}

double ClusteringValidation::ClusterCounter::GetCleanliness(ClusterID clusID)
{
  return (double)numSignalHitsPostClustering[clusID] / (double)(GetNumberHitsInCluster(clusID));
}

double ClusteringValidation::ClusterCounter::GetEfficiency(TrackID trackID)
{
  return 1 / (double)trackToClusterIDs.at(trackID).size();
}

int ClusteringValidation::ClusterCounter::GetNumberHitsFromTrack(TrackID trackID)
{
  return numHitsPreClustering[trackID];
}

int ClusteringValidation::ClusterCounter::GetNumberHitsInCluster(ClusterID clusID)
{
  return numSignalHitsPostClustering[clusID] + numNoiseHitsPostClustering[clusID];
}

int ClusteringValidation::ClusterCounter::GetNumberHitsInPlane()
{
  int nHits = 0;
  for (auto& trackHits : numHitsPreClustering)
    nHits += trackHits.second;
  return nHits;
}

ClusterIDs ClusteringValidation::ClusterCounter::GetListOfClusterIDs()
{
  ClusterIDs v;
  for (std::map<ClusterID, TrackID>::iterator i = clusterToTrackID.begin();
       i != clusterToTrackID.end();
       i++)
    v.push_back(i->first);
  return v;
}

TrackIDs ClusteringValidation::ClusterCounter::GetListOfTrackIDs()
{
  TrackIDs v;
  for (std::map<TrackID, ClusterIDs>::iterator i = trackToClusterIDs.begin();
       i != trackToClusterIDs.end();
       i++)
    v.push_back(i->first);
  return v;
}

std::vector<std::pair<TrackID, ClusterIDs>> ClusteringValidation::ClusterCounter::GetPhotons()
{
  std::vector<std::pair<TrackID, ClusterIDs>> photonVector;
  for (unsigned int track = 0; track < GetListOfTrackIDs().size(); ++track)
    if (!IsNoise(GetListOfTrackIDs().at(track)) &&
        pi_serv->TrackIdToParticle_P((int)GetListOfTrackIDs().at(track))->PdgCode() == 22)
      photonVector.push_back(std::pair<TrackID, ClusterIDs>(
        GetListOfTrackIDs().at(track), trackToClusterIDs.at(GetListOfTrackIDs().at(track))));
  return photonVector;
}

TrackID ClusteringValidation::ClusterCounter::GetTrack(ClusterID id)
{
  return clusterToTrackID.at(id);
}

bool ClusteringValidation::ClusterCounter::IsNoise(ClusterID clusID)
{
  return IsNoise(clusterToTrackID.at(clusID));
}

bool ClusteringValidation::ClusterCounter::IsNoise(TrackID trackID)
{
  return (int)trackID == 0 ? true : false;
}

bool ClusteringValidation::ClusterCounter::PassesCut()
{
  if (GetPhotons().size() > 2 || GetPhotons().size() == 0) return false;
  TrackIDs goodPhotons;
  for (unsigned int photon = 0; photon < GetPhotons().size(); ++photon)
    for (unsigned int cluster = 0; cluster < GetPhotons().at(photon).second.size(); ++cluster)
      if (GetCompleteness(GetPhotons().at(photon).second.at(cluster)) > 0.5)
        goodPhotons.push_back(GetPhotons().at(photon).first);
  if ((GetPhotons().size() == 2 && goodPhotons.size() > 2) ||
      (GetPhotons().size() == 1 && goodPhotons.size() > 1))
    std::cout << "More than 2 with >50%?!" << std::endl;
  bool pass = ((GetPhotons().size() == 2 && goodPhotons.size() == 2) ||
               (GetPhotons().size() == 1 && goodPhotons.size() == 1));
  return pass;
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class ClusteringValidation::ClusterAnalyser {
public:
  explicit ClusterAnalyser(std::string& label);

  void Analyse(detinfo::DetectorClocksData const& clockData,
               std::vector<art::Ptr<recob::Hit>>& hits,
               std::vector<art::Ptr<recob::Cluster>>& clusters,
               const art::FindManyP<recob::Hit>& fmh,
               int numHits);
  TrackID FindTrackID(detinfo::DetectorClocksData const& clockData, art::Ptr<recob::Hit>& hit);
  TrackID FindTrueTrack(detinfo::DetectorClocksData const& clockData,
                        std::vector<art::Ptr<recob::Hit>>& clusterHits);
  double FindPhotonAngle();
  double GetEndTrackDistance(TrackID id1, TrackID id2);
  const simb::MCParticle* GetPi0();
  TObjArray GetHistograms();
  void MakeHistograms();
  void WriteFile();

private:
  // Clustering properties
  std::string fClusterLabel;

  // hists
  TH1 *hCompleteness, *hCleanliness, *hComplCleanl;
  TH1 *hPi0Angle, *hPi0Energy, *hPi0ConversionDistance, *hPi0ConversionSeparation, *hPi0AngleCut,
    *hPi0EnergyCut, *hPi0ConversionDistanceCut, *hPi0ConversionSeparationCut;
  TH2 *hNumHitsCompleteness, *hNumHitsEnergy;
  TProfile *hCompletenessEnergy, *hCompletenessAngle, *hCompletenessConversionDistance,
    *hCompletenessConversionSeparation;
  TProfile *hCleanlinessEnergy, *hCleanlinessAngle, *hCleanlinessConversionDistance,
    *hCleanlinessConversionSeparation;
  TProfile *hComplCleanlEnergy, *hComplCleanlAngle, *hComplCleanlConversionDistance,
    *hComplCleanlConversionSeparation;
  TEfficiency *hEfficiencyAngle, *hEfficiencyEnergy, *hEfficiencyConversionDistance,
    *hEfficiencyConversionSeparation;
  TObjArray fHistArray;

  std::map<unsigned int, std::map<unsigned int, std::unique_ptr<ClusterCounter>>> clusterMap;
  std::map<TrackID, const simb::MCParticle*> trueParticles;

  // Services
  art::ServiceHandle<geo::Geometry const> geometry;
  art::ServiceHandle<cheat::ParticleInventoryService const> pi_serv;
  art::ServiceHandle<cheat::BackTrackerService const> bt_serv;
};

ClusteringValidation::ClusterAnalyser::ClusterAnalyser(std::string& clusterLabel)
{

  fClusterLabel = clusterLabel;

  // Make the histograms
  hCompleteness = new TH1D("Completeness", ";Completeness;", 101, 0, 1.01);
  hCompletenessEnergy =
    new TProfile("CompletenessEnergy", ";True Energy (GeV);Completeness", 100, 0, 2);
  hCompletenessAngle =
    new TProfile("CompletenessAngle", ";True Angle (deg);Completeness;", 100, 0, 180);
  hCompletenessConversionDistance = new TProfile(
    "CompletenessConversionDistance", ";True Distance from Vertex (cm);Completeness", 100, 0, 200);
  hCompletenessConversionSeparation = new TProfile("CompletenessConversionSeparation",
                                                   ";True Conversion Separation (cm);Completeness",
                                                   100,
                                                   0,
                                                   200);
  hCleanliness = new TH1D("Cleanliness", ";Cleanliness;", 101, 0, 1.01);
  hCleanlinessEnergy =
    new TProfile("CleanlinessEnergy", ";True Energy (GeV);Cleanliness", 100, 0, 2);
  hCleanlinessAngle =
    new TProfile("CleanlinessAngle", ";True Angle (deg);Cleanliness;", 100, 0, 180);
  hCleanlinessConversionDistance = new TProfile(
    "CleanlinessConversionDistance", ";True Distance from Vertex (cm);Cleanliness", 100, 0, 200);
  hCleanlinessConversionSeparation = new TProfile(
    "CleanlinessConversionSeparation", ";True Conversion Separation (cm);Cleanliness", 100, 0, 200);
  hComplCleanl = new TH1D("CompletenessCleanliness", ";Completeness * Cleanliness;", 101, 0, 1.01);
  hComplCleanlEnergy = new TProfile(
    "CompletenessCleanlinessEnergy", ";True Energy (GeV);Completeness * Cleanliness", 100, 0, 2);
  hComplCleanlAngle = new TProfile(
    "CompletenessCleanlinessAngle", ";True Angle (deg);Completeness * Cleanliness;", 100, 0, 180);
  hComplCleanlConversionDistance =
    new TProfile("CompletenessCleanlinessConversionDistance",
                 ";True Distance from Vertex (cm);Completeness * Cleanliness",
                 100,
                 0,
                 200);
  hComplCleanlConversionSeparation =
    new TProfile("CompletenessCleanlinessConversionSeparation",
                 ";True Conversion Separation (cm);Completeness * Cleanliness",
                 100,
                 0,
                 200);
  hPi0Energy = new TH1D("Pi0EnergyCut", ";True Energy (GeV);", 25, 0, 2);
  hPi0Energy->Sumw2();
  hPi0Angle = new TH1D("Pi0AngleCut", ";True Angle (deg);", 25, 0, 180);
  hPi0Angle->Sumw2();
  hPi0ConversionDistance =
    new TH1D("Pi0ConversionDistanceCut", ";True Distance from Vertex (cm);", 25, 0, 200);
  hPi0ConversionDistance->Sumw2();
  hPi0ConversionSeparation =
    new TH1D("Pi0ConversionSeparationCut", ";True Separation from Vertex (cm);", 25, 0, 200);
  hPi0ConversionSeparation->Sumw2();
  hPi0EnergyCut = new TH1D("Pi0EnergyCut", ";True Energy (GeV);Efficiency", 25, 0, 2);
  hPi0EnergyCut->Sumw2();
  hPi0AngleCut = new TH1D("Pi0AngleCut", ";True Angle (deg);Efficiency", 25, 0, 180);
  hPi0AngleCut->Sumw2();
  hPi0ConversionDistanceCut =
    new TH1D("Pi0ConversionDistanceCut", ";True Distance from Vertex (cm);Efficiency", 25, 0, 200);
  hPi0ConversionDistanceCut->Sumw2();
  hPi0ConversionSeparationCut = new TH1D(
    "Pi0ConversionSeparationCut", ";True Separation from Vertex (cm);Efficiency", 25, 0, 200);
  hPi0ConversionSeparationCut->Sumw2();
  hNumHitsCompleteness =
    new TH2D("NumHitsCompleteness", ";Completeness;Size of Cluster", 101, 0, 1.01, 100, 0, 100);
  hNumHitsEnergy =
    new TH2D("NumHitsEnergy", ";True Energy (GeV);Size of Cluster", 100, 0, 2, 100, 0, 100);
}

void ClusteringValidation::ClusterAnalyser::Analyse(detinfo::DetectorClocksData const& clockData,
                                                    std::vector<art::Ptr<recob::Hit>>& hits,
                                                    std::vector<art::Ptr<recob::Cluster>>& clusters,
                                                    const art::FindManyP<recob::Hit>& fmh,
                                                    int minHits)
{

  // Make a map of cluster counters in TPC/plane space
  for (unsigned int tpc = 0; tpc < geometry->NTPC(0); ++tpc) {
    for (unsigned int plane = 0; plane < geometry->Nplanes(tpc, 0); ++plane) {
      clusterMap[tpc][plane] = (std::unique_ptr<ClusterCounter>)new ClusterCounter(tpc, plane);
    }
  }

  // Save preclustered hits
  for (size_t hitIt = 0; hitIt < hits.size(); ++hitIt) {
    art::Ptr<recob::Hit> hit = hits.at(hitIt);
    TrackID trackID = FindTrackID(clockData, hit);
    clusterMap[hit->WireID().TPC % 2][hit->WireID().Plane]->AddHitPreClustering(trackID);
  }

  // Save true tracks
  trueParticles.clear();
  const sim::ParticleList& particles = pi_serv->ParticleList();
  for (sim::ParticleList::const_iterator particleIt = particles.begin();
       particleIt != particles.end();
       ++particleIt) {
    const simb::MCParticle* particle = particleIt->second;
    trueParticles[(TrackID)particle->TrackId()] = particle;
  }

  // Save the clustered hits
  for (size_t clusIt = 0; clusIt < clusters.size(); ++clusIt) {

    // Get cluster information
    unsigned int tpc = clusters.at(clusIt)->Plane().TPC % 2;
    unsigned int plane = clusters.at(clusIt)->Plane().Plane;
    ClusterID id = (ClusterID)clusters.at(clusIt)->ID();

    // Only analyse planes with enough hits in to fairly represent the clustering
    if (clusterMap[tpc][plane]->GetNumberHitsInPlane() < minHits) continue;

    // Get the hits from the cluster
    std::vector<art::Ptr<recob::Hit>> clusterHits = fmh.at(clusIt);

    if (clusterHits.size() < 10) continue;

    // Find which track this cluster belongs to
    TrackID trueTrackID = FindTrueTrack(clockData, clusterHits);

    // Save the info for this cluster
    clusterMap[tpc][plane]->AssociateClusterAndTrack(id, trueTrackID);
    for (std::vector<art::Ptr<recob::Hit>>::iterator clusHitIt = clusterHits.begin();
         clusHitIt != clusterHits.end();
         ++clusHitIt) {
      art::Ptr<recob::Hit> hit = *clusHitIt;
      TrackID trackID = FindTrackID(clockData, hit);
      if (trackID == trueTrackID)
        clusterMap[tpc][plane]->AddSignalHitPostClustering(id);
      else
        clusterMap[tpc][plane]->AddNoiseHitPostClustering(id);
    }
  }

  this->MakeHistograms();
}

TrackID ClusteringValidation::ClusterAnalyser::FindTrackID(
  detinfo::DetectorClocksData const& clockData,
  art::Ptr<recob::Hit>& hit)
{
  double particleEnergy = 0;
  TrackID likelyTrackID = (TrackID)0;
  std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
  for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
    if (trackIDs.at(idIt).energy > particleEnergy) {
      particleEnergy = trackIDs.at(idIt).energy;
      likelyTrackID = (TrackID)TMath::Abs(trackIDs.at(idIt).trackID);
    }
  }
  return likelyTrackID;
}

TrackID ClusteringValidation::ClusterAnalyser::FindTrueTrack(
  detinfo::DetectorClocksData const& clockData,
  std::vector<art::Ptr<recob::Hit>>& clusterHits)
{
  std::map<TrackID, double> trackMap;
  for (std::vector<art::Ptr<recob::Hit>>::iterator clusHitIt = clusterHits.begin();
       clusHitIt != clusterHits.end();
       ++clusHitIt) {
    art::Ptr<recob::Hit> hit = *clusHitIt;
    TrackID trackID = FindTrackID(clockData, hit);
    trackMap[trackID] += hit->Integral();
  }
  //return std::max_element(trackMap.begin(), trackMap.end(), [](const std::pair<int,double>& p1, const std::pair<int,double>& p2) {return p1.second < p2.second;} )->first;
  double highestCharge = 0;
  TrackID clusterTrack = (TrackID)0;
  for (std::map<TrackID, double>::iterator trackIt = trackMap.begin(); trackIt != trackMap.end();
       ++trackIt)
    if (trackIt->second > highestCharge) {
      highestCharge = trackIt->second;
      clusterTrack = trackIt->first;
    }
  return clusterTrack;
}

double ClusteringValidation::ClusterAnalyser::FindPhotonAngle()
{
  const simb::MCParticle* pi0 = GetPi0();
  if (pi0->NumberDaughters() != 2) return -999;
  double angle = (trueParticles.at((TrackID)pi0->Daughter(0))
                    ->Momentum()
                    .Angle(trueParticles.at((TrackID)pi0->Daughter(1))->Momentum().Vect()) *
                  (180 / TMath::Pi()));
  return angle;
}

const simb::MCParticle* ClusteringValidation::ClusterAnalyser::GetPi0()
{
  const simb::MCParticle* pi0 = nullptr;
  for (std::map<TrackID, const simb::MCParticle*>::iterator particleIt = trueParticles.begin();
       particleIt != trueParticles.end();
       ++particleIt)
    if (particleIt->second->PdgCode() == 111) pi0 = particleIt->second;
  return pi0;
}

double ClusteringValidation::ClusterAnalyser::GetEndTrackDistance(TrackID id1, TrackID id2)
{
  return TMath::Sqrt(
    TMath::Power(
      trueParticles.at(id1)->EndPosition().X() - trueParticles.at(id2)->EndPosition().X(), 2) +
    TMath::Power(
      trueParticles.at(id1)->EndPosition().Y() - trueParticles.at(id2)->EndPosition().Y(), 2) +
    TMath::Power(
      trueParticles.at(id1)->EndPosition().Z() - trueParticles.at(id2)->EndPosition().Z(), 2));
}

TObjArray ClusteringValidation::ClusterAnalyser::GetHistograms()
{

  // Make efficiency histograms
  hEfficiencyEnergy = new TEfficiency(*hPi0EnergyCut, *hPi0Energy);
  hEfficiencyAngle = new TEfficiency(*hPi0AngleCut, *hPi0Angle);
  hEfficiencyConversionDistance =
    new TEfficiency(*hPi0ConversionDistanceCut, *hPi0ConversionDistance);
  hEfficiencyConversionSeparation =
    new TEfficiency(*hPi0ConversionSeparationCut, *hPi0ConversionSeparation);

  hEfficiencyEnergy->SetName("EfficiencyEnergy");
  hEfficiencyAngle->SetName("EnergyAngle");
  hEfficiencyConversionDistance->SetName("EfficiencyConversionDistance");
  hEfficiencyConversionSeparation->SetName("EfficiencyConversionSeparation");

  // Add all the histograms to the object array
  fHistArray.Add(hCompleteness);
  fHistArray.Add(hCompletenessEnergy);
  fHistArray.Add(hCompletenessAngle);
  fHistArray.Add(hCompletenessConversionDistance);
  fHistArray.Add(hCompletenessConversionSeparation);
  fHistArray.Add(hCleanliness);
  fHistArray.Add(hCleanlinessEnergy);
  fHistArray.Add(hCleanlinessAngle);
  fHistArray.Add(hCleanlinessConversionDistance);
  fHistArray.Add(hCleanlinessConversionSeparation);
  fHistArray.Add(hComplCleanl);
  fHistArray.Add(hComplCleanlEnergy);
  fHistArray.Add(hComplCleanlAngle);
  fHistArray.Add(hComplCleanlConversionDistance);
  fHistArray.Add(hComplCleanlConversionSeparation);
  fHistArray.Add(hEfficiencyEnergy);
  fHistArray.Add(hEfficiencyAngle);
  fHistArray.Add(hEfficiencyConversionDistance);
  fHistArray.Add(hEfficiencyConversionSeparation);
  fHistArray.Add(hNumHitsCompleteness);
  fHistArray.Add(hNumHitsEnergy);

  return fHistArray;
}

void ClusteringValidation::ClusterAnalyser::MakeHistograms()
{

  // Loop over the tpcs and planes in the geometry
  for (unsigned int tpc = 0; tpc < geometry->NTPC(0); ++tpc) {
    for (unsigned int plane = 0; plane < geometry->Nplanes(tpc, 0); ++plane) {

      ClusterIDs clusterIDs = clusterMap[tpc][plane]->GetListOfClusterIDs();

      // Fill histograms for the efficiency
      if (clusterMap[tpc][plane]->GetPhotons().size() == 2) {
        hPi0Angle->Fill(FindPhotonAngle());
        hPi0Energy->Fill(GetPi0()->Momentum().E());
        hPi0ConversionDistance->Fill(
          std::min(GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(0).first,
                                       (TrackID)GetPi0()->TrackId()),
                   GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(1).first,
                                       (TrackID)GetPi0()->TrackId())));
        hPi0ConversionSeparation->Fill(
          GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(0).first,
                              clusterMap[tpc][plane]->GetPhotons().at(1).first));
        if (clusterMap[tpc][plane]->PassesCut()) {
          hPi0AngleCut->Fill(FindPhotonAngle());
          hPi0EnergyCut->Fill(GetPi0()->Momentum().E());
          hPi0ConversionDistanceCut->Fill(
            std::min(GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(0).first,
                                         (TrackID)GetPi0()->TrackId()),
                     GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(1).first,
                                         (TrackID)GetPi0()->TrackId())));
          hPi0ConversionSeparationCut->Fill(
            GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(0).first,
                                clusterMap[tpc][plane]->GetPhotons().at(1).first));
        }
        else
          std::cout << "TPC " << tpc << ", Plane " << plane << " fails the cut" << std::endl;
      }

      // Look at all the clusters
      for (unsigned int cluster = 0; cluster < clusterIDs.size(); ++cluster) {

        ClusterID clusID = clusterIDs.at(cluster);
        double completeness = clusterMap[tpc][plane]->GetCompleteness(clusID);
        double cleanliness = clusterMap[tpc][plane]->GetCleanliness(clusID);
        int numClusterHits = clusterMap[tpc][plane]->GetNumberHitsInCluster(clusID);

        // Fill histograms for this cluster
        hCompleteness->Fill(completeness, numClusterHits);
        hCleanliness->Fill(cleanliness, numClusterHits);
        hComplCleanl->Fill(completeness * cleanliness, numClusterHits);
        hNumHitsCompleteness->Fill(completeness, numClusterHits);

        // Is this cluster doesn't correspond to a true particle continue
        if (clusterMap[tpc][plane]->IsNoise(clusID)) continue;

        double pi0Energy = GetPi0()->Momentum().E();
        double pi0DecayAngle = FindPhotonAngle();
        double conversionDistance = GetEndTrackDistance(clusterMap[tpc][plane]->GetTrack(clusID),
                                                        (TrackID)GetPi0()->TrackId());

        hCompletenessEnergy->Fill(pi0Energy, completeness, numClusterHits);
        hCompletenessAngle->Fill(pi0DecayAngle, completeness, numClusterHits);
        hCompletenessConversionDistance->Fill(conversionDistance, completeness, numClusterHits);
        hCleanlinessEnergy->Fill(pi0Energy, cleanliness, numClusterHits);
        hCleanlinessAngle->Fill(pi0DecayAngle, cleanliness, numClusterHits);
        hCleanlinessConversionDistance->Fill(conversionDistance, cleanliness, numClusterHits);
        hComplCleanlEnergy->Fill(pi0Energy, cleanliness * completeness, numClusterHits);
        hComplCleanlAngle->Fill(pi0DecayAngle, cleanliness * completeness, numClusterHits);
        hComplCleanlConversionDistance->Fill(
          conversionDistance, cleanliness * completeness, numClusterHits);
        hNumHitsEnergy->Fill(pi0Energy, numClusterHits);

        // Continue if there are not two photons in the view
        if (clusterMap[tpc][plane]->GetPhotons().size() != 2) continue;

        double conversionSeparation =
          GetEndTrackDistance(clusterMap[tpc][plane]->GetPhotons().at(0).first,
                              clusterMap[tpc][plane]->GetPhotons().at(1).first);

        hCompletenessConversionSeparation->Fill(conversionSeparation, completeness, numClusterHits);
        hCleanlinessConversionSeparation->Fill(conversionSeparation, cleanliness, numClusterHits);
      }
    }
  }
}

void ClusteringValidation::ClusterAnalyser::WriteFile()
{

  // Average completeness/cleanliness
  double avCompleteness = hCompleteness->GetMean();
  double avCleanliness = hCleanliness->GetMean();

  // Write file
  std::ofstream outFile("effpur");
  outFile << avCompleteness << " " << avCleanliness;
  outFile.close();
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
class ClusteringValidation::ClusteringValidation : public art::EDAnalyzer {
public:
  explicit ClusteringValidation(fhicl::ParameterSet const& p);

private:
  void analyze(art::Event const& evt);
  void beginJob();
  void endJob();

  // Clusterings to compare and the hits which the clustering was run over
  std::vector<std::string> fClusterModuleLabels;
  std::string fHitsModuleLabel;

  // Minimum hits needed to analyse a plane
  int fMinHitsInPlane;

  // Canvas on which to save histograms
  TCanvas* fCanvas;

  // The cluster analysers
  // -- one for each clustering type being compared
  std::map<std::string, std::unique_ptr<ClusterAnalyser>> clusterAnalysis;
};

ClusteringValidation::ClusteringValidation::ClusteringValidation(fhicl::ParameterSet const& pset)
  : EDAnalyzer(pset)
{
  fMinHitsInPlane = pset.get<int>("MinHitsInPlane");
  fClusterModuleLabels = pset.get<std::vector<std::string>>("ClusterModuleLabels");
  fHitsModuleLabel = pset.get<std::string>("HitsModuleLabel");

  fCanvas = new TCanvas("fCanvas", "", 800, 600);
  gStyle->SetOptStat(0);
}

void ClusteringValidation::ClusteringValidation::analyze(art::Event const& evt)
{
  // Get the hits from the event
  art::Handle<std::vector<recob::Hit>> hitHandle;
  std::vector<art::Ptr<recob::Hit>> hits;
  if (evt.getByLabel(fHitsModuleLabel, hitHandle)) art::fill_ptr_vector(hits, hitHandle);

  // Get clustering information from event
  // and give to the ClusterAnalyser to analyse
  for (auto clustering : fClusterModuleLabels) {

    // Get the clusters from the event
    art::Handle<std::vector<recob::Cluster>> clusterHandle;
    std::vector<art::Ptr<recob::Cluster>> clusters;
    if (evt.getByLabel(clustering, clusterHandle)) art::fill_ptr_vector(clusters, clusterHandle);

    // Find the associations (the hits) for the clusters
    art::FindManyP<recob::Hit> fmh(clusterHandle, evt, clustering);

    // Analyse this particular clustering
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    clusterAnalysis.at(clustering)->Analyse(clockData, hits, clusters, fmh, fMinHitsInPlane);
  }
}

void ClusteringValidation::ClusteringValidation::beginJob()
{

  // Construct the cluster analysers here
  // -- one for each of the clustering types to compare
  for (auto clustering : fClusterModuleLabels)
    clusterAnalysis[clustering] = (std::unique_ptr<ClusterAnalyser>)new ClusterAnalyser(clustering);
}

void ClusteringValidation::ClusteringValidation::endJob()
{

  // Make a map of all the histograms for each clustering
  std::map<std::string, TObjArray> allHistograms;
  for (auto clustering : fClusterModuleLabels)
    allHistograms[clustering] = clusterAnalysis.at(clustering)->GetHistograms();

  // Write histograms to file
  TFile* file = TFile::Open("validationHistograms.root", "RECREATE");
  for (auto clustering : fClusterModuleLabels) {
    file->cd();
    file->mkdir(clustering.c_str());
    file->cd(clustering.c_str());
    for (int histIt = 0; histIt < allHistograms.begin()->second.GetEntriesFast(); ++histIt)
      allHistograms.at(clustering).At(histIt)->Write();
  }

  // Write images of overlaid histograms
  for (int histIt = 0; histIt < allHistograms.begin()->second.GetEntriesFast(); ++histIt) {
    fCanvas->cd();
    fCanvas->Clear();
    const char* name = allHistograms.begin()->second.At(histIt)->GetName();
    TLegend* l = new TLegend(0.6, 0.8, 0.8, 0.9, name, "brNDC");
    int clusterings = 1;
    for (std::map<std::string, TObjArray>::iterator clusteringIt = allHistograms.begin();
         clusteringIt != allHistograms.end();
         ++clusteringIt, ++clusterings) {
      TH1* h = (TH1*)allHistograms.at(clusteringIt->first).At(histIt);
      h->SetLineColor(clusterings);
      h->SetMarkerColor(clusterings);
      if (clusterings == 1)
        h->Draw();
      else
        h->Draw("same");
      l->AddEntry(h, clusteringIt->first.c_str(), "lp");
    }
    l->Draw("same");
    //fCanvas->SaveAs(name+TString(".png"));
    file->cd();
    fCanvas->Write(name);
  }

  file->Close();
  delete file;

  if (clusterAnalysis.find("blurredclusteringdc") != clusterAnalysis.end())
    clusterAnalysis.at("blurredclusteringdc")->WriteFile();
}

DEFINE_ART_MODULE(ClusteringValidation::ClusteringValidation)
