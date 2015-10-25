////////////////////////////////////////////////////////////////////////////
// Class:       EMShower
// Module Type: producer
// File:        EMShower_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), September 2015
//
// Module to make EM showers.
// Takes the output from cluster finding and track finding and combines
// information to make complete 3D shower.
////////////////////////////////////////////////////////////////////////////

// Framework includes:
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
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/FindManyP.h"

// LArSoft includes
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "AnalysisAlg/CalorimetryAlg.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/Shower.h"
#include "RecoAlg/EMShowerAlg.h"

// ROOT includes
#include "TPrincipal.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TF2.h"
#include "TTree.h"

namespace shower {
  class EMShower;
}

class shower::EMShower : public art::EDProducer {
public:

  EMShower(fhicl::ParameterSet const& pset);

  void produce(art::Event& evt);
  void reconfigure(fhicl::ParameterSet const& p);

private:

  std::string fHitsModuleLabel, fClusterModuleLabel, fTrackModuleLabel;
  int fMinTrackLength;

  EMShowerAlg fEMShowerAlg;
  calo::CalorimetryAlg fCalorimetryAlg;

  art::ServiceHandle<geo::Geometry> fGeom;
  art::ServiceHandle<util::DetectorProperties> fDetProp;

};

shower::EMShower::EMShower(fhicl::ParameterSet const& pset) : fEMShowerAlg(),
							      fCalorimetryAlg(pset.get<fhicl::ParameterSet>("CalorimetryAlg")) {
  this->reconfigure(pset);
  produces<std::vector<recob::Shower> >();
  produces<art::Assns<recob::Shower, recob::Hit> >();
  produces<art::Assns<recob::Shower, recob::Cluster> >();
  produces<art::Assns<recob::Shower, recob::Track> >();
  produces<art::Assns<recob::Shower, recob::SpacePoint> >();
}

void shower::EMShower::reconfigure(fhicl::ParameterSet const& p) {
  fHitsModuleLabel    = p.get<std::string>("HitsModuleLabel");
  fClusterModuleLabel = p.get<std::string>("ClusterModuleLabel");
  fTrackModuleLabel   = p.get<std::string>("TrackModuleLabel");
  fMinTrackLength     = p.get<int>        ("MinTrackLength");
}

void shower::EMShower::produce(art::Event& evt) {

  // Output -- showers and associations with hits and clusters
  std::unique_ptr<std::vector<recob::Shower> > showers(new std::vector<recob::Shower>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Cluster> > clusterAssociations(new art::Assns<recob::Shower, recob::Cluster>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Hit> > hitAssociations(new art::Assns<recob::Shower, recob::Hit>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Track> > trackAssociations(new art::Assns<recob::Shower, recob::Track>);
  std::unique_ptr<art::Assns<recob::Shower, recob::SpacePoint> > spacePointAssociations(new art::Assns<recob::Shower, recob::SpacePoint>);

  // Event has hits, tracks and clusters found already

  // Hits
  art::Handle<std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if (evt.getByLabel(fHitsModuleLabel,hitHandle))
    art::fill_ptr_vector(hits, hitHandle);

  // Tracks
  art::Handle<std::vector<recob::Track> > trackHandle;
  std::vector<art::Ptr<recob::Track> > tracks;
  if (evt.getByLabel(fTrackModuleLabel,trackHandle))
    art::fill_ptr_vector(tracks, trackHandle);

  // Clusters
  art::Handle<std::vector<recob::Cluster> > clusterHandle;
  std::vector<art::Ptr<recob::Cluster> > clusters;
  if (evt.getByLabel(fClusterModuleLabel,clusterHandle))
    art::fill_ptr_vector(clusters, clusterHandle);

  // Associations
  art::FindManyP<recob::Hit> fmh(clusterHandle, evt, fClusterModuleLabel);
  art::FindManyP<recob::Track> fmt(hitHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::SpacePoint> fmsp(trackHandle, evt, fTrackModuleLabel);

  // Map between tracks and clusters
  std::map<int,std::vector<int> > trackToClusters;
  std::map<int,std::vector<int> > clusterToTracks;

  // Look through all the clusters
  for (std::vector<art::Ptr<recob::Cluster> >::iterator clusterIt = clusters.begin(); clusterIt != clusters.end(); ++clusterIt) {

    // Get the hits in this cluster
    std::vector<art::Ptr<recob::Hit> > clusterHits = fmh.at(clusterIt->key());

    // Look at all these hits and find the associated tracks
    for (std::vector<art::Ptr<recob::Hit> >::iterator clusterHitIt = clusterHits.begin(); clusterHitIt != clusterHits.end(); ++clusterHitIt) {

      // Get the tracks associated with this hit
      std::vector<art::Ptr<recob::Track> > clusterHitTracks = fmt.at(clusterHitIt->key());
      if (clusterHitTracks.size() > 1) { std::cout << "More than one track associated with this hit!" << std::endl; continue; }
      if (clusterHitTracks.size() < 1) continue;
      if (clusterHitTracks.at(0)->Length() < fMinTrackLength) continue;

      // Add this cluster to the track map
      int track = clusterHitTracks.at(0).key();
      int cluster = (*clusterIt).key();
      if (std::find(trackToClusters[track].begin(), trackToClusters[track].end(), cluster) == trackToClusters[track].end())
	trackToClusters[track].push_back(cluster);
      if (std::find(clusterToTracks[cluster].begin(), clusterToTracks[cluster].end(), track) == clusterToTracks[cluster].end())
	clusterToTracks[cluster].push_back(track);

    }

  }

  // Make showers
  std::vector<std::vector<int> > newShowers;
  fEMShowerAlg.MakeShowers(trackToClusters, newShowers);

  // Make output larsoft products
  int showerNum = 0;
  for (std::vector<std::vector<int> >::iterator newShower = newShowers.begin(); newShower != newShowers.end(); ++newShower, ++showerNum) {

    // New shower
    std::cout << std::endl << "Start shower " << showerNum << std::endl;

    // New associations
    art::PtrVector<recob::Hit> showerHits;
    art::PtrVector<recob::Cluster> showerClusters;
    art::PtrVector<recob::Track> showerTracks;
    art::PtrVector<recob::SpacePoint> showerSpacePoints;

    std::vector<int> associatedTracks;

    // Make showers and associations
    for (std::vector<int>::iterator showerCluster = (*newShower).begin(); showerCluster != (*newShower).end(); ++showerCluster) {

      // Clusters
      art::Ptr<recob::Cluster> cluster = clusters.at(*showerCluster);
      showerClusters.push_back(cluster);

      // Hits
      std::vector<art::Ptr<recob::Hit> > showerClusterHits = fmh.at(cluster.key());
      for (std::vector<art::Ptr<recob::Hit> >::iterator showerClusterHit = showerClusterHits.begin(); showerClusterHit != showerClusterHits.end(); ++showerClusterHit)
  	showerHits.push_back(*showerClusterHit);

      // Tracks
      std::vector<int> clusterTracks = clusterToTracks.at(*showerCluster);
      for (std::vector<int>::iterator clusterTracksIt = clusterTracks.begin(); clusterTracksIt != clusterTracks.end(); ++clusterTracksIt)
  	if (std::find(associatedTracks.begin(), associatedTracks.end(), *clusterTracksIt) == associatedTracks.end())
  	  associatedTracks.push_back(*clusterTracksIt);

    }

    // Tracks and space points
    for (std::vector<int>::iterator associatedTracksIt = associatedTracks.begin(); associatedTracksIt != associatedTracks.end(); ++associatedTracksIt) {
      art::Ptr<recob::Track> showerTrack = tracks.at(*associatedTracksIt);
      showerTracks.push_back(showerTrack);
      std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmsp.at(showerTrack.key());
      for (std::vector<art::Ptr<recob::SpacePoint> >::iterator spacePointsIt = spacePoints.begin(); spacePointsIt != spacePoints.end(); ++spacePointsIt)
  	showerSpacePoints.push_back(*spacePointsIt);
    }

    // Find the properties of this shower
    TVector3 direction, directionError, vertex, vertexError;
    std::vector<double> totalEnergy, totalEnergyError, dEdx, dEdxError;
    int bestPlane;
    fEMShowerAlg.FindShowerProperties(showerHits, fmt, fCalorimetryAlg, direction, directionError, vertex, vertexError, totalEnergy, totalEnergyError, dEdx, dEdxError, bestPlane);

    // Make shower object and associations
    showers->emplace_back(direction, directionError, vertex, vertexError, totalEnergy, totalEnergyError, dEdx, dEdxError, bestPlane, showerNum);
    util::CreateAssn(*this, evt, *(showers.get()), showerHits,        *(hitAssociations.get()));
    util::CreateAssn(*this, evt, *(showers.get()), showerClusters,    *(clusterAssociations.get()));
    util::CreateAssn(*this, evt, *(showers.get()), showerTracks,      *(trackAssociations.get()));
    util::CreateAssn(*this, evt, *(showers.get()), showerSpacePoints, *(spacePointAssociations.get()));

   }

  // Put in event
  evt.put(std::move(showers));
  evt.put(std::move(hitAssociations));
  evt.put(std::move(clusterAssociations));
  evt.put(std::move(trackAssociations));
  evt.put(std::move(spacePointAssociations));

}

DEFINE_ART_MODULE(shower::EMShower)
