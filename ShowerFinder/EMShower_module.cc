////////////////////////////////////////////////////////////////////////////
// Class:       EMShower
// Module Type: producer
// File:        EMShower_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), September 2015
//
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
#include "MCCheater/BackTracker.h"

// LArSoft includes
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
#include "Utilities/AssociationUtil.h"

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

  explicit EMShower(fhicl::ParameterSet const& pset);

  void produce(art::Event &evt);
  void reconfigure(fhicl::ParameterSet const& p);

private:

};

shower::EMShower::EMShower(fhicl::ParameterSet const& pset) {
  this->reconfigure(pset);
  produces<std::vector<recob::Shower> >();
  produces<art::Assns<recob::Shower, recob::Hit> >();
  produces<art::Assns<recob::Shower, recob::Cluster> >();
}

void shower::EMShower::reconfigure(fhicl::ParameterSet const& p) { }

void shower::EMShower::produce(art::Event &evt) {

  // Output -- showers and associations with hits and clusters
  std::unique_ptr<std::vector<recob::Shower> > showers(new std::vector<recob::Shower>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Cluster> > clusterAssociations(new art::Assns<recob::Shower, recob::Cluster>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Hit> > hitAssociations(new art::Assns<recob::Shower, recob::Hit>);

  // Event has hits, tracks and clusters found already

  // Hits
  art::Handle<std::vector<recob::Hit> > hitHandle;
  std::vector<art::Ptr<recob::Hit> > hits;
  if (evt.getByLabel("lineclusterdc",hitHandle))
    art::fill_ptr_vector(hits, hitHandle);

  // Tracks
  art::Handle<std::vector<recob::Track> > trackHandle;
  std::vector<art::Ptr<recob::Track> > tracks;
  if (evt.getByLabel("pmtrackdc",trackHandle))
    art::fill_ptr_vector(tracks, trackHandle);

  // Clusters
  art::Handle<std::vector<recob::Cluster> > clusterHandle;
  std::vector<art::Ptr<recob::Cluster> > clusters;
  if (evt.getByLabel("blurredclusterdc",clusterHandle))
    art::fill_ptr_vector(clusters, clusterHandle);

  // Associations
  art::FindManyP<recob::Hit> fmh(clusterHandle, evt, "blurredclusterdc");
  art::FindManyP<recob::Track> fmt(hitHandle, evt, "pmtrackdc");

  // Map between tracks and clusters
  std::map<int,std::vector<int> > trackToClusters;

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

      // Add this cluster to the track map
      int track = clusterHitTracks.at(0)->ID();
      if (std::find(trackToClusters[track].begin(), trackToClusters[track].end(), (*clusterIt)->ID()) == trackToClusters[track].end())
	trackToClusters[track].push_back((*clusterIt)->ID());

    }     

  }

  std::cout << "End of event " << evt.event() << "; here are the clusters associated with each track..." << std::endl;
  for (std::map<int,std::vector<int> >::iterator trackToClusterIt = trackToClusters.begin(); trackToClusterIt != trackToClusters.end(); ++trackToClusterIt) {
    std::cout << "Track " << trackToClusterIt->first << " has the following clusters: " << std::endl;
    for (std::vector<int>::iterator clusterIt = trackToClusterIt->second.begin(); clusterIt != trackToClusterIt->second.end(); ++clusterIt)
      std::cout << *clusterIt << ", ";
    std::cout << std::endl << std::endl;
  }

  // Make showers

  // Shower vector
  std::vector<std::vector<int> > newShowers;

  // Loop over all tracks 
  for (std::map<int,std::vector<int> >::iterator trackToClusterIt = trackToClusters.begin(); trackToClusterIt != trackToClusters.end(); ++ trackToClusterIt) {

    // Find which showers already made are associated with this track
    std::vector<int> matchingShowers;
    for (unsigned int shower = 0; shower < newShowers.size(); ++shower)
      for (std::vector<int>::iterator cluster = trackToClusterIt->second.begin(); cluster != trackToClusterIt->second.end(); ++cluster)
	if ( (std::find(newShowers.at(shower).begin(), newShowers.at(shower).end(), *cluster) != newShowers.at(shower).end()) and
	     (std::find(matchingShowers.begin(), matchingShowers.end(), shower)) == matchingShowers.end() )
	  matchingShowers.push_back(shower);

    // Shouldn't be more than one
    if (matchingShowers.size() > 1)
      std::cout << "WARNING! Number of showers this track matches is " << matchingShowers.size() << std::endl;

    // New shower
    if (matchingShowers.size() < 1)
      newShowers.push_back(trackToClusterIt->second);

    // Add to existing shower
    else {
      for (std::vector<int>::iterator cluster = trackToClusterIt->second.begin(); cluster != trackToClusterIt->second.end(); ++cluster)
	if (std::find(newShowers.at(matchingShowers.at(0)).begin(), newShowers.at(matchingShowers.at(0)).end(), *cluster) == newShowers.at(matchingShowers.at(0)).end())
	  newShowers.at(matchingShowers.at(0)).push_back(*cluster);
    }

  }

  std::cout << "There are " << newShowers.size() << " showers: " << std::endl;
  for (unsigned int shower = 0; shower < newShowers.size(); ++shower) {
    std::cout << "Shower " << shower << " is composed of the following clusters..." << std::endl;
    for (std::vector<int>::iterator clusterIt = newShowers.at(shower).begin(); clusterIt != newShowers.at(shower).end(); ++clusterIt)
      std::cout << *clusterIt << ", ";
    std::cout << std::endl << std::endl;
  }

}


DEFINE_ART_MODULE(shower::EMShower)
