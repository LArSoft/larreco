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
#include "RecoBase/PFParticle.h"
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

  std::string fHitsModuleLabel, fClusterModuleLabel, fTrackModuleLabel, fPFParticleModuleLabel;
  EMShowerAlg fEMShowerAlg;

  art::ServiceHandle<geo::Geometry> fGeom;
  art::ServiceHandle<util::DetectorProperties> fDetProp;

  int fShower;
  int fPlane;

};

shower::EMShower::EMShower(fhicl::ParameterSet const& pset) : fEMShowerAlg(pset.get<fhicl::ParameterSet>("EMShowerAlg")) {
  this->reconfigure(pset);
  produces<std::vector<recob::Shower> >();
  produces<art::Assns<recob::Shower, recob::Hit> >();
  produces<art::Assns<recob::Shower, recob::Cluster> >();
  produces<art::Assns<recob::Shower, recob::Track> >();
  produces<art::Assns<recob::Shower, recob::SpacePoint> >();
}

void shower::EMShower::reconfigure(fhicl::ParameterSet const& p) {
  fHitsModuleLabel       = p.get<std::string>("HitsModuleLabel");
  fClusterModuleLabel    = p.get<std::string>("ClusterModuleLabel");
  fTrackModuleLabel      = p.get<std::string>("TrackModuleLabel");
  fPFParticleModuleLabel = p.get<std::string>("PFParticleModuleLabel","");
  fShower = p.get<int>("Shower",-1);
  fPlane = p.get<int>("Plane",-1);
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

  // PFParticles
  art::Handle<std::vector<recob::PFParticle> >pfpHandle;
  std::vector<art::Ptr<recob::PFParticle> >pfps;
  if (evt.getByLabel(fPFParticleModuleLabel, pfpHandle))
    art::fill_ptr_vector(pfps, pfpHandle);

  // Associations
  art::FindManyP<recob::Hit> fmh(clusterHandle, evt, fClusterModuleLabel);
  art::FindManyP<recob::Track> fmt(hitHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::SpacePoint> fmsp(trackHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Cluster> fmc(hitHandle, evt, fHitsModuleLabel);

  // Make showers
  std::vector<std::vector<int> > newShowers;

  std::map<int,std::vector<int> > clusterToTracks;
  std::map<int,std::vector<int> > trackToClusters;

  if (!pfpHandle.isValid()){

    // Map between tracks and clusters
    fEMShowerAlg.AssociateClustersAndTracks(clusters, fmh, fmt, clusterToTracks, trackToClusters);

    // Fix issues where one view screws things up
    std::vector<std::vector<int> > initialShowers;
    fEMShowerAlg.FindShowers(trackToClusters, initialShowers);
    if (fGeom->MaxPlanes() > 2) {
      std::vector<int> clustersToIgnore;
      fEMShowerAlg.CheckShowerPlanes(initialShowers, clustersToIgnore, clusters, fmh);
      if (clustersToIgnore.size() > 0) {
	clusterToTracks.clear();
	trackToClusters.clear();
	fEMShowerAlg.AssociateClustersAndTracks(clusters, fmh, fmt, clustersToIgnore, clusterToTracks, trackToClusters);
	fEMShowerAlg.FindShowers(trackToClusters, newShowers);
      }
      else
	newShowers = initialShowers;
    }
    else
      newShowers = initialShowers;
  }
  else{
    //use pfparticle information
    art::FindManyP<recob::Cluster> fmcp(pfpHandle, evt, fPFParticleModuleLabel);
    for (size_t ipfp = 0; ipfp<pfps.size(); ++ipfp){
      art::Ptr<recob::PFParticle> pfp = pfps[ipfp];
      if (pfp->PdgCode()==11){//shower particle
	if (fmcp.isValid()){
	  std::vector<int> clusters;
	  std::vector<art::Ptr<recob::Cluster> > clus = fmcp.at(ipfp);
	  for (size_t iclu = 0; iclu<clus.size(); ++iclu){
	    clusters.push_back(clus[iclu].key());
	  }
	  if (clusters.size()) newShowers.push_back(clusters);
	}
      }
    }
  }
  
  // Make output larsoft products
  int showerNum = 0;
  for (std::vector<std::vector<int> >::iterator newShower = newShowers.begin(); newShower != newShowers.end(); ++newShower, ++showerNum) {

    if (showerNum != fShower and fShower != -1) continue;

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

      if (!pfpHandle.isValid()){//only do this for non-pfparticle mode
	// Tracks
	std::vector<int> clusterTracks = clusterToTracks.at(*showerCluster);
	for (std::vector<int>::iterator clusterTracksIt = clusterTracks.begin(); clusterTracksIt != clusterTracks.end(); ++clusterTracksIt)
	  if (std::find(associatedTracks.begin(), associatedTracks.end(), *clusterTracksIt) == associatedTracks.end())
	    associatedTracks.push_back(*clusterTracksIt);
      }
    }

    if (!pfpHandle.isValid()){
      //for non-pfparticles, get space points from tracks
      // Tracks and space points
      for (std::vector<int>::iterator associatedTracksIt = associatedTracks.begin(); associatedTracksIt != associatedTracks.end(); ++associatedTracksIt) {
	art::Ptr<recob::Track> showerTrack = tracks.at(*associatedTracksIt);
	showerTracks.push_back(showerTrack);
	std::vector<art::Ptr<recob::SpacePoint> > spacePoints = fmsp.at(showerTrack.key());
	for (std::vector<art::Ptr<recob::SpacePoint> >::iterator spacePointsIt = spacePoints.begin(); spacePointsIt != spacePoints.end(); ++spacePointsIt)
	  showerSpacePoints.push_back(*spacePointsIt);
      }
    }
    else{
      //for pfparticles, get space points from hits
      art::FindManyP<recob::SpacePoint> fmspp(showerHits, evt, fPFParticleModuleLabel);
      for (size_t ihit = 0; ihit<showerHits.size(); ++ihit){
	if (fmspp.isValid()){
	  std::vector<art::Ptr<recob::SpacePoint> >spacePoints = fmspp.at(ihit);
	  for (std::vector<art::Ptr<recob::SpacePoint> >::iterator spacePointsIt = spacePoints.begin(); spacePointsIt != spacePoints.end(); ++spacePointsIt)
	    showerSpacePoints.push_back(*spacePointsIt);
	}
      }
    }
    
    if (!pfpHandle.isValid()){
      // Find the track at the start of the shower
      std::unique_ptr<recob::Track> initialTrack;
      std::map<int,std::vector<art::Ptr<recob::Hit> > > initialTrackHits;
      fEMShowerAlg.FindInitialTrack(showerHits, initialTrack, initialTrackHits, fmc, fPlane);
      
      // Make shower object and associations
      recob::Shower shower = fEMShowerAlg.MakeShower(showerHits, initialTrack, initialTrackHits);
      shower.set_id(showerNum);
      showers->push_back(shower);
    }
    else{
      art::FindManyP<recob::Vertex> fmv(pfpHandle, evt, fPFParticleModuleLabel);
      std::vector<art::Ptr<recob::Vertex> > vertices = fmv.at(showerNum);
      if (vertices.size()){
	recob::Shower shower = fEMShowerAlg.MakeShower(showerHits, vertices[0]);
	shower.set_id(showerNum);
	showers->push_back(shower);
      }
    }

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
