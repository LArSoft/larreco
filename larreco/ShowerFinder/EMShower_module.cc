////////////////////////////////////////////////////////////////////////////
// Class:       EMShower
// Module Type: producer
// File:        EMShower_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), September 2015
//
// Module to make EM showers.
// Takes the output from cluster finding and track finding and combines
// information to make complete 3D shower.
//
// See DUNE-DocDB 1369 (public) for a detailed description.
////////////////////////////////////////////////////////////////////////////

// Framework includes:
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larreco/RecoAlg/EMShowerAlg.h"

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
  bool fSaveNonCompleteShowers;
  bool fFindBadPlanes;

  art::ServiceHandle<geo::Geometry> fGeom;
  detinfo::DetectorProperties const* fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();

  int fShower;
  int fPlane;

  int fDebug;

};

shower::EMShower::EMShower(fhicl::ParameterSet const& pset) : fEMShowerAlg(pset.get<fhicl::ParameterSet>("EMShowerAlg")) {
  this->reconfigure(pset);
  produces<std::vector<recob::Shower> >();
  produces<std::vector<recob::SpacePoint> >();
  produces<art::Assns<recob::Shower, recob::Hit> >();
  produces<art::Assns<recob::Shower, recob::Cluster> >();
  produces<art::Assns<recob::Shower, recob::Track> >();
  produces<art::Assns<recob::Shower, recob::SpacePoint> >();
  produces<art::Assns<recob::SpacePoint, recob::Hit> >();
}

void shower::EMShower::reconfigure(fhicl::ParameterSet const& p) {
  fHitsModuleLabel        = p.get<std::string>("HitsModuleLabel");
  fClusterModuleLabel     = p.get<std::string>("ClusterModuleLabel");
  fTrackModuleLabel       = p.get<std::string>("TrackModuleLabel");
  fPFParticleModuleLabel  = p.get<std::string>("PFParticleModuleLabel","");
  fFindBadPlanes          = p.get<bool>       ("FindBadPlanes","true");
  fSaveNonCompleteShowers = p.get<bool>       ("SaveNonCompleteShowers","true");
  fShower = p.get<int>("Shower",-1);
  fPlane = p.get<int>("Plane",-1);
  fDebug = p.get<int>("Debug",0);
  fEMShowerAlg.fDebug = fDebug;
}

void shower::EMShower::produce(art::Event& evt) {

  // Output -- showers and associations with hits and clusters
  std::unique_ptr<std::vector<recob::Shower> > showers(new std::vector<recob::Shower>);
  std::unique_ptr<std::vector<recob::SpacePoint> > spacePoints(new std::vector<recob::SpacePoint>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Cluster> > clusterAssociations(new art::Assns<recob::Shower, recob::Cluster>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Hit> > hitShowerAssociations(new art::Assns<recob::Shower, recob::Hit>);
  std::unique_ptr<art::Assns<recob::Shower, recob::Track> > trackAssociations(new art::Assns<recob::Shower, recob::Track>);
  std::unique_ptr<art::Assns<recob::Shower, recob::SpacePoint> > spShowerAssociations(new art::Assns<recob::Shower, recob::SpacePoint>);
  std::unique_ptr<art::Assns<recob::SpacePoint, recob::Hit> > hitSpAssociations(new art::Assns<recob::SpacePoint, recob::Hit>);

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
  art::Handle<std::vector<recob::PFParticle> > pfpHandle;
  std::vector<art::Ptr<recob::PFParticle> > pfps;
  if (evt.getByLabel(fPFParticleModuleLabel, pfpHandle))
    art::fill_ptr_vector(pfps, pfpHandle);

  // Associations
  art::FindManyP<recob::Hit> fmh(clusterHandle, evt, fClusterModuleLabel);
  art::FindManyP<recob::Track> fmt(hitHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::SpacePoint> fmsp(trackHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Cluster> fmc(hitHandle, evt, fHitsModuleLabel);

  // Make showers
  std::vector<std::vector<int> > newShowers;
  std::vector<unsigned int> pfParticles;

  std::map<int,std::vector<int> > clusterToTracks;
  std::map<int,std::vector<int> > trackToClusters;

  if (!pfpHandle.isValid()) {

    // Map between tracks and clusters
    fEMShowerAlg.AssociateClustersAndTracks(clusters, fmh, fmt, clusterToTracks, trackToClusters);

    // Make initial showers
    std::vector<std::vector<int> > initialShowers = fEMShowerAlg.FindShowers(trackToClusters);

    // Deal with views in which 2D reconstruction failed
    std::vector<int> clustersToIgnore;
    if (fFindBadPlanes)
      clustersToIgnore = fEMShowerAlg.CheckShowerPlanes(initialShowers, clusters, fmh);
    if (clustersToIgnore.size() > 0) {
      clusterToTracks.clear();
      trackToClusters.clear();
      fEMShowerAlg.AssociateClustersAndTracks(clusters, fmh, fmt, clustersToIgnore, clusterToTracks, trackToClusters);
      newShowers = fEMShowerAlg.FindShowers(trackToClusters);
    }
    else
      newShowers = initialShowers;

  }

  else {

    // Use pfparticle information
    art::FindManyP<recob::Cluster> fmcp(pfpHandle, evt, fPFParticleModuleLabel);
    for (size_t ipfp = 0; ipfp<pfps.size(); ++ipfp){
      art::Ptr<recob::PFParticle> pfp = pfps[ipfp];
      if (pfp->PdgCode()==11){ //shower particle
	if (fmcp.isValid()){
	  std::vector<int> clusters;
	  std::vector<art::Ptr<recob::Cluster> > clus = fmcp.at(ipfp);
	  for (size_t iclu = 0; iclu<clus.size(); ++iclu){
	    clusters.push_back(clus[iclu].key());
	  }
	  if (clusters.size()){
	    newShowers.push_back(clusters);
	    pfParticles.push_back(ipfp);
	  }
	}
      }
    }
  }

  // Make output larsoft products
  int showerNum = 0;
  for (std::vector<std::vector<int> >::iterator newShower = newShowers.begin(); newShower != newShowers.end(); ++newShower, ++showerNum) {

    if (showerNum != fShower and fShower != -1) continue;

    // New shower
    if (fDebug > 0)
      std::cout << std::endl << std::endl << "Start shower " << showerNum << std::endl;

    // New associations
    art::PtrVector<recob::Hit> showerHits;
    art::PtrVector<recob::Cluster> showerClusters;
    art::PtrVector<recob::Track> showerTracks;
    art::PtrVector<recob::SpacePoint> showerSpacePoints_p;

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
      if (!pfpHandle.isValid()) { // Only do this for non-pfparticle mode
	std::vector<int> clusterTracks = clusterToTracks.at(*showerCluster);
	for (std::vector<int>::iterator clusterTracksIt = clusterTracks.begin(); clusterTracksIt != clusterTracks.end(); ++clusterTracksIt)
	  if (std::find(associatedTracks.begin(), associatedTracks.end(), *clusterTracksIt) == associatedTracks.end())
	    associatedTracks.push_back(*clusterTracksIt);
      }
    }

    if (!pfpHandle.isValid()) { // For non-pfparticles, get space points from tracks
      // Tracks and space points
      for (std::vector<int>::iterator associatedTracksIt = associatedTracks.begin(); associatedTracksIt != associatedTracks.end(); ++associatedTracksIt) {
	art::Ptr<recob::Track> showerTrack = tracks.at(*associatedTracksIt);
	showerTracks.push_back(showerTrack);
      }
    }

    else { // For pfparticles, get space points from hits
      art::FindManyP<recob::SpacePoint> fmspp(showerHits, evt, fPFParticleModuleLabel);
      for (size_t ihit = 0; ihit<showerHits.size(); ++ihit){
	if (fmspp.isValid()){
	  std::vector<art::Ptr<recob::SpacePoint> > spacePoints_pfp = fmspp.at(ihit);
	  for (std::vector<art::Ptr<recob::SpacePoint> >::iterator spacePointsIt = spacePoints_pfp.begin(); spacePointsIt != spacePoints_pfp.end(); ++spacePointsIt)
	    showerSpacePoints_p.push_back(*spacePointsIt);
	}
      }
    }

    if (!pfpHandle.isValid()) {

      // Find the track at the start of the shower
      std::unique_ptr<recob::Track> initialTrack;
      std::map<int,std::vector<art::Ptr<recob::Hit> > > initialTrackHits;
      fEMShowerAlg.FindInitialTrack(showerHits, initialTrack, initialTrackHits, fPlane);

      // Make space points
      std::vector<std::vector<art::Ptr<recob::Hit> > > hitAssns;
      std::vector<recob::SpacePoint> showerSpacePoints = fEMShowerAlg.MakeSpacePoints(showerHits, hitAssns);
      int firstSpacePoint = spacePoints->size(), nSpacePoint = 0;
      for (std::vector<recob::SpacePoint>::const_iterator spacePointIt = showerSpacePoints.begin(); spacePointIt != showerSpacePoints.end(); ++spacePointIt, ++nSpacePoint) {
	spacePoints->emplace_back(spacePointIt->XYZ(), spacePointIt->ErrXYZ(), spacePointIt->Chisq(), spacePoints->size());
	util::CreateAssn(*this, evt, *(spacePoints.get()), hitAssns.at(nSpacePoint), *(hitSpAssociations.get()));
      }
      int lastSpacePoint = spacePoints->size();

      // Make shower object and associations
      recob::Shower shower = fEMShowerAlg.MakeShower(showerHits, initialTrack, initialTrackHits);
      shower.set_id(showerNum);
      if ( fSaveNonCompleteShowers or (!fSaveNonCompleteShowers and shower.ShowerStart() != TVector3(0,0,0)) ) {
	showers->push_back(shower);
	util::CreateAssn(*this, evt, *(showers.get()), showerHits,           *(hitShowerAssociations.get()));
	util::CreateAssn(*this, evt, *(showers.get()), showerClusters,       *(clusterAssociations.get()));
	util::CreateAssn(*this, evt, *(showers.get()), showerTracks,         *(trackAssociations.get()));
	util::CreateAssn(*this, evt, *(showers.get()), *(spacePoints.get()), *(spShowerAssociations.get()), firstSpacePoint, lastSpacePoint);
      }
      else
	mf::LogInfo("EMShower") << "Discarding shower " << showerNum << " due to incompleteness (SaveNonCompleteShowers == false)";
    }

    else { // pfParticle
      art::FindManyP<recob::Vertex> fmv(pfpHandle, evt, fPFParticleModuleLabel);
      std::vector<art::Ptr<recob::Vertex> > vertices = fmv.at(pfParticles[newShower-newShowers.begin()]);
      if (vertices.size()) {
	int iok = 0;
	recob::Shower shower = fEMShowerAlg.MakeShower(showerHits, vertices[0], iok);
	//shower.set_id(showerNum);
	if (iok==0) {
	  showers->push_back(shower);
	  showers->back().set_id(showers->size()-1);
	  util::CreateAssn(*this, evt, *(showers.get()), showerHits,          *(hitShowerAssociations.get()));
	  util::CreateAssn(*this, evt, *(showers.get()), showerClusters,      *(clusterAssociations.get()));
	  util::CreateAssn(*this, evt, *(showers.get()), showerTracks,        *(trackAssociations.get()));
	  util::CreateAssn(*this, evt, *(showers.get()), showerSpacePoints_p, *(spShowerAssociations.get()));
	}
      }

    }

  }

  // Put in event
  evt.put(std::move(showers));
  evt.put(std::move(spacePoints));
  evt.put(std::move(hitShowerAssociations));
  evt.put(std::move(clusterAssociations));
  evt.put(std::move(trackAssociations));
  evt.put(std::move(spShowerAssociations));
  evt.put(std::move(hitSpAssociations));

}

DEFINE_ART_MODULE(shower::EMShower)
