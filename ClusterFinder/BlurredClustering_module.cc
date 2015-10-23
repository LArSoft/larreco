////////////////////////////////////////////////////////////////////////
// Class:       BlurredClustering
// Module Type: producer
// File:        BlurredClustering_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), May 2015
//
// Reconstructs showers by blurring the hit map image to introduce fake
// hits before clustering to make fuller and more complete clusters
////////////////////////////////////////////////////////////////////////

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

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Utilities/AssociationUtil.h"
#include "ClusterFinder/ClusterCreator.h"
#include "RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "RecoAlg/ClusterParamsImportWrapper.h"
#include "RecoAlg/BlurredClusteringAlg.h"
#include "RecoAlg/MergeClusterAlg.h"

// ROOT & C++ includes
#include <string>
#include <vector>
#include <map>

namespace cluster {
  class BlurredClustering;
}

class cluster::BlurredClustering: public art::EDProducer {
public:

  explicit BlurredClustering(fhicl::ParameterSet const& pset);
  virtual ~BlurredClustering();

  void cluster(std::vector<art::Ptr<recob::Hit> > const &hits, std::vector<art::PtrVector<recob::Hit> > &clusters, int tpc, int plane);
  void produce(art::Event &evt);
  void reconfigure(fhicl::ParameterSet const &p);
  void showerHits(std::vector<art::Ptr<recob::Hit> > const& hits, art::FindManyP<recob::Track> const& fmt, std::vector<art::Ptr<recob::Hit> >& hitsToCluster);

private:

  int fEvent, fRun, fSubrun;
  std::string fHitsModuleLabel, fTrackModuleLabel;
  bool fCreateDebugPDF, fMergeClusters, fGlobalTPCRecon, fShowerReconOnly;

  // Create instances of algorithm classes to perform the clustering
  cluster::BlurredClusteringAlg fBlurredClusteringAlg;
  cluster::MergeClusterAlg fMergeClusterAlg;

  // Output containers to place in event
  std::unique_ptr<std::vector<recob::Cluster> > clusters;
  std::unique_ptr<art::Assns<recob::Cluster,recob::Hit> > associations;

};

cluster::BlurredClustering::BlurredClustering(fhicl::ParameterSet const &pset) : fBlurredClusteringAlg(pset.get<fhicl::ParameterSet>("BlurredClusterAlg")),
                                                                                 fMergeClusterAlg(pset.get<fhicl::ParameterSet>("MergeClusterAlg")) {
  this->reconfigure(pset);
  produces<std::vector<recob::Cluster> >();
  produces<art::Assns<recob::Cluster,recob::Hit> >();
}

cluster::BlurredClustering::~BlurredClustering() { }

void cluster::BlurredClustering::reconfigure(fhicl::ParameterSet const& p) {
  fHitsModuleLabel  = p.get<std::string>("HitsModuleLabel");
  fTrackModuleLabel = p.get<std::string>("TrackModuleLabel");
  fCreateDebugPDF   = p.get<bool>       ("CreateDebugPDF");
  fMergeClusters    = p.get<bool>       ("MergeClusters");
  fGlobalTPCRecon   = p.get<bool>       ("GlobalTPCRecon");
  fShowerReconOnly  = p.get<bool>       ("ShowerReconOnly");
  fBlurredClusteringAlg.reconfigure(p.get<fhicl::ParameterSet>("BlurredClusterAlg"));
  fMergeClusterAlg.reconfigure(p.get<fhicl::ParameterSet>("MergeClusterAlg"));
}

void cluster::BlurredClustering::produce(art::Event &evt) {

  fEvent  = evt.event();
  fRun    = evt.run();
  fSubrun = evt.subRun();

  // Create debug pdf to illustrate the blurring process
  if (fCreateDebugPDF)
    fBlurredClusteringAlg.CreateDebugPDF(fRun, fSubrun, fEvent);

  // Output containers -- collection of clusters and associations
  clusters.reset(new std::vector<recob::Cluster>);
  associations.reset(new art::Assns<recob::Cluster,recob::Hit>);

  // Compute the cluster characteristics
  // Just use default for now, but configuration will go here
  ClusterParamsImportWrapper<StandardClusterParamsAlg> ClusterParamAlgo;

  // Create geometry handle
  art::ServiceHandle<geo::Geometry> geom;

  // Get the hit collection from the event
  art::Handle<std::vector<recob::Hit> > hitCollection;
  evt.getByLabel(fHitsModuleLabel,hitCollection);

  // Get the tracks from the event
  art::Handle<std::vector<recob::Track> > trackCollection;
  std::vector<art::Ptr<recob::Track> > tracks;
  if (evt.getByLabel(fTrackModuleLabel,trackCollection))
    art::fill_ptr_vector(tracks, trackCollection);

  // Global recon -- merged TPCs
  if (fGlobalTPCRecon) {

    // Make a map between the planes and the hits on each
    std::map<std::pair<int,int>,std::vector<art::Ptr<recob::Hit> > > planeToHits;
    for (size_t hitIt = 0; hitIt < hitCollection->size(); ++hitIt)
      planeToHits[std::make_pair(hitCollection->at(hitIt).WireID().Plane,hitCollection->at(hitIt).WireID().TPC%2)].push_back(art::Ptr<recob::Hit>(hitCollection,hitIt));

    // Loop over views
    for (std::map<std::pair<int,int>,std::vector<art::Ptr<recob::Hit> > >::iterator planeIt = planeToHits.begin(); planeIt != planeToHits.end(); ++planeIt) {

      //std::cout << "Clustering in plane " << planeIt->first.first << " in global TPC " << planeIt->first.second << std::endl;

      // Make the clusters
      std::vector<art::PtrVector<recob::Hit> > finalClusters;
      std::vector<art::Ptr<recob::Hit> > hitsToCluster;
      if (fShowerReconOnly and trackCollection.isValid()) {
	art::FindManyP<recob::Track> fmt(hitCollection, evt, fTrackModuleLabel);
	showerHits(planeIt->second, fmt, hitsToCluster);
      }
      else
	hitsToCluster = planeIt->second;

      cluster(hitsToCluster, finalClusters, planeIt->first.second, planeIt->first.first);

      for (std::vector<art::PtrVector<recob::Hit> >::iterator clusIt = finalClusters.begin(); clusIt != finalClusters.end(); ++clusIt) {

	art::PtrVector<recob::Hit> clusterHits = *clusIt;
	if (clusterHits.size() > 0) {

	  // Get the start and end wires of the cluster
	  unsigned int startWire = fBlurredClusteringAlg.FindGlobalWire(clusterHits.front()->WireID());
	  unsigned int endWire = fBlurredClusteringAlg.FindGlobalWire(clusterHits.back()->WireID());

	  // Put cluster hits in the algorithm
	  ClusterParamAlgo.ImportHits(clusterHits);

	  // Create the recob::Cluster and place in the vector of clusters
	  ClusterCreator cluster(
				 ClusterParamAlgo,                        // algo
				 float(startWire),                        // start_wire
				 0.,                                      // sigma_start_wire
				 clusterHits.front()->PeakTime(),         // start_tick
				 clusterHits.front()->SigmaPeakTime(),    // sigma_start_tick
				 float(endWire),                          // end_wire
				 0.,                                      // sigma_end_wire,
				 clusterHits.back()->PeakTime(),          // end_tick
				 clusterHits.back()->SigmaPeakTime(),     // sigma_end_tick
				 clusters->size(),                        // ID
				 clusterHits.front()->View(),             // view
				 clusterHits.front()->WireID().planeID(), // plane
				 recob::Cluster::Sentry                   // sentry
				 );

	  clusters->emplace_back(cluster.move());

	  // Associate the hits to this cluster
	  util::CreateAssn(*this, evt, *(clusters.get()), clusterHits, *(associations.get()));

	} // End this cluster

      } // End loop over all clusters

    }
  }

  // Classic recon -- separate TPCs
  else {

    // Make a map between the planes and TPCs and hits on each
    std::map<geo::PlaneID,std::vector<art::Ptr<recob::Hit> > > planeIDToHits;
    for (size_t hitIt = 0; hitIt < hitCollection->size(); ++hitIt)
      planeIDToHits[hitCollection->at(hitIt).WireID().planeID()].push_back(art::Ptr<recob::Hit>(hitCollection,hitIt));

    // Loop over views
    for (std::map<geo::PlaneID,std::vector<art::Ptr<recob::Hit> > >::iterator planeIt = planeIDToHits.begin(); planeIt != planeIDToHits.end(); ++planeIt) {

      // Make the clusters
      std::vector<art::PtrVector<recob::Hit> > finalClusters;
      std::vector<art::Ptr<recob::Hit> > hitsToCluster;
      if (fShowerReconOnly and trackCollection.isValid()) {
	art::FindManyP<recob::Track> fmt(hitCollection, evt, fTrackModuleLabel);
	showerHits(planeIt->second, fmt, hitsToCluster);
      }
      else
	hitsToCluster = planeIt->second;
      cluster(hitsToCluster, finalClusters, planeIt->first.TPC, planeIt->first.Plane);

      for (std::vector<art::PtrVector<recob::Hit> >::iterator clusIt = finalClusters.begin(); clusIt != finalClusters.end(); ++clusIt) {

	art::PtrVector<recob::Hit> clusterHits = *clusIt;
	if (clusterHits.size() > 0) {

	  // Get the start and end wires of the cluster
	  unsigned int startWire = fBlurredClusteringAlg.FindGlobalWire(clusterHits.front()->WireID());
	  unsigned int endWire = fBlurredClusteringAlg.FindGlobalWire(clusterHits.back()->WireID());

	  // Put cluster hits in the algorithm
	  ClusterParamAlgo.ImportHits(clusterHits);
	
	  // Create the recob::Cluster and place in the vector of clusters
	  ClusterCreator cluster(
				 ClusterParamAlgo,                        // algo
				 float(startWire),                        // start_wire
				 0.,                                      // sigma_start_wire
				 clusterHits.front()->PeakTime(),         // start_tick
				 clusterHits.front()->SigmaPeakTime(),    // sigma_start_tick
				 float(endWire),                          // end_wire
				 0.,                                      // sigma_end_wire,
				 clusterHits.back()->PeakTime(),          // end_tick
				 clusterHits.back()->SigmaPeakTime(),     // sigma_end_tick
				 clusters->size(),                        // ID
				 clusterHits.front()->View(),             // view
				 clusterHits.front()->WireID().planeID(), // plane
				 recob::Cluster::Sentry                   // sentry
				 );

	  clusters->emplace_back(cluster.move());

	  // Associate the hits to this cluster
	  util::CreateAssn(*this, evt, *(clusters.get()), clusterHits, *(associations.get()));

	} // End this cluster

      } // End loop over all clusters

    }
  }

  evt.put(std::move(clusters));
  evt.put(std::move(associations));

  return;
    
}

void cluster::BlurredClustering::cluster(std::vector<art::Ptr<recob::Hit> > const& allHits, std::vector<art::PtrVector<recob::Hit> >& finalClusters, int tpc, int plane) {

  /// Takes vector of recob::Hits and returns vector of clusters

  // Implement the algorithm
  if (allHits.size() >= fBlurredClusteringAlg.GetMinSize()) {

    // Convert hit map to TH2 histogram and blur it
    TH2F image = fBlurredClusteringAlg.ConvertRecobHitsToTH2(allHits);
    TH2F* blurred = fBlurredClusteringAlg.GaussianBlur(&image);

    // Find clusters in histogram
    std::vector<std::vector<int> > allClusterBins; // Vector of clusters (clusters are vectors of hits)
    int numClusters = fBlurredClusteringAlg.FindClusters(blurred, allClusterBins);
    mf::LogVerbatim("Blurred Clustering") << "Found " << numClusters << " clusters" << std::endl;

    // Create output clusters from the vector of clusters made in FindClusters
    std::vector<art::PtrVector<recob::Hit> > planeClusters;
    fBlurredClusteringAlg.ConvertBinsToClusters(&image, allClusterBins, planeClusters);

    // Use the cluster merging algorithm
    if (fMergeClusters) {
      int numMergedClusters = fMergeClusterAlg.MergeClusters(planeClusters, finalClusters);
      mf::LogVerbatim("Blurred Clustering") << "After merging, there are " << numMergedClusters << " clusters" << std::endl;
    }
    else finalClusters = planeClusters;

    // Make the debug PDF
    if (fCreateDebugPDF) {
      fBlurredClusteringAlg.SaveImage(&image, 1, tpc, plane);
      fBlurredClusteringAlg.SaveImage(blurred, 2, tpc, plane);
      fBlurredClusteringAlg.SaveImage(blurred, allClusterBins, 3, tpc, plane);
      fBlurredClusteringAlg.SaveImage(&image, finalClusters, 4, tpc, plane);
    }

    blurred->Delete();

  } // End min hits check

  fBlurredClusteringAlg.fHitMap.clear();

  return;

}

void cluster::BlurredClustering::showerHits(std::vector<art::Ptr<recob::Hit> > const& initialHits, art::FindManyP<recob::Track> const& fmt, std::vector<art::Ptr<recob::Hit> >& hitsToCluster) {

  /// Takes all hits and track associations and returns just hits which are not determined to be track-like

  for (std::vector<art::Ptr<recob::Hit> >::const_iterator initialHit = initialHits.begin(); initialHit != initialHits.end(); ++initialHit) {
    std::vector<art::Ptr<recob::Track> > tracks = fmt.at(initialHit->key());
    // Shower-like tracks have this bit set
    if (tracks.size() and (tracks.at(0)->ID() & 65536) == 0)
      continue;
    hitsToCluster.push_back(*initialHit);
  }

}

DEFINE_ART_MODULE(cluster::BlurredClustering)
