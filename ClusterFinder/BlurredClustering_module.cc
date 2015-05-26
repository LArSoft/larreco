////////////////////////////////////////////////////////////////////////
// Class:       BlurredClustering
// Module Type: producer
// File:        BlurredClustering_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), May 2015
//
// Reconstructs showers by blurred the hit map image to introduce fake
// hits before clustering to make fully and more complete clusters
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
#include "Filters/ChannelFilter.h"
#include "RecoAlg/BlurredClusteringAlg.h"
#include "ClusterFinder/ClusterCreator.h"
#include "RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "RecoAlg/ClusterParamsImportWrapper.h"

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

  void produce(art::Event &evt);
  void reconfigure(fhicl::ParameterSet const &p);

private:

  int fEvent, fRun, fSubrun;
  std::string fHitsModuleLabel;
  bool fCreateDebugPDF;

  // Create instance of algorithm class to perform the clustering
  cluster::BlurredClusteringAlg fBlurredClusteringAlg;

};

cluster::BlurredClustering::BlurredClustering(fhicl::ParameterSet const &pset) : fBlurredClusteringAlg(pset.get<fhicl::ParameterSet>("BlurredClusteringAlg")) {
  this->reconfigure(pset);
  produces<std::vector<recob::Cluster> >();
  produces<art::Assns<recob::Cluster,recob::Hit> >();
}

cluster::BlurredClustering::~BlurredClustering() { }

void cluster::BlurredClustering::reconfigure(fhicl::ParameterSet const& p) {
  fHitsModuleLabel = p.get<std::string>("HitsModuleLabel");
  fCreateDebugPDF  = p.get<bool>       ("CreateDebugPDF");
  fBlurredClusteringAlg.reconfigure(p.get<fhicl::ParameterSet>("BlurredClusteringAlg"));
}

void cluster::BlurredClustering::produce(art::Event &evt) {

  fEvent  = evt.event();
  fRun    = evt.run();
  fSubrun = evt.subRun();

  // Create debug pdf to illustrate the blurring process
  if (fCreateDebugPDF)
    fBlurredClusteringAlg.CreateDebugPDF(fEvent, fRun, fSubrun, fCreateDebugPDF);

  // Output containers -- collection of clusters and associations
  std::unique_ptr<std::vector<recob::Cluster> > clusters(new std::vector<recob::Cluster>);
  std::unique_ptr<art::Assns<recob::Cluster,recob::Hit> > associations(new art::Assns<recob::Cluster,recob::Hit>);

  // Compute the cluster characteristics
  // Just use default for now, but configuration will go here
  ClusterParamsImportWrapper<StandardClusterParamsAlg> ClusterParamAlgo;

  // Create geometry handle
  art::ServiceHandle<geo::Geometry> geom;

  // Get the hit collection from the event
  art::Handle<std::vector<recob::Hit> > hitCollection;
  evt.getByLabel(fHitsModuleLabel,hitCollection);

  // Get the channel filter
  filter::ChannelFilter channelFilter;

  // Make a map between the planes and the hits on each
  std::map<geo::PlaneID,std::vector<art::Ptr<recob::Hit> > > planeIDToHits;
  for(size_t hitIt = 0; hitIt < hitCollection->size(); ++hitIt)
    planeIDToHits[hitCollection->at(hitIt).WireID().planeID()].push_back(art::Ptr<recob::Hit>(hitCollection,hitIt));

  // Loop over the planes
  for (std::map<geo::PlaneID,std::vector<art::Ptr<recob::Hit> > >::iterator planeIt = planeIDToHits.begin(); planeIt != planeIDToHits.end(); ++planeIt) {

    // Find the plane type
    //geo::SigType_t sigType = geom->SignalType(planeIt->first);
    fBlurredClusteringAlg.fPlane = planeIt->first.Plane;
    fBlurredClusteringAlg.fTPC   = planeIt->first.TPC;

    // Vector to hold all hits in event for particular plane
    std::vector<art::Ptr<recob::Hit> > *allHits(&planeIt->second);

    // Implement the algorithm
    if (allHits->size() >= fBlurredClusteringAlg.fMinHits) {

      // Convert hit map to TH2 histogram
      TH2F image = fBlurredClusteringAlg.ConvertRecobHitsToTH2(allHits);

      /// Find clusters in histogram
      std::vector<std::vector<int> > allClusterBins; /// Vector of clusters (clusters are vectors of hits)
      int numClusters = fBlurredClusteringAlg.FindClusters(&image, allClusterBins);
      mf::LogVerbatim("Blurred Clustering") << "Found "<< numClusters << " clusters" << std::endl;

      /// Create output clusters from the vector of clusters made in FindClusters
      std::vector<art::PtrVector<recob::Hit> > planeClusters = fBlurredClusteringAlg.ConvertBinsToClusters(&image, allHits, allClusterBins); // returns vector of clusters

      for (std::vector<art::PtrVector<recob::Hit> >::iterator clusIt = planeClusters.begin(); clusIt != planeClusters.end(); ++clusIt) {

	art::PtrVector<recob::Hit> clusterHits = *clusIt;
	if (clusterHits.size() > 0) {

	  // Get the start and end wires of the cluster
	  unsigned int startWire = clusterHits.front()->WireID().Wire;
	  unsigned int endWire = clusterHits.back()->WireID().Wire;
	
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

    } // End min hits check

    fBlurredClusteringAlg.fHitMap.clear();
    allHits->clear();

  } // End loop over planes

  evt.put(std::move(clusters));
  evt.put(std::move(associations));

  return;
    
}

DEFINE_ART_MODULE(cluster::BlurredClustering)
