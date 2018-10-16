////////////////////////////////////////////////////////////////////////
// Class:       BlurredClustering
// Module Type: producer
// File:        BlurredClustering_module.cc
// Author:      Mike Wallbank (m.wallbank@sheffield.ac.uk), May 2015
//
// Reconstructs showers by blurring the hit map image to introduce fake
// hits before clustering to make fuller and more complete clusters.
//
// See DUNE-DocDB 54 (public) for details.
////////////////////////////////////////////////////////////////////////

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

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/ClusterFinder/ClusterCreator.h"
#include "larreco/RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "larreco/RecoAlg/ClusterParamsImportWrapper.h"
#include "larreco/RecoAlg/BlurredClusteringAlg.h"
#include "larreco/RecoAlg/MergeClusterAlg.h"
#include "larreco/RecoAlg/TrackShowerSeparationAlg.h"

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
  std::string fHitsModuleLabel, fTrackModuleLabel, fVertexModuleLabel, fPFParticleModuleLabel;
  bool fCreateDebugPDF, fMergeClusters, fGlobalTPCRecon, fShowerReconOnly;

  // Create instances of algorithm classes to perform the clustering
  cluster::BlurredClusteringAlg fBlurredClusteringAlg;
  cluster::MergeClusterAlg fMergeClusterAlg;
  shower::TrackShowerSeparationAlg fTrackShowerSeparationAlg;

  // Output containers to place in event
  std::unique_ptr<std::vector<recob::Cluster> > clusters;
  std::unique_ptr<art::Assns<recob::Cluster,recob::Hit> > associations;

};

cluster::BlurredClustering::BlurredClustering(fhicl::ParameterSet const &pset) : fBlurredClusteringAlg(pset.get<fhicl::ParameterSet>("BlurredClusterAlg")),
                                                                                 fMergeClusterAlg(pset.get<fhicl::ParameterSet>("MergeClusterAlg")),
										 fTrackShowerSeparationAlg(pset.get<fhicl::ParameterSet>("TrackShowerSeparationAlg")) {
  this->reconfigure(pset);
  produces<std::vector<recob::Cluster> >();
  produces<art::Assns<recob::Cluster,recob::Hit> >();
}

cluster::BlurredClustering::~BlurredClustering() { }

void cluster::BlurredClustering::reconfigure(fhicl::ParameterSet const& p) {
  fHitsModuleLabel       = p.get<std::string>("HitsModuleLabel");
  fTrackModuleLabel      = p.get<std::string>("TrackModuleLabel");
  fVertexModuleLabel     = p.get<std::string>("VertexModuleLabel");
  fPFParticleModuleLabel = p.get<std::string>("PFParticleModuleLabel");
  fCreateDebugPDF        = p.get<bool>       ("CreateDebugPDF");
  fMergeClusters         = p.get<bool>       ("MergeClusters");
  fGlobalTPCRecon        = p.get<bool>       ("GlobalTPCRecon");
  fShowerReconOnly       = p.get<bool>       ("ShowerReconOnly");
  fBlurredClusteringAlg.reconfigure(p.get<fhicl::ParameterSet>("BlurredClusterAlg"));
  fMergeClusterAlg.reconfigure(p.get<fhicl::ParameterSet>("MergeClusterAlg"));
  fTrackShowerSeparationAlg.reconfigure(p.get<fhicl::ParameterSet>("TrackShowerSeparationAlg"));
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

  // Get the hits from the event
  art::Handle<std::vector<recob::Hit> > hitCollection;
  std::vector<art::Ptr<recob::Hit> > hits;
  std::vector<art::Ptr<recob::Hit> > hitsToCluster;
  if (evt.getByLabel(fHitsModuleLabel,hitCollection))
    art::fill_ptr_vector(hits, hitCollection);

  if (fShowerReconOnly) {

    // Get the tracks from the event
    art::Handle<std::vector<recob::Track> > trackCollection;
    std::vector<art::Ptr<recob::Track> > tracks;
    if (evt.getByLabel(fTrackModuleLabel,trackCollection))
      art::fill_ptr_vector(tracks, trackCollection);

    // Get the space points from the event
    art::Handle<std::vector<recob::SpacePoint> > spacePointCollection;
    std::vector<art::Ptr<recob::SpacePoint> > spacePoints;
    if (evt.getByLabel(fTrackModuleLabel,spacePointCollection))
      art::fill_ptr_vector(spacePoints, spacePointCollection);

    // Get vertices from the event
    art::Handle<std::vector<recob::Vertex> > vertexCollection;
    std::vector<art::Ptr<recob::Vertex> > vertices;
    if (evt.getByLabel(fVertexModuleLabel, vertexCollection))
      art::fill_ptr_vector(vertices, vertexCollection);

    // Get pandora pfparticles and clusters from the event
    art::Handle<std::vector<recob::PFParticle> > pfParticleCollection;
    std::vector<art::Ptr<recob::PFParticle> > pfParticles;
    if (evt.getByLabel(fPFParticleModuleLabel, pfParticleCollection))
      art::fill_ptr_vector(pfParticles, pfParticleCollection);
    art::Handle<std::vector<recob::Cluster> > clusterCollection;
    evt.getByLabel(fPFParticleModuleLabel, clusterCollection);

    if (trackCollection.isValid()) {
      art::FindManyP<recob::Hit> fmht(trackCollection, evt, fTrackModuleLabel);
      art::FindManyP<recob::Track> fmth(hitCollection, evt, fTrackModuleLabel);
      art::FindManyP<recob::SpacePoint> fmspt(trackCollection, evt, fTrackModuleLabel);
      art::FindManyP<recob::Track> fmtsp(spacePointCollection, evt, fTrackModuleLabel);
      hitsToCluster = fTrackShowerSeparationAlg.SelectShowerHits(evt.event(), hits, tracks, spacePoints, fmht, fmth, fmspt, fmtsp);
    }

    // // Remove hits from tracks before performing any clustering
    // if (pfParticleCollection.isValid() and clusterCollection.isValid()) {
    //   mf::LogInfo("BlurredCluster") << "Removing track-like hits before clustering: will use information from PFParticles." << std::endl;
    //   art::FindManyP<recob::Cluster> fmcpfp(pfParticleCollection, evt, fPFParticleModuleLabel);
    //   art::FindManyP<recob::Hit> fmhpfp(clusterCollection, evt, fPFParticleModuleLabel);
    //   hitsToCluster = fTrackShowerSeparationAlg.RemoveTrackHits(hits, pfParticles, fmcpfp, fmhpfp);
    // }
    // else if (trackCollection.isValid() and trackCollection.isValid() and spacePointCollection.isValid()) {
    //   mf::LogInfo("BlurredCluster") << "Removing track-like hits before clustering: no PFParticle information available so will use tracks and vertices." << std::endl;
    //   art::FindManyP<recob::Track> fmth(hitCollection, evt, fTrackModuleLabel);
    //   art::FindManyP<recob::Track> fmtsp(spacePointCollection, evt, fTrackModuleLabel);
    //   art::FindManyP<recob::Hit> fmh(trackCollection, evt, fTrackModuleLabel);
    //   hitsToCluster = fTrackShowerSeparationAlg.RemoveTrackHits(hits, tracks, spacePoints, vertices, fmth, fmtsp, fmh, evt.event(), evt.run());
    // }
    // else
    //   throw art::Exception(art::errors::Configuration) << "Error: configuration is set to remove track-like hits before clustering but no prior reconstruction is provided... "
    // 						       << std::endl
    // 						       << "Require either a) Pandora or b) track, vertices and space points to have already been found." << std::endl;

  }

  else
    hitsToCluster = hits;

  // Make a map between the planes and the hits on each
  std::map<std::pair<int,int>,std::vector<art::Ptr<recob::Hit> > > planeToHits;
  for (std::vector<art::Ptr<recob::Hit> >::iterator hitToCluster = hitsToCluster.begin(); hitToCluster != hitsToCluster.end(); ++hitToCluster) {
    if (fGlobalTPCRecon)
      planeToHits[std::make_pair((*hitToCluster)->WireID().Plane,(*hitToCluster)->WireID().TPC%2)].push_back(*hitToCluster);
    else
      planeToHits[std::make_pair((*hitToCluster)->WireID().Plane,(*hitToCluster)->WireID().TPC)].push_back(*hitToCluster);
  }

  // Loop over views
  for (std::map<std::pair<int,int>,std::vector<art::Ptr<recob::Hit> > >::iterator planeIt = planeToHits.begin(); planeIt != planeToHits.end(); ++planeIt) {

    //std::cout << "Clustering in plane " << planeIt->first.first << " in global TPC " << planeIt->first.second << std::endl;
    // if (!(planeIt->first.first == 1 and planeIt->first.second == 1))
    //   continue;

    std::vector<art::PtrVector<recob::Hit> > finalClusters;

    // Implement the algorithm
    if (planeIt->second.size() >= fBlurredClusteringAlg.GetMinSize()) {

      // Convert hit map to TH2 histogram and blur it
      std::vector<std::vector<double> > image = fBlurredClusteringAlg.ConvertRecobHitsToVector(planeIt->second);
      std::vector<std::vector<double> > blurred = fBlurredClusteringAlg.GaussianBlur(image);

       // Find clusters in histogram
      std::vector<std::vector<int> > allClusterBins; // Vector of clusters (clusters are vectors of hits)
      int numClusters = fBlurredClusteringAlg.FindClusters(blurred, allClusterBins);
      mf::LogVerbatim("Blurred Clustering") << "Found " << numClusters << " clusters" << std::endl;

      // Create output clusters from the vector of clusters made in FindClusters
      std::vector<art::PtrVector<recob::Hit> > planeClusters;
      fBlurredClusteringAlg.ConvertBinsToClusters(image, allClusterBins, planeClusters);

      // Use the cluster merging algorithm
      if (fMergeClusters) {
	int numMergedClusters = fMergeClusterAlg.MergeClusters(planeClusters, finalClusters);
	mf::LogVerbatim("Blurred Clustering") << "After merging, there are " << numMergedClusters << " clusters" << std::endl;
      }
      else finalClusters = planeClusters;

      // Make the debug PDF
      if (fCreateDebugPDF) {
	std::stringstream name;
	name << "blurred_image";
	TH2F* imageHist = fBlurredClusteringAlg.MakeHistogram(image, TString(name.str()));
	name << "_convolved";
	TH2F* blurredHist = fBlurredClusteringAlg.MakeHistogram(blurred, TString(name.str()));
      	fBlurredClusteringAlg.SaveImage(imageHist, 1, planeIt->first.second, planeIt->first.first);
      	fBlurredClusteringAlg.SaveImage(blurredHist, 2, planeIt->first.second, planeIt->first.first);
      	fBlurredClusteringAlg.SaveImage(blurredHist, allClusterBins, 3, planeIt->first.second, planeIt->first.first);
      	fBlurredClusteringAlg.SaveImage(imageHist, finalClusters, 4, planeIt->first.second, planeIt->first.first);
	imageHist->Delete();
	blurredHist->Delete();
      }

    } // End min hits check

    //fBlurredClusteringAlg.fHitMap.clear();

    // Make the output cluster objects
    for (std::vector<art::PtrVector<recob::Hit> >::iterator clusIt = finalClusters.begin(); clusIt != finalClusters.end(); ++clusIt) {

      art::PtrVector<recob::Hit> clusterHits = *clusIt;
      if (clusterHits.size() > 0) {

	// Get the start and end wires of the cluster
	unsigned int startWire = fBlurredClusteringAlg.GlobalWire(clusterHits.front()->WireID());
	unsigned int endWire = fBlurredClusteringAlg.GlobalWire(clusterHits.back()->WireID());

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

  evt.put(std::move(clusters));
  evt.put(std::move(associations));

  return;
    
}

DEFINE_ART_MODULE(cluster::BlurredClustering)
