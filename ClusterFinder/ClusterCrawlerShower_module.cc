////////////////////////////////////////////////////////////////////////
// Class:       ClusterCrawlerShower
// Module Type: producer
// File:        ClusterCrawlerShower_module.cc
//
// Generated at Tue May 27 14:15:41 2014 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "Geometry/Geometry.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/PxHitConverter.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoAlg/CMTool/CMToolApp/CMergeHelper.h"
#include "RecoAlg/CMTool/CMTAlgMerge/CBAlgoMergeAll.h"
#include "RecoAlg/CMTool/CMTAlgMerge/CBAlgoArray.h"
#include "RecoAlg/CMTool/CMTAlgMerge/CBAlgoShortestDist.h"
#include "RecoAlg/CMTool/CMTAlgMerge/CBAlgoShortestDistSmallCluster.h"

#include <memory>

namespace cluster {

  class ClusterCrawlerShower : public art::EDProducer {
    
  public:

    explicit ClusterCrawlerShower(fhicl::ParameterSet const & p);

    virtual ~ClusterCrawlerShower();
    
    void produce(art::Event & evt) override;
    
  private:

    /// ClusterMergeHelper
    ::cmtool::CMergeHelper fCMerge;

    /// Input cluster data product producer name label
    std::string fClusterModuleLabel;

    bool fDistAlg0Verbose;        // 
    int fDistAlg0MinHits;           // 
    int fDistAlg0MaxHits;           // 
    double fDistAlg0SquaredDistanceCut; // 

    bool fDistAlg1Verbose;        // 
    int fDistAlg1MinHits;           // 
    double fDistAlg1SquaredDistanceCut; // 

    /// GeometryUtilities instance
    ::util::GeometryUtilities fGeoU;

    //--- CMAlgo instances ---//
    /*
      CMergeManager takes pointer to CBoolAlgoBase inherit class instance.
      So the simplest way is to create them on heap and keep the pointers.
      But that has a concern to computation speed (though I think that would
      never be a slow enough to be a concern)... so here we take an example
      of using CMAlgo instances created on stack.
    */

    /// Algo 0: Merge clusters that are small in order to do shower recognition
    ::cmtool::CBAlgoArray fMergeAlg0;
    /// Algo 0: Distance algorithm for small cluster merging
    ::cmtool::CBAlgoShortestDistSmallCluster fDistAlg0;

    /// Algo 1: Merge shower-like clusters that intercept
    ::cmtool::CBAlgoArray fMergeAlg1;
    /// Algo 0: Distance algorithm for shower-like clusters that intercept, it should be 0
    ::cmtool::CBAlgoShortestDist	fDistAlg1;
  
  };
}

namespace cluster {

  ClusterCrawlerShower::ClusterCrawlerShower(fhicl::ParameterSet const & p)
  {
    // Declare output data products
    produces< std::vector<recob::Cluster> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    
    // Fill fcl parameter
    fClusterModuleLabel = p.get<std::string>("InputClusterLabel");
    
    fDistAlg0Verbose = p.get<bool>("DistAlg0Verbose");
    fDistAlg0MinHits = p.get<int>("DistAlg0MinHits");
    fDistAlg0MaxHits = p.get<int>("DistAlg0MaxHits");
    fDistAlg0SquaredDistanceCut = p.get<double>("DistAlg0SquaredDistanceCut");
    
    fDistAlg1Verbose = p.get<bool>("DistAlg1Verbose");
    fDistAlg1MinHits = p.get<int>("DistAlg1MinHits");
    fDistAlg1SquaredDistanceCut = p.get<double>("DistAlg1SquaredDistanceCut");

    //--- Configure Merging Algorithm ---//
    /*
      This is where we should pass fcl parameters to configure various algorithms
      we defined as class member. Here, for simplicity, I don't pass any configuration
      parameters but you should in your custom merging module.
    */

    /// Configure algorithm 0, Merges small clusters in order to perform shower recognition
    fDistAlg0.SetVerbose(fDistAlg0Verbose);        // No verbous mode ... annoying
    fDistAlg0.SetMinHits(fDistAlg0MinHits);           // Set minimum # hits 
    fDistAlg0.SetMaxHits(fDistAlg0MaxHits);           // Set minimum # hits 
    fDistAlg0.SetSquaredDistanceCut(fDistAlg0SquaredDistanceCut); // Set distance-squared cut, is in cm^2 
    fMergeAlg0.AddAlgo(&fDistAlg0,true); // attach to CMAlgoArray in AND condition

    /// Configure algorithm 0, Merge shower like clusters that intercept
    fDistAlg1.SetVerbose(fDistAlg1Verbose);        // No verbous mode ... annoying
    fDistAlg1.SetMinHits(fDistAlg1MinHits);           // Set minimum # hits 
    fDistAlg1.SetSquaredDistanceCut(fDistAlg1SquaredDistanceCut); // Set distance-squared cut, it should be 0 and is in cm^2 
    fMergeAlg1.AddAlgo(&fDistAlg1,true); // attach to CMAlgoArray in AND condition


    //--- Configure Merger ---//
    fCMerge.GetManager(0).AddMergeAlgo(&fMergeAlg0);        // Attach merging  algorithm
    //fCMerge.GetManager(0).AddMergeAlgo(&fMergeAlg1);        // Attach merging  algorithm

    fCMerge.GetManager(0).DebugMode(::cmtool::CMergeManager::kPerIteration); // Set verbosity level to be per-merging-iteration report
    //fCMerge.GetManager(0).DebugMode(::cluster::CMergeManager::kPerMerging); // Set verbosity level to be per-merging-iteration report
    fCMerge.GetManager(0).MergeTillConverge(true);         // Set to iterate over till it finds no more newly merged clusters

    //
    // FYI there's an algorithm to just-merge-everything if you want to do a simple test (line commented out below)
    //
    //fCMerge.GetManager(0).AddMergeAlgo( new CMAlgoMergeAll );

    
  }
  
  ClusterCrawlerShower::~ClusterCrawlerShower()
  {
    // Clean up dynamic memory and other resources here.
  }
  
  void ClusterCrawlerShower::produce(art::Event & evt)
  {
    std::unique_ptr<std::vector<recob::Cluster> > out_clusters(new std::vector<recob::Cluster>);
    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit> > out_assn(new art::Assns<recob::Cluster, recob::Hit>);

    art::ServiceHandle<geo::Geometry> geo;

    //
    // Preparation
    //

    // Retrieve input clusters
    art::Handle<std::vector<recob::Cluster> > cHandle;
    evt.getByLabel(fClusterModuleLabel,cHandle);

    if(!cHandle.isValid())
      throw cet::exception(__FUNCTION__) << "Invalid input cluster label!" << std::endl;

    // Cluster type conversion: recob::Hit => util::PxHit
    std::vector<std::vector< ::util::PxHit> > local_clusters;
    art::FindManyP<recob::Hit> hit_m(cHandle, evt, fClusterModuleLabel);
    ::util::PxHitConverter conv;
    for(size_t i=0; i<cHandle->size(); ++i) {

      local_clusters.push_back(std::vector< ::util::PxHit>());

      const std::vector<art::Ptr<recob::Hit> >& hits = hit_m.at(i);

      conv.GeneratePxHit(hits, local_clusters.back());
    }

    //--- Process merging ---//
    fCMerge.Process(local_clusters);

    // Store output
    auto merged_clusters = fCMerge.GetResult().GetResult();

    auto const& cpan_v = fCMerge.GetClusters();
    if(merged_clusters.size()!=cpan_v.size())

      throw cet::exception(__FUNCTION__) << "LOGIC ERROR: merged cluster id length != output cluster counts..." << std::endl;

    for(size_t out_index=0; out_index < merged_clusters.size(); ++out_index) {

      // To save typing let's just retrieve const cluster_params instance
      const cluster_params &res = cpan_v[out_index].GetParams();

      // View_t needed but not a part of cluster_params, so retrieve it here
      geo::View_t view_id = geo->Plane(cpan_v[out_index].Plane()).View();

      // Push back a new cluster data product with parameters copied from cluster_params
      out_clusters->push_back( recob::Cluster( res.start_point.w / fGeoU.WireToCm(), 0,  // start wire & error
                                               res.start_point.t / fGeoU.TimeToCm(), 0,  // start time & error
                                               res.end_point.w   / fGeoU.WireToCm(), 0,  // end   wire & error
                                               res.end_point.t   / fGeoU.TimeToCm(), 0,  // end   time & error
                                               res.cluster_angle_2d,                 0,  // dT/dW (slope)
                                               0,                                    0,  // dQ/dW (what is that?)
                                               res.sum_charge,                           // charge sum
                                               view_id,                                  // geo::View_t
                                               out_clusters->size()                      // Cluster ID
                                               )
                               );


      std::vector<art::Ptr<recob::Hit> > merged_hits;

      for(auto const& c_index : merged_clusters[out_index]) {

	const std::vector<art::Ptr<recob::Hit> >& hits = hit_m.at(c_index);

	merged_hits.reserve(merged_hits.size() + hits.size());

	for(auto const& ptr : hits) merged_hits.push_back(ptr);

      } 

      util::CreateAssn(*this,
                       evt,
                       *(out_clusters.get()),
                       merged_hits,
		       *(out_assn.get())
                       );

    }

    evt.put(std::move(out_clusters));
    evt.put(std::move(out_assn));
  }
}

DEFINE_ART_MODULE(cluster::ClusterCrawlerShower)
