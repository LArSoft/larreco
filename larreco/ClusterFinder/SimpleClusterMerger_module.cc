////////////////////////////////////////////////////////////////////////
// Class:       SimpleClusterMerger
// Module Type: producer
// File:        SimpleClusterMerger_module.cc
//
// Generated at Tue May 27 14:15:41 2014 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/PxHitConverter.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larreco/RecoAlg/CMTool/CMToolApp/CMergeHelper.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoMergeAll.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoArray.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoShortestDist.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoAngleCompat.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoTrackSeparate.h"
#include "larreco/RecoAlg/ClusterRecoUtil/LazyClusterParamsAlg.h"

#include <memory>

namespace cluster {

  class SimpleClusterMerger : public art::EDProducer {
    
  public:

    explicit SimpleClusterMerger(fhicl::ParameterSet const & p);

    virtual ~SimpleClusterMerger();
    
    void produce(art::Event & evt) override;
    
  private:

    /// ClusterMergeHelper
    ::cmtool::CMergeHelper fCMerge;

    /// Input cluster data product producer name label
    std::string fClusterModuleLabel;

    /// GeometryUtilities instance
    ::util::GeometryUtilities fGeoU;

    //--- CBAlgo instances ---//
    /*
      CMergeManager takes pointer to CBoolAlgoBase inherit class instance.
      So the simplest way is to create them on heap and keep the pointers.
      But that has a concern to computation speed (though I think that would
      never be a slow enough to be a concern)... so here we take an example
      of using CBAlgo instances created on stack.
    */

    /// Example merging algorithm: algorithm array container
    ::cmtool::CBAlgoArray fMergeAlg;

    /// Merging algorithm 1
    ::cmtool::CBAlgoShortestDist fDistAlg;

    /// Merging algorithm 2
    ::cmtool::CBAlgoAngleCompat fAngleAlg;

    /// Example prohibit algorithm
    ::cmtool::CBAlgoTrackSeparate fProhibitAlg;
  
  };
}

namespace cluster {

  SimpleClusterMerger::SimpleClusterMerger(fhicl::ParameterSet const & p)
  {
    // Declare output data products
    produces< std::vector<recob::Cluster> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    
    // Fill fcl parameter
    fClusterModuleLabel = p.get<std::string>("InputClusterLabel");

    //--- Configure Merging Algorithm ---//
    /*
      This is where we should pass fcl parameters to configure various algorithms
      we defined as class member. Here, for simplicity, I don't pass any configuration
      parameters but you should in your custom merging module.
    */
    
    // Configure angle algorithm
    fAngleAlg.SetVerbose(false);       // no verbous mode... annoying
    fAngleAlg.SetMinHits(10);          // Set minimum # hits to be 10
    fAngleAlg.SetAngleCut(180);        // Set angle-diff cut < 180 degree (i.e. anything :))
    fAngleAlg.SetAllow180Ambig(true);  // Allow 180-degree ambiguity (direction mis-reco)

    // Configure distance algorithm
    fDistAlg.SetVerbose(false);        // No verbous mode ... annoying
    fDistAlg.SetMinHits(10);           // Set minimum # hits to be 10 
    fDistAlg.SetSquaredDistanceCut(9); // Set distance-squared cut to be 9 cm^2

    // Attach them to CBAlgoArray to make it into one merging algorithm
    fMergeAlg.AddAlgo(&fAngleAlg,true); // attach to CBAlgoArray in AND condition
    fMergeAlg.AddAlgo(&fDistAlg,true);  // attach to CBAlgoArray in AND condition

    //--- Configure Prohibit Algorithm ---//

    // I configure this using totally arbitrary numbers + I do not configure all parameters... 
    // This is just for example.
    fProhibitAlg.SetVerbose(false);
    fProhibitAlg.SetDebug(false);
    fProhibitAlg.SetMinNumHits(10);       
    fProhibitAlg.SetMinAngleDiff(5);
    fProhibitAlg.SetMaxOpeningAngle(10);
    fProhibitAlg.SetMinLength(5);
    

    //--- Configure Merger ---//

    fCMerge.GetManager(0).AddMergeAlgo(&fMergeAlg);        // Attach merging  algorithm
    fCMerge.GetManager(0).AddSeparateAlgo(&fProhibitAlg);  // Attach prohibit algorithm
    fCMerge.GetManager(0).DebugMode(::cmtool::CMergeManager::kPerIteration); // Set verbosity level to be per-merging-iteration report
    fCMerge.GetManager(0).MergeTillConverge(true);         // Set to iterate over till it finds no more newly merged clusters

    //
    // FYI there's an algorithm to just-merge-everything if you want to do a simple test (line commented out below)
    //
    //fCMerge.GetManager(0).AddMergeAlgo( new CBAlgoMergeAll );

    
  }
  
  SimpleClusterMerger::~SimpleClusterMerger()
  {
    // Clean up dynamic memory and other resources here.
  }
  
  void SimpleClusterMerger::produce(art::Event & evt)
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
      
      // this "algo" is actually parroting its cluster_params
      LazyClusterParamsAlg algo(res);
      
      std::vector<art::Ptr<recob::Hit> > merged_hits;
      for(auto const& c_index : merged_clusters[out_index]) {
        const std::vector<art::Ptr<recob::Hit> >& hits = hit_m.at(c_index);
        merged_hits.reserve(merged_hits.size() + hits.size());
        for(auto const& ptr : hits) merged_hits.push_back(ptr);
      }
      
      // the full plane needed but not a part of cluster_params...
      // get the one from the first hit
      geo::PlaneID plane; // invalid by default
      if (!merged_hits.empty()) plane = merged_hits.front()->WireID().planeID();
      
      // View_t needed but not a part of cluster_params, so retrieve it here
      geo::View_t view_id = geo->Plane(plane).View();
      
      // Push back a new cluster data product with parameters copied from cluster_params
      out_clusters->emplace_back(
        res.start_point.w / fGeoU.WireToCm(), // start_wire
        0.,                                   // sigma_start_wire
        res.start_point.t / fGeoU.TimeToCm(), // start_tick
        0.,                                   // sigma_start_tick
        algo.StartCharge().value(),           // start_charge
        algo.StartAngle().value(),            // start_angle
        algo.StartOpeningAngle().value(),     // start_opening
        res.end_point.w   / fGeoU.WireToCm(), // end_wire
        0.,                                   // sigma_end_wire
        res.end_point.t   / fGeoU.TimeToCm(), // end_tick
        0.,                                   // sigma_end_tick
        algo.EndCharge().value(),             // end_charge
        algo.EndAngle().value(),              // end_angle
        algo.EndOpeningAngle().value(),       // end_opening
        algo.Integral().value(),              // integral
        algo.IntegralStdDev().value(),        // integral_stddev
        algo.SummedADC().value(),             // summedADC
        algo.SummedADCStdDev().value(),       // summedADC_stddev
        algo.NHits(),                         // n_hits
        algo.MultipleHitDensity(),              // multiple_hit_density
        algo.Width(),                         // width
        out_clusters->size(),                 // ID
        view_id,                              // view
        plane,                                // plane
        recob::Cluster::Sentry                // sentry
        );
      
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

DEFINE_ART_MODULE(cluster::SimpleClusterMerger)
