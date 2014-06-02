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
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "Geometry/Geometry.h"
#include "Utilities/AssociationUtil.h"
#include "RecoAlg/ClusterMergeHelper.h"
#include "RecoAlg/CMTool/CMAlgoMergeAll.h"
#include "RecoAlg/CMTool/CMAlgoArray.h"
#include "RecoAlg/CMTool/CMAlgoShortestDist.h"
#include "RecoAlg/CMTool/CMAlgoAngleCompat.h"
#include "RecoAlg/CMTool/CMAlgoTrackSeparate.h"

#include <memory>

namespace cluster {

  class SimpleClusterMerger : public art::EDProducer {
    
  public:

    explicit SimpleClusterMerger(fhicl::ParameterSet const & p);

    virtual ~SimpleClusterMerger();
    
    void produce(art::Event & evt) override;
    
  private:

    /// ClusterMergeHelper
    ClusterMergeHelper fCMerge;

    /// Input cluster data product producer name label
    std::string fClusterModuleLabel;

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

    /// Example merging algorithm: algorithm array container
    CMAlgoArray fMergeAlg;

    /// Merging algorithm 1
    CMAlgoShortestDist fDistAlg;

    /// Merging algorithm 2
    CMAlgoAngleCompat fAngleAlg;

    /// Example prohibit algorithm
    CMAlgoTrackSeparate fProhibitAlg;
  
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

    // Attach them to CMAlgoArray to make it into one merging algorithm
    fMergeAlg.AddAlgo(&fAngleAlg,true); // attach to CMAlgoArray in AND condition
    fMergeAlg.AddAlgo(&fDistAlg,true);  // attach to CMAlgoArray in AND condition

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

    fCMerge.GetManager().AddMergeAlgo(&fMergeAlg);        // Attach merging  algorithm
    fCMerge.GetManager().AddSeparateAlgo(&fProhibitAlg);  // Attach prohibit algorithm
    fCMerge.GetManager().SetMergePriority(::cluster::CMergeManager::kIndex); // Choose input cluster vector ordering as merging priority order
    fCMerge.GetManager().DebugMode(::cluster::CMergeManager::kPerIteration); // Set verbosity level to be per-merging-iteration report
    fCMerge.GetManager().MergeTillConverge(true);         // Set to iterate over till it finds no more newly merged clusters

    //
    // FYI there's an algorithm to just-merge-everything if you want to do a simple test (line commented out below)
    //
    //fCMerge.GetManager().AddMergeAlgo( new CMAlgoMergeAll );

    
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
    
    // Implementation of required member function here.
    
    //--- Set input cluster data ---//
    fCMerge.SetClusters(evt, fClusterModuleLabel);
    
    //--- Process merging ---//
    fCMerge.Process();
    
    //--- Retrieving Output ---//
    /*
      There are 2 kinds of output (that I can think of) you may want to retrieve:
      
      (0) A vector of merged cluster set (i.e. output) in terms of a collection of input cluster vector index
          Most basic information (which may be too basic to be useful for most users).
          You can retrieve through CBookKeeper instance:
      
	  fCMerge.GetManager().GetBookKeeper().GetResult();

      (1) A vector of CPAN
          This is practically the most useful thing. CPAN = a set of reconstructed cluster parameters.
	  You get a vector of CPAN corresponding a vector of merged clusters.
	  Note ClusterMergeHelper does not create recob::Cluster by itself after running Process() function.
	  To retrieve, try:

	  fCMerge.GetMergedCPAN();

      (2) A vector of clusters in terms of hit list
          which means std::vector<std::vector<art::Ptr<recob::Hit> > >, each content representing a
	  cluster by its hit list. 
	  To retrieve, try:

	  fCMerge.GetMergedClusterHits();

      (3) Fill std::vector<recob::Cluster> and art::Assns<recob::Cluster,recob::Hit> through pass-by-reference
          This may be useful if you want to just store the output in LArSoft data product (and if the output
	  of this merging process is the final result). This basically uses a return vector from (1) and (2)
	  above to create output recob::Cluster and associations to a hit list. ClusterMergeHelper can help
	  you this by a single function call.
	  For this, try:

	  fCMerge.AppendResult( this,
	                        evt,
			        *(out_clusters.get()), 
			        *(out_assn.get())
			       );

      In the following, as an example, I show how (3) can be done using (1) and (2). Hopefully this gives an
      idea of what kind of (detail) information you can retrieve. 
    */

    // Store output
    for(size_t out_index=0; out_index < fCMerge.GetMergedCPAN().size(); ++out_index) {
      
      // To save typing let's just retrieve const cluster_params instance
      const cluster_params &res = fCMerge.GetMergedCPAN().at(out_index).GetParams();
      
      // View_t needed but not a part of cluster_params, so retrieve it here
      geo::View_t view_id = geo->Plane(fCMerge.GetMergedCPAN().at(out_index).Plane()).View();
      
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
      
      util::CreateAssn(*this, 
		       evt, 
		       *(out_clusters.get()), 
		       fCMerge.GetMergedClusterHits().at(out_index), 
		       *(out_assn.get())
		       );

    }
    
    // As advertized, the following (commented-out) line does the same thing
    /*
      
    fCMerge.AppendResult( *this, evt, *(out_clusters.get()), *(out_assn.get()) );
    
  */

  evt.put(std::move(out_clusters));
  evt.put(std::move(out_assn));
    
  }
}

DEFINE_ART_MODULE(cluster::SimpleClusterMerger)
