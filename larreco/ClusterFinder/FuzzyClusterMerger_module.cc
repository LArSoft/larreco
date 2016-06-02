////////////////////////////////////////////////////////////////////////
// Class:       FuzzyClusterMerger
// Module Type: producer
// File:        FuzzyClusterMerger_module.cc
//
// Generated at Tue May 27 14:15:41 2014 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/Utilities/PxHitConverter.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"
#include "larreco/RecoAlg/CMTool/CMToolApp/CMergeHelper.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoMergeAll.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoArray.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoShortestDist.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoShortestDistSmallCluster.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoTrackSeparate.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoOutOfConeSeparate.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoAngleIncompat.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoCenterOfMass.h"

#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoStartTrack.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoPolyContain.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoPolyOverlap.h"
#include "larreco/RecoAlg/CMTool/CMTAlgMerge/CBAlgoPolyShortestDist.h"
#include "larreco/RecoAlg/CMTool/CMTAlgPriority/CPAlgoIgnoreTracks.h"
#include "larreco/RecoAlg/ClusterRecoUtil/LazyClusterParamsAlg.h"



#include <memory>

namespace cluster {

  class FuzzyClusterMerger : public art::EDProducer {
    
  public:

    explicit FuzzyClusterMerger(fhicl::ParameterSet const & p);

    virtual ~FuzzyClusterMerger();
    
    void produce(art::Event & evt) override;
    
  private:

    /// ClusterMergeHelper
    ::cmtool::CMergeHelper fCMerge;

    /// Input cluster data product producer name label
    std::string fClusterModuleLabel;

    

    /// GeometryUtilities instance
    ::util::GeometryUtilities fGeoU;

    bool fTSep1UseEP;
    
    double fOOCS1MaxAngle; 
    
    int fSDMinHits;
    double fSDSqDistCut;
    
    bool fTSep2UseEP;
    
    double fOOCS2MaxAngle;
    
    int fAI2MinHits;
    bool fAI2Allow180Ambig;
    bool fAI2UseOpeningAngle;
    double fAI2AngleCut;
    double fAI2MinLength;
    
    bool fCOM2UseCOMInPoly;
    bool fCOM2UseCOMInCone;
    bool fCOM2UseCOMNearClus;
    double fCOM2SetLengthReach;

    int fPO2MinHits;
 
    int fPSD2MinHits;
    int fPSD2MaxHits; 
    double fPSD2MinDistSqd; 
    
    
//     TSep1UseEP:  true
//  OOCS1MaxAngle:   20.     # OutOfConeSeparate 1st stage
//     
//  SDMinHits:      10.     # CBAlgoShortestDist # SetMinHits; SetSquaredDistanceCut
//  SDSqDistCut:     5.  
//     
//  # 2nd stage merge
// #TSep2Debug:         false      # CBAlgoTrackSeparate // SetDebug(false); SetVerbose(false); SetUseEP(true);
// # TSep2Verbose:          false
//  TSep2UseEP:        true
//  
// # OOCS2Debug:         false        # CBAlgoOutOfConeSeparate // SetDebug(false); SetVerbose(false); SetMaxAngleSep(20.);
// # OOCS2Verbose:       false
//  OOCS2MaxAngle:      20.
//    
//  AI2MinHits:         50            # CBAlgoAngleIncompat // SetMinHits(50); SetAllow180Ambig(true); SetUseOpeningAngle(false); 
//  AI2Allow180Ambig:   true         #                     // SetAngleCut(10.); SetMinLength(20.); SetDebug(false);
//  AI2UseOpeningAngle: false
//  AI2AngleCut:        10.
//  AI2MinLength:       20.
// # AI2Debug:           false
//  
// # COM2Debug:         false      # CBAlgoCenterOfMass // SetDebug(false); SetVerbose(false); UseCOMInPoly(true);
// # COM2Verbose:       false      # 		 	   //	UseCOMInCone(true); UseCOMNearClus(true); SetLengthReach(3.);
//  COM2UseCOMInPoly:  true
//  COM2UseCOMInCone:  true
//  COM2UseCOMNearClus: true 
//  COM2SetLengthReach:  3.
// 
//  PO2MinHits:        0          # CBAlgoPolyOverlap  // SetMinNumHits(0);
//  
//  PSD2MinHits:        30        # CCBAlgoPolyShortestDist; // SetMinNumHits(30); SetMaxNumHits(9999); SetMinDistSquared(1.); SetDebug(false);
//  PSD2MaxHits:        9999
//  PSD2MinDistSqd:    1
  
   
    
  
  };
}

namespace cluster {

  FuzzyClusterMerger::FuzzyClusterMerger(fhicl::ParameterSet const & p)
  {
    // Declare output data products
    produces< std::vector<recob::Cluster> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    
    // Fill fcl parameter
    fClusterModuleLabel = p.get<std::string>("InputClusterLabel");
    
    fTSep1UseEP = p.get<bool>("TSep1UseEP");
    
    fOOCS1MaxAngle = p.get<double>("OOCS1MaxAngle");
    
    fSDMinHits = p.get<int>("SDMinHits");
    fSDSqDistCut = p.get<double>("SDSqDistCut");
    
    fTSep2UseEP = p.get<bool>("TSep2UseEP");
    
    fOOCS2MaxAngle = p.get<double>("OOCS2MaxAngle");
    
    fAI2MinHits = p.get<int>("AI2MinHits");
    fAI2Allow180Ambig = p.get<bool>("AI2Allow180Ambig");
    fAI2UseOpeningAngle = p.get<bool>("AI2UseOpeningAngle");
    fAI2AngleCut = p.get<double>("AI2AngleCut");
    fAI2MinLength = p.get<double>("AI2MinLength");
    
    
    fCOM2UseCOMInPoly = p.get<bool>("COM2UseCOMInPoly");
    fCOM2UseCOMInCone = p.get<bool>("COM2UseCOMInCone");
    fCOM2UseCOMNearClus = p.get<bool>("COM2UseCOMNearClus");
    fCOM2SetLengthReach = p.get<double>("COM2SetLengthReach");
    
   
    fPO2MinHits = p.get<int>("PO2MinHits");
    
    fPSD2MinHits = p.get<int>("PSD2MinHits");
    fPSD2MaxHits = p.get<int>("PSD2MaxHits");
    fPSD2MinDistSqd = p.get<double>("PSD2MinDistSqd");
 
 
    
    //--- Configure Merging Algorithm ---//
   

    fCMerge.GetManager(0).MergeTillConverge(true);
    //fCMerge.GetManager(0).DebugMode(::cmtool::CMergeManager::kPerIteration);

    // Prohibit algorithms    
    auto prohib_algo_1_1 = new ::cmtool::CBAlgoTrackSeparate;
    prohib_algo_1_1->SetUseEP(fTSep1UseEP);

    auto prohib_algo_1_2 = new ::cmtool::CBAlgoOutOfConeSeparate;
    prohib_algo_1_2->SetMaxAngleSep(fOOCS1MaxAngle);

    auto prohib_algo_1 = new ::cmtool::CBAlgoArray;
    prohib_algo_1->AddAlgo(prohib_algo_1_1,false);
    prohib_algo_1->AddAlgo(prohib_algo_1_2,false);

    fCMerge.GetManager(0).AddSeparateAlgo( prohib_algo_1 );

    // Merge algorithms
    auto merge_algo_1_1 = new ::cmtool::CBAlgoShortestDist;
    merge_algo_1_1->SetMinHits(fSDMinHits);
    merge_algo_1_1->SetSquaredDistanceCut(fSDSqDistCut);

    auto merge_algo_1_2 = new ::cmtool::CBAlgoStartTrack;

    auto merge_algo_1_3 = new ::cmtool::CBAlgoPolyContain;

    auto merge_algo_1 = new ::cmtool::CBAlgoArray;
    merge_algo_1->AddAlgo(merge_algo_1_1,false);
    merge_algo_1->AddAlgo(merge_algo_1_2,false);
    merge_algo_1->AddAlgo(merge_algo_1_3,false);

    fCMerge.GetManager(0).AddMergeAlgo( merge_algo_1 );

    // Configure 2nd stage merging
    //auto& fCMerge.GetManager(1) = GetManager(1);
    fCMerge.GetManager(1).MergeTillConverge(true);
    //fCMerge.GetManager(0).DebugMode(::cmtool::CMergeManager::kPerIteration);
    
    // Prohibit algorithms
    auto prohib_algo_2_1 = new ::cmtool::CBAlgoTrackSeparate;
  //  prohib_algo_2_1->SetDebug(false);
  //  prohib_algo_2_1->SetVerbose(false);
    prohib_algo_2_1->SetUseEP(fTSep2UseEP);

    auto prohib_algo_2_2 = new ::cmtool::CBAlgoOutOfConeSeparate;
  //  prohib_algo_2_2->SetDebug(false);
  //  prohib_algo_2_2->SetVerbose(false);
    prohib_algo_2_2->SetMaxAngleSep(fOOCS2MaxAngle);

    auto prohib_algo_2_3 = new ::cmtool::CBAlgoAngleIncompat;
    prohib_algo_2_3->SetMinHits(fAI2MinHits);
    prohib_algo_2_3->SetAllow180Ambig(fAI2Allow180Ambig);
    prohib_algo_2_3->SetUseOpeningAngle(fAI2UseOpeningAngle);
    prohib_algo_2_3->SetAngleCut(fAI2AngleCut);
    prohib_algo_2_3->SetMinLength(fAI2MinLength);
   // prohib_algo_2_3->SetDebug(false);
    
       
    
    auto prohib_algo_2 = new ::cmtool::CBAlgoArray;
    prohib_algo_2->AddAlgo(prohib_algo_2_1,false);
    prohib_algo_2->AddAlgo(prohib_algo_2_2,false);
    prohib_algo_2->AddAlgo(prohib_algo_2_3,false);

    fCMerge.GetManager(1).AddSeparateAlgo(prohib_algo_2);

    // Merge algorithms
    auto merge_algo_2_1 = new ::cmtool::CBAlgoCenterOfMass;
  //  merge_algo_2_1->SetDebug(false);
  //  merge_algo_2_1->SetVerbose(false);
    merge_algo_2_1->UseCOMInPoly(fCOM2UseCOMInPoly);
    merge_algo_2_1->UseCOMInCone(fCOM2UseCOMInCone);
    merge_algo_2_1->UseCOMNearClus(fCOM2UseCOMNearClus);
    merge_algo_2_1->SetLengthReach(fCOM2SetLengthReach);

    
    auto merge_algo_2_2 = new ::cmtool::CBAlgoPolyOverlap;
    merge_algo_2_2->SetMinNumHits(fPO2MinHits);

    auto merge_algo_2_3 = new ::cmtool::CBAlgoPolyShortestDist;
    merge_algo_2_3->SetMinNumHits(fPSD2MinHits);
    merge_algo_2_3->SetMaxNumHits(fPSD2MaxHits);
    merge_algo_2_3->SetMinDistSquared(fPSD2MinDistSqd);
    merge_algo_2_3->SetDebug(false);

    auto merge_algo_2 = new ::cmtool::CBAlgoArray;
    merge_algo_2->AddAlgo(merge_algo_2_1,false);
    merge_algo_2->AddAlgo(merge_algo_2_2,false);
    merge_algo_2->AddAlgo(merge_algo_2_3,false);
    
    fCMerge.GetManager(1).AddMergeAlgo(merge_algo_2);

    // Prohibit algorithms
    auto priority_algo_2 = new ::cmtool::CPAlgoIgnoreTracks;

    fCMerge.GetManager(1).AddPriorityAlgo(priority_algo_2);

    //
    // FYI there's an algorithm to just-merge-everything if you want to do a simple test (line commented out below)
    //
    //fCMerge.GetManager(0).AddMergeAlgo( new CMAlgoMergeAll );

    
  }
  
  FuzzyClusterMerger::~FuzzyClusterMerger()
  {
    // Clean up dynamic memory and other resources here.
  }
  
  void FuzzyClusterMerger::produce(art::Event & evt)
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
      
      
      
      
      
      
      
  /*    
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

    }*/

    evt.put(std::move(out_clusters));
    evt.put(std::move(out_assn));
  }
}

DEFINE_ART_MODULE(cluster::FuzzyClusterMerger)
