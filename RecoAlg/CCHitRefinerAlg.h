////////////////////////////////////////////////////////////////////////
// CCHitRefiner.h
//
// Cluster Crawler Hit Refiner class
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef CCHITREFINERALG_H
#define CCHITREFINERALG_H

#include "TMath.h"

#include <vector>

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 

#include "Geometry/Geometry.h"
#include "RecoBase/Hit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "RecoAlg/ClusterCrawlerAlg.h"
#include "RecoAlg/CCHitFinderAlg.h"

struct MinuitStruct{
  
    std::vector< std::vector<float> > WireSignals;
    std::vector< std::vector<float> > HitSignals;
    float MatchChisq;

    // Parameters for extrapolating clusters into the vertex region
    // This struct is defined for separate RAT ranges
    struct VtxCluster {
      unsigned short clIndx; // cluster index
      short Wire0;  // The wire index for cluster Time0
      float Time0;  // Time index
      float Slope;  // 
      float RMS;    // average hit RMS
      float Amp; // average hit Amplitude
      float AmpErr; // error on the amplitude
      bool isDS; // cluster is DS of the vertex
    };
    std::vector< VtxCluster > vcl;

};

namespace cluster {

  class CCHitRefinerAlg {
  
  public:

    CCHitRefinerAlg(fhicl::ParameterSet const& pset);
    virtual ~CCHitRefinerAlg();

    void reconfigure(fhicl::ParameterSet const& pset);

    void RunCCHitRefiner(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      CCHitFinderAlg::HitCuts& hitcuts,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl,
      std::vector<ClusterCrawlerAlg::VtxStore>& vtx,
      ClusterCrawlerAlg& fCCAlg);

    bool fRefineHits;   ///< refine hits?
    float fBEChgRat;    ///< Begin/End Charge ratio to determine cluster direction

//    const std::vector< std::vector<float> >& WireSignals() const;

  private:
    
    
    bool prt;
    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;
    
    unsigned short plane;
    std::vector< std::pair<short, short> > WireHitRange;
    unsigned short fFirstWire;
    unsigned short fLastWire;

    unsigned short theVtx;
    // list of clusters associated with theVtx
    std::vector<unsigned short> clBeg;
    std::vector<unsigned short> clEnd;

    // boundaries of the RAT range
    unsigned short loWire;
    unsigned short hiWire;
    unsigned short loTime;
    unsigned short hiTime;
    
    void RefineHits(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl,
      std::vector<ClusterCrawlerAlg::VtxStore>& vtx);

    // define a Region Above Threhold (RAT) range of wires loWire - hiWire
    // and times loTime - hiTime surrounding a vertex in which hit
    // refining will be done
    void FindRATRange(
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl, 
      std::vector<ClusterCrawlerAlg::VtxStore>& vtx);

    // Fills the WireSignals vector in the RAT range
    void FillWireSignals(std::vector<CCHitFinderAlg::CCHit>& allhits);

    // Fits clusters at the boundary of the RAT range
    void RefitClusters(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl, 
      std::vector<ClusterCrawlerAlg::VtxStore>& vtx,
      ClusterCrawlerAlg& fCCAlg);
    
    void PrintSignals();
    
    void FitHitSignals(std::vector<ClusterCrawlerAlg::VtxStore>& vtx);

    // Sets the beginning and end direction of the cluster
    void SetClusterBeginEnd(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl);

  }; // class CCHitRefinerAlg

} // cluster
/*
    inline const std::vector< std::vector<float> >& cluster::CCHitRefinerAlg::WireSignals()
          const { return fWireSignals; };
*/
#endif // ifndef CCHITREFINERALG_H
