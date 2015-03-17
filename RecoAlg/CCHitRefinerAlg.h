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

struct MinuitStruct{
  
    std::vector< std::vector<float> > WireSignals; // wire signals - data
    std::vector< std::vector<float> > HitSignals; // hit signals - fit
    std::vector< double > WireWght; // fcnA weighting by wire
    unsigned short WireIndex; // index of the wire used in fcnW
    double fcnVal;

    // Parameters for extrapolating clusters into the vertex region
    // This struct is defined for separate RAT ranges
    struct VtxCluster {
      short tclID; // Cluster ID in the tcl struct
      short Wire0;  // The wire index for cluster Time0
      float Time0;  // Time on Wire0
      float Slope;  // 
      float SlopeErr;  // slope error from FitClusterMid
      float RMS;    // average hit RMS
      std::vector<float> Time; // Hit time. < 0 means no hit expected
      std::vector<float> Amp; // Hit amplitude on each wire
      std::vector<float> AmpErr; // error on the amplitude
      std::vector<float> HitChiDOF; // Chisq of the fit to the hit amplitudes
      bool isDS; // cluster is DS of the vertex
    };
    std::vector< VtxCluster > vcl;

};

namespace hit {

  class CCHitRefinerAlg {
  
  public:
    using ClusterCrawlerAlg = cluster::ClusterCrawlerAlg;
    using CCHitFinderAlg = hit::CCHitFinderAlg;

    CCHitRefinerAlg(fhicl::ParameterSet const& pset);
    virtual ~CCHitRefinerAlg();

    void reconfigure(fhicl::ParameterSet const& pset);

    void RunCCHitRefiner(
      std::vector<recob::Hit>& allhits,
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
    
    /// association between clusters and hits
    ClusterCrawlerAlg::HitInCluster_t HitInCluster;
    
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
      std::vector<recob::Hit>& allhits,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl,
      std::vector<ClusterCrawlerAlg::VtxStore>& vtx);

    // define a Region Above Threhold (RAT) range of wires loWire - hiWire
    // and times loTime - hiTime surrounding a vertex in which hit
    // refining will be done
    void FindRATRange(
      std::vector<recob::Hit>& allhits,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl, 
      std::vector<ClusterCrawlerAlg::VtxStore>& vtx,
      bool& SkipIt);

    // Fills the WireSignals vector in the RAT range
    void FillWireSignals(std::vector<recob::Hit>& allhits);

    // Fits clusters at the boundary of the RAT range
    void FillVcl(std::vector<recob::Hit>& allhits,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl, 
      std::vector<ClusterCrawlerAlg::VtxStore>& vtx);
    
    void PrintSignals();
    
    // fits the vertex position and average hit amplitudes in the RAT range
    void FitVtxPos(std::vector<ClusterCrawlerAlg::VtxStore>& vtx);
    // fits the hit amplitudes on each wire in the RAT range
    void FitHitAmplitudes(std::vector<ClusterCrawlerAlg::VtxStore>& vtx);

    // Sets the beginning and end direction of the cluster
    void SetClusterBeginEnd(std::vector<recob::Hit>& allhits,
      std::vector<ClusterCrawlerAlg::ClusterStore>& tcl);

    void Printvcl(std::vector<ClusterCrawlerAlg::VtxStore>& vtx);


  }; // class CCHitRefinerAlg

} // hit
/*
    inline const std::vector< std::vector<float> >& cluster::CCHitRefinerAlg::WireSignals()
          const { return fWireSignals; };
*/
#endif // ifndef CCHITREFINERALG_H
