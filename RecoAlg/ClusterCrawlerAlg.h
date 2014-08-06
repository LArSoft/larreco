////////////////////////////////////////////////////////////////////////
// ClusterCrawlerAlg.h
//
// ClusterCrawlerAlg class
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef CLUSTERCRAWLERALG_H
#define CLUSTERCRAWLERALG_H

#include "TMath.h"

#include <vector>

#include "fhiclcpp/ParameterSet.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 

#include "SimpleTypesAndConstants/geo_types.h"
#include "Geometry/Geometry.h"
#include "RecoBase/Hit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RecoAlg/CCHitFinderAlg.h"

namespace recob { class Hit; }

namespace cluster {


  class ClusterCrawlerAlg {
    public:

    // some functions to handle the CTP_t type
    typedef unsigned int CTP_t;
    
    static constexpr unsigned int CTPpad = 1000; // alignment for CTP sub-items
    static CTP_t EncodeCTP(unsigned int cryo, unsigned int tpc, unsigned int plane)
      { return cryo * CTPpad*CTPpad + tpc * CTPpad + plane; }
    static CTP_t EncodeCTP(const geo::PlaneID& planeID)
      { return EncodeCTP(planeID.Cryostat, planeID.TPC, planeID.Plane); }
    static geo::PlaneID DecodeCTP(CTP_t CTP)
      { return { CTP / (CTPpad*CTPpad), CTP / CTPpad % CTPpad, CTP % CTPpad }; }
    
    // struct of temporary clusters
    struct ClusterStore {
      short ID;         // Cluster ID. ID < 0 = abandoned cluster
      short ProcCode;   // Processor code for debugging
      short StopCode;   // code for the reason for stopping cluster tracking
      CTP_t CTP;        // Cryostat/TPC/Plane code
      float BeginSlp;   // beginning slope (= DS end = high wire number)
      float BeginSlpErr; // error
      unsigned short BeginWir;   // begin wire
      float BeginTim;   // begin time
      float BeginChg;   // beginning average charge
      float EndSlp;     // end slope (= US end = low  wire number)
      float EndSlpErr;  // error
      unsigned short EndWir;     // end wire
      float EndTim;     // ending time
      float EndChg;     // ending average charge
      short BeginVtx;   // ID of the begin vertex
      short EndVtx;     // ID of the end vertex
      std::vector<unsigned short> tclhits; // hits on the cluster
    }; // ClusterStore
    std::vector< ClusterStore > tcl;

    // struct of temporary 2D vertices
    struct VtxStore {
      float Wire;
      float Time;
      float Wght;       // Wght < 0 for abandoned vertices 
      short Topo;
      CTP_t CTP;
    };
    std::vector< VtxStore > vtx;
    
    // struct of temporary 3D vertices
    struct Vtx3Store {
      std::array<short, 3> Ptr2D; // pointers to 2D vertices in each plane
      float X;                    // x position
      float Y;                    // y position
      float Z;                    // z position
      short Wire;                 // wire number for an incomplete 3D vertex
      unsigned short CStat;
      unsigned short TPC;
    };
    std::vector< Vtx3Store > vtx3;

    ClusterCrawlerAlg(fhicl::ParameterSet const& pset);
    virtual ~ClusterCrawlerAlg();

    void reconfigure(fhicl::ParameterSet const& pset);
    void RunCrawler(std::vector<CCHitFinderAlg::CCHit>& allhits);
    
    unsigned short fNumPass;                 ///< number of passes over the hit collection
    std::vector<unsigned short> fMaxHitsFit; ///< Max number of hits fitted
    std::vector<unsigned short> fMinHits;    ///< Min number of hits to make a cluster
    std::vector<unsigned short> fNHitsAve;   ///< number of US hits used to compute fAveChg
                                    ///< set to > 2 to do a charge fit using fNHitsAve hits
    std::vector<float> fChiCut;     ///< stop adding hits to clusters if chisq too high
    std::vector<float> fKinkChiRat; ///< Max consecutive chisq increase for the last 
                                    ///< 3 hits on the cluster
    std::vector<float> fKinkAngCut; ///< kink angle cut made after fKinkChiRat
    std::vector<float> fChgCut;     ///< charge difference cut for adding a hit to a cluster
    std::vector<unsigned short> fMaxWirSkip; ///< max number of wires that can be skipped while crawling
    std::vector<unsigned short> fMinWirAfterSkip; ///< minimum number of hits on consecutive wires
                                    ///< after skipping
    std::vector<bool> fDoMerge;     ///< try to merge clusters?
    std::vector<float> fTimeDelta;  ///< max time difference for matching
    std::vector<float> fMergeChgCut;  ///< max charge ratio for matching
    std::vector<bool> fFindVertices;    ///< run vertexing code after clustering?
    std::vector<bool> fLACrawl;    ///< Crawl Large Angle clusters on pass?

    // global cuts and parameters 
    float fHitErrFac;   ///< hit time error = fHitErrFac * hit RMS
                        ///< used for cluster fit
    float fLAClusAngleCut;  ///< call Large Angle Clustering code if > 0
    float fLAClusSlopeCut;
    float fHitMergeChiCut; ///< Merge cluster hit-multiplets if the separation chisq
                             ///< is < cut. Set < 0 for no merging
    bool fMergeGhostClusters; ///< Merge clusters if they share hits in a muliplet
                           ///< in the same hit multiplet. Set < 0 for no merging
    unsigned short fAllowNoHitWire; 
    float fVertex3DCut;   ///< 2D vtx -> 3D vtx matching cut (cm)
    short fDebugPlane;
    short fDebugWire;  ///< set to the Begin Wire and Hit of a cluster to print
    short fDebugHit;   ///< out detailed information while crawling

    // fills a wirehitrange vector for the supplied Cryostat/TPC/Plane code
    void GetHitRange(std::vector<CCHitFinderAlg::CCHit>& allhits,
      CTP_t CTP, 
      std::vector< std::pair<short, short> >& wirehitrange,
      unsigned short& firstwire, unsigned short& lastwire);

    // Fits the middle of a temporary cluster it1 using hits iht to iht + nhit
    void FitClusterMid(std::vector<CCHitFinderAlg::CCHit>& allhits, 
      std::vector<ClusterStore>& tcl, unsigned short it1, unsigned short iht,
      short nhit);

    // these variables define the cluster used during crawling
    float clpar[2];     ///< cluster parameters for the current fit with
                        ///< origin at the US wire on the cluster
    float clparerr[2];  ///< cluster parameter errors
    float clChisq;     ///< chisq of the current fit
    float fAveChg;  ///< average charge at leading edge of cluster
    float fChgSlp;  ///< slope of the  charge vs wire 

////////////////////////////////////

    private:
    
    bool prt;
    bool vtxprt;
    unsigned short NClusters;
    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;

    float clBeginSlp;  ///< begin slope (= DS end = high wire number)
    float clBeginSlpErr; 
    unsigned short clBeginWir;  ///< begin wire
    float clBeginTim;  ///< begin time
    float clBeginChg;  ///< begin average charge
    float clEndSlp;    ///< slope at the end   (= US end = low  wire number)
    float clEndSlpErr; 
    unsigned short clEndWir;    ///< begin wire
    float clEndTim;    ///< begin time
    float clEndChg;    ///< end average charge
    short clStopCode;     ///< code for the reason for stopping cluster tracking
                        ///< 0 = no signal on the next wire
                        ///< 1 = skipped too many occupied/dead wires
                        ///< 2 = failed the fMinWirAfterSkip cut
                        ///< 3 = ended on a kink. Fails fKinkChiRat
                        ///< 4 = failed the fChiCut cut
                        ///< 5 = cluster split by VtxClusterSplit
                        ///< 6 = stop at a vertex
    short clProcCode;     ///< Processor code = pass number
                        ///< +   10 ChkMerge
                        ///< +   20 ChkMerge with overlapping hits
                        ///< +  100 ChkMerge12
                        ///< +  200 ClusterFix
                        ///< +  300 LACrawlUS
                        ///< + 1000 VtxClusterSplit
                        ///< + 2000 failed pass N cuts but passes pass N=1 cuts
                        ///< + 3000 Cluster hits merged
                        ///< +10000 modified by CCHitRefiner
    CTP_t clCTP;        ///< Cryostat/TPC/Plane code
    
    unsigned short fFirstWire;    ///< the first wire with a hit
    unsigned short fFirstHit;     ///< first hit used
    unsigned short fLastWire;      ///< the last wire with a hit
    float fAveWid;  ///< average hit width at the leading edge of the cluster
    unsigned int cstat;         // the current cryostat
    unsigned int tpc;         // the current TPC
    unsigned int plane;         // the current plane
    unsigned short fNumWires;   // number of wires in the current plane
    unsigned short fMaxTime;    // number of time samples in the current plane
    float fScaleF;     ///< scale factor from Tick/Wire to dx/du
    
    unsigned short pass;

    // vector of pairs of first (.first) and last+1 (.second) hit on each wire
    // in the range fFirstWire to fLastWire. A value of -2 indicates that there
    // are no hits on the wire. A value of -1 indicates that the wire is dead
    std::vector< std::pair<short, short> > WireHitRange;
    
    std::vector<unsigned short> fcl2hits;  ///< vector of hits used in the cluster
    std::vector<float> chifits;   ///< fit chisq for monitoring kinks, etc

    std::string fhitsModuleLabel;
    

    // ******** crawling routines *****************

    // Loops over wires looking for seed clusters
    void ClusterLoop(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx);
    // Finds a hit on wire kwire, adds it to the cluster and re-fits it
    void AddHit(std::vector<CCHitFinderAlg::CCHit>& allhits, 
      unsigned short kwire, bool& HitOK, bool& SigOK);
    // Finds a hit on wire kwire, adds it to a LargeAngle cluster and re-fits it
    void AddLAHit(std::vector<CCHitFinderAlg::CCHit>& allhits, 
      std::vector<VtxStore>& vtx,
      unsigned short kwire, bool& ChkCharge, bool& HitOK, bool& SigOK);
    // Fits the cluster hits in fcl2hits to a straight line
    void FitCluster(std::vector<CCHitFinderAlg::CCHit>& allhits);
    // Fits the charge of the cluster hits in fcl2hits
    void FitClusterChg(std::vector<CCHitFinderAlg::CCHit>& allhits);
    // Crawls along a trail of hits UpStream
    void CrawlUS(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<VtxStore>& vtx);
    // Crawls along a trail of hits UpStream - Large Angle version
    void LACrawlUS(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<VtxStore>& vtx);
    // returns true if the cluster is nearing a vertex on the next wire

    // ************** cluster merging routines *******************

    // Compares two cluster combinations to see if they should be merged
    void ChkMerge(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl);
    // Checks merge for cluster cl2 within the bounds of cl1
    void ChkMerge12(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl, unsigned short it1, unsigned short it2, bool& didit);
    // Merges clusters cl1 and cl2
    void DoMerge(std::vector<CCHitFinderAlg::CCHit>& allhits, 
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
      unsigned short it1, unsigned short it2, 
      short ProcCode);
    // Merges clusters that share a significant fraction of hits in
    // a multiplet
    void MergeGhostClusters(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl);

    // ************** vertex routines *******************

    // Find 2D vertices
    void FindVertices(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx);
    // check a vertex (vw, fvt) made with clusters it1, and it2 against the
    // vector of existing clusters
    void ChkVertex(std::vector<CCHitFinderAlg::CCHit>& allhits,
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
        float fvw, float fvt, unsigned short it1, unsigned short it2, short topo);
    // try to attach a cluster to an existing vertex
    void ClusterVertex(std::vector<CCHitFinderAlg::CCHit>& allhits, 
        std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
        unsigned short it2);
    // Split clusters that cross a vertex
    void VtxClusterSplit(std::vector<CCHitFinderAlg::CCHit>& allhits,
       std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx);
    // returns true if a vertex is encountered while crawling
    bool CrawlVtxChk( 
      std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<VtxStore>& vtx, unsigned short kwire);
    // use a vertex constraint to start a cluster
    void VtxConstraint(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<VtxStore>& vtx, unsigned short iwire, unsigned short ihit,
      unsigned short jwire, unsigned short& useHit, bool& doConstrain);
    // fit the vertex position
    void FitVtx(std::vector<ClusterStore>& tcl,
        std::vector<VtxStore>& vtx, unsigned short iv, float& ChiDOF);
    // match vertices between planes
    void VtxMatch(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
      std::vector<Vtx3Store>& vtx3, unsigned int cstat, unsigned int tpc);
    // split clusters using 3D vertex information
    void Vtx3ClusterSplit(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx,
      std::vector<Vtx3Store>& vtx3, unsigned int cstat, unsigned int tpc);

    // ************** hit routines *******************

    // merge the hits in a multiplet into one hit
    void MergeHits(std::vector<CCHitFinderAlg::CCHit>& allhits,
      unsigned short theHit);
    // merge all unused cluster hits in a multiplet into one hit
    void MergeClusterHits(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl);   

    // ************** utility routines *******************

    // inits everything
    void CrawlInit();
    // Sorts clusters in tcl by decreasing number of hits, while ignoring
    // abandoned clusters
    void cl2SortByLength(std::vector<ClusterStore>& tcl,
        std::map<unsigned short, unsigned short>& sortindex);
    // Stores cluster information in a temporary vector
    void TmpStore(std::vector<CCHitFinderAlg::CCHit>& allhits, 
      std::vector<ClusterStore>& tcl);
    // Gets a temp cluster and puts it into the working cluster variables
    void TmpGet(std::vector<ClusterStore>& tcl, unsigned short it1);
    // Splits a cluster into two clusters at position pos. Associates the
    // new clusters with a vertex
    void SplitCluster(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl, unsigned short icl, unsigned short pos,
      unsigned short ivx);
    // Prints cluster information to the screen
    void PrintClusters(std::vector<CCHitFinderAlg::CCHit>& allhits,
      std::vector<ClusterStore>& tcl, std::vector<VtxStore>& vtx);
    // check for a signal on all wires between two points
    void ChkSignal(std::vector<CCHitFinderAlg::CCHit>& allhits,
      unsigned short wire1, float time1, unsigned short wire2, float time2, bool& SigOK);
    // line fitter
    void LinFit(std::vector<float>& x, std::vector<float>& y, 
      std::vector<float>& ey2, float& Intercept, float& Slope, 
      float& InterceptError, float& SlopeError, float& ChiDOF);
    float AngleFactor(float slope);

  }; // class ClusterCrawlerAlg


} // namespace cluster

#endif // ifndef CLUSTERCRAWLERALG_H
