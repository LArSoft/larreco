////////////////////////////////////////////////////////////////////////
//
//
// TrajClusterAlg
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALG_H
#define TRAJCLUSTERALG_H


// C/C++ standard libraries
#include <array>
#include <map>
#include <vector>
#include <memory> // std::move()
#include <utility> // std::pair<>

// framework libraries
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 

// LArSoft libraries
#include "SimpleTypesAndConstants/geo_types.h"
#include "Geometry/Geometry.h"
#include "RecoBase/Hit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "RecoAlg/LinFitAlg.h"

namespace cluster {
  
  class TrajClusterAlg {
    public:
    
    // some functions to handle the CTP_t type
    typedef unsigned int CTP_t;
    
    static constexpr unsigned int CTPpad = 1000; // alignment for CTP sub-items
    static CTP_t EncodeCTP
      (unsigned int cryo, unsigned int tpc, unsigned int plane)
      { return cryo * CTPpad*CTPpad + tpc * CTPpad + plane; }
    static CTP_t EncodeCTP(const geo::PlaneID& planeID)
      { return EncodeCTP(planeID.Cryostat, planeID.TPC, planeID.Plane); }
    static geo::PlaneID DecodeCTP(CTP_t CTP)
      { return { CTP / (CTPpad*CTPpad), CTP / CTPpad % CTPpad, CTP % CTPpad }; }

    /// @{
    /// @name Data structures for the reconstruction results

    /// struct of temporary clusters
    struct ClusterStore {
      short ID;         // Cluster ID. ID < 0 = abandoned cluster
//      short ProcCode;   // Processor code for debugging
//      short StopCode;   // code for the reason for stopping cluster tracking
      CTP_t CTP;        // Cryostat/TPC/Plane code
      float BeginWir;   // begin wire
      float BeginTim;   // begin tick
      float BeginAng;   // begin angle
      float BeginChg;   // beginning average charge
			short BeginVtx; 	// ID of the Begin vertex
      float EndWir;     // end wire
      float EndTim;     // end tick
      float EndAng;     // end angle
      float EndChg;     // ending average charge
      short EndVtx;     // ID of the end vertex
      std::vector<unsigned int> tclhits; // hits on the cluster
    }; // ClusterStore
    std::vector<short> inClus;

    /// struct of temporary 2D vertices (end points)
    struct VtxStore {
      float Wire;
			float WireErr;
      float Time;
			float TimeErr;
      unsigned short NClusters;  // = 0 for abandoned vertices
			float ChiDOF;
      short Topo; 			// 1 = US-US, 2 = US-DS, 3 = DS-US, 4 = DS-DS, 5 = Star,
												// 6 = hammer, 7 = vtx3clustermatch, 8 = vtx3clustersplit, 9 = FindTrajVertices
      CTP_t CTP;
      bool Fixed;                 // Vertex position fixed (should not be re-fit)
    };
    
    /// struct of temporary 3D vertices
    struct Vtx3Store {
      std::array<short, 3> Ptr2D; // pointers to 2D vertices in each plane
      float X;                    // x position
      float XErr;                 // x position error
      float Y;                    // y position
      float YErr;                 // y position error
      float Z;                    // z position
      float ZErr;                 // z position error
      short Wire;                 // wire number for an incomplete 3D vertex
      unsigned short CStat;
      unsigned short TPC;
      unsigned short ProcCode;
    };

    TrajClusterAlg(fhicl::ParameterSet const& pset);

    virtual void reconfigure(fhicl::ParameterSet const& pset);
//    void RunTrajClusterAlg(std::vector<art::Ptr<recob::Hit>> const& allhits);
    void RunTrajClusterAlg(art::Event & evt);
    
    std::vector<short> const& GetinClus() const {return inClus; }
    
    /// Returns a constant reference to the clusters found
    std::vector<ClusterStore> const& GetClusters() const { return tcl; }
    
    /// Returns a constant reference to the 2D end points found
    std::vector<VtxStore> const& GetEndPoints() const { return vtx; }
    
    /// Returns a constant reference to the 3D vertices found
    std::vector<Vtx3Store> const& GetVertices() const { return vtx3; }
    
    /// Deletes all the results
    void ClearResults();
    
    
////////////////////////////////////

    private:
    
    art::InputTag fHitFinderModuleLabel; ///< label of module producing input hits
    std::vector<art::Ptr<recob::Hit>> fHits;
    
    short fMode;            ///  StepCrawl mode (0 = turn off)
    short fStepDir;             /// US->DS (1), DS->US (-1)
    float fChgDiffCut;      ///  Max charge difference (Q - Q_ave) / Q_rms
    float fKinkAngCut;     ///  kink angle cut
    float fMaxWireSkip;    ///< max number of wires to skip w/o a signal on them
    float fMaxDeltaJump;   /// Ignore hits have a Delta larger than this value
    float fAveHitResCut;
    float fFitChiCut;
    bool fTagClusters;              ///< tag clusters as shower-like or track-like
    bool fStudyMode;       ///< study cuts
		
//		float fMinAmp;									///< expected minimum signal

    bool fFindTrajVertices;

    float fHitErrFac;   ///< hit time error = fHitErrFac * hit RMS used for cluster fit
//    float fMinHitFrac;
    float fLAClusAngleCut;  ///< call Large Angle Clustering code if > 0
    float fLAClusSlopeCut;
//    float fHitCloseChi;     ///< Hits are considered close if in same multiplet, etc
    float fMergeOverlapAngCut;   ///< angle cut for merging overlapping clusters
    unsigned short fAllowNoHitWire;
		float fVertex2DCut; 	///< 2D vtx -> cluster matching cut (chisq/dof)
    float fVertex2DWireErrCut;
    float fVertex3DCut;   ///< 2D vtx -> 3D vtx matching cut (chisq/dof)

    int fDebugPlane;
    int fDebugWire;  ///< set to the Begin Wire and Hit of a cluster to print
    int fDebugHit;   ///< out detailed information while crawling
    
    bool prt;
    
    art::ServiceHandle<geo::Geometry> geom;
    art::ServiceHandle<util::LArProperties> larprop;
    art::ServiceHandle<util::DetectorProperties> detprop;
    
    std::vector< ClusterStore > tcl; ///< the clusters we are creating
    std::vector< VtxStore > vtx; ///< the endpoints we are reconstructing
    std::vector< Vtx3Store > vtx3; ///< the 3D vertices we are reconstructing
    
    trkf::LinFitAlg fLinFitAlg;

    unsigned int fFirstWire;    ///< the first wire with a hit
    unsigned int fFirstHit;     ///< first hit used
    unsigned int fLastWire;      ///< the last wire with a hit
    unsigned int fCstat;         // the current cryostat
    unsigned int fTpc;         // the current TPC
    CTP_t fCTP;        ///< Cryostat/TPC/Plane code
    unsigned int fPlane;         // the current plane
    unsigned int fNumWires;   // number of wires in the current plane
    unsigned int fMaxTime;    // number of time samples in the current plane
    float fScaleF;     ///< scale factor from Tick/Wire to dx/du

    // vector of pairs of first (.first) and last+1 (.second) hit on each wire
    // in the range fFirstWire to fLastWire. A value of -2 indicates that there
    // are no hits on the wire. A value of -1 indicates that the wire is dead
    std::vector< std::pair<int, int> > WireHitRange;

    std::string fhitsModuleLabel;
    
    float fAveChg;
    float fChgRMS;
    float fAveHitRMS;
    
    // struct for step crawling
    struct TrajPoint {
      std::array<float, 2> HitPos; // Charge weighted position of hits in wire equivalent units
      std::array<float, 2> Pos; // Trajectory position in wire equivalent units
      std::array<float, 2> Dir; // Direction
      float Ang;                // Trajectory angle (-pi, +pi)
      float AngErr;             // Trajectory angle error
      float AveHitRes;          // Average normalized hit residual (< 1 => Good)
      float Chg;                // Charge
      float AveChg;             // Average charge of last ~20 TPs
      float ChgDiff;            //  = (Chg - fAveChg) / fChgRMS
      float Delta;              // Deviation between trajectory and hits (WSE)
      float DeltaRMS;           // RMS deviation of all TP Deltas
      unsigned short NumTPsFit; // Number of trajectory points fitted to make this point
      unsigned short Step;      // Step number at which this TP was created
      float FitChi;             // Chi/DOF of the fit
      unsigned short NumClose;  // Number of hits not used but close to the TP
      bool UsedNotCloseHit;        // true if a hit was was used but failed the initial acceptance cut
      std::vector<unsigned int> Hits; // vector of fHits indices
      // default constructor
      TrajPoint() {
        Ang = 0; AngErr = 0.5; AveHitRes = 0; Chg = 0; ChgDiff = 0; AveChg = 0; Delta = 0; DeltaRMS = 0.2;
        NumTPsFit = 2; Step = 0; FitChi = 0; NumClose = 0; UsedNotCloseHit = false; Hits.clear();
      }
    };
    
    // associated information for the trajectory
    struct Trajectory {
      CTP_t CTP;                      ///< Cryostat, TPC, Plane code
      short StepDir;                 /// -1 = going US (CC proper order), 1 = going DS
      unsigned short ClusterIndex;   ///< Not needed?
      unsigned short ProcCode;       ///< USHRT_MAX = abandoned trajectory
      int TruPDG;                    ///< MC truth
      int TruKE;                     ///< MeV
      bool IsPrimary;                ///< MC truth
      bool IsSecondary;              ///< MC truth
      std::array<short, 2> Vtx;      ///< Index of 2D vertex
      std::vector<TrajPoint> Pts;    ///< Trajectory points
      unsigned short NHits;          ///< Number of hits in all TPs
      Trajectory() {
        CTP = 0; StepDir = 0; ClusterIndex = USHRT_MAX; ProcCode = 900; Pts.clear();
        Vtx[0] = -1; Vtx[1] = -1; TruPDG = 0; TruKE = 0; IsPrimary = false; IsSecondary = false; NHits = 0;
      }
    };
    Trajectory work;      ///< trajectory under construction
    //  trajectories
    std::vector<Trajectory> allTraj; ///< vector of all trajectories in each plane
    std::vector<short> inTraj;
    //  trial    trajectories
    std::vector<std::vector<Trajectory>> trial; ///< vector of all trajectories for all trials in one plane
    std::vector<std::vector<short>> inTrialTraj;
    
    struct TjPairHitShare {
      // Trajectories in two different trials that share hits
      unsigned short iTrial;
      unsigned short iTj;
      unsigned short jTrial;
      unsigned short jTj;
      unsigned short nSameHits;
    };
    std::vector<TjPairHitShare> tjphs;
    
    
    // runs the TrajCluster algorithm on one plane specified by the calling routine
    // (which should also have called GetHitRange)
    void RunStepCrawl();
    void GetHitRange();
    void HitSanityCheck();
    bool TrajHitsOK(unsigned int iht, unsigned int jht);
    // Find all trajectories in the plane using the pre-defined step direction, cuts, etc.
    // Fills the allTraj vector
    void ReconstructAllTraj();

    
/////// old stuff below
    // Main stepping/crawling routine
    void StepCrawl(bool& success);
    // Calculate the number of missed steps that should be expected for the given
    // trajectory point direction
    unsigned short SetMissedStepCut(TrajPoint const& tp);
    // Start a trajectory going from fromHit to (toWire, toTick)
    void StartTraj(float fromWire, float fromTick, float toWire, float toTick);
    // Returns the charge weighted wire, time position of all hits in the multiplet
    // of which hit is a member
    void HitMultipletPosition(unsigned int hit, float& hitTick);
    // Add hits to the trajectory point
    void AddTrajHits(TrajPoint& tp, bool& SignalPresent);
    // Set inTraj = -3 for all hits in the work trajectory
    void FlagWorkHits(short flag);
    // Checks to see if any hit has inTraj = -3 when it shouldn't
    bool CheckAllHitsFlag();
    void FindHit(std::string someText, unsigned int iht); // temp routine
    void FlagSet(std::string someText, unsigned int iht); // temp routine
    // analyze the sat vector to construct a vector of trajectories that is
    // the best
    void AnalyzeTrials();
    void CountSameHits(std::vector<unsigned int>& iHitVec, std::vector<unsigned int>& jHitVec, unsigned short& nSameHits);
    void AdjudicateTrials(bool& reAnalyze);
    void MergeTrajPair(unsigned short ipr, bool& reAnalyze);
    // Split the allTraj trajectory itj at position pos into two trajectories
    // with an optional vertex assignment
    bool SplitAllTraj(unsigned short itj, unsigned short pos, unsigned short ivx);
    // Make clusters from all trajectories in allTraj
    void MakeAllTrajClusters();
    // Push the work trajectory into allTraj
    void StoreWork();
    // Reverse the work trajectory
    void ReverseTraj(Trajectory& tj);
    void UpdateWork(bool& success);
    void UpdateWorkChgDiff();
    void UpdateWorkDelta();
    void FitTraj();
    float TrajLength(Trajectory& tj);
    bool GottaKink();
    void FitTrajMid(unsigned short fromIndex, unsigned short toIndex, TrajPoint& tp);
    void CheckTrajEnd();
    void TrajClosestApproach(Trajectory const& tj, float x, float y, unsigned short& iClosePt, float& Distance);
    float PointTrajDOCA(float wire, float time, TrajPoint const& tp);
    void TrajSeparation(Trajectory& iTj, Trajectory& jTj, std::vector<float>& tSep);
    void FindTrajVertices();
    unsigned short VtxClusterEnd(unsigned short ivx, unsigned short icl);
    void PutTrajHitsInVector(Trajectory const& tj, std::vector<unsigned int>& hitVec);
    bool TrajIntersection(TrajPoint& tp1, TrajPoint tp2, float& x, float& y);
    float TrajPointSeparation(TrajPoint& tp1, TrajPoint& tp2);
    void PrintTrajectory(Trajectory& tj ,unsigned short tPoint);
    void PrintAllTraj(unsigned short itj, unsigned short ipt);
    void PrintHeader();
    void PrintTrajPoint(unsigned short ipt, short dir, TrajPoint tp);
    // Tag as shower-like or track-like
    void TagClusters();
    void FillTrajTruth();


    // ************** 2D vertex routines *******************
/*
    // Find 2D vertices
    void FindVertices();
    // Find 2D star topology vertices
    void FindStarVertices();
    // check a vertex (vw, fvt) made with clusters it1, and it2 against the
    // vector of existing clusters
    void ChkVertex(float fvw, float fvt, unsigned short it1, unsigned short it2, short topo);
    // try to attach a cluster to an existing vertex
    void ClusterVertex(unsigned short it2);
    // try to attach a cluster to a specified vertex
    void VertexCluster(unsigned short ivx);
    // Refine cluster ends near vertices
    void RefineVertexClusters(unsigned short ivx);
    // Split clusters that cross a vertex
    bool VtxClusterSplit();
    // returns true if a vertex is encountered while crawling
    bool CrawlVtxChk(unsigned short kwire);
    // returns true if this cluster is between a vertex and another
    // cluster that is associated with the vertex
    bool CrawlVtxChk2();
    // use a vertex constraint to start a cluster
    void VtxConstraint(unsigned short iwire, unsigned int ihit,
      unsigned short jwire, unsigned int& useHit, bool& doConstrain);
    // fit the vertex position
    void FitVtx(unsigned short iv);
    // weight and fit all vertices
    void FitAllVtx(CTP_t inCTP);

 // ************** 3D vertex routines *******************

    // match vertices between planes
    void VtxMatch(geo::TPCID const& tpcid);
    // Match clusters to endpoints using 3D vertex information
    void Vtx3ClusterMatch(geo::TPCID const& tpcid);
    // split clusters using 3D vertex information
    void Vtx3ClusterSplit(geo::TPCID const& tpcid);
		// look for a long cluster that stops at a short cluster in two views
		void FindHammerClusters();
*/
    // ************** utility routines *******************
/*
    // inits everything
    void CrawlInit();
    // inits the cluster stuff
    void ClusterInit();
    // fills the wirehitrange vector for the supplied Cryostat/TPC/Plane code
    void GetHitRange(CTP_t CTP);
    // Stores cluster information in a temporary vector
    bool TmpStore();
    // Gets a temp cluster and puts it into the working cluster variables
    void TmpGet(unsigned short it1);
    // Does just what it says
    void CalculateAveHitWidth();
    // Shortens the fcl2hits, chifits, etc vectors by the specified amount
    void FclTrimUS(unsigned short nTrim);
    // Splits a cluster into two clusters at position pos. Associates the
    // new clusters with a vertex
    bool SplitCluster(unsigned short icl, unsigned short pos, unsigned short ivx);
    // Prints cluster information to the screen
    void PrintClusters();
    // check for a signal on all wires between two points
    bool ChkSignal(unsigned short wire1, float time1, unsigned short wire2, float time2);
    // returns an angle-dependent scale factor for weighting fits, etc
    float AngleFactor(float slope);
    // calculate the kink angle between hits 0-2 and 3 - 5 on the leading edge of
    // the cluster under construction
    float EndKinkAngle();
    /// Returns true if there are no duplicates in the hit list for next cluster
    bool CheckHitDuplicates
      (std::string location, std::string marker = "") const;
    // Find the distance of closest approach between the end of a cluster and a (wire,tick) position
    float DoCA(short icl, unsigned short end, float vwire, float vtick);
    // Find the Chisq/DOF between the end of a cluster and a (wire,tick) position
    float ClusterVertexChi(short icl, unsigned short end, unsigned short ivx);
    // Find the Chisq/DOF between a point and a vertex
    float PointVertexChi(float wire, float tick, unsigned short ivx);
    
    /// Returns a pair of first and past-the-last index
    /// of all the contiguous hits belonging to the same multiplet
    std::pair<size_t, size_t> FindHitMultiplet(size_t iHit) const;

    void CheckHitClusterAssociations();
    
    /// Returns whether the two hits belong to the same multiplet
    static bool areInSameMultiplet(recob::Hit const& first_hit, recob::Hit const& second_hit);
*/
    
  }; // class TrajClusterAlg


} // namespace cluster

#endif // ifndef TRAJCLUSTERALG_H
