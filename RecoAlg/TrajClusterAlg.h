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
#include "art/Framework/Services/Optional/TFileService.h"

// LArSoft libraries
#include "SimpleTypesAndConstants/geo_types.h"
#include "Geometry/Geometry.h"
#include "RecoBase/Hit.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "RecoAlg/LinFitAlg.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>

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
      unsigned short PDG; // PDG-like code shower-like or line-like
      unsigned short ParentCluster;
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
    
    /// Returns (and loses) the collection of reconstructed hits
//    std::vector<recob::Hit>&& YieldHits() { return std::move(fHits); }
    std::vector<art::Ptr<recob::Hit>> const& YieldHits() const { return fHits; }
    
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
    short fNPtsAve;         /// number of points to find AveChg
    std::vector<unsigned short> fMinNPtsFit; ///< Reconstruct in two passes
    float fMultHitSep;      ///< preferentially "merge" hits with < this separation
    float fTP3ChiCut;       ///<
    float fChgDiffCut;      ///  Max charge difference (Q - Q_ave) / Q_rms
    float fKinkAngCut;     ///  kink angle cut
    float fMaxWireSkip;    ///< max number of wires to skip w/o a signal on them
    float fMaxDeltaJump;   /// Ignore hits have a Delta larger than this value
    float fProjectionErrFactor;
    bool fTagAllTraj;              ///< tag clusters as shower-like or track-like
    float fMaxTrajSep;     ///< max trajectory point separation for making showers
    bool fStudyMode;       ///< study cuts
    bool fShowerStudy;    ///< study shower identification cuts
		
//		float fMinAmp;									///< expected minimum signal

    bool fFindTrajVertices;

    float fHitErrFac;   ///< hit time error = fHitErrFac * hit RMS used for cluster fit
    float fLargeAngle;  ///< (degrees) call Large Angle Clustering code if TP angle > this value
    float fLAClusSlopeCut;
    float fMergeOverlapAngCut;   ///< angle cut for merging overlapping clusters
    unsigned short fAllowNoHitWire;
		float fVertex2DCut; 	///< 2D vtx -> cluster matching cut (chisq/dof)
    float fVertex2DWireErrCut;
    float fVertex3DCut;   ///< 2D vtx -> 3D vtx matching cut (chisq/dof)

    int fDebugPlane;
    int fDebugWire;  ///< set to the Begin Wire and Hit of a cluster to print
    int fDebugHit;   ///< out detailed information while crawling
    
    TH2F *fnHitsPerTP_Angle[3];
    TProfile *fnHitsPerTP_AngleP[3];
    
    TH1F *fShowerNumTrjint;
    TH2F *fShowerTheta_Sep;
    TH1F *fShowerDVtx;
    TH2F *fShowerDVtx_Sep;
    
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
    unsigned short fPass;
    CTP_t fCTP;        ///< Cryostat/TPC/Plane code
    unsigned int fPlane;         // the current plane
    unsigned int fNumWires;   // number of wires in the current plane
    float fMaxTime;    // number of time samples in the current plane
    float fScaleF;     ///< scale factor from Tick/Wire to dx/du

    // vector of pairs of first (.first) and last+1 (.second) hit on each wire
    // in the range fFirstWire to fLastWire. A value of -2 indicates that there
    // are no hits on the wire. A value of -1 indicates that the wire is dead
    std::vector< std::pair<int, int> > WireHitRange;

    std::string fhitsModuleLabel;
    
    // variables for step crawling - updated at each TP
    float fHitChgRMS, fDefaultHitChgRMS;
    float fAveHitRMS, fDefaultAveHitRMS;
    
    // struct for step crawling
    struct TrajPoint {
      std::array<float, 2> HitPos; // Charge weighted position of hits in wire equivalent units
      std::array<float, 2> Pos; // Trajectory position in wire equivalent units
      std::array<float, 2> Dir; // Direction
      float HitsTimeErr2;       // Uncertainty^2 in the time component using all hits
      float Ang;                // Trajectory angle (-pi, +pi)
      float AngErr;             // Trajectory angle error
      float Chg;                // Charge
      float AveChg;             // Average charge of last ~20 TPs
      float ChgDiff;            //  = (Chg - fAveChg) / fChgRMS
      float Delta;              // Deviation between trajectory and hits (WSE)
      float TP3Chi;             // Chisq fit of TP, TP-1 and TP+1
      float DeltaRMS;           // RMS deviation of all TP Deltas
      unsigned short NTPsFit; // Number of trajectory points fitted to make this point
      unsigned short Step;      // Step number at which this TP was created
      float FitChi;             // Chi/DOF of the fit
      unsigned short NCloseNotUsed;   //NUmber of close hits that weren't used in the TP
      std::vector<unsigned int> Hits; // vector of fHits indices
      // default constructor
      TrajPoint() {
        Ang = 0; AngErr = 0.5; TP3Chi = 0; Chg = 0; ChgDiff = 0.1; AveChg = 0; Delta = 0;
        DeltaRMS = 0.02; NTPsFit = 2; Step = 0; FitChi = 0; NCloseNotUsed = 0; Hits.clear();
      }
    };
    
    // associated information for the trajectory
    struct Trajectory {
      CTP_t CTP;                      ///< Cryostat, TPC, Plane code
      unsigned short Pass;            ///< the pass on which it was created
      short StepDir;                 /// -1 = going US (CC proper order), 1 = going DS
      unsigned short ClusterIndex;   ///< Index not the ID...
      unsigned short ProcCode;       ///< USHRT_MAX = abandoned trajectory
      float AveTP3Chi;               ///< average of all TP
      unsigned short PDG;            ///< shower-like or line-like
      unsigned short ParentTraj;     ///< index of the parent (if PDG = 12)
      int TruPDG;                    ///< MC truth
      int TruKE;                     ///< MeV
      bool IsPrimary;                ///< MC truth
      std::array<short, 2> Vtx;      ///< Index of 2D vertex
      std::vector<TrajPoint> Pts;    ///< Trajectory points
      unsigned short NHits;          ///< Number of hits in all TPs
      Trajectory() {
        CTP = 0; Pass = 0; PDG = 0; StepDir = 0; ClusterIndex = USHRT_MAX; ProcCode = 900; Pts.clear();
        Vtx[0] = -1; Vtx[1] = -1; TruPDG = 0; TruKE = 0; IsPrimary = false; NHits = 0; ParentTraj = USHRT_MAX;
      }
    };
    Trajectory work;      ///< trajectory under construction
    //  trajectories
    std::vector<Trajectory> allTraj; ///< vector of all trajectories in each plane
    std::vector<short> inTraj;
    //  trial    trajectories
    std::vector<std::vector<Trajectory>> trial; ///< vector of all trajectories for all trials in one plane
    std::vector<std::vector<short>> inTrialTraj;

    // Trajectory "intersections" used to search for superclusters (aka showers)
    struct TrjInt {
      unsigned short itj1;
      unsigned short ipt1;
      unsigned short itj2;
      unsigned short ipt2;
      float sep2;   // separation^2 at closest point
      float dang;   // opening angle at closest point
      float vw;     // intersection wire
      float vt;     // intersection time
    };
    std::vector<TrjInt> trjint;
    // vectors of clusters of trajectories => a set of somewhat intersecting trajectories
    std::vector<std::vector<unsigned short>> ClsOfTrj;

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
    void FindDefaults();
    void HitSanityCheck();
    bool TrajHitsOK(unsigned int iht, unsigned int jht);
    // Find all trajectories in the plane using the pre-defined step direction, cuts, etc.
    // Fills the allTraj vector
    void ReconstructAllTraj();

    // Main stepping/crawling routine
    void StepCrawl(bool& success);
    // Calculate the number of missed steps that should be expected for the given
    // trajectory point direction
    unsigned short SetMissedStepCut(TrajPoint const& tp);
    // Start a trajectory going from fromHit to (toWire, toTick)
    void StartTraj(float fromWire, float fromTick, float toWire, float toTick);
    // Make a bare trajectory point that only has position and direction defined
    void MakeBareTrajPoint(unsigned int fromHit, unsigned int toHit, TrajPoint& tp);
    void MakeBareTrajPoint(float fromWire, float fromTick, float toWire, float toTick, TrajPoint& tp);
    // Returns the charge weighted wire, time position of all hits in the multiplet
    // of which hit is a member
    void HitMultipletPosition(unsigned int hit, float& hitTick, float& deltaRms, float& qtot);
    bool LargeHitSep(unsigned int iht, unsigned int jht);
    float HitsTimeErr2(std::vector<unsigned int> const& hitVec);
    // Add hits to the trajectory point
    void AddTrajHits(TrajPoint& tp, bool& SignalPresent);
    // Set inTraj = -3 for all hits in the work trajectory
    void FlagWorkHits(short flag);
    bool IsLargeAngle(TrajPoint const& tp);
    // Checks to see if any hit has inTraj = -3 when it shouldn't
    bool CheckAllHitsFlag();
    void FindHit(std::string someText, unsigned int iht); // temp routine
    void IsHitFlagSet(std::string someText, unsigned int iht); // temp routine
    bool SignalAtTp(TrajPoint const& tp);
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
    // returns the DOCA between a point and a trajectory
    float PointTrajDOCA(float wire, float time, TrajPoint const& tp);
    // returns the DOCA^2 between a point and a trajectory
    float PointTrajDOCA2(float wire, float time, TrajPoint const& tp);
    void TrajTrajDOCA(Trajectory const& tp1, Trajectory const& tp2, unsigned short& ipt1, unsigned short& ipt2, float& minSep);
//    void TrajSeparation(Trajectory& iTj, Trajectory& jTj, std::vector<float>& tSep);
    void FindTrajVertices();
    void PutTrajHitsInVector(Trajectory const& tj, std::vector<unsigned int>& hitVec);
    void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, float& x, float& y);
    float TrajPointSeparation(TrajPoint& tp1, TrajPoint& tp2);
    void PrintTrajectory(Trajectory& tj ,unsigned short tPoint);
    void PrintAllTraj(unsigned short itj, unsigned short ipt);
    void PrintHeader();
    void PrintTrajPoint(unsigned short ipt, short dir, TrajPoint tp);
    // Tag as shower-like or track-like
    void TagAllTraj();
    void FindClustersOfTrajectories(std::vector<std::vector<unsigned short>>& trjintIndices);
    void DefineShowerTraj(unsigned short icot, std::vector<std::vector<unsigned short>> trjintIndices);
    // Make a track-like cluster using primTraj and a shower-like cluster consisting
    // of all other trajectories in ClsOfTrj[icot]
    void TagShowerTraj(unsigned short icot, unsigned short primTraj, unsigned short primTrajEnd, float showerAngle);
    // make a vertex from a trajectory intersection
    void MakeTrajVertex(TrjInt const& aTrjInt, unsigned short bin1, unsigned short bin2, bool& madeVtx);
    void PrintClusters();
    void FillTrajTruth();
    
  }; // class TrajClusterAlg


} // namespace cluster

#endif // ifndef TRAJCLUSTERALG_H
