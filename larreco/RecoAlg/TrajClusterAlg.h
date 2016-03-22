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
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Utilities/Exception.h"

// LArSoft libraries
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/RecoAlg/LinFitAlg.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

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
    static CTP_t EncodeCTP(const geo::WireID& wireID)
    { return EncodeCTP(wireID.Cryostat, wireID.TPC, wireID.Plane); }
    
    static geo::PlaneID DecodeCTP(CTP_t CTP)
      { return { CTP / (CTPpad*CTPpad), CTP / CTPpad % CTPpad, CTP % CTPpad }; }

    /// @{
    /// @name Data structures for the reconstruction results

    /// struct of temporary clusters
    struct ClusterStore {
      short ID;         // Cluster ID. ID < 0 = abandoned cluster
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
      float Wire {0};
      float WireErr {2};
      float Time {0};
      float TimeErr {0.1};
      unsigned short NTraj {0};  // = 0 for abandoned vertices
      float ChiDOF {0};
      short Topo {0}; 			// 1 = US-US, 2 = US-DS, 3 = DS-US, 4 = DS-DS, 5 = Star, 6 = hammer
      CTP_t CTP {0};
      bool Fixed {false};                 // Vertex position fixed (should not be re-fit)
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
    std::vector<unsigned short> fMinPtsFit; ///< Reconstruct in two passes
    std::vector<unsigned short> fMinPts;    ///< min number of Pts required to make a cluster
    std::vector<bool> fLAStep;              ///< Allow LA stepping on pass?
    float fMultHitSep;      ///< preferentially "merge" hits with < this separation
    float fMaxChi;
    float fHitFOMCut;      ///  Combined charge & delta difference Figure of Merit cut
    float fChgRatCut;      ///  Relative charge to delta difference weighting
    float fKinkAngCut;     ///  kink angle cut
    float fMaxWireSkipNoSignal;    ///< max number of wires to skip w/o a signal on them
    float fMaxWireSkipWithSignal;  ///< max number of wires to skip with a signal on them
    float fProjectionErrFactor;
    short fReversePropagate;            /// reverse the trajectory and propogate back to the beginning
    short fRecoveryAlgs;             ///< try to recover poor trajectories with different algs
    
    float fJTMaxHitSep2;  /// Max hit separation for making junk trajectories. < 0 to turn off
    
    bool fTagAllTraj;              ///< tag clusters as shower-like or track-like
    float fMaxTrajSep;     ///< max trajectory point separation for making showers
    bool fStudyMode;       ///< study cuts
    bool fShowerStudy;    ///< study shower identification cuts
    short fShowerPrtPlane; ///< set to plane number to print out

    bool fFind2DVertices;
    std::vector<float> fMaxVertexTrajSep;

    float fHitErrFac;   ///< hit time error = fHitErrFac * hit RMS used for cluster fit
    float fMinAmp;      ///< min amplitude required for declaring a wire signal is present
    float fLargeAngle;  ///< (degrees) call Large Angle Clustering code if TP angle > this value
    float fLAClusSlopeCut;
    float fMergeOverlapAngCut;   ///< angle cut for merging overlapping clusters
    unsigned short fAllowNoHitWire;
		float fVertex2DCut; 	///< 2D vtx -> cluster matching cut (chisq/dof)
    float fVertex2DWireErrCut;
    float fVertex3DChiCut;   ///< 2D vtx -> 3D vtx matching cut (chisq/dof)

    int fDebugPlane;
    int fDebugWire;  ///< set to the Begin Wire and Hit of a cluster to print
    int fDebugHit;   ///< out detailed information while crawling
    
    bool fIsRealData;
    TH1F *fnHitsPerTP[3];
    TH1F *fnHitsFitPerTP[3];
    TH1F *fDelta[3];
    TH1F *fDeltaN[3];
    TH1F *fCharge[3];
    TH2F *fnHitsPerTP_Angle[3];
    TProfile *fnHitsPerTP_AngleP[3];
    
    TH1F *fShowerNumTrjint;
    TH2F *fShowerTheta_Sep;
    TH1F *fShowerDVtx;
    TH2F *fShowerDVtx_Sep;
    
    bool prt;
    bool vtxPrt;
    bool didPrt; // Set true if a print condition was met
    bool shPrt; /// print shower info
    
    art::ServiceHandle<geo::Geometry> geom;
    const detinfo::LArProperties* larprop;
    const detinfo::DetectorProperties* detprop;
    
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
    float fMaxWire;     // max wire in WSE units
    float fMaxTime;    // max time in WSE units
    float fScaleF;     ///< scale factor from Tick/Wire to dx/du

    // vector of pairs of first (.first) and last+1 (.second) hit on each wire
    // in the range fFirstWire to fLastWire. A value of -2 indicates that there
    // are no hits on the wire. A value of -1 indicates that the wire is dead
    std::vector< std::pair<int, int> > WireHitRange;

    std::string fhitsModuleLabel;
    
    // variables for step crawling - updated at each TP
    float fAveChg;
    bool fHitDoublet; // Set true if there are a lot of hit doublets on this trajectory
    // tracking flags
    bool fQuitAlg;          // A significant error occurred. Delete everything and return
    bool fGoodWork;         // the work trajectory is good and should be stored
    bool fAddedBigDeltaHit;
    bool fTryWithNextPass;  // fGoodWork false, try with next pass settings
    bool fCheckWorkModified;
    bool fUpdateTrajOK;     // update
    bool fMaskedLastTP;
    bool fRevProp;          // in reverse propagation phase

    bool fSplitTrajOK;
    
    struct TrajPoint {
      CTP_t CTP {0};                   ///< Cryostat, TPC, Plane code
      std::array<float, 2> HitPos {{0,0}}; // Charge weighted position of hits in wire equivalent units
      std::array<float, 2> Pos {{0,0}}; // Trajectory position in wire equivalent units
      std::array<float, 2> Dir {{0,0}}; // Direction
      float HitPosErr2 {0};         // Uncertainty^2 of the hit position perpendiclar to the direction
                                // HitPosErr2 < 0 = HitPos not defined because no hits used
      float Ang {0};                // Trajectory angle (-pi, +pi)
      float AngErr {0.1};             // Trajectory angle error
      float KinkAng {-1};            // Just what it says
      float Chg {0};                // Charge
      float AveChg {0};             // Average charge of last ~20 TPs
      float ChgRat {0.1};            //  = (Chg - fAveChg) / fChgRMS
      float Delta {0};              // Deviation between trajectory and hits (WSE)
      float DeltaRMS {0};           // RMS of Deviation between trajectory and hits (WSE)
      unsigned short NTPsFit {2}; // Number of trajectory points fitted to make this point
      unsigned short Step {0};      // Step number at which this TP was created
      float FitChi {0};             // Chi/DOF of the fit
      std::vector<unsigned int> Hits; // vector of fHits indices
      std::vector<bool> UseHit; // set true if the hit is used in the fit
    };
    
    // Global information for the trajectory
    struct Trajectory {
      short ID;
      CTP_t CTP {0};                      ///< Cryostat, TPC, Plane code
      unsigned short Pass {0};            ///< the pass on which it was created
      short StepDir {0};                 /// -1 = going US (CC proper order), 1 = going DS
      unsigned short ClusterIndex {USHRT_MAX};   ///< Index not the ID...
      unsigned short ProcCode {1};       ///< USHRT_MAX = abandoned trajectory
      std::bitset<32> AlgMod;        ///< Bit set if algorithm AlgBit_t modifed the trajectory
      // ProcCode = 1 (normal step), 2 (junk traj), 3 (split trajectory)
      unsigned short PDG {13};            ///< shower-like or track-like {default is track-like}
      unsigned short ParentTraj {USHRT_MAX};     ///< index of the parent (if PDG = 12)
      int TruPDG {0};                    ///< MC truth
      int TruKE {0};                     ///< MeV
      bool IsPrimary {false};                ///< MC truth
      std::array<short, 2> Vtx {{-1,-1}};      ///< Index of 2D vertex
      std::array<unsigned short, 2> EndPt {{0,0}}; ///< First and last point in the trajectory that has a hit
      std::vector<TrajPoint> Pts;    ///< Trajectory points
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
    
    // Algorithm modification bits
    typedef enum {
      kMaskedWorkHits,
      kGottaKink,     ///< GottaKink found a kink
      kCWKink,        ///< kink found in CheckWork
      kCWStepChk,
      kMaybeDeltaRay,
      kModifyShortTraj,
      kTryWithNextPass,
      kRevProp,
      kRecovery1,
      kRecovery2,
      kManyHitsAdded,
      kSplitTraj,
      kComp3DVx,
      kHiDelta,
      kHammer2DVx,
      kJunkTj,
      kAlgBitSize
    } AlgBit_t;
    
    std::vector<std::string> AlgBitNames {
      "MaskedHits",
      "Kink",
      "CWKink",
      "CWStepChk",
      "DeltaRay?",
      "ModifyShortTj",
      "TryNextPass",
      "RevProp",
      "Recovery1",
      "Recovery2",
      "ManyHitsAdded",
      "SplitTraj",
      "Comp3DVx",
      "HiDelta",
      "Hammer2DVx",
      "JunkTj"
    };
    
    // runs the TrajCluster algorithm on one plane specified by the calling routine
    // (which should also have called GetHitRange)
    void RunStepCrawl();
    void ReconstructAllTraj();
    // Check the work trajectory and try to recover it if it has poor quality
    // Main stepping/crawling routine
    void StepCrawl();
    // Add hits on the trajectory point ipt that are close to the trajectory point Pos
    void AddHits(Trajectory& tj, unsigned short ipt, bool& sigOK);
    void GetHitRange();
    float DeadWireCount(float inWirePos1, float inWirePos2);
    void HitSanityCheck();
    // Hits on two adjacent wires have an acceptable signal overlap
    bool TrajHitsOK(unsigned int iht, unsigned int jht);
    // Find all trajectories in the plane using the pre-defined step direction, cuts, etc.
    void TryRecoveryAlgs();
    // Finds junk trajectories using unassigned hits
    void FindJunkTraj();
    // Finds junk trajectories using unassigned hits
    void MakeJunkTraj(std::vector<unsigned int> tHits);
    // Step through TPs starting at the end and moving to the beginning
    void ReversePropagate(Trajectory& tj);
    // Calculate the number of missed steps that should be expected for the given
    // trajectory point direction
    unsigned short SetMissedStepCut(TrajPoint const& tp);
    // Start a trajectory going from fromHit to (toWire, toTick)
    void StartWork(float fromWire, float fromTick, float toWire, float toTick, CTP_t tCTP);
    void StartWork(unsigned int fromHit, unsigned int toHit);
    // Make a bare trajectory point that only has position and direction defined
    void MakeBareTrajPoint(unsigned int fromHit, unsigned int toHit, TrajPoint& tp);
    void MakeBareTrajPoint(float fromWire, float fromTick, float toWire, float toTick, CTP_t tCTP, TrajPoint& tp);
    // Returns the charge weighted wire, time position of all hits in the multiplet
    // of which hit is a member
    void HitMultipletPosition(unsigned int hit, float& hitTick, float& deltaRms, float& qtot);
    // Return true if iht and jht are both in a multiplet but have the wrong local index to start a trajectory
    bool SkipHighMultHitCombo(unsigned int iht, unsigned int jht);
    // Returns fHits[iht]->RMS() * fScaleF * fHitErrFac * fHits[iht]->Multiplicity();
    float HitTimeErr(unsigned int iht);
    // Estimates the error^2 of the time using all hits in hitVec
    float HitsTimeErr2(std::vector<unsigned int> const& hitVec);
    // defines HitPos, HitPosErr2 and Chg for the used hits in the trajectory point
    void DefineHitPos(TrajPoint& tp);
    // Decide which hits to use to determine the trajectory point
    // fit, charge, etc. This is done by setting UseHit true and
    // setting inTraj < 0.
    void FindUseHits(Trajectory& tj, unsigned short ipt);
    // Try to use the hits on this TP by reducing the number of points fitted. This
    // should only be done for reasonably long TJ's
    void SetPoorUsedHits(Trajectory& tj, unsigned short ipt);
    // Returns true if the charge of the hit pointed to by ipt and iht is consistent
    // with the average charge fAveChg
    bool HitChargeOK(Trajectory& tj, unsigned short ipt, unsigned short iht);
    // Sets inTraj[] = 0 and UseHit false for all used hits in tp
    void UnsetUsedHits(TrajPoint& tp);
    void SetAllHitsUsed(TrajPoint& tp);
    // Returns  true if there is a signal on the line between (wire1, time1) and (wire2, time2).
    bool SignalPresent(float wire1, float time1, TrajPoint const& tp);
    bool SignalPresent(unsigned int wire1, float time1, unsigned int wire2, float time2);
    bool SignalPresent(float wire1, float time1, float wire2, float time2);
    bool SignalPresent(TrajPoint const& tp);
    // Counts the number of used hits in tp
    unsigned short NumUsedHits(TrajPoint& tp);
    // Counts the number of TPs in the trajectory that have charge
    unsigned short NumPtsWithCharge(Trajectory& tj);
    // Sets inTraj[] = 0 and UseHit false for all TPs in work. Called when abandoning work
    void ReleaseWorkHits();
    // Returns true if the TP angle exceeds user cut fLargeAngle
    bool IsLargeAngle(TrajPoint const& tp);
    // Checks to see if any hit has inTraj = -1 when it shouldn't
    bool CheckAllHitsFlag();
    // Print debug output if hit iht exists in work or allTraj
    void FindHit(std::string someText, unsigned int iht);
    // Check allTraj -> inTraj associations
    void CheckInTraj(std::string someText);
    // Returns true if there a is wire signal at tp
    bool SignalAtTp(TrajPoint const& tp);
    // analyze the sat vector to construct a vector of trajectories that is the best
    void AnalyzeTrials();
    // Counts the number of hits that are used in two different vectors of hits
    void CountSameHits(std::vector<unsigned int>& iHitVec, std::vector<unsigned int>& jHitVec, unsigned short& nSameHits);
    void AdjudicateTrials(bool& reAnalyze);
    // merge the trajectories and put the results into allTraj. Returns
    // reAnalyze true if AnalyzeTrials needs to be called again
    void MergeTrajPair(unsigned short ipr, bool& reAnalyze);
    // Split the allTraj trajectory itj at position pos into two trajectories
    // with an optional vertex assignment
    void SplitAllTraj(unsigned short itj, unsigned short pos, unsigned short ivx);
    // Make clusters from all trajectories in allTraj
    void MakeAllTrajClusters();
    void CheckHitClusterAssociations();
    // Push the work trajectory into allTraj
    void StoreWork();
    // Check the quality of the work trajectory and possibly trim it
    void CheckWork();
    // Check for many unused hits in work and try to use them
    void CheckHiMultUnusedHits();
    // Check for high values of Delta at the beginning of the trajectory
    void CheckHiDeltas();
    // Reverse a trajectory
    void ReverseTraj(Trajectory& tj);
    // Find the first (last) TPs, EndPt[0] (EndPt[1], that have charge
    void SetEndPoints(Trajectory& tj);
    // Updates the last added trajectory point fit, average hit rms, etc.
    void UpdateWork();
    void UpdateTraj(Trajectory& tj);
    // Project TP to a "wire position" Pos[0] and update Pos[1]
    void MoveTPToWire(TrajPoint& tp, float wire);
    // Find the average charge using fNPtsAve values of TP Chg. Result put in fAveChg and Pts[updatePt].Chg
    void UpdateAveChg(Trajectory& tj, unsigned short updatePt, short dir);
   // Estimate the Delta RMS of the TPs on the end of work.
    void UpdateDeltaRMS(Trajectory& tj);
    // The hits in the TP at the end of the trajectory were masked off. Decide whether to continue stepping with the
    // current configuration or whether to stop and possibly try with the next pass settings
    bool MaskedWorkHitsOK();
    // Any re-sizing should have been done by the calling routine. This code updates the Pass and adjusts the number of
    // fitted points to get FitCHi < 2
    void PrepareWorkForNextPass();
    // Fit the supplied trajectory using HitPos positions with the origin at originPt.
    void FitTraj(Trajectory& tj, unsigned short originPt, unsigned short npts, short fitDir, TrajPoint& tpFit);
    // A version just like FitWork (that will be replaced with this one)
    void FitTraj(Trajectory& tj);
   // Jacket around FitTraj to fit the work trajectory
    void FitWork();
    float TrajLength(Trajectory& tj);
    // Does a local fit of just-added TPs on work to identify a kink while stepping.
    // Truncates the work vector and returns true if one is found.
    void GottaKink(Trajectory& tj, unsigned short& killPts);
    // Combines hit multiplets in short trajectories
    void ModifyShortTraj(Trajectory& tj);
    // See if the trajectory appears to be a delta ray. This is characterized by a significant fraction of hits
    // in the trajectory belonging to an existing trajectory. This may also flag ghost trajectories...
    bool MaybeDeltaRay(Trajectory& tj);
    void CheckTrajEnd();
    void TrajClosestApproach(Trajectory const& tj, float x, float y, unsigned short& iClosePt, float& Distance);
    // returns the DOCA between a hit and a trajectory
    float PointTrajDOCA(unsigned int iht, TrajPoint const& tp);
    // returns the DOCA between a (W,T) point and a trajectory
    float PointTrajDOCA(float wire, float time, TrajPoint const& tp);
    // returns the DOCA^2 between a point and a trajectory
    float PointTrajDOCA2(float wire, float time, TrajPoint const& tp);
    // returns the separation^2 between two TPs
    float TrajPointHitSep2(TrajPoint const& tp1, TrajPoint const& tp2);
    // returns the separation^2 between a point and a TP
    float PointTrajSep2(float wire, float time, TrajPoint const& tp);
    // finds the point on trajectory tj that is closest to trajpoint tp
    void TrajPointTrajDOCA(TrajPoint const& tp, Trajectory const& tj, unsigned short& closePt, float& minSep);
   // returns the intersection position, intPos, of two trajectory points
    void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, float& x, float& y);
    // Returns the separation distance between two trajectory points
    float TrajPointSeparation(TrajPoint& tp1, TrajPoint& tp2);
    // returns the separation^2 between two hits in WSE units
    float HitSep2(unsigned int iht, unsigned int jht);
    // Returns true if the time separation of two hits exceeds fMultHitSep
    bool LargeHitSep(unsigned int iht, unsigned int jht);
    // Find the Distance Of Closest Approach between two trajectories, exceeding minSep
    void TrajTrajDOCA(Trajectory const& tp1, Trajectory const& tp2, unsigned short& ipt1, unsigned short& ipt2, float& minSep);
    // Calculates the angle between two TPs
    float TwoTPAngle(TrajPoint& tp1, TrajPoint& tp2);
    // Put hits in each trajectory point into a flat vector. Only hits with UseHit if onlyUsedHits == true
    void PutTrajHitsInVector(Trajectory const& tj, bool onlyUsedHits, std::vector<unsigned int>& hitVec);
    // Returns a vector with the size of the number of trajectory points on iTj
    // with the separation distance between the closest traj points
//    void TrajSeparation(Trajectory& iTj, Trajectory& jTj, std::vector<float>& tSep);
    // ****************************** Vertex code  ******************************
    void Find2DVertices();
    void AttachAnyTrajToVertex(unsigned short iv, float docaCut2, bool requireSignal);
    // make a vertex from a trajectory intersection
    void MakeTrajVertex(TrjInt const& aTrjInt, unsigned short bin1, unsigned short bin2, bool& madeVtx);
    void SplitTrajCrossingVertices();
    void FindHammerVertices();
    void Find3DVertices(geo::TPCID const& tpcid);
    void CompleteIncomplete3DVertices(geo::TPCID const& tpcid);
    // ****************************** Printing code  ******************************
    // Print trajectories, TPs, etc to mf::LogVerbatim
    void PrintTrajectory(Trajectory const& tj ,unsigned short tPoint);
    void PrintAllTraj(unsigned short itj, unsigned short ipt);
    void PrintHeader();
    void PrintTrajPoint(unsigned short ipt, short dir, unsigned short pass, TrajPoint const& tp);
    // Print clusters after calling MakeAllTrajClusters
    void PrintClusters();
    // Print a single hit in the standard format
    std::string PrintHit(unsigned int iht);
    // Print Trajectory position in the standard format
    std::string PrintPos(TrajPoint const& tp);
    // ****************************** Shower/Track ID code  ******************************
    // Tag as shower-like or track-like
    void TagAllTraj();
    void KillAllTrajInCTP(CTP_t tCTP);
    // Associate trajectories that are close to each into Clusters Of Trajectories (COTs)
    void FindClustersOfTrajectories(std::vector<std::vector<unsigned short>>& trjintIndices);
    // Try to define a shower trajectory consisting of a single track-like trajectory (the electron/photon)
    // plus a group of shower-like trajectories, using the vector of trjintIndices
    void DefineShowerTraj(unsigned short icot, std::vector<std::vector<unsigned short>> trjintIndices);
    // Make a track-like cluster using primTraj and a shower-like cluster consisting
    // of all other trajectories in ClsOfTrj[icot]
    void TagShowerTraj(unsigned short icot, unsigned short primTraj, unsigned short primTrajEnd, float showerAngle);
    void FillTrajTruth();
    
  }; // class TrajClusterAlg


} // namespace cluster

#endif // ifndef TRAJCLUSTERALG_H
