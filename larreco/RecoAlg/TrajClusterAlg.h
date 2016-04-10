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

#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"

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

namespace tca {

  
  class TrajClusterAlg {
    public:
    
 
    /// @{
    /// @name Data structures for the reconstruction results


    TrajClusterAlg(fhicl::ParameterSet const& pset);

    virtual void reconfigure(fhicl::ParameterSet const& pset);

    void RunTrajClusterAlg(art::Event & evt);
    
    std::vector<short> const& GetinClus() const {return tjs.inClus; }
    
    /// Returns (and loses) the collection of reconstructed hits
    std::vector<art::Ptr<recob::Hit>> const& YieldHits() const { return tjs.fHits; }
    
    /// Returns a constant reference to the clusters found
    std::vector<ClusterStore> const& GetClusters() const { return tjs.tcl; }
    
    /// Returns a constant reference to the 2D end points found
    std::vector<VtxStore> const& GetEndPoints() const { return tjs.vtx; }
    
    /// Returns a constant reference to the 3D vertices found
    std::vector<Vtx3Store> const& GetVertices() const { return tjs.vtx3; }
    
    /// Deletes all the results
    void ClearResults();
    
    TjStuff tjs;
    DebugStuff Debug;
    
    private:
    
    art::InputTag fHitFinderModuleLabel; ///< label of module producing input hits
    
    short fMode;            ///  StepCrawl mode (0 = turn off)
    short fStepDir;             /// US->DS (1), DS->US (-1)
    short fNPtsAve;         /// number of points to find AveChg
    std::vector<unsigned short> fMinPtsFit; ///< Reconstruct in two passes
    std::vector<unsigned short> fMinPts;    ///< min number of Pts required to make a cluster
    std::vector<bool> fLAStep;              ///< Allow LA stepping on pass?
    float fMultHitSep;      ///< preferentially "merge" hits with < this separation
    float fMaxChi;
    float fKinkAngCut;     ///  kink angle cut
    float fChgPullCut;
    float fMaxWireSkipNoSignal;    ///< max number of wires to skip w/o a signal on them
    float fMaxWireSkipWithSignal;  ///< max number of wires to skip with a signal on them
    float fProjectionErrFactor;
    short fReversePropagate;            /// reverse the trajectory and propogate back to the beginning
    short fRecoveryAlgs;             ///< try to recover poor trajectories with different algs
    
    float fJTMaxHitSep2;  /// Max hit separation for making junk trajectories. < 0 to turn off
    
    bool fTagAllTraj;              ///< tag clusters as shower-like or track-like
    float fMaxTrajSep;     ///< max trajectory point separation for making showers
    bool fMerge;
    bool fStudyMode;       ///< study cuts
    bool fShowerStudy;    ///< study shower identification cuts
    short fShowerPrtPlane; ///< set to plane number to print out

     std::vector<float> fMaxVertexTrajSep;

    float fHitErrFac;   ///< hit time error = fHitErrFac * hit RMS used for cluster fit
    float fMinAmp;      ///< min amplitude required for declaring a wire signal is present
    float fLargeAngle;  ///< (degrees) call Large Angle Clustering code if TP angle > this value
    float fLAClusSlopeCut;
    unsigned short fAllowNoHitWire;
		float fVertex2DIPCut; 	///< 2D vtx -> cluster Impact Parameter cut (WSE)
    float fVertex3DChiCut;   ///< 2D vtx -> 3D vtx matching cut (chisq/dof)
    // TEMP variables for summing Eff*Pur
    float PiPrSum, MuSum;
    unsigned short nPiPr, nMu;

    
    bool fIsRealData;
    TH1F *fnHitsPerTP[3];
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
    bool mrgPrt;
    bool vtxPrt;
    bool didPrt; // Set true if a print condition was met
    bool shPrt; /// print shower info
    
    art::ServiceHandle<geo::Geometry> geom;
    const detinfo::LArProperties* larprop;
    const detinfo::DetectorProperties* detprop;
    
    
    trkf::LinFitAlg fLinFitAlg;

    unsigned int fFirstWire;    ///< the first wire with a hit
    unsigned int fFirstHit;     ///< first hit used
    unsigned int fLastWire;      ///< the last wire with a hit
    unsigned int fCstat;         // the current cryostat
    unsigned int fTpc;         // the current TPC
    unsigned short fPass;
    unsigned int fEvent;
    unsigned int fEventsProcessed;
    CTP_t fCTP;        ///< Cryostat/TPC/Plane code
    unsigned int fPlane;         // the current plane
    unsigned int fNumWires;   // number of wires in the current plane
    float fMaxWire;     // max wire in WSE units
    float fMaxTime;    // max time in WSE units

    // vector of pairs of first (.first) and last+1 (.second) hit on each wire
    // in the range fFirstWire to fLastWire. A value of -2 indicates that there
    // are no hits on the wire. A value of -1 indicates that the wire is dead
    std::vector< std::pair<int, int> > WireHitRange;

    std::string fhitsModuleLabel;
    
    // variables for step crawling - updated at each TP
    bool fHitDoublet; // Set true if there are a lot of hit doublets on this trajectory
    // tracking flags
    bool fGoodWork;         // the work trajectory is good and should be stored
    bool fAddedBigDeltaHit;
    bool fTryWithNextPass;  // fGoodWork false, try with next pass settings
    bool fCheckWorkModified;
    bool fUpdateTrajOK;     // update
    bool fMaskedLastTP;
    bool fRevProp;          // in reverse propagation phase
    bool fQuitAlg;          // A significant error occurred. Delete everything and return

    //  trajectories
    Trajectory work;      ///< trajectory under construction
    //  trial    trajectories
    std::vector<std::vector<Trajectory>> trial; ///< vector of all trajectories for all trials in one plane
    std::vector<std::vector<short>> inTrialTraj;

    std::vector<TrjInt> trjint;
    // vectors of clusters of trajectories => a set of somewhat intersecting trajectories
    std::vector<std::vector<unsigned short>> ClsOfTrj;
    std::vector<TjPairHitShare> tjphs;
    
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
    void MakeJunkTraj(std::vector<unsigned int> tHits, unsigned short& newTjIndex);
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
    void GetHitMultiplet(unsigned int theHit, std::vector<unsigned int>& hitsInMultiplet);
    void GetHitMultiplet(unsigned int theHit, std::vector<unsigned int>& hitsInMultiplet, unsigned short& localIndex);
    // Return true if iht and jht are both in a multiplet but have the wrong local index to start a trajectory
    bool SkipHighMultHitCombo(unsigned int iht, unsigned int jht);
    // Returns fHits[iht]->RMS() * fScaleF * fHitErrFac * fHits[iht]->Multiplicity();
    float HitTimeErr(unsigned int iht);
    // Estimates the error^2 of the time using all hits in hitVec
    float HitsTimeErr2(std::vector<unsigned int> const& hitVec);
    // estimate the number of hits expected for the provided angle
    unsigned short NumHitsExpected(float angle);
    // defines HitPos, HitPosErr2 and Chg for the used hits in the trajectory point
    void DefineHitPos(TrajPoint& tp);
    // Decide which hits to use to determine the trajectory point
    // fit, charge, etc. This is done by setting UseHit true and
    // setting inTraj < 0.
    void FindUseHits(Trajectory& tj, unsigned short ipt);
    // Try to use the hits on this TP by reducing the number of points fitted. This
    // should only be done for reasonably long TJ's
    void SetPoorUsedHits(Trajectory& tj, unsigned short ipt);
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
    // Append the allTraj trajectory to work
    void AppendToWork(unsigned short itj);
    void MakeTrajectoryObsolete(unsigned short itj);
    // Make clusters from all trajectories in allTraj
    void MakeAllTrajClusters();
    void CheckHitClusterAssociations();
    // Push the work trajectory into allTraj
    void StoreWork();
    // Check the quality of the work trajectory and possibly trim it
    void CheckWork();
    // See if another trajectory can be appended to work
    void CheckAppend();
    // Check for many unused hits in work and try to use them
    void CheckHiMultUnusedHits();
    void UseHiMultEndHits(unsigned short lastMult1Pt);
    // Check for high values of Delta at the beginning of the trajectory
    void CheckHiDeltas();
    // Check for a TJ that is close to the Large Angle cut
    void CheckNearLA();
    // Reverse a trajectory
    void ReverseTraj(Trajectory& tj);
    // Updates the last added trajectory point fit, average hit rms, etc.
    void UpdateWork();
    void UpdateTraj(Trajectory& tj);
    // Find the average charge using fNPtsAve values of TP Chg.
    void UpdateAveChg(Trajectory& tj);
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
    void EndMerge();
    void ChainMerge();
    // ****************************** Vertex code  ******************************
    void Find2DVertices();
    void AttachAnyTrajToVertex(unsigned short iv, float docaCut2, bool requireSignal);
    // make a vertex from a trajectory intersection
    void MakeTrajVertex(TrjInt const& aTrjInt, unsigned short bin1, unsigned short bin2, bool& madeVtx);
    void SplitTrajCrossingVertices();
    void FindHammerVertices();
    void Find3DVertices(geo::TPCID const& tpcid);
    void CompleteIncomplete3DVertices(geo::TPCID const& tpcid);
    short TPNearVertex(const TrajPoint& tp);
    // ****************************** Shower/Track ID code  ******************************
    // Tag as shower-like or track-like
    void TagAllTraj();
    void TagPhotons();
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
    void KillVerticesInShowers();

    
  }; // class TrajClusterAlg


} // namespace cluster

#endif // ifndef TRAJCLUSTERALG_H
