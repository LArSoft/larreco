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
#include <vector>
#include <utility> // std::pair<>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

// framework libraries
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Utilities/Exception.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
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
    
    /// Returns (and loses) the art::Ptr collection of previously reconstructed hits (e.g. gaushit)
    std::vector<recob::Hit> YieldHits();
    art::InputTag const& GetHitFinderModuleLabel() { return fHitFinderModuleLabel; }
    
    /// Returns a constant reference to the clusters found
    std::vector<ClusterStore> const& GetClusters() const { return tjs.tcl; }
    
    /// Returns a constant reference to the 2D end points found
    std::vector<VtxStore> const& GetEndPoints() const { return tjs.vtx; }
    
    /// Returns a constant reference to the 3D vertices found
    std::vector<Vtx3Store> const& GetVertices() const { return tjs.vtx3; }
    
    /// Returns a constant reference to the 3D vertices found
    std::vector<std::vector<unsigned short>> const& Get3DMatchedClusters() const { return tjs.MatchedClusters; }
    
    std::vector<unsigned int> const& GetAlgModCount() const {return fAlgModCount; }
    std::vector<std::string> const& GetAlgBitNames() const {return AlgBitNames; }
    
    static bool SortByMultiplet(TCHit const& a, TCHit const& b);
    
    /// Deletes all the results
    void ClearResults();
    
    TjStuff tjs;
    DebugStuff Debug;
    
    private:
    
    art::InputTag fHitFinderModuleLabel; ///< label of module producing input hits
    art::InputTag fCalDataModuleLabel; ///< label of module producing ROIs on wires
    
    short fMode;            ///  StepCrawl mode (0 = turn off)
    short fStepDir;             /// US->DS (1), DS->US (-1)
    short fNPtsAve;         /// number of points to find AveChg
    std::vector<unsigned short> fMinPtsFit; ///< Reconstruct in two passes
    std::vector<unsigned short> fMinPts;    ///< min number of Pts required to make a trajectory
    std::vector<unsigned short> fMaxAngleCode;   ///< max allowed angle code for each pass
    float fMultHitSep;      ///< preferentially "merge" hits with < this separation
    float fMaxChi;
    std::vector<float> fKinkCuts; ///< kink angle, nPts fit, (alternate) kink angle significance
    std::vector<float> fQualityCuts; ///< Min points/wire, min consecutive pts after a gap
    std::vector<float> fChargeCuts;
    float fMaxWireSkipNoSignal;    ///< max number of wires to skip w/o a signal on them
    float fMaxWireSkipWithSignal;  ///< max number of wires to skip with a signal on them
    float fProjectionErrFactor;
    bool fMakeNewHits;
   
    float fJTMaxHitSep2;  /// Max hit separation for making junk trajectories. < 0 to turn off
    
    bool fTagAllTraj;              ///< tag clusters as shower-like or track-like
    float fMaxTrajSep;     ///< max trajectory point separation for making showers
    bool fStudyMode;       ///< study cuts
    std::vector<float> fMatchTruth;     ///< Match to MC truth
 
    std::vector<float> fMaxVertexTrajSep;
    std::bitset<32> fUseAlg;  ///< Allow user to mask off specific algorithms

    float fHitErrFac;   ///< hit time error = fHitErrFac * hit RMS used for cluster fit
    float fMinAmp;      ///< min amplitude required for declaring a wire signal is present
    std::vector<float> fAngleRanges; ///< list of max angles for each angle range
    std::vector<float> fAngleRangesMaxHitsRMS;
    float fVLAStepSize;
    
    float fLAClusSlopeCut;
    unsigned short fAllowNoHitWire;
		float VertexPullCut; 	///< maximum 2D vtx - trajectory significance
    std::vector<short> fDeltaRayTag; ///< min length, min MCSMom and min separation (WSE) for a delta ray tag
    std::vector<short> fMuonTag; ///< min length and min MCSMom for a muon tag
    std::vector<short> fShowerTag; ///< [min MCSMom, max separation, min # Tj < separation] for a shower tag
    std::vector<float> fChkStopCuts; ///< [Min Chg ratio, Chg slope pull cut, Chg fit chi cut]

    std::vector<float> fVertex2DCuts; ///< Max position pull, max Position error rms
    float fVertex3DChiCut;   ///< 2D vtx -> 3D vtx matching cut (chisq/dof)
    std::vector<float> fMatch3DCuts;  ///< Max dX separation
    
    // Variables for summing Eff*Pur for electrons, muons, pions, kaons and protons
    std::array<short, 5> EPCnts;
    std::array<float, 5> EPSums;
    std::array<float, 5> EPTSums;
    bool fIsRealData;
/*
    TH2F *fMCSMom_KE_e;
    TH2F *fMCSMom_KE_mu;
    TH2F *fMCSMom_KE_pi;
    TH2F *fMCSMom_KE_p;
    TH1F *fnHitsPerTP[3];
    TH1F *fDelta[3];
    TH1F *fDeltaN[3];
    TH1F *fCharge[3];
    TH1F *fPrEP;
    TH1F *fMuPiEP;
    TH2F *fnHitsPerTP_Angle[3];
    TProfile *fnHitsPerTP_AngleP[3];
*/
    // Reco-MC vertex position difference
    TH1F* fNuVtx_dx;
    TH1F* fNuVtx_dy;
    TH1F* fNuVtx_dz;
    
    // Reco-MC stopping wire difference for different MC Particles
    TH1F* fdWire[5];
    // EP vs KE for different MC Particles
    TProfile* fEP_T[5];

    TH1F *fDeltaN[3];
    TH1F *fHitRMS[3];
    TH2F *fTPWidth_Angle[3];
    TProfile *fTPWidth_AngleP[3];
    TProfile *fExpect_Angle[3];
    
    // number of primary particles in the event
    unsigned short nTruPrimary;
    // number of reconstructable primary particles in the event
    unsigned short nTruPrimaryOK;
    // number of reconstructable neutrino vertices in ALL events
    unsigned short nTruPrimaryVtxOK;
    // number of reconstructable neutrino vertices in ALL events that were reconstructed
    unsigned short nTruPrimaryVtxReco;
 
    bool prt;
    bool mrgPrt;
    bool vtxPrt;
    bool didPrt;
    short TJPrt; // Set to the WorkID of a trajectory that is being debugged
    bool shPrt; /// print shower info
    
    art::ServiceHandle<geo::Geometry> geom;
    const detinfo::LArProperties* larprop;
    const detinfo::DetectorProperties* detprop;
    // TEMP for writing event filter selection
//    std::ofstream outFile;
    
    trkf::LinFitAlg fLinFitAlg;

    unsigned int fCstat;         // the current cryostat
    unsigned int fTpc;         // the current TPC
    unsigned int fRun, fSubRun;
    unsigned int fEvent;
    unsigned int fEventsProcessed;
    CTP_t fCTP;        ///< Cryostat/TPC/Plane code
    unsigned int fPlane;         // the current plane
    int fWorkID;


    std::string fhitsModuleLabel;
    
    // variables for step crawling - updated at each TP
    // tracking flags
    bool fGoodTraj;         // the working trajectory is good and should be stored
    bool fTryWithNextPass;  // Try with next pass settings
    bool fUpdateTrajOK;     // update
    bool fMaskedLastTP;
    bool fQuitAlg;          // A significant error occurred. Delete everything and return
    
    std::vector<float> fAveHitRMS;      ///< average RMS of an isolated hit
    
    std::vector<unsigned int> fAlgModCount;
    
//    void WatchHit(std::string someText, int wid);
//    short watchInTraj;
    // runs the TrajCluster algorithm on one plane specified by the calling routine
    void RunStepCrawl();
    void InitializeAllTraj();
    void ReconstructAllTraj();
    // Main stepping/crawling routine
    void StepCrawl(Trajectory& tj);
    // Add hits on the trajectory point ipt that are close to the trajectory point Pos
    void AddHits(Trajectory& tj, unsigned short ipt, bool& sigOK);
    // Large Angle version
    void AddLAHits(Trajectory& tj, unsigned short ipt, bool& sigOK);
    // Hits on two adjacent wires have an acceptable signal overlap
    bool TrajHitsOK(const std::vector<unsigned int>& iHitsInMultiplet, const std::vector<unsigned int>& jHitsInMultiplet);
    // Try to use unused nearby hits in all trajectories after stepping is done
    void UseUnusedHits();
    // Finds junk trajectories using unassigned hits
    void FindJunkTraj();
    // Finds junk trajectories using unassigned hits
    void MakeJunkTraj(std::vector<unsigned int> tHits, unsigned short& newTjIndex);
    // Step through TPs starting at the end and moving to the beginning
    void ReversePropagate(Trajectory& tj);
    // Start a trajectory going from fromHit to (toWire, toTick)
    bool StartTraj(Trajectory& tj, const unsigned int& fromHit, const unsigned int& toHit, const unsigned short& pass);
    bool StartTraj(Trajectory& tj, const float& fromWire, const float& fromTick, const float& toWire, const float& toTick, const CTP_t& tCTP, const unsigned short& pass);
    void GetHitMultiplet(unsigned int theHit, std::vector<unsigned int>& hitsInMultiplet);
    void GetHitMultiplet(unsigned int theHit, std::vector<unsigned int>& hitsInMultiplet, unsigned short& localIndex);
    // Returns fHits[iht]->RMS() * fScaleF * fHitErrFac * fHits[iht]->Multiplicity();
    float HitTimeErr(const unsigned int iht);
    // Estimates the error^2 of the time using all hits in hitVec
    float HitsTimeErr2(const std::vector<unsigned int>& hitVec);
    // defines HitPos, HitPosErr2 and Chg for the used hits in the trajectory point
    void DefineHitPos(TrajPoint& tp);
    // Decide which hits to use to determine the trajectory point
    // fit, charge, etc. This is done by setting UseHit true and
    // setting inTraj < 0
     void FindUseHits(Trajectory& tj, unsigned short ipt, float maxDelta, bool useChg);
    // Sets the tp AngCode = 0 (small angle), 1 (large angle) or 2 (very large angle)
    void SetAngleCode(TrajPoint& tp);
    // returns the index of the angle range that tp is in
    unsigned short AngleRange(TrajPoint const& tp);
    unsigned short AngleRange(float angle);
    // Print debug output if hit iht exists in work or allTraj
    void FindHit(std::string someText, unsigned int iht);
    // Check allTraj -> inTraj associations
    void ChkInTraj(std::string someText);
    // Counts the number of hits that are used in two different vectors of hits
    void CountSameHits(std::vector<unsigned int>& iHitVec, std::vector<unsigned int>& jHitVec, unsigned short& nSameHits);
    // Merge and store the two trajectories in allTraj
    bool MergeAndStore(unsigned short tj1,  unsigned short tj2);
    // Make clusters from all trajectories in allTraj
    void MakeAllTrajClusters();
    // Push the trajectory into allTraj
    void StoreTraj(Trajectory& tj);
    // Check the quality of the trajectory and possibly trim it
    void CheckTraj(Trajectory& tj);
     // Truncates the trajectory if a soft kink is found in it
    void FindSoftKink(Trajectory& tj);
    // Fill gaps in the trajectory
    void FillGaps(Trajectory& tj);
    // Check for many unused hits and try to use them
    void CheckHiMultUnusedHits(Trajectory& tj);
    void CheckHiMultEndHits(Trajectory& tj);
    // Check for high values of Delta at the beginning of the trajectory
    void CheckHiDeltas(Trajectory& tj);
    // Check for a TJ that is close to the Large Angle cut
    void CheckNearLA();
    // Updates the last added trajectory point fit, average hit rms, etc.
    void UpdateTraj(Trajectory& tj);
    // Find the average charge using fNPtsAve values of TP Chg.
    void UpdateAveChg(Trajectory& tj);
   // Estimate the Delta RMS of the TPs on the end of tj.
    void UpdateDeltaRMS(Trajectory& tj);
    void MaskBadTPs(Trajectory& tj, float const& maxChi);
    // The hits in the TP at the end of the trajectory were masked off. Decide whether to continue stepping with the
    // current configuration or whether to stop and possibly try with the next pass settings
    bool MaskedHitsOK(Trajectory& tj);
    // Any re-sizing should have been done by the calling routine. This code updates the Pass and adjusts the number of
    // fitted points to get FitCHi < 2
    void PrepareForNextPass(Trajectory& tj);
    // Fit the supplied trajectory using HitPos positions with the origin at originPt.
    void FitTraj(Trajectory& tj, unsigned short originPt, unsigned short npts, short fitDir, TrajPoint& tpFit);
    // Fit the leading edge of the trajectory
    void FitTraj(Trajectory& tj);
    // Does a local fit of just-added TPs to identify a kink while stepping.
    // Truncates the vector and returns true if one is found.
    void GottaKink(Trajectory& tj, unsigned short& killPts);
    // Update the parameters at the beginning of the trajectory
    void FixTrajBegin(Trajectory& tj);
    void FixTrajBegin(Trajectory& tj, unsigned short atPt);
    void FixTrajEnd(Trajectory& tj, unsigned short atPt);
    bool IsGhost(std::vector<unsigned int>& tHits, unsigned short& ofTraj);
    bool IsGhost(Trajectory& tj);
    void CheckTrajEnd();
    void EndMerge();
    void EndMerge2();
    void FillWireHitRange(const geo::TPCID& tpcid);
    float ExpectedHitsRMS(TrajPoint const& tp);
    /// sets fQuitAlg true if WireHitRange has a problem
    bool CheckWireHitRange();
    // Erases delHit and makes corrections to inTraj, allTraj and WireHitRange
    bool EraseHit(const unsigned int& delHit);
    // Creates a hit in tjs.fHits using the supplied information. Returns UINT_MAX if there is failure.
    // Returns the index of the newly created hit
    void DefineHit(TCHit& tcHit, CTP_t& hitCTP, unsigned int& hitWire);
    unsigned int CreateHit(TCHit tcHit);
    // Merges all of the hits used in each TP into one hit
    void MergeTPHits();
    void MaskTrajEndPoints(Trajectory& tj, unsigned short nPts);
    // Sets the StopsAtEnd bits for all trajectories in the current CTP
    void ChkAllStop();
    // Sets the StopsAtEnd bits for the trajectory
    void ChkStop(Trajectory& tj);
    // Check the Michel electron topology, lastGoodPt is the last point of muon
    bool ChkMichel(Trajectory& tj, unsigned short& lastGoodPt);
    void MatchTruth();
    // ****************************** Vertex code  ******************************
    void Find2DVertices();
    void FindVtxTraj(unsigned short ivx);
//    void Refine2DVertices();
    void SplitTrajCrossingVertices();
    void FindHammerVertices();
    void FindHammerVertices2();
    void Find3DVertices(const geo::TPCID& tpcid);
    void CompleteIncomplete3DVertices(const geo::TPCID& tpcid);
    void CompleteIncomplete3DVerticesInGaps(const geo::TPCID& tpcid);
    // Improve hit assignments near vertex 
    void VtxHitsSwap();
    // ****************************** 3D Tj matching code  ******************************
    void Match3D(const geo::TPCID& tpcid);
    void Match3D2Views(const geo::TPCID& tpcid, const std::vector<float>& xx);
    
  }; // class TrajClusterAlg


} // namespace cluster

#endif // ifndef TRAJCLUSTERALG_H
