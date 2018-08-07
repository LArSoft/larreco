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

#include "larreco/RecoAlg/TCAlg/Utils.h"
#include "larreco/RecoAlg/TCAlg/TCTruth.h"
#include "larreco/RecoAlg/TCAlg/TCCR.h"

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
#include "larreco/RecoAlg/LinFitAlg.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

//#include "TH1F.h"
//#include "TH2F.h"
//#include "TProfile.h"
#include "TTree.h"

namespace tca {

  class TrajClusterAlg {
    public:
    
 
    /// @{
    /// @name Data structures for the reconstruction results
    
    TrajClusterAlg(fhicl::ParameterSet const& pset);

    virtual void reconfigure(fhicl::ParameterSet const& pset);

    bool SetInputHits(std::vector<recob::Hit> const& inputHits);
//    void SetMCParticles(std::vector<simb::MCParticle> const& mcps) { fmcpList = &mcps;}
    void RunTrajClusterAlg(std::vector<unsigned int>& hitsInSlice);
    bool CreateSlice(std::vector<unsigned int>& hitsInSlice);
    

    void DefineShTree(TTree* t);
    
//    void DefineCRTree(TTree* t);

    unsigned short GetSlicesSize() { return slices.size(); }
    TCSlice const& GetSlice(unsigned short sliceIndex) const {return slices[sliceIndex]; }
    recob::Hit MergeTPHits(std::vector<unsigned int>& tpHits);
    
    std::vector<unsigned int> const& GetAlgModCount() const {return fAlgModCount; }
    std::vector<std::string> const& GetAlgBitNames() const {return AlgBitNames; }
    
    /// Deletes all the results
    void ClearResults() { slices.resize(0); }
    TruthMatcher fTM;
    
    private:
    
    bool fUseOldBackTracker {false};

    // SHOWER VARIABLE TREE
    TTree* showertree;

    // Cosmic Removal Variable Tree
    TTree* crtree;
    
    trkf::LinFitAlg fLinFitAlg;
    calo::CalorimetryAlg fCaloAlg;
    TMVA::Reader fMVAReader;

    int fWorkID;
    
    // variables for step crawling - updated at each TP
    // tracking flags
    bool fGoodTraj;         // the working trajectory is good and should be stored
    bool fTryWithNextPass;  // Try with next pass settings
    bool fMaskedLastTP;
    
    std::vector<unsigned int> fAlgModCount;

    void ReconstructAllTraj(TCSlice& slc, CTP_t inCTP);
    // Main stepping/crawling routine
    void StepCrawl(TCSlice& slc, Trajectory& tj);
    // Add hits on the trajectory point ipt that are close to the trajectory point Pos
    void AddHits(TCSlice& slc, Trajectory& tj, unsigned short ipt, bool& sigOK);
    // Large Angle version
    void AddLAHits(TCSlice& slc, Trajectory& tj, unsigned short ipt, bool& sigOK);
    // Try to use unused nearby hits in all trajectories after stepping is done
    void UseUnusedHits(TCSlice& slc);
    // Finds junk trajectories using unassigned hits
    void FindJunkTraj(TCSlice& slc, CTP_t inCTP);
    // Finds junk trajectories using unassigned hits
    bool MakeJunkTraj(TCSlice& slc, std::vector<unsigned int> tHits);
    // Step through TPs starting at the end and moving to the beginning
    void ReversePropagate(TCSlice& slc, Trajectory& tj);
    // Start a trajectory going from fromHit to (toWire, toTick)
    bool StartTraj(TCSlice& slc, Trajectory& tj, unsigned int fromhit, unsigned int tohit, unsigned short pass);
    bool StartTraj(TCSlice& slc, Trajectory& tj, float fromWire, float fromTick, float toWire, float toTick, CTP_t& tCTP, unsigned short pass);
    void GetHitMultiplet(TCSlice& slc, unsigned int theHit, std::vector<unsigned int>& hitsInMultiplet);
    void GetHitMultiplet(TCSlice& slc, unsigned int theHit, std::vector<unsigned int>& hitsInMultiplet, unsigned short& localIndex);
    // Returns fHits[iht]->RMS() * fScaleF * fHitErrFac * fHits[iht]->Multiplicity();
    float HitTimeErr(TCSlice& slc, const unsigned int iht);
    // Estimates the error^2 of the time using all hits in hitVec
    float HitsTimeErr2(TCSlice& slc, const std::vector<unsigned int>& hitVec);
    // defines HitPos, HitPosErr2 and Chg for the used hits in the trajectory point
    void ChkStopEndPts(TCSlice& slc, Trajectory& tj, bool prt);
    void DefineHitPos(TCSlice& slc, TrajPoint& tp);
    // Decide which hits to use to determine the trajectory point
    // fit, charge, etc. This is done by setting UseHit true and
    // setting inTraj < 0
     void FindUseHits(TCSlice& slc, Trajectory& tj, unsigned short ipt, float maxDelta, bool useChg);
    // Print debug output if hit iht exists in work or allTraj
    void FindHit(TCSlice& slc, std::string someText, unsigned int iht);
    // Check allTraj -> inTraj associations
    void ChkInTraj(std::string someText, TCSlice& slc);
    void FindMissedVxTjs(TCSlice& slc);
    void FindVtxTjs(TCSlice& slc, CTP_t inCTP);
    void FindVtxTraj(TCSlice& slc, VtxStore& theVtx);
    // Check the quality of the trajectory and possibly trim it
    void CheckTraj(TCSlice& slc, Trajectory& tj);
     // Truncates the trajectory if a soft kink is found in it
    void FindSoftKink(TCSlice& slc, Trajectory& tj);
    // Fill gaps in the trajectory
    void FillGaps(TCSlice& slc, Trajectory& tj);
    // Check for many unused hits and try to use them
    void CheckHiMultUnusedHits(TCSlice& slc, Trajectory& tj);
    void CheckHiMultEndHits(TCSlice& slc, Trajectory& tj);
    // Check for high values of Delta at the beginning of the trajectory
    void HiEndDelta(TCSlice& slc, Trajectory& tj);
    // Updates the last added trajectory point fit, average hit rms, etc.
    void UpdateTraj(TCSlice& slc, Trajectory& tj);
   // Estimate the Delta RMS of the TPs on the end of tj.
    void UpdateDeltaRMS(TCSlice& slc, Trajectory& tj);
    void MaskBadTPs(TCSlice& slc, Trajectory& tj, float const& maxChi);
    // The hits in the TP at the end of the trajectory were masked off. Decide whether to continue stepping with the
    // current configuration or whether to stop and possibly try with the next pass settings
    bool MaskedHitsOK(TCSlice& slc, Trajectory& tj);
    // Any re-sizing should have been done by the calling routine. This code updates the Pass and adjusts the number of
    // fitted points to get FitCHi < 2
    bool StopIfBadFits(TCSlice& slc, Trajectory& tj);
    void PrepareForNextPass(TCSlice& slc, Trajectory& tj);
    // Does a local fit of just-added TPs to identify a kink while stepping.
    // Truncates the vector and returns true if one is found.
    void GottaKink(TCSlice& slc, Trajectory& tj, unsigned short& killPts);
    // Update the parameters at the beginning of the trajectory
    void FixTrajBegin(TCSlice& slc, Trajectory& tj);
    void FixTrajBegin(TCSlice& slc, Trajectory& tj, unsigned short atPt);
    void FixTrajEnd(TCSlice& slc, Trajectory& tj, unsigned short atPt);
    bool IsGhost(TCSlice& slc, std::vector<unsigned int>& tHits);
    bool IsGhost(TCSlice& slc, Trajectory& tj);
    void CheckTrajEnd(TCSlice& slc);
    void EndMerge(TCSlice& slc, CTP_t inCTP, bool lastPass);
    void MaskTrajEndPoints(TCSlice& slc, Trajectory& tj, unsigned short nPts);
    // Sets the StopsAtEnd bits for the trajectory
    void ChkStop(TCSlice& slc, Trajectory& tj);
    // Check the Michel electron topology, lastGoodPt is the last point of muon
    bool ChkMichel(TCSlice& slc, Trajectory& tj, unsigned short& lastGoodPt);
    // TY: Split high charge hits near the trajectory end
    void ChkHiChgHits(TCSlice& slc, CTP_t inCTP);
    void SplitHiChgHits(TCSlice& slc, Trajectory& tj);

  }; // class TrajClusterAlg

} // namespace cluster

#endif // ifndef TRAJCLUSTERALG_H
