////////////////////////////////////////////////////////////////////////
//
//
// StepUtils - Core stepping algorithms
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef STEPUTILS_H
#define STEPUTILS_H

// LArSoft libraries
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"

// C/C++ standard libraries
#include <vector>

namespace tca {

  // Main stepping/crawling routine
  void StepAway(TCSlice& slc, Trajectory& tj);
  bool StopShort(TCSlice& slc, Trajectory& tj, bool prt);
  void SetStrategy(TCSlice& slc, Trajectory& tj);
  void Forecast(TCSlice& slc, const Trajectory& tj);
  // Updates the last added trajectory point fit, average hit rms, etc.
  void UpdateTraj(TCSlice& slc, Trajectory& tj);
  // Version with a different strategy for tracking high energy electrons
  void UpdateStiffEl(TCSlice& slc, Trajectory& tj);
  // Add hits on the trajectory point ipt that are close to the trajectory point Pos
  void AddHits(TCSlice& slc, Trajectory& tj, unsigned short ipt, bool& sigOK);
  // Large Angle version
  void AddLAHits(TCSlice& slc, Trajectory& tj, unsigned short ipt, bool& sigOK);
  // Step through TPs starting at the end and moving to the beginning
  void ReversePropagate(TCSlice& slc, Trajectory& tj);
  void GetHitMultiplet(const TCSlice& slc, unsigned int theHit, std::vector<unsigned int>& hitsInMultiplet,
                       bool useLongPulseHits);
  // Returns fHits[iht]->RMS() * fScaleF * fHitErrFac * fHits[iht]->Multiplicity();
  float HitTimeErr(const TCSlice& slc, const unsigned int iht);
  // Estimates the error^2 of the time using all hits in hitVec
  float HitsTimeErr2(const TCSlice& slc, const std::vector<unsigned int>& hitVec);
  void DefineHitPos(TCSlice& slc, TrajPoint& tp);
  // Decide which hits to use to determine the trajectory point
  // fit, charge, etc. This is done by setting UseHit true and
  // setting inTraj < 0
  void FindUseHits(TCSlice& slc, Trajectory& tj, unsigned short ipt, float maxDelta, bool useChg);
  // Estimate the Delta RMS of the TPs on the end of tj.
  void UpdateDeltaRMS(TCSlice& slc, Trajectory& tj);
  void MaskBadTPs(TCSlice& slc, Trajectory& tj, float const& maxChi);
  // The hits in the TP at the end of the trajectory were masked off. Decide whether to continue stepping with the
  // current configuration or whether to stop and possibly try with the next pass settings
  bool MaskedHitsOK(TCSlice& slc, Trajectory& tj);
  // Any re-sizing should have been done by the calling routine. This code updates the Pass and adjusts the number of
  // fitted points to get FitCHi < 2
  bool StopIfBadFits(TCSlice& slc, Trajectory& tj);
  // Does a local fit of just-added TPs to identify a kink while stepping.
  // Truncates the vector and returns true if one is found.
  bool GottaKink(TCSlice& slc, Trajectory& tj, bool doTrim);
  TrajPoint CreateTPFromTj(TCSlice& slc, const Trajectory& tj);
  void MaskTrajEndPoints(TCSlice& slc, Trajectory& tj, unsigned short nPts);
} // namespace tca

#endif // ifndef STEPUTILS_H
