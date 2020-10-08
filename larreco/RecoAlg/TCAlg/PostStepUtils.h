////////////////////////////////////////////////////////////////////////
//
//
// PostStepUtils - Utilities called after stepping is done
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef POSTSTEPUTILS_H
#define POSTSTEPUTILS_H

#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"
#include "larreco/RecoAlg/TCAlg/StepUtils.h"

namespace tca {

  // Check the quality of the trajectory and possibly trim it
  void CheckTraj(TCSlice& slc, Trajectory& tj);
  void CheckStiffEl(TCSlice& slc, Trajectory& tj);
  // defines HitPos, HitPosErr2 and Chg for the used hits in the trajectory point
  void ChkStopEnd1(TCSlice& slc, Trajectory& tj, bool prt);
  // Check the parameters at the start of the trajectory
  void ChkBegin(TCSlice& slc, Trajectory& tj);
  // Fix the parameters at the start of the trajectory
  void FixBegin(TCSlice& slc, Trajectory& tj, unsigned short atPt);
  bool IsGhost(TCSlice& slc, std::vector<unsigned int>& tHits);
  bool IsGhost(TCSlice& slc, Trajectory& tj);
  // Fill gaps in the trajectory
  void FillGaps(TCSlice& slc, Trajectory& tj);
  void TrimHiChgEndPts(TCSlice& slc, Trajectory& tj);
  // Trim high multiplicity TPs at the end
  void TrimHiMultEndPts(TCSlice& slc, Trajectory& tj);
  void ChkHiMultEndHits(TCSlice& slc, Trajectory& tj);
  void LastEndMerge(TCSlice& slc, CTP_t inCTP);
  void EndMerge(TCSlice& slc, CTP_t inCTP, bool lastPass);
  // Sets the StopsAtEnd bits for the trajectory
  void ChkStop(TCSlice& slc, Trajectory& tj);
  // Check the Michel electron topology, lastGoodPt is the last point of muon
  bool ChkMichel(TCSlice& slc, Trajectory& tj, unsigned short& lastGoodPt);
  // Make a junk trajectory using the list of hits in tHits
  bool MakeJunkTraj(TCSlice& slc, std::vector<unsigned int> tHits);
  void MergeShortWithJunk(TCSlice& slc, CTP_t inCTP);
  bool BraggSplit(TCSlice& slc, unsigned short itj);
  void ChkChgAsymmetry(TCSlice& slc, Trajectory& tj, bool prt);
} // namespace

#endif // ifndef POSTSTEPUTILS_H
