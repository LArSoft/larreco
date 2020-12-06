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

namespace detinfo {
  class DetectorPropertiesData;
}

namespace tca {

  // Check the quality of the trajectory and possibly trim it
  void CheckTraj(TCSlice& slc, Trajectory& tj);
  // Check the parameters at the start of the trajectory
  void ChkBegin(TCSlice& slc, Trajectory& tj);
  void ChkBeginOld(TCSlice& slc, Trajectory& tj);
  void ChkBeginChg(TCSlice& slc, unsigned short itj);
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
  void MergeJunk(TCSlice& slc, CTP_t inCTP);
  void ChkChgAsymmetry(TCSlice& slc, Trajectory& tj, bool prt);
  void TagShowerLike(TCSlice& slc, const CTP_t& inCTP);
} // namespace

#endif // ifndef POSTSTEPUTILS_H
