////////////////////////////////////////////////////////////////////////
//
//
// TCAlg utilities
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALGUTILS_H
#define TRAJCLUSTERALGUTILS_H


// C/C++ standard libraries
#include <array>
#include <vector>
#include <bitset>
#include <utility> // std::pair<>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/RecoBase/Hit.h"
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"

namespace tca {

  // ****************************** General purpose  ******************************
  // Returns  true if there is a signal on the line between (wire1, time1) and (wire2, time2).
  bool SignalPresent(TjStuff& tjs, float wire1, float time1, TrajPoint const& tp, float minAmp);
  bool SignalPresent(TjStuff& tjs, unsigned int wire1, float time1, unsigned int wire2, float time2, CTP_t pCTP, float minAmp);
  bool SignalPresent(TjStuff& tjs, float wire1, float time1, float wire2, float time2, CTP_t pCTP, float minAmp);
  bool SignalPresent(TrajPoint const& tp, float minAmp);
  void MakeTrajectoryObsolete(TjStuff& tjs, unsigned short itj);
  void RestoreObsoleteTrajectory(TjStuff& tjs, unsigned short itj);
  // Split the allTraj trajectory itj at position pos into two trajectories
  // with an optional vertex assignment
  bool SplitAllTraj(TjStuff& tjs, unsigned short itj, unsigned short pos, unsigned short ivx, bool prt);
  void TrajClosestApproach(Trajectory const& tj, float x, float y, unsigned short& iClosePt, float& Distance);
  // returns the DOCA between a hit and a trajectory
  float PointTrajDOCA(TjStuff& tjs, unsigned int iht, TrajPoint const& tp);
  // returns the DOCA between a (W,T) point and a trajectory
  float PointTrajDOCA(TjStuff& tjs, float wire, float time, TrajPoint const& tp);
  // returns the DOCA^2 between a point and a trajectory
  float PointTrajDOCA2(TjStuff& tjs, float wire, float time, TrajPoint const& tp);
  // returns the separation^2 between two TPs
  float TrajPointHitSep2(TrajPoint const& tp1, TrajPoint const& tp2);
  // returns the separation^2 between a point and a TP
  float PointTrajSep2(float wire, float time, TrajPoint const& tp);
  // finds the point on trajectory tj that is closest to trajpoint tp
  void TrajPointTrajDOCA(TjStuff& tjs, TrajPoint const& tp, Trajectory const& tj, unsigned short& closePt, float& minSep);
  // returns the intersection position, intPos, of two trajectory points
  void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, float& x, float& y);
  // Returns the separation distance between two trajectory points
  float TrajPointSeparation(TrajPoint& tp1, TrajPoint& tp2);
  float TrajLength(Trajectory& tj);
  // returns the separation^2 between two hits in WSE units
  float HitSep2(TjStuff& tjs, unsigned int iht, unsigned int jht);
  // Find the Distance Of Closest Approach between two trajectories, exceeding minSep
  void TrajTrajDOCA(Trajectory const& tp1, Trajectory const& tp2, unsigned short& ipt1, unsigned short& ipt2, float& minSep);
  // Calculates the angle between two TPs
  float TwoTPAngle(TrajPoint& tp1, TrajPoint& tp2);
  // Put hits in each trajectory point into a flat vector. Only hits with UseHit if onlyUsedHits == true
  void PutTrajHitsInVector(Trajectory const& tj, bool onlyUsedHits, std::vector<unsigned int>& hitVec);
  // Project TP to a "wire position" Pos[0] and update Pos[1]
  void MoveTPToWire(TrajPoint& tp, float wire);
  void MaskTrajEndPoints(TjStuff& tjs, Trajectory& tj, unsigned short nPts);
  float DeltaAngle(float Ang1, float Ang2);
  // Find the first (last) TPs, EndPt[0] (EndPt[1], that have charge
  void SetEndPoints(TjStuff& tjs, Trajectory& tj);
  // ****************************** Printing  ******************************
  // Print trajectories, TPs, etc to mf::LogVerbatim
  void PrintTrajectory(TjStuff& tjs, Trajectory const& tj ,unsigned short tPoint);
  void PrintAllTraj(TjStuff& tjs, DebugStuff& Debug, unsigned short itj, unsigned short ipt);
  void PrintHeader();
  void PrintTrajPoint(TjStuff& tjs, unsigned short ipt, short dir, unsigned short pass, TrajPoint const& tp);
  // Print clusters after calling MakeAllTrajClusters
  void PrintClusters();
  // Print a single hit in the standard format
  std::string PrintHit(const art::Ptr<recob::Hit>& hit);
  // Print Trajectory position in the standard format
  std::string PrintPos(TjStuff& tjs, TrajPoint const& tp);
} // namespace tca

#endif // ifndef TRAJCLUSTERALGUTILS_H
