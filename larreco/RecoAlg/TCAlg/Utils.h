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
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"

namespace tca {

  typedef enum {
    kAllHits,
    kUsedHits,
    kUnusedHits,
  } HitStatus_t ;

  // ****************************** General purpose  ******************************
  void WatchHit(std::string someText, TjStuff& tjs, const unsigned int& watchHit, short& watchInTraj, const unsigned short& tjID);
  // Return true if the 3D matched trajectories in tjs.matchVecPFPList are in the wrong order in terms of
  // physics standpoint, e.g. dQ/dx, muon delta-ray tag, cosmic rays entering the detector, etc
  bool Reverse3DMatchTjs(TjStuff& tjs, unsigned short im, bool prt);
  unsigned short Matched3DVtx(TjStuff& tjs, unsigned short im);
  void ReleaseHits(TjStuff& tjs, Trajectory& tj);
  void UnsetUsedHits(TjStuff& tjs, TrajPoint& tp);
  void TrimEndPts(TjStuff& tjs, Trajectory& tj, const std::vector<float>& fQualityCuts, bool prt);
  bool SignalBetween(TjStuff& tjs, const TrajPoint& tp1, const TrajPoint& tp2, const float& MinWireSignalFraction, bool prt);
  bool SignalBetween(TjStuff& tjs, TrajPoint tp, float toPos0, const float& MinWireSignalFraction, bool prt);
  bool SignalAtTp(TjStuff& tjs, TrajPoint const& tp);
  bool SignalAtPos(TjStuff& tjs, const float& pos0, const float& pos1, CTP_t tCTP);
  float TpSumHitChg(TjStuff& tjs, TrajPoint const& tp);
  bool CheckHitClusterAssociations(TjStuff& tjs);
  unsigned short NumPtsWithCharge(TjStuff& tjs, const Trajectory& tj, bool includeDeadWires);
  unsigned short NumPtsWithCharge(TjStuff& tjs, const Trajectory& tj, bool includeDeadWires, unsigned short firstPt, unsigned short lastPt);
  float DeadWireCount(TjStuff& tjs, const TrajPoint& tp1, const TrajPoint& tp2);
  float DeadWireCount(TjStuff& tjs, const float& inWirePos1, const float& inWirePos2, CTP_t tCTP);
  unsigned short PDGCodeIndex(TjStuff& tjs, int PDGCode);
  bool WireHitRangeOK(const TjStuff& tjs, const CTP_t& inCTP);
  // Returns  true if there is a signal on the line between (wire1, time1) and (wire2, time2).
  bool SignalPresent(TjStuff& tjs, float wire1, float time1, TrajPoint const& tp, float minAmp);
  bool SignalPresent(TjStuff& tjs, unsigned int wire1, float time1, unsigned int wire2, float time2, CTP_t pCTP, float minAmp);
  bool SignalPresent(TjStuff& tjs, float wire1, float time1, float wire2, float time2, CTP_t pCTP, float minAmp);
  bool SignalPresent(TrajPoint const& tp, float minAmp);
  void MakeTrajectoryObsolete(TjStuff& tjs, unsigned short itj);
  void MakeVertexObsolete(TjStuff& tjs, unsigned short ivx);
  void RestoreObsoleteTrajectory(TjStuff& tjs, unsigned short itj);
  void CheckVtxAssociations(TjStuff& tjs, const CTP_t& inCTP);
  bool TjHasNiceVtx(TjStuff& tjs, const Trajectory& tj);
  // Split the allTraj trajectory itj at position pos into two trajectories
  // with an optional vertex assignment
  bool SplitAllTraj(TjStuff& tjs, unsigned short itj, unsigned short pos, unsigned short ivx, bool prt);
  bool SplitAllTraj(TjStuff& tjs, Trajectory& tj, unsigned short pos, unsigned short ivx, bool prt);
  void TrajClosestApproach(Trajectory const& tj, float x, float y, unsigned short& iClosePt, float& Distance);
  // returns the DOCA between a hit and a trajectory
  float PointTrajDOCA(TjStuff const& tjs, unsigned int iht, TrajPoint const& tp);
  // returns the DOCA between a (W,T) point and a trajectory
  float PointTrajDOCA(TjStuff const& tjs, float wire, float time, TrajPoint const& tp);
  // returns the DOCA^2 between a point and a trajectory
  float PointTrajDOCA2(TjStuff const& tjs, float wire, float time, TrajPoint const& tp);
  // Fills tp.Hits sets tp.UseHit true for hits that are close to tp.Pos. Returns true if there are
  // close hits OR if the wire at this position is dead
  bool FindCloseHits(TjStuff const& tjs, TrajPoint& tp, float const& maxDelta, HitStatus_t hitRequest);
  std::vector<unsigned int> FindCloseHits(TjStuff const& tjs, std::array<int, 2> const& wireWindow, std::array<float, 2> const& timeWindow, const unsigned short plane, HitStatus_t hitRequest, bool usePeakTime, bool& hitsNear);
  float MaxHitDelta(TjStuff& tjs, Trajectory& tj);
  void ReverseTraj(TjStuff& tjs, Trajectory& tj);
  // returns the separation^2 between a point and a TP
  float PointTrajSep2(float wire, float time, TrajPoint const& tp);
  float PosSep(const std::array<float, 2>& pos1, const std::array<float, 2>& pos2);
  float PosSep2(const std::array<float, 2>& pos1, const std::array<float, 2>& pos2);
  float PosSep2(const std::array<float, 3>& pos1, const std::array<float, 3>& pos2);
  // finds the point on trajectory tj that is closest to trajpoint tp
  void TrajPointTrajDOCA(TjStuff& tjs, TrajPoint const& tp, Trajectory const& tj, unsigned short& closePt, float& minSep);
  // returns the intersection position, intPos, of two trajectory points
  void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, std::array<float, 2>& pos);
  void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, float& x, float& y);
  // Returns the separation distance between two trajectory points
  float TrajPointSeparation(TrajPoint& tp1, TrajPoint& tp2);
  float TrajLength(Trajectory& tj);
  // returns the separation^2 between two hits in WSE units
  float HitSep2(TjStuff& tjs, unsigned int iht, unsigned int jht);
  // Find the Distance Of Closest Approach between two trajectories, exceeding minSep
  void TrajTrajDOCA(TjStuff& tjs, Trajectory const& tp1, Trajectory const& tp2, unsigned short& ipt1, unsigned short& ipt2, float& minSep);
  void TrajTrajDOCA(TjStuff& tjs, Trajectory const& tp1, Trajectory const& tp2, unsigned short& ipt1, unsigned short& ipt2, float& minSep, bool considerDeadWires);
  // Calculates the angle between two TPs
  float TwoTPAngle(TrajPoint& tp1, TrajPoint& tp2);
  // Put hits in each trajectory point into a flat vector.
  std::vector<unsigned int> PutTrajHitsInVector(Trajectory const& tj, HitStatus_t hitRequest);
  // returns true if hit iht appears in trajectory tj. The last nPtsToCheck points are checked
  bool HitIsInTj(Trajectory const& tj, const unsigned int& iht, short nPtsToCheck);
  // returns true if a hit is associated with more than one point
  bool HasDuplicateHits(TjStuff const& tjs, Trajectory const& tj, bool prt);
  // Project TP to a "wire position" Pos[0] and update Pos[1]
  void MoveTPToWire(TrajPoint& tp, float wire);
  bool PointInsideEnvelope(const std::array<float, 2>& Point, const std::vector<std::array<float, 2>>& Envelope);
  double DeltaAngle(double Ang1, double Ang2);
  double DeltaAngle2(double Ang1, double Ang2);
  double DeltaAngle(const std::array<float,2>& p1, const std::array<float,2>& p2);
  // Find the first (last) TPs, EndPt[0] (EndPt[1], that have charge
  void SetEndPoints(TjStuff& tjs, Trajectory& tj);
  // Returns the hit width using StartTick() and EndTick()
  float TPHitsRMSTick(TjStuff& tjs, TrajPoint& tp, HitStatus_t hitRequest);
  float TPHitsRMSTime(TjStuff& tjs, TrajPoint& tp, HitStatus_t hitRequest);
  float HitsRMSTick(TjStuff& tjs, const std::vector<unsigned int>& hitsInMultiplet, HitStatus_t hitRequest);
  float HitsRMSTime(TjStuff& tjs, const std::vector<unsigned int>& hitsInMultiplet, HitStatus_t hitRequest);
  float HitsPosTick(TjStuff& tjs, const std::vector<unsigned int>& hitsInMultiplet, float& chg, HitStatus_t hitRequest);
  float HitsPosTime(TjStuff& tjs, const std::vector<unsigned int>& hitsInMultiplet, float& chg, HitStatus_t hitRequest);
  unsigned short NumHitsInTP(const TrajPoint& tp, HitStatus_t hitRequest);
  unsigned short NearestPtWithChg(TjStuff& tjs, Trajectory& tj, unsigned short thePt);
  // Calculate MCS momentum
  short MCSMom(TjStuff& tjs, Trajectory& tj);
  short MCSMom(TjStuff& tjs, Trajectory& tj, unsigned short FirstPt, unsigned short lastPt);
  // Calculate MCS theta RMS over the points specified. Returns MCS angle for the full length
  double MCSThetaRMS(TjStuff& tjs, Trajectory& tj, unsigned short firstPt, unsigned short lastPt);
  // Calculate MCS theta RMS over the entire length. Returns MCS angle for 1 WSE unit
  float MCSThetaRMS(TjStuff& tjs, Trajectory& tj);
  // Returns true if the trajectory has low hit multiplicity and is in a clean environment
  bool TrajIsClean(TjStuff& tjs, Trajectory& tj, bool prt);
  // Flag delta ray trajectories in allTraj
  void TagDeltaRays(TjStuff& tjs, const CTP_t& inCTP, const std::vector<short>& fDeltaRayTag, short debugWorkID);
  // Tag muon directions using delta proximity
  void TagMuonDirections(TjStuff& tjs, const short& minDeltaRayLength, short debugWorkID);
  // Make a bare trajectory point that only has position and direction defined
  bool MakeBareTrajPoint(TjStuff& tjs, unsigned int fromHit, unsigned int toHit, TrajPoint& tp);
  bool MakeBareTrajPoint(TjStuff& tjs, float fromWire, float fromTick, float toWire, float toTick, CTP_t tCTP, TrajPoint& tp);
  bool MakeBareTrajPoint(TjStuff& tjs, const TrajPoint& tpIn1, const TrajPoint& tpIn2, TrajPoint& tpOut);
  // ****************************** Shower finding  ******************************
  // Create showers (aka clusters of trajectories, tjs.cots)
  void FindShowers(TjStuff& tjs, const CTP_t& inCTP, const std::vector<float>& fShowerTag);
  void TagShowerTjs(TjStuff& tjs, const CTP_t& inCTP, const std::vector<float>& fShowerTag, std::vector<std::vector<unsigned short>>& tjList);
  void DefineShowerTj(TjStuff& tjs, const CTP_t& inCTP, const unsigned short& showerIndex, const std::vector<float>& fShowerTag, const CTP_t& printCTP);
  void FindShowerParent(TjStuff& tjs, const CTP_t& inCTP, const unsigned short& showerIndex, const std::vector<float>& fShowerTag, const CTP_t& printCTP);
  void CollectHits(TjStuff& tjs, const CTP_t& inCTP, const unsigned short& showerIndex, const CTP_t& printCTP);
  // ****************************** Vertex finding  ******************************
  unsigned short TPNearVertex(TjStuff& tjs, const TrajPoint& tp);
  bool AttachAnyTrajToVertex(TjStuff& tjs, unsigned short iv, const std::vector<float>& fVertex2DCuts, bool prt);
  bool AttachTrajToAnyVertex(TjStuff& tjs, unsigned short itj, const std::vector<float>& fVertex2DCuts, bool prt);
  bool AttachTrajToVertex(TjStuff& tjs, Trajectory& tj, VtxStore& vx, const std::vector<float>& fVertex2DCuts, bool prt);
  float TrajPointVertexPull(TjStuff& tjs, const TrajPoint& tp, const VtxStore& vx);
  float VertexVertexPull(TjStuff& tjs, const VtxStore& vx1, const VtxStore& vx2);
  bool FitVertex(TjStuff& tjs, VtxStore& vx, const std::vector<float>& fVertex2DCuts, bool prt);
  // ****************************** Printing  ******************************
  // Print trajectories, TPs, etc to mf::LogVerbatim
  void PrintTrajectory(std::string someText, TjStuff& tjs, Trajectory const& tj ,unsigned short tPoint);
  void PrintAllTraj(std::string someText, TjStuff& tjs, DebugStuff& Debug, unsigned short itj, unsigned short ipt);
  void PrintHeader(std::string someText);
  void PrintTrajPoint(std::string someText, TjStuff& tjs, unsigned short ipt, short dir, unsigned short pass, TrajPoint const& tp);
  // Print clusters after calling MakeAllTrajClusters
  void PrintClusters();
  // Print a single hit in the standard format
  std::string PrintHit(const TCHit& hit);
  std::string PrintHitShort(const TCHit& hit);
  // Print Trajectory position in the standard format
  std::string PrintPos(TjStuff& tjs, const TrajPoint& tp);
  std::string PrintPos(TjStuff& tjs, const std::array<float, 2>& pos);
  std::string PrintStopFlag(TjStuff& tjs, const Trajectory& tj, unsigned short end);
} // namespace tca

#endif // ifndef TRAJCLUSTERALGUTILS_H
