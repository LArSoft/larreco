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
#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/TCShower.h"
#include "larreco/RecoAlg/TCAlg/TCVertex.h"
#include "larreco/RecoAlg/TCAlg/PFPUtils.h"

namespace tca {

  typedef enum {
    kAllHits,
    kUsedHits,
    kUnusedHits,
  } HitStatus_t ;

  // ****************************** General purpose  ******************************
  
  void DefineTjParents(TCSlice& slc, bool prt);
  float MaxChargeAsymmetry(TCSlice& slc, std::vector<int>& tjIDs);
  int PDGCodeVote(TCSlice& slc, std::vector<int>& tjIDs, bool prt);
  unsigned short NumDeltaRays(TCSlice& slc, const Trajectory& tj);
  unsigned short NumDeltaRays(TCSlice& slc, std::vector<int>& tjIDs);
  int NeutrinoPrimaryTjID(TCSlice& slc, const Trajectory& tj);
  int PrimaryID(TCSlice& slc, const Trajectory& tj);
  int PrimaryUID(TCSlice& slc, const PFPStruct& pfp);
  bool MergeTjIntoPFP(TCSlice& slc, int mtjid, PFPStruct& pfp, bool prt);
  bool CompatibleMerge(TCSlice& slc, std::vector<int>& tjIDs, bool prt);
  bool CompatibleMerge(TCSlice& slc, const Trajectory& tj1, const Trajectory& tj2, bool prt);
  float OverlapFraction(TCSlice& slc, const Trajectory& tj1, const Trajectory& tj2);
  unsigned short AngleRange(TrajPoint const& tp);
  void SetAngleCode(TrajPoint& tp);
  unsigned short AngleRange(float angle);
  void FitTraj(TCSlice& slc, Trajectory& tj);
  void FitTraj(TCSlice& slc, Trajectory& tj, unsigned short originPt, unsigned short npts, short fitDir, TrajPoint& tpFit);
  float TjDirFOM(TCSlice& slc, const Trajectory& tj, bool prt);
//  void WatchHit(std::string someText, TCSlice& slc);
  unsigned short GetPFPIndex(TCSlice& slc, int tjID);
  unsigned short MatchVecIndex(TCSlice& slc, int tjID);
  void ReleaseHits(TCSlice& slc, Trajectory& tj);
  void UnsetUsedHits(TCSlice& slc, TrajPoint& tp);
  bool StoreTraj(TCSlice& slc, Trajectory& tj);
  void ChgSlope(TCSlice& slc, Trajectory& tj, float& slope, float& slopeErr, float& chiDOF);
  void ChgSlope(TCSlice& slc, Trajectory& tj, unsigned short fromPt, unsigned short toPt, float& slope, float& slopeErr, float& chiDOF);
  bool InTrajOK(TCSlice& slc, std::string someText);
  void CheckTrajBeginChg(TCSlice& slc, unsigned short itj);
  void TrimEndPts(std::string fcnLabel, TCSlice& slc, Trajectory& tj, const std::vector<float>& fQualityCuts, bool prt);
  void ChkChgAsymmetry(TCSlice& slc, Trajectory& tj, bool prt);
  bool SignalBetween(TCSlice& slc, const TrajPoint& tp1, const TrajPoint& tp2, const float& MinWireSignalFraction);
  bool SignalBetween(TCSlice& slc, TrajPoint tp, float toPos0, const float& MinWireSignalFraction);
  float ChgFracBetween(TCSlice& slc, TrajPoint tp, float toPos0);
  bool TrajHitsOK(TCSlice& slc, const std::vector<unsigned int>& iHitsInMultiplet, const std::vector<unsigned int>& jHitsInMultiplet);
  bool TrajHitsOK(TCSlice& slc, const unsigned int iht, const unsigned int jht);
  float ExpectedHitsRMS(TCSlice& slc, const TrajPoint& tp);
  bool SignalAtTp(TCSlice& slc, TrajPoint const& tp);
  float TpSumHitChg(TCSlice& slc, TrajPoint const& tp);
  unsigned short NumPtsWithCharge(TCSlice& slc, const Trajectory& tj, bool includeDeadWires);
  unsigned short NumPtsWithCharge(TCSlice& slc, const Trajectory& tj, bool includeDeadWires, unsigned short firstPt, unsigned short lastPt);
  float DeadWireCount(TCSlice& slc, const TrajPoint& tp1, const TrajPoint& tp2);
  float DeadWireCount(TCSlice& slc, const float& inWirePos1, const float& inWirePos2, CTP_t tCTP);
  unsigned short PDGCodeIndex(int PDGCode);
  void MakeTrajectoryObsolete(TCSlice& slc, unsigned int itj);
  void RestoreObsoleteTrajectory(TCSlice& slc, unsigned int itj);
  void MergeGhostTjs(TCSlice& slc, CTP_t inCTP);
  // Split the allTraj trajectory itj at position pos into two trajectories
  // with an optional vertex assignment
  bool SplitTraj(TCSlice& slc, unsigned short itj, unsigned short pos, unsigned short ivx, bool prt);
  bool SplitTraj(TCSlice& slc, unsigned short itj, float XPos, bool makeVx2, bool prt);
  bool TrajClosestApproach(Trajectory const& tj, float x, float y, unsigned short& closePt, float& DOCA);
  // returns the DOCA between a hit and a trajectory
  float PointTrajDOCA(TCSlice& slc, unsigned int iht, TrajPoint const& tp);
  // returns the DOCA between a (W,T) point and a trajectory
  float PointTrajDOCA(TCSlice& slc, float wire, float time, TrajPoint const& tp);
  // returns the DOCA^2 between a point and a trajectory
  float PointTrajDOCA2(TCSlice& slc, float wire, float time, TrajPoint const& tp);
  // Fills tp.Hits sets tp.UseHit true for hits that are close to tp.Pos. Returns true if there are
  // close hits OR if the wire at this position is dead
  bool FindCloseHits(TCSlice& slc, TrajPoint& tp, float const& maxDelta, HitStatus_t hitRequest);
  std::vector<unsigned int> FindCloseHits(TCSlice& slc, std::array<int, 2> const& wireWindow, Point2_t const& timeWindow, const unsigned short plane, HitStatus_t hitRequest, bool usePeakTime, bool& hitsNear);
  std::vector<int> FindCloseTjs(TCSlice& slc, const TrajPoint& fromTp, const TrajPoint& toTp, const float& maxDelta);
//  void PrimaryElectronLikelihood(TCSlice& slc, Trajectory& tj, float& likelihood, bool& flipDirection, bool prt);
  float ChgFracNearPos(TCSlice& slc, const Point2_t& pos, const std::vector<int>& tjIDs);
  float MaxHitDelta(TCSlice& slc, Trajectory& tj);
  void ReverseTraj(TCSlice& slc, Trajectory& tj);
  // returns the end of a trajectory that is closest to a point
  unsigned short CloseEnd(TCSlice& slc, const Trajectory& tj, const Point2_t& pos);
  // returns the separation^2 between a point and a TP
  float PointTrajSep2(float wire, float time, TrajPoint const& tp);
  float PosSep(const Point2_t& pos1, const Point2_t& pos2);
  float PosSep2(const Point2_t& pos1, const Point2_t& pos2);
  // finds the point on trajectory tj that is closest to trajpoint tp
  void TrajPointTrajDOCA(TCSlice& slc, TrajPoint const& tp, Trajectory const& tj, unsigned short& closePt, float& minSep);
  // returns the intersection position, intPos, of two trajectory points
  void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, Point2_t& pos);
  void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, float& x, float& y);
  float MaxTjLen(TCSlice& slc, std::vector<int>& tjIDs);
  // Returns the separation distance between two trajectory points
  float TrajPointSeparation(TrajPoint& tp1, TrajPoint& tp2);
  float TrajLength(Trajectory& tj);
  // returns the separation^2 between two hits in WSE units
  float HitSep2(TCSlice& slc, unsigned int iht, unsigned int jht);
  // Find the Distance Of Closest Approach between two trajectories, exceeding minSep
  bool TrajTrajDOCA(TCSlice& slc, const Trajectory& tp1, const Trajectory& tp2, unsigned short& ipt1, unsigned short& ipt2, float& minSep);
  bool TrajTrajDOCA(TCSlice& slc, const Trajectory& tp1, const Trajectory& tp2, unsigned short& ipt1, unsigned short& ipt2, float& minSep, bool considerDeadWires);
  // Calculates the angle between two TPs
  float TwoTPAngle(TrajPoint& tp1, TrajPoint& tp2);
  void TagJunkTj(TCSlice& slc, Trajectory& tj, bool prt);
  // Put hits in each trajectory point into a flat vector.
  std::vector<unsigned int> PutTrajHitsInVector(Trajectory const& tj, HitStatus_t hitRequest);
  // returns true if a hit is associated with more than one point
  bool HasDuplicateHits(TCSlice& slc, Trajectory const& tj, bool prt);
  // Project TP to a "wire position" Pos[0] and update Pos[1]
  void MoveTPToWire(TrajPoint& tp, float wire);
  bool PointInsideEnvelope(const Point2_t& Point, const std::vector<Point2_t>& Envelope);
  bool SetMag(Vector2_t& v1, double mag);
  void FindAlongTrans(Point2_t pos1, Vector2_t dir1, Point2_t pos2, Point2_t& alongTrans);
  inline double DotProd(const Vector2_t& v1, const Vector2_t& v2) {return v1[0]*v2[0] + v1[1]*v2[1]; }
  double DeltaAngle(double Ang1, double Ang2);
  double DeltaAngle2(double Ang1, double Ang2);
  double DeltaAngle(const Point2_t& p1, const Point2_t& p2);
  // Find the first (last) TPs, EndPt[0] (EndPt[1], that have charge
  void SetEndPoints(Trajectory& tj);
  // Returns the hit width using StartTick() and EndTick()
  float TPHitsRMSTick(TCSlice& slc, TrajPoint& tp, HitStatus_t hitRequest);
  float TPHitsRMSTime(TCSlice& slc, TrajPoint& tp, HitStatus_t hitRequest);
  float HitsRMSTick(TCSlice& slc, const std::vector<unsigned int>& hitsInMultiplet, HitStatus_t hitRequest);
  float HitsRMSTime(TCSlice& slc, const std::vector<unsigned int>& hitsInMultiplet, HitStatus_t hitRequest);
  float HitsPosTick(TCSlice& slc, const std::vector<unsigned int>& hitsInMultiplet, float& chg, HitStatus_t hitRequest);
  float HitsPosTime(TCSlice& slc, const std::vector<unsigned int>& hitsInMultiplet, float& chg, HitStatus_t hitRequest);
  unsigned short NumHitsInTP(const TrajPoint& tp, HitStatus_t hitRequest);
  unsigned short NumUsedHitsInTj(TCSlice& slc, const Trajectory& tj);
  unsigned short NearestPtWithChg(TCSlice& slc, Trajectory& tj, unsigned short thePt);
  // Calculate MCS momentum
  short MCSMom(TCSlice& slc, const std::vector<int>& tjIDs);
  short MCSMom(TCSlice& slc, Trajectory& tj);
  short MCSMom(TCSlice& slc, Trajectory& tj, unsigned short FirstPt, unsigned short lastPt);
  // Calculate MCS theta RMS over the points specified. Returns MCS angle for the full length
  double MCSThetaRMS(TCSlice& slc, Trajectory& tj, unsigned short firstPt, unsigned short lastPt);
  // Calculate MCS theta RMS over the entire length. Returns MCS angle for 1 WSE unit
  float MCSThetaRMS(TCSlice& slc, Trajectory& tj);
  void TjDeltaRMS(TCSlice& slc, Trajectory& tj, unsigned short firstPt, unsigned short lastPt, double& rms, unsigned short& cnt);
  // Returns true if the trajectory has low hit multiplicity and is in a clean environment
  bool TrajIsClean(TCSlice& slc, Trajectory& tj, bool prt);
  // Flag delta ray trajectories in allTraj
  void TagDeltaRays(TCSlice& slc, const CTP_t& inCTP);
  // Tag muon directions using delta proximity
//  void TagMuonDirections(TCSlice& slc, short debugWorkID);
  void UpdateTjChgProperties(std::string inFcnLabel, TCSlice& slc, Trajectory& tj, bool prt);
  void UpdateVxEnvironment(std::string inFcnLabel, TCSlice& slc, VtxStore& vx2, bool prt);
  // Make a bare trajectory point that only has position and direction defined
  TrajPoint MakeBareTP(TCSlice& slc, Point3_t& pos, Vector3_t& dir, CTP_t inCTP);
  bool MakeBareTrajPoint(TCSlice& slc, unsigned int fromHit, unsigned int toHit, TrajPoint& tp);
  bool MakeBareTrajPoint(TCSlice& slc, float fromWire, float fromTick, float toWire, float toTick, CTP_t tCTP, TrajPoint& tp);
  bool MakeBareTrajPoint(const Point2_t& fromPos, const Point2_t& toPos, TrajPoint& tpOut);
  bool MakeBareTrajPoint(TCSlice& slc, const TrajPoint& tpIn1, const TrajPoint& tpIn2, TrajPoint& tpOut);
  unsigned short FarEnd(TCSlice& slc, const Trajectory& tj, const Point2_t& pos);
  Vector2_t PointDirection(const Point2_t p1, const Point2_t p2);
  void SetPDGCode(TCSlice& slc, Trajectory& tj, bool tjDone);
  void SetPDGCode(TCSlice& slc, unsigned short itj, bool tjDone);
//  void AnalyzeHits(TCSlice& slc);
  bool AnalyzeHits();
  bool LongPulseHit(const recob::Hit& hit);
  bool FillWireHitRange(TCSlice& slc);
//  bool CheckWireHitRange(TCSlice& slc);
  bool WireHitRangeOK(TCSlice& slc, const CTP_t& inCTP);
  bool MergeAndStore(TCSlice& slc, unsigned int itj1, unsigned int itj2, bool doPrt);
  std::vector<int> GetAssns(TCSlice& slc, std::string type1Name, int id, std::string type2Name);
  // Start a trajectory going from fromHit to (toWire, toTick)
  bool StartTraj(TCSlice& slc, Trajectory& tj, unsigned int fromhit, unsigned int tohit, unsigned short pass);
  bool StartTraj(TCSlice& slc, Trajectory& tj, float fromWire, float fromTick, float toWire, float toTick, CTP_t& tCTP, unsigned short pass);
  bool Fit2D(short mode, Point2_t inPt, float& inPtErr, Vector2_t& outVec, Vector2_t& outVecErr, float& chiDOF);
  std::pair<unsigned short, unsigned short> GetSliceIndex(std::string typeName, int uID);
  template <typename T>
  std::vector<T> SetIntersection(const std::vector<T>& set1, const std::vector<T>& set2);
  template <typename T>
  std::vector<T> SetDifference(const std::vector<T>& set1, const std::vector<T>& set2);
  bool DecodeDebugString(std::string ctpwt);
  // ****************************** Printing  ******************************
  void DumpTj();
  void PrintAll(std::string someText, const std::vector<simb::MCParticle*>& mcpList);
  void PrintP(std::string someText, mf::LogVerbatim& myprt, PFPStruct& pfp, bool& printHeader);
  void Print3V(std::string someText, mf::LogVerbatim& myprt, Vtx3Store& vx3);
  void Print2V(std::string someText, mf::LogVerbatim& myprt, VtxStore& vx2);
  void Print3S(std::string someText, mf::LogVerbatim& myprt, ShowerStruct3D& ss3);
  void PrintT(std::string someText, mf::LogVerbatim& myprt, Trajectory& tj, bool& printHeader);
  void PrintTrajectory(std::string someText, TCSlice& slc, const Trajectory& tj ,unsigned short tPoint);
  void PrintAllTraj(std::string someText, TCSlice& slc, unsigned short itj, unsigned short ipt, bool printVtx = true);
  void PrintHeader(std::string someText);
  void PrintTrajPoint(std::string someText, TCSlice& slc, unsigned short ipt, short dir, unsigned short pass, TrajPoint const& tp);
  void PrintPFP(std::string someText, TCSlice& slc, const PFPStruct& pfp, bool printHeader);
  void PrintPFPs(std::string someText, TCSlice& slc);
  // Print clusters after calling MakeAllTrajClusters
  void PrintClusters();
  // Print a single hit in the standard format
  std::string PrintHit(const TCHit& hit);
  std::string PrintHitShort(const TCHit& hit);
  // Print Trajectory position in the standard format
  std::string PrintPos(TCSlice& slc, const TrajPoint& tp);
  std::string PrintPos(TCSlice& slc, const Point2_t& pos);
  std::string PrintStopFlag(const Trajectory& tj, unsigned short end);
  
  ////////////////////////////////////////////////
  template <typename T>
  std::vector<T> SetIntersection(const std::vector<T>& set1, const std::vector<T>& set2)
  {
    // returns a vector containing the elements of set1 and set2 that are common. This function
    // is a replacement for std::set_intersection which fails in the following situation:
    // set1 = {11 12 17 18} and set2 = {6 12 18}
    // There is no requirement that the elements be sorted, unlike std::set_intersection
    std::vector<T> shared;
    
    if(set1.empty()) return shared;
    if(set2.empty()) return shared;
    for(auto element1 : set1) {
      // check for a common element
      if(std::find(set2.begin(), set2.end(), element1) == set2.end()) continue;
      // check for a duplicate
      if(std::find(shared.begin(), shared.end(), element1) != shared.end()) continue;
      shared.push_back(element1);
    } // element1
    return shared;
  } // SetIntersection
  
  ////////////////////////////////////////////////
  template <typename T>
  std::vector<T> SetDifference(const std::vector<T>& set1, const std::vector<T>& set2)
  {
    // returns the elements of set1 and set2 that are different
    std::vector<T> different;
    if(set1.empty() && set2.empty()) return different;
    if(!set1.empty() && set2.empty()) return set1;
    if(set1.empty() && !set2.empty()) return set2;
    for(auto element1 : set1) {
      // check for a common element
      if(std::find(set2.begin(), set2.end(), element1) != set2.end()) continue;
      // check for a duplicate
      if(std::find(different.begin(), different.end(), element1) != different.end()) continue;
      different.push_back(element1);
    } // element1
    for(auto element2 : set2) {
      // check for a common element
      if(std::find(set1.begin(), set1.end(), element2) != set1.end()) continue;
      // check for a duplicate
      if(std::find(different.begin(), different.end(), element2) != different.end()) continue;
      different.push_back(element2);
    } // element1
    return different;
  } // SetDifference

} // namespace tca

#endif // ifndef TRAJCLUSTERALGUTILS_H
