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

// C++ standard libraries
#include <algorithm>
#include <array>
#include <string>
#include <utility>
#include <vector>

// Framework libraries
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
namespace detinfo {
  class DetectorClocksData;
  class DetectorPropertiesData;
}
namespace geo {
  struct TPCID;
}
namespace recob {
  class Hit;
}

bool valsDecreasing(const SortEntry& c1, const SortEntry& c2);
bool valsIncreasing(const SortEntry& c1, const SortEntry& c2);

namespace tca {

  typedef enum {
    kAllHits,
    kUsedHits,
    kUnusedHits,
  } HitStatus_t;

  // ****************************** General purpose  ******************************
  // dressed muons
  void MakeHaloTj(TCSlice& slc, Trajectory& muTj, bool prt);
  void DefineTjParents(TCSlice& slc, bool prt);
  float MaxChargeAsymmetry(TCSlice& slc, std::vector<int>& tjIDs);
  int PDGCodeVote(const TCSlice& slc, const std::vector<int>& tjIDs);
  int NeutrinoPrimaryTjID(const TCSlice& slc, const Trajectory& tj);
  int PrimaryID(const TCSlice& slc, const Trajectory& tj);
  int PrimaryUID(const TCSlice& slc, const PFPStruct& pfp);
  bool MergeTjIntoPFP(TCSlice& slc, int mtjid, PFPStruct& pfp, bool prt);
  float PointPull(TCSlice& slc, Point2_t pos, float chg, const Trajectory& tj);
  bool CompatibleMerge(const TCSlice& slc, std::vector<int>& tjIDs, bool prt);
  bool CompatibleMerge(const TCSlice& slc, const Trajectory& tj1, const Trajectory& tj2, bool prt);
  float OverlapFraction(const TCSlice& slc, const Trajectory& tj1, const Trajectory& tj2);
  unsigned short AngleRange(TrajPoint const& tp);
  void SetAngleCode(TrajPoint& tp);
  unsigned short AngleRange(float angle);
  void FitTraj(TCSlice& slc, Trajectory& tj);
  void FitTraj(TCSlice& slc,
               Trajectory& tj,
               unsigned short originPt,
               unsigned short npts,
               short fitDir,
               TrajPoint& tpFit);
  unsigned short GetPFPIndex(const TCSlice& slc, int tjID);
  void ReleaseHits(TCSlice& slc, Trajectory& tj);
  void UnsetUsedHits(TCSlice& slc, TrajPoint& tp);
  bool StoreTraj(TCSlice& slc, Trajectory& tj);
  void FitPar(const TCSlice& slc,
              const Trajectory& tj,
              unsigned short originPt,
              unsigned short npts,
              short fitDir,
              ParFit& pFit,
              unsigned short usePar);
  bool InTrajOK(TCSlice& slc, std::string someText);
  void CheckTrajBeginChg(TCSlice& slc, unsigned short itj);
  bool BraggSplit(TCSlice& slc, unsigned short itj);
  void ChkEndKink(TCSlice& slc, Trajectory& tj, bool prt);
  void TrimHiChgEndPts(TCSlice& slc, Trajectory& tj, bool prt);
  void TrimEndPts(std::string fcnLabel,
                  TCSlice& slc,
                  Trajectory& tj,
                  const std::vector<float>& fQualityCuts,
                  bool prt);
  void ChkMissedKink(TCSlice& slc, Trajectory& tj, bool prt);
  void ChkChgAsymmetry(TCSlice& slc, Trajectory& tj, bool prt);
  bool SignalBetween(const TCSlice& slc,
                     const TrajPoint& tp1,
                     const TrajPoint& tp2,
                     const float& MinWireSignalFraction);
  bool SignalBetween(const TCSlice& slc,
                     TrajPoint tp,
                     float toPos0,
                     const float& MinWireSignalFraction);
  float ChgFracBetween(const TCSlice& slc, TrajPoint tp, float toPos0);
  bool TrajHitsOK(TCSlice& slc,
                  const std::vector<unsigned int>& iHitsInMultiplet,
                  const std::vector<unsigned int>& jHitsInMultiplet);
  bool TrajHitsOK(TCSlice& slc, const unsigned int iht, const unsigned int jht);
  float ExpectedHitsRMS(TCSlice& slc, const TrajPoint& tp);
  bool SignalAtTpInSlc(const TCSlice& slc, const TrajPoint& tp);
  bool SignalAtTp(TrajPoint& tp);
  bool NearbySrcHit(geo::PlaneID plnID, unsigned int wire, float loTick, float hiTick);
  float TpSumHitChg(const TCSlice& slc, TrajPoint const& tp);
  unsigned short NumPtsWithCharge(const TCSlice& slc, const Trajectory& tj, bool includeDeadWires);
  unsigned short NumPtsWithCharge(const TCSlice& slc,
                                  const Trajectory& tj,
                                  bool includeDeadWires,
                                  unsigned short firstPt,
                                  unsigned short lastPt);
  float DeadWireCount(const TCSlice& slc, const TrajPoint& tp1, const TrajPoint& tp2);
  float DeadWireCount(const TCSlice& slc,
                      const float& inWirePos1,
                      const float& inWirePos2,
                      CTP_t tCTP);
  unsigned short PDGCodeIndex(int PDGCode);
  void MakeTrajectoryObsolete(TCSlice& slc, unsigned int itj);
  void RestoreObsoleteTrajectory(TCSlice& slc, unsigned int itj);
  void MergeGhostTjs(TCSlice& slc, CTP_t inCTP);
  // Split the allTraj trajectory itj at position pos into two trajectories
  // with an optional vertex assignment
  bool SplitTraj(TCSlice& slc,
                 unsigned short itj,
                 unsigned short pos,
                 unsigned short ivx,
                 bool prt);
  bool SplitTraj(detinfo::DetectorPropertiesData const& detProp,
                 TCSlice& slc,
                 unsigned short itj,
                 float XPos,
                 bool makeVx2,
                 bool prt);
  bool TrajClosestApproach(Trajectory const& tj,
                           float x,
                           float y,
                           unsigned short& closePt,
                           float& DOCA);
  // returns the DOCA between a hit and a trajectory
  float PointTrajDOCA(const TCSlice& slc, unsigned int iht, TrajPoint const& tp);
  // returns the DOCA between a (W,T) point and a trajectory
  float PointTrajDOCA(const TCSlice& slc, float wire, float time, TrajPoint const& tp);
  // returns the DOCA^2 between a point and a trajectory
  float PointTrajDOCA2(const TCSlice& slc, float wire, float time, TrajPoint const& tp);
  // Fills tp.Hits sets tp.UseHit true for hits that are close to tp.Pos. Returns true if there are
  // close hits OR if the wire at this position is dead
  bool FindCloseHits(TCSlice& slc, TrajPoint& tp, float const& maxDelta, HitStatus_t hitRequest);
  std::vector<unsigned int> FindCloseHits(const TCSlice& slc,
                                          std::array<int, 2> const& wireWindow,
                                          Point2_t const& timeWindow,
                                          const unsigned short plane,
                                          HitStatus_t hitRequest,
                                          bool usePeakTime,
                                          bool& hitsNear);
  unsigned short NearbyCleanPt(const TCSlice& slc, const Trajectory& tj, unsigned short nearPt);
  std::vector<int> FindCloseTjs(const TCSlice& slc,
                                const TrajPoint& fromTp,
                                const TrajPoint& toTp,
                                const float& maxDelta);
  float ElectronLikelihood(const TCSlice& slc, const Trajectory& tj);
  float KinkSignificance(TCSlice& slc,
                         Trajectory& tj1,
                         unsigned short end1,
                         Trajectory& tj2,
                         unsigned short end2,
                         unsigned short nPtsFit,
                         bool useChg,
                         bool prt);
  float KinkSignificance(TCSlice& slc,
                         Trajectory& tj,
                         unsigned short kinkPt,
                         unsigned short nPtsFit,
                         bool useChg,
                         bool prt);
  float ChgFracNearPos(const TCSlice& slc, const Point2_t& pos, const std::vector<int>& tjIDs);
  float MaxHitDelta(TCSlice& slc, Trajectory& tj);
  void ReverseTraj(TCSlice& slc, Trajectory& tj);
  // returns the end of a trajectory that is closest to a point
  unsigned short CloseEnd(const TCSlice& slc, const Trajectory& tj, const Point2_t& pos);
  // returns the separation^2 between a point and a TP
  float PointTrajSep2(float wire, float time, TrajPoint const& tp);
  float PosSep(const Point2_t& pos1, const Point2_t& pos2);
  float PosSep2(const Point2_t& pos1, const Point2_t& pos2);
  // finds the point on trajectory tj that is closest to trajpoint tp
  void TrajPointTrajDOCA(const TCSlice& slc,
                         TrajPoint const& tp,
                         Trajectory const& tj,
                         unsigned short& closePt,
                         float& minSep);
  // returns the intersection position, intPos, of two trajectory points
  void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, Point2_t& pos);
  void TrajIntersection(TrajPoint const& tp1, TrajPoint const& tp2, float& x, float& y);
  float MaxTjLen(const TCSlice& slc, std::vector<int>& tjIDs);
  // Returns the separation distance between two trajectory points
  float TrajPointSeparation(const TrajPoint& tp1, const TrajPoint& tp2);
  float TrajLength(const Trajectory& tj);
  // returns the separation^2 between two hits in WSE units
  float HitSep2(const TCSlice& slc, unsigned int iht, unsigned int jht);
  // Find the Distance Of Closest Approach between two trajectories, exceeding minSep
  bool TrajTrajDOCA(const TCSlice& slc,
                    const Trajectory& tp1,
                    const Trajectory& tp2,
                    unsigned short& ipt1,
                    unsigned short& ipt2,
                    float& minSep);
  bool TrajTrajDOCA(const TCSlice& slc,
                    const Trajectory& tp1,
                    const Trajectory& tp2,
                    unsigned short& ipt1,
                    unsigned short& ipt2,
                    float& minSep,
                    bool considerDeadWires);
  // Calculates the angle between two TPs
  float TwoTPAngle(const TrajPoint& tp1, const TrajPoint& tp2);
  void TagJunkTj(TCSlice& slc, Trajectory& tj, bool prt);
  std::vector<unsigned int> PutHitsInVector(const TCSlice& slc,
                                            PFPStruct const& pfp,
                                            HitStatus_t hitRequest);
  // Put hits in each trajectory point into a flat vector.
  std::vector<unsigned int> PutTrajHitsInVector(const Trajectory& tj, HitStatus_t hitRequest);
  // returns true if a hit is associated with more than one point
  bool HasDuplicateHits(const TCSlice& slc, Trajectory const& tj, bool prt);
  // Project TP to a "wire position" Pos[0] and update Pos[1]
  void MoveTPToWire(TrajPoint& tp, float wire);
  bool PointInsideEnvelope(const Point2_t& Point, const std::vector<Point2_t>& Envelope);
  bool SetMag(Vector2_t& v1, double mag);
  void FindAlongTrans(Point2_t pos1, Vector2_t dir1, Point2_t pos2, Point2_t& alongTrans);
  inline double
  DotProd(const Vector2_t& v1, const Vector2_t& v2)
  {
    return v1[0] * v2[0] + v1[1] * v2[1];
  }
  double DeltaAngle(double Ang1, double Ang2);
  double DeltaAngle2(double Ang1, double Ang2);
  double DeltaAngle(const Point2_t& p1, const Point2_t& p2);
  // Find the first (last) TPs, EndPt[0] (EndPt[1], that have charge
  void SetEndPoints(Trajectory& tj);
  // Returns the hit width using StartTick() and EndTick()
  float TPHitsRMSTick(const TCSlice& slc, const TrajPoint& tp, HitStatus_t hitRequest);
  float TPHitsRMSTime(const TCSlice& slc, const TrajPoint& tp, HitStatus_t hitRequest);
  float HitsRMSTick(const TCSlice& slc,
                    const std::vector<unsigned int>& hitsInMultiplet,
                    HitStatus_t hitRequest);
  float HitsRMSTime(const TCSlice& slc,
                    const std::vector<unsigned int>& hitsInMultiplet,
                    HitStatus_t hitRequest);
  float HitsPosTick(const TCSlice& slc,
                    const std::vector<unsigned int>& hitsInMultiplet,
                    float& chg,
                    HitStatus_t hitRequest);
  float HitsPosTime(const TCSlice& slc,
                    const std::vector<unsigned int>& hitsInMultiplet,
                    float& chg,
                    HitStatus_t hitRequest);
  unsigned short NumHitsInTP(const TrajPoint& tp, HitStatus_t hitRequest);
  unsigned short NumUsedHitsInTj(const TCSlice& slc, const Trajectory& tj);
  unsigned short NearestPtWithChg(const TCSlice& slc, const Trajectory& tj, unsigned short thePt);
  // Calculate MCS momentum
  short MCSMom(const TCSlice& slc, const std::vector<int>& tjIDs);
  short MCSMom(const TCSlice& slc, const Trajectory& tj);
  short MCSMom(const TCSlice& slc,
               const Trajectory& tj,
               unsigned short FirstPt,
               unsigned short lastPt);
  // Calculate MCS theta RMS over the points specified. Returns MCS angle for the full length
  double MCSThetaRMS(const TCSlice& slc,
                     const Trajectory& tj,
                     unsigned short firstPt,
                     unsigned short lastPt);
  // Calculate MCS theta RMS over the entire length. Returns MCS angle for 1 WSE unit
  float MCSThetaRMS(const TCSlice& slc, const Trajectory& tj);
  void TjDeltaRMS(const TCSlice& slc,
                  const Trajectory& tj,
                  unsigned short firstPt,
                  unsigned short lastPt,
                  double& rms,
                  unsigned short& cnt);
  void SetTPEnvironment(TCSlice& slc, CTP_t inCTP);
  // Returns true if the trajectory has low hit multiplicity and is in a clean environment
  bool TrajIsClean(TCSlice& slc, Trajectory& tj, bool prt);
  void UpdateTjChgProperties(std::string inFcnLabel, TCSlice& slc, Trajectory& tj, bool prt);
  void UpdateVxEnvironment(TCSlice& slc);
  void UpdateVxEnvironment(TCSlice& slc, VtxStore& vx2, bool prt);
  TrajPoint MakeBareTP(detinfo::DetectorPropertiesData const& detProp,
                       const TCSlice& slc,
                       const Point3_t& pos,
                       CTP_t inCTP);
  // Make a bare trajectory point that only has position and direction defined
  TrajPoint MakeBareTP(detinfo::DetectorPropertiesData const& detProp,
                       const TCSlice& slc,
                       const Point3_t& pos,
                       const Vector3_t& dir,
                       CTP_t inCTP);
  bool MakeBareTrajPoint(const TCSlice& slc,
                         unsigned int fromHit,
                         unsigned int toHit,
                         TrajPoint& tp);
  bool MakeBareTrajPoint(const TCSlice& slc,
                         float fromWire,
                         float fromTick,
                         float toWire,
                         float toTick,
                         CTP_t tCTP,
                         TrajPoint& tp);
  bool MakeBareTrajPoint(const Point2_t& fromPos, const Point2_t& toPos, TrajPoint& tpOut);
  bool MakeBareTrajPoint(const TCSlice& slc,
                         const TrajPoint& tpIn1,
                         const TrajPoint& tpIn2,
                         TrajPoint& tpOut);
  unsigned short FarEnd(TCSlice& slc, const Trajectory& tj, const Point2_t& pos);
  Vector2_t PointDirection(const Point2_t p1, const Point2_t p2);
  void SetPDGCode(TCSlice& slc, Trajectory& tj);
  void SetPDGCode(TCSlice& slc, unsigned short itj);
  bool AnalyzeHits();
  bool LongPulseHit(const recob::Hit& hit);
  void FillWireHitRange(geo::TPCID inTPCID);
  bool FillWireHitRange(detinfo::DetectorClocksData const& clockData,
                        detinfo::DetectorPropertiesData const& detProp,
                        TCSlice& slc);
  //  bool CheckWireHitRange(TCSlice& slc);
  bool WireHitRangeOK(TCSlice& slc, const CTP_t& inCTP);
  bool MergeAndStore(TCSlice& slc, unsigned int itj1, unsigned int itj2, bool doPrt);
  std::vector<int> GetAssns(TCSlice& slc, std::string type1Name, int id, std::string type2Name);
  // Start a trajectory going from fromHit to (toWire, toTick)
  bool StartTraj(TCSlice& slc,
                 Trajectory& tj,
                 unsigned int fromhit,
                 unsigned int tohit,
                 unsigned short pass);
  bool StartTraj(TCSlice& slc,
                 Trajectory& tj,
                 float fromWire,
                 float fromTick,
                 float toWire,
                 float toTick,
                 CTP_t& tCTP,
                 unsigned short pass);
  bool Fit2D(short mode,
             Point2_t inPt,
             float& inPtErr,
             Vector2_t& outVec,
             Vector2_t& outVecErr,
             float& chiDOF);
  std::pair<unsigned short, unsigned short> GetSliceIndex(std::string typeName, int uID);
  template <typename T>
  std::vector<T> SetIntersection(const std::vector<T>& set1, const std::vector<T>& set2);
  template <typename T>
  std::vector<T> SetDifference(const std::vector<T>& set1, const std::vector<T>& set2);
  bool DecodeDebugString(std::string ctpwt);
  // ****************************** Printing  ******************************
  void DumpTj();
  void PrintDebugMode();
  void PrintAll(detinfo::DetectorPropertiesData const& detProp, std::string someText);
  void PrintP(std::string someText, mf::LogVerbatim& myprt, PFPStruct& pfp, bool& printHeader);
  void Print3V(detinfo::DetectorPropertiesData const& detProp,
               std::string someText,
               mf::LogVerbatim& myprt,
               Vtx3Store& vx3,
               bool& printHeader);
  void Print2V(std::string someText, mf::LogVerbatim& myprt, VtxStore& vx2, bool& printHeader);
  void Print3S(detinfo::DetectorPropertiesData const& detProp,
               std::string someText,
               mf::LogVerbatim& myprt,
               ShowerStruct3D& ss3);
  void PrintT(std::string someText, mf::LogVerbatim& myprt, Trajectory& tj, bool& printHeader);
  void PrintTrajectory(std::string someText,
                       const TCSlice& slc,
                       const Trajectory& tj,
                       unsigned short tPoint);
  void PrintAllTraj(detinfo::DetectorPropertiesData const& detProp,
                    std::string someText,
                    TCSlice& slc,
                    unsigned short itj,
                    unsigned short ipt,
                    bool printVtx = true);
  void PrintTPHeader(std::string someText);
  void PrintTP(std::string someText,
               const TCSlice& slc,
               unsigned short ipt,
               short dir,
               unsigned short pass,
               const TrajPoint& tp);
  std::string TPEnvString(const TrajPoint& tp);
  void PrintPFP(std::string someText, TCSlice& slc, const PFPStruct& pfp, bool printHeader);
  void PrintPFPs(std::string someText, TCSlice& slc);
  // Print clusters after calling MakeAllTrajClusters
  void PrintClusters();
  // Print a single hit in the standard format
  std::string PrintHit(const TCHit& hit);
  std::string PrintHitShort(const TCHit& hit);
  // Print Trajectory position in the standard format
  std::string PrintPos(const TCSlice& slc, const TrajPoint& tp);
  std::string PrintPos(const TCSlice& slc, const Point2_t& pos);
  std::string PrintEndFlag(const Trajectory& tj, unsigned short end);
  std::string PrintEndFlag(const PFPStruct& pfp, unsigned short end);

  ////////////////////////////////////////////////
  template <typename T>
  std::vector<T>
  SetIntersection(const std::vector<T>& set1, const std::vector<T>& set2)
  {
    // returns a vector containing the elements of set1 and set2 that are common. This function
    // is a replacement for std::set_intersection which fails in the following situation:
    // set1 = {11 12 17 18} and set2 = {6 12 18}
    // There is no requirement that the elements be sorted, unlike std::set_intersection
    std::vector<T> shared;

    if (set1.empty()) return shared;
    if (set2.empty()) return shared;
    for (auto element1 : set1) {
      // check for a common element
      if (std::find(set2.begin(), set2.end(), element1) == set2.end()) continue;
      // check for a duplicate
      if (std::find(shared.begin(), shared.end(), element1) != shared.end()) continue;
      shared.push_back(element1);
    } // element1
    return shared;
  } // SetIntersection

  ////////////////////////////////////////////////
  template <typename T>
  std::vector<T>
  SetDifference(const std::vector<T>& set1, const std::vector<T>& set2)
  {
    // returns the elements of set1 and set2 that are different
    std::vector<T> different;
    if (set1.empty() && set2.empty()) return different;
    if (!set1.empty() && set2.empty()) return set1;
    if (set1.empty() && !set2.empty()) return set2;
    for (auto element1 : set1) {
      // check for a common element
      if (std::find(set2.begin(), set2.end(), element1) != set2.end()) continue;
      // check for a duplicate
      if (std::find(different.begin(), different.end(), element1) != different.end()) continue;
      different.push_back(element1);
    } // element1
    for (auto element2 : set2) {
      // check for a common element
      if (std::find(set1.begin(), set1.end(), element2) != set1.end()) continue;
      // check for a duplicate
      if (std::find(different.begin(), different.end(), element2) != different.end()) continue;
      different.push_back(element2);
    } // element1
    return different;
  } // SetDifference

} // namespace tca

#endif // ifndef TRAJCLUSTERALGUTILS_H
