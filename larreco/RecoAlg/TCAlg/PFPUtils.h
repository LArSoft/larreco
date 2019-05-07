
////////////////////////////////////////////////////////////////////////
//
//
// PFParticle  utilities
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALGSPTUTILS_H
#define TRAJCLUSTERALGSPTUTILS_H

// C/C++ standard libraries
#include <array>
#include <vector>
#include <bitset>
#include <utility> // std::pair<>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TMatrixD.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"

namespace tca {

  void StitchPFPs();
  void FindPFParticles(TCSlice& slc);
  void ChkPFPMC(TCSlice& slc, PFPStruct& pfp);
  void ReconcileTPs(TCSlice& slc);
  void MakePFPTjs(TCSlice& slc);
  void MatchPlanes(TCSlice& slc, unsigned short numPlanes, std::vector<MatchStruct>& matVec, bool prt);
  bool Define(TCSlice& slc, PFPStruct& pfp, bool prt);
  void FixOrder(TCSlice& slc, PFPStruct& pfp, bool prt);
  bool Update(TCSlice& slc, PFPStruct& pfp, bool prt);
  bool ReSection(TCSlice& slc, PFPStruct& pfp, bool prt);
  void CountBadPoints(TCSlice& slc, PFPStruct& pfp, unsigned short fromPt, unsigned short toPt, unsigned short& nBadPts, unsigned short& firstBadPt);
  bool KillBadPoint(TCSlice& slc, PFPStruct& pfp, unsigned short fromPt, unsigned short toPt, bool prt);
  unsigned short Find3DRecoRange(TCSlice& slc, PFPStruct& pfp, unsigned short fromPt, unsigned short min2DPts, short dir);
  bool FitSection(TCSlice& slc, PFPStruct& pfp, unsigned short sfIndex);
  void ReconcileVertices(TCSlice& slc, PFPStruct& pfp, bool prt);
  void TrimEndPts(TCSlice& slc, PFPStruct& pfp, bool prt);
  unsigned short FirstPointInPlane(PFPStruct& pfp, unsigned short plane, unsigned short end);
  void FillGaps3D(TCSlice& slc, PFPStruct& pfp, bool prt);
  void AddPointsInRange(TCSlice& slc, PFPStruct& pfp, unsigned short fromPt, unsigned short toPt, CTP_t inCTP, float maxPull, float maxDelta, bool prt);
  bool InsertTP3D(PFPStruct& pfp, TP3D& tp3d);
  bool SortSection(PFPStruct& pfp, unsigned short sectionFitIndex);
  void CountOrder(TCSlice& slc, int tid, const std::vector<TP3D>& tp3ds, unsigned short& nNeg, unsigned short& nPos);
  void MakeTP3Ds(TCSlice& slc, PFPStruct& pfp);
  void Reverse(TCSlice& slc, PFPStruct& pfp);
  void FillmAllTraj(TCSlice& slc);
  bool SharesHighScoreVx(TCSlice& slc, const PFPStruct& pfp, const Trajectory& tj);
  bool MakeTp3(TCSlice& slc, const TrajPoint& itp, const TrajPoint& jtp, TrajPoint3& tp3, bool findDirection);
  double DeltaAngle(const Vector3_t v1, const Vector3_t v2);
  inline double DotProd(const Vector3_t& v1, const Vector3_t& v2) {return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; }
  Vector3_t PointDirection(const Point3_t p1, const Point3_t p2);
  double PosSep(const Point3_t& pos1, const Point3_t& pos2);
  double PosSep2(const Point3_t& pos1, const Point3_t& pos2);
  bool SetMag(Vector3_t& v1, double mag);
  void FilldEdx(TCSlice& slc, PFPStruct& pfp);
  float dEdx(TCSlice& slc, TP3D& tp3d);
  TP3D CreateTP3D(TCSlice& slc, unsigned short slHitsIndex, bool vlaInduction);
  TP3D CreateTP3D(TCSlice& slc, int tjID, unsigned short tjPt, bool vlaInduction);
  bool SetSection(TCSlice& slc, PFPStruct& pfp, TP3D& tp3d);
  float PointPull(const PFPStruct& pfp, TP3D& tp3d);
  PFPStruct CreatePFP(TCSlice& slc);
  void PFPVertexCheck(TCSlice& tcs);
  void DefinePFPParents(TCSlice& slc, bool prt);
//  void DefinePFPParentsTestBeam(TCSlice& slc, bool prt);
  bool StorePFP(TCSlice& slc, PFPStruct& pfp);
  bool InsideFV(TCSlice& slc, PFPStruct& pfp, unsigned short end);
  bool InsideTPC(const Point3_t& pos, geo::TPCID& inTPCID);
  void FindAlongTrans(Point3_t pos1, Vector3_t dir1, Point3_t pos2, Point2_t& alongTrans);
  bool PointDirIntersect(Point3_t p1, Vector3_t p1Dir, Point3_t p2, Vector3_t p2Dir, Point3_t& intersect, float& doca);
  bool LineLineIntersect(Point3_t p1, Point3_t p2, Point3_t p3, Point3_t p4, Point3_t& intersect, float& doca);
  float ChgFracBetween(TCSlice& slc, Point3_t pos1, Point3_t pos2);
  float ChgFracNearEnd(TCSlice& slc, PFPStruct& pfp, unsigned short end);
  Point3_t PosAtEnd(const PFPStruct& pfp, unsigned short end);
  Vector3_t DirAtEnd(const PFPStruct& pfp, unsigned short end);
  float Length(const PFPStruct& pfp);
  bool SectionStartEnd(const PFPStruct& pfp, unsigned short sfIndex, unsigned short& startPt, unsigned short& endPt);
  unsigned short FarEnd(TCSlice& slc, const PFPStruct& pfp, const Point3_t& pos);
  void PrintTP3Ds(std::string someText, TCSlice& slc, const PFPStruct& pfp, short printPts);
  
} // namespace tca

#endif // ifndef TRAJCLUSTERALGSPTUTILS_H
