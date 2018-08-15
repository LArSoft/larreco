
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
  void UpdateMatchStructs(TCSlice& slc, int oldTj, int newTj);
  void UpdateTp3s(TCSlice& slc, PFPStruct& pfp, int oldTj, int newTj);
  void FillmAllTraj(TCSlice& slc);
  bool SetStart(TCSlice& slc, PFPStruct& pfp, bool prt);
  void FollowTp3s(TCSlice& slc, PFPStruct& pfp, bool prt);
  bool FitTp3s(TCSlice& slc, const std::vector<TrajPoint3>& tp3s, Point3_t& pos, Vector3_t& dir, float& rCorr);
  bool FitTp3s(TCSlice& slc, const std::vector<TrajPoint3>& tp3s, unsigned short fromPt, unsigned short toPt, Point3_t& pos, Vector3_t& dir, float& rCorr);
  bool FitTp3(TCSlice& slc, TrajPoint3& tp3, const std::vector<Tj2Pt>& tj2pts);
  void FindCompleteness(TCSlice& slc, PFPStruct& pfp, bool doFit, bool fillTp3s, bool prt);
  void FindMissedTjsInTp3s(TCSlice& slc, PFPStruct& pfp, std::vector<int>& missTjs, std::vector<float>& missFrac);
  bool SharesHighScoreVx(TCSlice& slc, const PFPStruct& pfp, const Trajectory& tj);
  void Fit3D(unsigned short mode, Point3_t point, Vector3_t dir, Point3_t& fitPos, Vector3_t& fitDir);
//  bool CheckAndMerge(TCSlice& slc, PFPStruct& pfp, bool prt);
  float AspectRatio(TCSlice& slc, std::vector<int>& tjids, CTP_t inCTP);
  unsigned short WiresSkippedInCTP(TCSlice& slc, std::vector<int>& tjids, CTP_t inCTP);
  float LengthInCTP(TCSlice& slc, std::vector<int>& tjids, CTP_t inCTP);
  bool AddMissedTj(TCSlice& slc, PFPStruct& pfp, unsigned short itj, bool looseCuts, bool prt);
  void CleanTjs(TCSlice& slc, PFPStruct& pfp, bool prt);
  bool MergePFPTjs(TCSlice& slc, PFPStruct& pfp, bool prt);
  void FindXMatches(TCSlice& slc, unsigned short numPlanes, short maxScore, std::vector<MatchStruct>& matVec, bool prt);
  bool MakeTp3(TCSlice& slc, const TrajPoint& itp, const TrajPoint& jtp, TrajPoint3& tp3, bool findDirection);
  double DeltaAngle(const Vector3_t v1, const Vector3_t v2);
  inline double DotProd(const Vector3_t& v1, const Vector3_t& v2) {return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; }
  Vector3_t PointDirection(const Point3_t p1, const Point3_t p2);
  double PosSep(const Point3_t& pos1, const Point3_t& pos2);
  double PosSep2(const Point3_t& pos1, const Point3_t& pos2);
  bool SetMag(Vector3_t& v1, double mag);
  void FilldEdx(TCSlice& slc, PFPStruct& pfp);
  void FilldEdx(TCSlice& slc, TrajPoint3& tp3);
  float PFPDOCA(const PFPStruct& pfp1,  const PFPStruct& pfp2, unsigned short& close1, unsigned short& close2);
  bool Split3DKink(TCSlice& slc, PFPStruct& pfp, double sep, bool prt);
  std::vector<unsigned short> FindKinks(TCSlice& slc, PFPStruct& pfp, double sep, bool prt);
  double KinkAngle(TCSlice& slc, const std::vector<TrajPoint3>& tp3s, unsigned short atPt, double sep);
  PFPStruct CreatePFP(TCSlice& slc);
  void FindPFParticles(TCSlice& slc);
  bool DefinePFP(std::string inFcnLabel, TCSlice& slc, PFPStruct& pfp, bool prt);
  bool PFPVxTjOK(TCSlice& slc, PFPStruct& pfp, bool prt);
  void PFPVertexCheck(TCSlice& tcs);
  bool AnalyzePFP(TCSlice& slc, PFPStruct& pfp, bool prt);
  void DefinePFPParents(TCSlice& slc, bool prt);
  void DefinePFPParentsTestBeam(TCSlice& slc, bool prt);
  bool StorePFP(TCSlice& slc, PFPStruct& pfp);
  bool InsideFV(TCSlice& slc, PFPStruct& pfp, unsigned short end);
  bool InsideTPC(const Point3_t& pos, geo::TPCID& inTPCID);
  void FindAlongTrans(Point3_t pos1, Vector3_t dir1, Point3_t pos2, Point2_t& alongTrans);
  bool PointDirIntersect(Point3_t p1, Vector3_t p1Dir, Point3_t p2, Vector3_t p2Dir, Point3_t& intersect, float& doca);
  bool LineLineIntersect(Point3_t p1, Point3_t p2, Point3_t p3, Point3_t p4, Point3_t& intersect, float& doca);
  void ReversePFP(TCSlice& slc, PFPStruct& pfp);
  float ChgFracBetween(TCSlice& slc, Point3_t pos1, Point3_t pos2);
  float ChgFracNearEnd(TCSlice& slc, PFPStruct& pfp, unsigned short end);
  unsigned short FarEnd(TCSlice& slc, const PFPStruct& pfp, const Point3_t& pos);
  void PrintTp3(std::string fcnLabel, TCSlice& slc, const TrajPoint3& tp3);
  void PrintTp3s(std::string someText, TCSlice& slc, const PFPStruct& pfp, short printPts);
  
} // namespace tca

#endif // ifndef TRAJCLUSTERALGSPTUTILS_H
