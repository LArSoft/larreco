
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

#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TCAlg/DataStructs.h"
#include "larreco/RecoAlg/TCAlg/DebugStruct.h"
#include "larreco/RecoAlg/TCAlg/Utils.h"

namespace tca {

  bool Repair(TjStuff& tjs, PFPStruct& pfp, int tjNotInVx, bool prt);
  void UpdateMatchStructs(TjStuff& tjs, std::vector<int> oldTjs, int newTj);
  void FillMatchVectors(TjStuff& tjs, const geo::TPCID& tpcid, bool prt);
  void FillmAllTraj(TjStuff& tjs, const geo::TPCID& tpcid);
  void MakePFPTp3s(TjStuff& tjs, PFPStruct& pfp, bool anyTj);
  bool SetNewStart(TjStuff& tjs, PFPStruct& pfp, bool prt);
  void SortByDistanceFromStart(TjStuff& tjs, PFPStruct& pfp);
  bool CheckTp3Validity(TjStuff& tjs, PFPStruct& pfp, Vector3_t generalDirection, bool prt);
  bool FitTp3(TjStuff& tjs, std::vector<TrajPoint3> tp3s, unsigned short originPt, unsigned short npts, short fitDir, TrajPoint3& outTp3);
  void FindXMatches(TjStuff& tjs, unsigned short numPlanes, short maxScore, std::vector<MatchStruct>& matVec, bool prt);
  bool MakeTp3(TjStuff& tjs, const TrajPoint& itp, const TrajPoint& jtp, TrajPoint3& tp3);
  double DeltaAngle(const Vector3_t v1, const Vector3_t v2);
  inline double DotProd(const Vector3_t& v1, const Vector3_t& v2) {return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; }
  Vector3_t PointDirection(const Point3_t p1, const Point3_t p2);
  double PosSep(const Point3_t& pos1, const Point3_t& pos2);
  double PosSep2(const Point3_t& pos1, const Point3_t& pos2);
  bool SetMag(Vector3_t& v1, double mag);
  void FilldEdx(TjStuff& tjs, TrajPoint3& tp3);
  void SplitAtKink(TjStuff& tjs, PFPStruct& pfp, double sep, bool prt);
  std::vector<unsigned short> FindKinks(const TjStuff& tjs, PFPStruct& pfp, double sep, bool prt);
  double KinkAngle(const TjStuff& tjs, const std::vector<TrajPoint3>& tp3s, unsigned short atPt, double sep);
  void PrintTp3(std::string fcnLabel, const TjStuff& tjs, const TrajPoint3& tp3);
  PFPStruct CreatePFP(const TjStuff& tjs, const geo::TPCID& tpcid);
  bool DefinePFP(TjStuff& tjs, PFPStruct& pfp, bool prt);
  void DefinePFPParents(TjStuff& tjs, const geo::TPCID& tpcid, bool prt);
  bool StorePFP(TjStuff& tjs, PFPStruct& pfp);
  void SetStopFlags(TjStuff& tjs, PFPStruct& pfp, bool prt);
  bool InsideTPC(const TjStuff& tjs, Point3_t& pos, geo::TPCID& inTPCID);
  
} // namespace tca

#endif // ifndef TRAJCLUSTERALGSPTUTILS_H
