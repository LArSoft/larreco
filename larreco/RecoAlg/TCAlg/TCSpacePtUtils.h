
////////////////////////////////////////////////////////////////////////
//
//
// TCAlg SpacePoint utilities
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
  
  void Match3DSpts(TjStuff& tjs, const art::Event& evt, const geo::TPCID& tpcid, const art::InputTag& fSpacePointModuleLabel);
  bool SetPFPEndPoints(TjStuff& tjs, PFPStruct& pfp, std::vector<std::vector<unsigned int>>& sptLists, int tjID, bool prt);
  bool MergeBrokenTjs(TjStuff& tjs, std::vector<int>& tjInPln, bool prt);
  std::vector<int> TjsNearSpacePts(TjStuff& tjs, std::vector<unsigned int> sptlist);
  std::vector<unsigned int> SpacePtsAssociatedWith(TjStuff& tjs, const TrajPoint& tp);
  std::vector<unsigned int> SpacePtsAssociatedWith(TjStuff& tjs, const Trajectory& tj);
  std::vector<unsigned int> SpacePtsAtHit(TjStuff& tjs, unsigned int iht);
  void UpdateMatchStructs(TjStuff& tjs, std::vector<int> oldTjs, int newTj);
  void FillmAllTraj(TjStuff& tjs, const geo::TPCID& tpcid);
  void MakePFPSpacePts(TjStuff& tjs, PFPStruct& pfp, bool anyTj);
  bool SetNewStart(TjStuff& tjs, PFPStruct& pfp, bool prt);
  void SortByDistanceFromStart(TjStuff& tjs, PFPStruct& pfp);
  void CheckSptValidity(TjStuff& tjs, PFPStruct& pfp, bool prt);
  void FindXMatches(TjStuff& tjs, unsigned short numPlanes, short maxScore, std::vector<MatchStruct>& matVec, bool prt);
  bool MakeSpt(TjStuff& tjs, const TrajPoint& itp, const TrajPoint& jtp, Point3_t& pos, Vector3_t& dir);
  double DeltaAngle(const Vector3_t v1, const Vector3_t v2);
  inline double DotProd(const Vector3_t& v1, const Vector3_t& v2) {return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; }
  Vector3_t PointDirection(const Point3_t p1, const Point3_t p2);
  bool TrajPoint3D(TjStuff& tjs, const TrajPoint& itp, const TrajPoint& jtp, TVector3& pos, TVector3& dir, bool prt);
  
} // namespace tca

#endif // ifndef TRAJCLUSTERALGSPTUTILS_H
