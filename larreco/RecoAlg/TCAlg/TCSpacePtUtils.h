
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
  
} // namespace tca

#endif // ifndef TRAJCLUSTERALGSPTUTILS_H
