////////////////////////////////////////////////////////////////////////
//
//
// TCAlg shower code
//
// Bruce Baller
//
///////////////////////////////////////////////////////////////////////
#ifndef TRAJCLUSTERALGSHOWERS_H
#define TRAJCLUSTERALGSHOWERS_H


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
#include "larreco/RecoAlg/TCAlg/TCTruth.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

namespace tca {

  void Find3DShowerEndPoints(TjStuff& tjs, const geo::TPCID& tpcid);
  void MakeShowers(TjStuff& tjs);
  void FindShowers(TjStuff& tjs, const CTP_t& inCTP);
  bool FindShowers3D(TjStuff& tjs, const geo::TPCID& tpcid);
  void FindMatchingTjs(TjStuff& tjs, unsigned short cotIndex, bool prt);
  void Match2DShowers(TjStuff& tjs, const geo::TPCID& tpcid, bool prt);
  void FillPts(TjStuff& tjs, unsigned short cotIndex, bool prt);
  bool DefineShower(TjStuff& tjs, unsigned short cotIndex, bool prt);
  void MakeShowerObsolete(TjStuff& tjs, unsigned short cotIndex, bool prt);
  bool AddTj(TjStuff& tjs, int TjID, unsigned short cotIndex, bool doUpdate, bool prt);
  bool RemoveTj(TjStuff& tjs, int TjID, unsigned short cotIndex, bool doUpdate, bool prt);
  bool FindChargeCenter(TjStuff& tjs, unsigned short cotIndex, bool prt);
  void FindAngle(TjStuff& tjs, unsigned short cotIndex, bool prt);
  void FillRotPos(TjStuff& tjs, unsigned short cotIndex, bool prt);
  bool AnalyzeRotPos(TjStuff& tjs, unsigned short cotIndex, bool prt);
  bool DefineShowerTj(TjStuff& tjs, unsigned short cotIndex, bool prt);
  void FindExternalParent(TjStuff& tjs, unsigned short cotIndex, bool prt);
  void UpdateShowerWithParent(TjStuff& tjs, unsigned short cotIndex, bool prt);
  bool WrongSplitTj(TjStuff& tjs, Trajectory& tj, unsigned short tjEnd, ShowerStruct& ss, bool prt);
  float ParentFOM(TjStuff& tjs, Trajectory& tj, unsigned short& tjEnd, ShowerStruct& ss, bool prt);
  void DefineEnvelope(TjStuff& tjs, unsigned short cotIndex, bool prt);
  void AddTjsInsideEnvelope(TjStuff& tjs, unsigned short cotIndex, bool prt);
  void RefineShowerTj(TjStuff& tjs, unsigned short cotIndex, bool prt);
  bool AddLooseHits(TjStuff& tjs, unsigned short cotIndex, bool prt);
  void FindStartChg(TjStuff& tjs, unsigned short cotIndex, bool prt);
  std::vector<float> StartChgVec(TjStuff& tjs, unsigned short cotIndex, bool prt);
  void DumpShowerPts(TjStuff& tjs, unsigned short cotIndex);
  void CheckQuality(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  
  void TagShowerTjs(TjStuff& tjs, const CTP_t& inCTP, std::vector<std::vector<int>>& tjList);
  void AddCloseTjsToList(TjStuff& tjs, unsigned short itj, std::vector<int> list);
  void MergeTjList(std::vector<std::vector<int>>& tjList);
  void MergeTjList2(TjStuff& tjs, std::vector<std::vector<int>>& tjList, bool prt);
  void MergeNearby2DShowers(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeOverlap(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  void MergeSubShowers(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  bool MergeShowersAndStore(TjStuff& tjs, unsigned short icotIndex, unsigned short jcotIndex, bool prt);
  bool MergeShowerTjsAndStore(TjStuff& tjs, unsigned short istj, unsigned short jstj, bool prt);
  bool TransferTjHits(TjStuff& tjs, const CTP_t& inCTP, bool prt);
  unsigned short GetCotsIndex(TjStuff& tjs, unsigned short ShowerTjID);
  float ShowerEnergy(const TjStuff& tjs, const ShowerStruct& ss);
  unsigned short FarEnd(TjStuff& tjs, const Trajectory& tj, ShowerStruct& ss);
  unsigned short Create2DShower(TjStuff& tjs, const std::vector<int>& tjl);
  void FindNearbyTjs(TjStuff& tjs, unsigned short cotIndex, bool prt);

}


#endif // ifndef TRAJCLUSTERALGSHOWERS_H
